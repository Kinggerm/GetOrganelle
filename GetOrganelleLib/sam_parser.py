# from GetOrganelleLib.statistical_func import *
import re
import sys
import time


DIGIT_REG = "([0-9]{,})"
CIGAR_ALPHA_REG = "([MIDNSHPX=])"
CIGAR_ALPHA = "MIDNSHPX="
# for counting length
CIGAR_LEN = {"M", "I", "S", "=", "X"}
# for counting coverages/query
# both in ref and read;
CIGAR_GO_BOTH = {"M", "=", "X"}
# in ref and not in read;
CIGAR_GO_REF_PAUSE_READ = {"D", "N"}
# not in ref
CIGAR_PAUSE_REF = {"I", "S", "H", "P"}
# in ref and not in read;
CIGAR_GO_READ_PAUSE_REF = {"I", "S"}
# not in ref
CIGAR_PAUSE_READ = {"D", "N", "H", "P"}
#
CIGAR_GO_REF = {"M", "=", "X", "D", "N"}
#
CIGAR_GO_READ = {"M", "=", "X", "I", "S"}


class MapRecord(object):
    def __init__(self, ref, position, map_quality, cigar, ref_next, pos_next, lib_len, seq, seq_quality, is_paired,
                 is_all_mapped, is_mapped, is_next_mapped, is_reversed, is_next_reversed, is_first, optional_fields,
                 ref_len):
        self.ref = ref
        self.ref_left_by_alignment = position
        self.map_quality = map_quality
        self.cigar_str = cigar
        self.cigar = split_cigar_str(self.cigar_str)
        self.ref_next = ref_next
        self.pos_next = pos_next
        self.lib_len = lib_len
        self.seq = seq
        self.seq_len = len(seq)
        self.seq_quality = seq_quality
        self.is_paired = is_paired
        self.is_all_mapped = is_all_mapped
        self.is_mapped = is_mapped
        self.is_mapped_next = is_next_mapped
        self.is_reverse_complementary = is_reversed
        self.is_reverse_complementary_next = is_next_reversed
        self.is_first = is_first
        self.optional_fields = optional_fields
        self.ref_len = ref_len

        self.query_start_align = 1
        self.query_end_align = self.seq_len
        self.ref_right_by_alignment = None
        self.left_to_end = 0
        self.right_to_end = 0
        self.min_dist_to_end = 0


class MapRecords(object):
    def __init__(self, sam_file, ref_real_len_dict=None, log_handler=None):
        self.sam_file = sam_file
        self.references = {}
        self.queries = {}
        self.library_len = 0
        self.lib_mean = 0
        self.lib_std = 0
        self.remains = 0
        self.dropped = 0
        self.parse_sam(sam_file, ref_real_len_dict=ref_real_len_dict, log_handler=log_handler)
        self.coverages = {}

    def parse_sam(self, sam_file, ref_real_len_dict=None, log_handler=None):

        candidate_library_lengths = []
        if ref_real_len_dict:
            for ref_name in ref_real_len_dict:
                self.references[ref_name] = {"index_len": ref_real_len_dict[ref_name],
                                             "real_len": ref_real_len_dict[ref_name]}
        with open(sam_file, "r") as open_sam:
            for aligned_line in open_sam:
                if aligned_line.startswith("@SQ\tSN"):
                    ref_name, ref_len = aligned_line.strip().split("\t")[1:3]
                    ref_name = ":".join(ref_name.split(":")[1:])
                    if ref_real_len_dict:
                        if ref_name in ref_real_len_dict:
                            self.references[ref_name] = {"index_len": int(ref_len.split(":")[-1]),
                                                         "real_len": ref_real_len_dict[ref_name]}
                        else:
                            if log_handler:
                                log_handler.warning(ref_name + " not found in given fasta!")
                            else:
                                sys.stdout.write("Warning: " + ref_name + " not found in given fasta!")
                            self.references[ref_name] = {"index_len": int(ref_len.split(":")[-1]),
                                                         "real_len": int(ref_len.split(":")[-1])}
                    else:
                        self.references[ref_name] = {"index_len": int(ref_len.split(":")[-1]),
                                                     "real_len": int(ref_len.split(":")[-1])}
                elif not aligned_line.startswith("@"):
                    aligned_split = aligned_line.strip().split("\t")
                    q_name, flag, ref_name, pos, map_quality, cigar, ref_next, pos_next, lib_len, seq, seq_quality = \
                        aligned_split[:11]
                    optional_fields = {}
                    for flag_type_val in aligned_split[11:]:
                        op_flag, op_type, op_val = flag_type_val.split(":")
                        if op_type == "i":
                            optional_fields[op_flag] = int(op_val)
                        elif op_type == "Z":
                            optional_fields[op_flag] = op_val
                    flag_fields = [bool(int(bit)) for bit in '{:012b}'.format(int(flag))][::-1]
                    is_paired, is_all_mp, not_mp, not_next_mp, is_reversed, is_next_reversed, is_first, is_last = \
                        flag_fields[:8]
                    if (q_name, is_first) not in self.queries:
                        self.queries[(q_name, is_first)] = []
                    if ref_real_len_dict:
                        pos = (int(pos) - 1) % ref_real_len_dict[ref_name] + 1
                        pos_next = (int(pos_next) - 1) % ref_real_len_dict[ref_name] + 1
                    else:
                        pos = int(pos)
                        pos_next = int(pos_next)
                    self.queries[(q_name, is_first)].append(
                        MapRecord(ref=ref_name,
                                  position=pos,
                                  map_quality=int(map_quality),
                                  cigar=cigar,
                                  ref_next=ref_next,
                                  pos_next=pos_next,
                                  lib_len=int(lib_len),
                                  seq=seq,
                                  seq_quality=seq_quality,
                                  is_paired=is_paired,
                                  is_all_mapped=is_all_mp,
                                  is_mapped=not not_mp,
                                  is_next_mapped=not not_next_mp,
                                  is_reversed=is_reversed,
                                  is_next_reversed=is_next_reversed,
                                  is_first=is_first,
                                  optional_fields=optional_fields,
                                  ref_len=self.references[ref_name]["real_len"]))

                    # all segments are mapped to the same reference
                    if int(lib_len) > 0:
                        candidate_library_lengths.append(int(lib_len))

    def update_coverages(self, multiple_hits_mode="best"):
        self.coverages = {ref: [0 for foo in range(self.references[ref]["index_len"])] for ref in self.references}
        if multiple_hits_mode == "best":
            for query_tuple in self.queries:
                record_of_a_read = sorted(self.queries[query_tuple], key=lambda x: -x.map_quality)[0]
                if record_of_a_read.ref_left_by_alignment > 0:  # only mapped reads are considered
                    reference = record_of_a_read.ref
                    go_to_base = record_of_a_read.ref_left_by_alignment
                    for op_len, operation in record_of_a_read.cigar:
                        if operation in CIGAR_GO_BOTH:
                            for add_go in range(op_len):
                                # list is zero-based
                                # cigar is one-based
                                position_1based = go_to_base + add_go
                                self.coverages[reference][position_1based - 1] += 1
                            go_to_base += op_len
                        elif operation in CIGAR_GO_REF_PAUSE_READ:
                            go_to_base += op_len
        elif multiple_hits_mode == "all":
            for query_tuple in self.queries:
                for record_of_a_read in self.queries[query_tuple]:
                    if record_of_a_read.ref_left_by_alignment > 0:
                        reference = record_of_a_read.ref
                        go_to_base = record_of_a_read.ref_left_by_alignment
                        for op_len, operation in record_of_a_read.cigar:
                            if operation in CIGAR_GO_BOTH:
                                for add_go in range(op_len):
                                    # list is zero-based
                                    # cigar is one-based
                                    position_1based = go_to_base + add_go
                                    self.coverages[reference][position_1based - 1] += 1
                                go_to_base += op_len
                            elif operation in CIGAR_GO_REF_PAUSE_READ:
                                go_to_base += op_len
        # circularize
        for ref in self.references:
            if self.references[ref]["index_len"] > self.references[ref]["real_len"]:
                this_real_len = self.references[ref]["real_len"]
                for merge_pos_0based in range(this_real_len, self.references[ref]["index_len"]):
                    self.coverages[ref][merge_pos_0based % this_real_len] += self.coverages[ref][merge_pos_0based]
                del self.coverages[ref][this_real_len:]

    def get_number_of_mapped_reads(self):
        paired_reads_mapped = 0
        single_end_mapped = 0
        for query_key in self.queries:
            is_this_mapped = False
            is_next_mapped = False
            for map_info in self.queries[query_key]:
                is_this_mapped |= map_info.is_mapped
                is_next_mapped |= map_info.is_mapped_next
            if is_this_mapped and is_next_mapped:
                paired_reads_mapped += 1
            elif is_this_mapped:
                single_end_mapped += 1
        return {"paired": paired_reads_mapped, "single": single_end_mapped}

    def get_customized_mapping_characteristics(self, cigar_char_list=CIGAR_ALPHA, multiple_hits_mode="best"):
        mapping_statistics = {cigar_char: {ref: [0 for foo in range(self.references[ref]["index_len"])]
                                           for ref in self.references}
                              for cigar_char in CIGAR_ALPHA}
        if multiple_hits_mode == "best":
            for query_tuple in self.queries:
                record_of_a_read = sorted(self.queries[query_tuple], key=lambda x: -x.map_quality)[0]
                if record_of_a_read.ref_left_by_alignment > 0:
                    reference = record_of_a_read.ref
                    go_to_base = record_of_a_read.ref_left_by_alignment
                    for op_len, operation in record_of_a_read.cigar:
                        if operation in CIGAR_GO_REF:
                            for add_go in range(op_len):
                                # list is zero-based
                                # cigar is one-based
                                position_1based = go_to_base + add_go
                                mapping_statistics[operation][reference][position_1based - 1] += 1
                            go_to_base += op_len
                        else:
                            mapping_statistics[operation][reference][go_to_base] += op_len
        elif multiple_hits_mode == "all":
            for query_tuple in self.queries:
                for record_of_a_read in self.queries[query_tuple]:
                    if record_of_a_read.ref_left_by_alignment > 0:
                        reference = record_of_a_read.ref
                        go_to_base = record_of_a_read.ref_left_by_alignment
                        for op_len, operation in record_of_a_read.cigar:
                            if operation in CIGAR_GO_REF:
                                for add_go in range(op_len):
                                    # list is zero-based
                                    # cigar is one-based
                                    position_1based = go_to_base + add_go
                                    mapping_statistics[operation][reference][position_1based - 1] += 1
                                go_to_base += op_len
                            else:
                                mapping_statistics[operation][reference][go_to_base] += op_len
        # PART1
        # remove redundant chars
        # for cigar_char in cigar_char_list:
        #     if sum([sum(mapping_statistics[cigar_char][ref]) for ref in mapping_statistics[cigar_char]]) == 0:
        #         del mapping_statistics[cigar_char]
        # for cigar_char in list(mapping_statistics):
        #     if cigar_char not in cigar_char_list:
        #         del mapping_statistics[cigar_char]

        # read mismatches from MD
        # activated only when PART1 is activated
        # if "X" in cigar_char_list and "X" not in mapping_statistics:
        #     mapping_statistics["X"] = {ref: [0 for foo in range(self.references[ref]["index_len"])]
        #                                for ref in self.references}
        if multiple_hits_mode == "best":
            for query_tuple in self.queries:
                record_of_a_read = sorted(self.queries[query_tuple], key=lambda x: -x.map_quality)[0]
                if record_of_a_read.ref_left_by_alignment > 0:
                    reference = record_of_a_read.ref
                    go_to_base = record_of_a_read.ref_left_by_alignment
                    for md_char in record_of_a_read.optional_fields["MD"]:
                        if md_char.isdigit():
                            go_to_base += int(md_char)
                        elif md_char.startswith("^"):
                            go_to_base += len(md_char) - 1
                        else:
                            for add_go in range(len(md_char)):
                                # list is zero-based
                                # cigar is one-based
                                position_1based = go_to_base + add_go
                                mapping_statistics["X"][reference][position_1based - 1] += 1
        elif multiple_hits_mode == "all":
            for query_tuple in self.queries:
                for record_of_a_read in self.queries[query_tuple]:
                    if record_of_a_read.ref_left_by_alignment > 0:
                        reference = record_of_a_read.ref
                        go_to_base = record_of_a_read.ref_left_by_alignment
                        for md_char in record_of_a_read.optional_fields["MD"]:
                            if md_char.isdigit():
                                go_to_base += int(md_char)
                            elif md_char.startswith("^"):
                                go_to_base += len(md_char) - 1
                            else:
                                for add_go in range(len(md_char)):
                                    # list is zero-based
                                    # cigar is one-based
                                    position_1based = go_to_base + add_go
                                    mapping_statistics["X"][reference][position_1based - 1] += 1
        # circularize
        for cigar_char in mapping_statistics:
            for ref in self.references:
                if self.references[ref]["index_len"] > self.references[ref]["real_len"]:
                    this_real_len = self.references[ref]["real_len"]
                    for merge_pos_0based in range(this_real_len, self.references[ref]["index_len"]):
                        mapping_statistics[cigar_char][ref][merge_pos_0based % this_real_len] += \
                            mapping_statistics[cigar_char][ref][merge_pos_0based]
                    del mapping_statistics[cigar_char][ref][this_real_len:]
        return mapping_statistics


def split_md_str(md_str):
    return re.split(DIGIT_REG, md_str)[1:-1]  # empty start and end


def split_cigar_str(cigar_str):
    cigar_split = re.split(CIGAR_ALPHA_REG, cigar_str)[:-1]  # empty end
    cigar_list = []
    for go_part in range(0, len(cigar_split), 2):
        cigar_list.append((int(cigar_split[go_part]), cigar_split[go_part + 1]))
    return cigar_list


def get_heads_from_sam_fast(*bowtie_sam_files):
    hit_heads = set()
    for bt_s_f in bowtie_sam_files:
        for line in open(bt_s_f, 'r'):
            if line.strip() and not line.startswith('@'):
                line_split = line.strip().split('\t')
                this_name = line_split[0]
                flag = int(line_split[1])
                direction = 1 if flag % 32 < 16 else 2
                hit_heads.add((this_name, direction))
    return hit_heads


def mapping_gap_percent_gap_lengths_from_list(coverage_list, circular=False):
    gaps = []
    gap_continued = False
    for this_cov in coverage_list:
        if not this_cov:
            if gap_continued:
                gaps[-1] += 1
            else:
                gaps.append(1)
                gap_continued = True
        else:
            gap_continued = False
    if circular:
        first_is_gap = coverage_list[0] == 0
        last_is_gap = coverage_list[1] == 0
        if first_is_gap and last_is_gap and len(gaps) > 1:
            gaps[0] += gaps.pop(-1)
    gaps.sort(reverse=True)
    gap_percentage = sum(gaps) / float(len(coverage_list))
    return gap_percentage, gaps


def mapping_gap_info_from_coverage_dict(coverage_info, circular=False):
    gap_per_ref = []
    for ref_seq in coverage_info:
        gap_percent, gaps = mapping_gap_percent_gap_lengths_from_list(coverage_info[ref_seq], circular=circular)
        gap_per_ref.append([ref_seq, gap_percent, gaps, sum(coverage_info[ref_seq])])
    gap_per_ref.sort(key=lambda x: (-x[-1], x[1], -sum(x[2][:5])))
    return gap_per_ref[0][:3]


def get_cover_range(all_coverages, guessing_percent=0.07):
    # empirical value
    center_p = 1 - guessing_percent
    all_coverages.sort()
    total_len = len(all_coverages)
    start_cov_pos = max(int((center_p - guessing_percent / 2) * total_len), 0)
    end_cov_pos = min(int((center_p + guessing_percent / 2) * total_len), total_len)
    used_len = end_cov_pos-start_cov_pos
    used_points = all_coverages[start_cov_pos:end_cov_pos]
    mean = sum(used_points) / used_len
    up_dev = (sum([(single_p-mean)**2 for single_p in used_points if single_p >= mean])/used_len) ** 0.5
    down_dev = (sum([(single_p-mean)**2 for single_p in used_points if single_p <= mean])/used_len) ** 0.5
    return round(max(mean - down_dev, min(used_points)), 2), round(mean, 2), round(min(mean + up_dev, max(used_points)), 2)
