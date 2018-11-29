from Library.statistical_func import *
import re
import time


CIGAR_ALPHA_REG = "([MIDNSHPX=])"
CIGAR_EQUAL = {"M", "="}
CIGAR_CONSTANT = {"X"}
# not in read
CIGAR_DEL = {"D", "N", "H"}
# not in ref
CIGAR_ADD = {"I", "S", "P"}


class MapRecord:
    def __init__(self, ref, position, cigar, ref_next, pos_next, lib_len, seq, seq_quality, is_paired, is_all_mapped,
                 is_mapped, is_next_mapped, is_reversed, is_next_reversed, is_first, optional_fields, ref_len):
        self.ref = ref
        self.ref_left_by_alignment = position
        self.cigar_str = cigar
        self.cigar = []
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

        self.ref_right_by_alignment = position + self.seq_len - 1
        for op_len, operation in self.cigar:
            if operation in CIGAR_ADD:
                self.ref_right_by_alignment -= op_len
            elif operation in CIGAR_DEL:
                self.ref_right_by_alignment += op_len

        self.query_start_align = 1
        self.query_end_align = self.seq_len

        self.left_to_end = 0
        self.right_to_end = 0
        self.min_dist_to_end = 0
        # time_test[0] -= time.time()
        self.parse_cigar()
        self.get_ref_end_info(ref_len=ref_len)
        # time_test[0] += time.time()
        
    def parse_cigar(self):
        cigar_split = re.split(CIGAR_ALPHA_REG, self.cigar_str)[:-1]  # ref order
        self.cigar = []
        for go_part in range(0, len(cigar_split), 2):
            self.cigar.append((int(cigar_split[go_part]), cigar_split[go_part + 1]))
        
    def get_ref_end_info(self, ref_len):

        # possible extension caused by false-insertion
        left_potential_extension = 0
        right_potential_extension = 0
        for op_len, operation in self.cigar:
            if operation in CIGAR_EQUAL:
                # longer than 12 are thought to be conserved ?
                if op_len > 12:
                    break
            elif operation in CIGAR_DEL:
                left_potential_extension -= op_len
            elif operation in CIGAR_ADD:
                left_potential_extension += op_len
        for op_len, operation in self.cigar[::-1]:
            if operation in CIGAR_EQUAL:
                if op_len > 12:
                    break
            elif operation in CIGAR_DEL:
                right_potential_extension -= op_len
            elif operation in CIGAR_ADD:
                right_potential_extension += op_len
        
        if self.lib_len != 0:
            if self.is_reverse_complementary:
                left_end_case_1 = (self.ref_left_by_alignment + self.seq_len + self.lib_len - left_potential_extension) - 1
                left_end_case_2 = (self.pos_next - left_potential_extension) - 1
                right_end_case_1 = ref_len - (self.ref_left_by_alignment + self.seq_len - 1 + right_potential_extension)
                right_end_case_2 = ref_len - (self.pos_next - (self.lib_len + 1) + right_potential_extension)
            else:
                left_end_case_1 = (self.ref_left_by_alignment - left_potential_extension) - 1
                left_end_case_2 = (self.pos_next + self.seq_len - self.lib_len - left_potential_extension) - 1
                right_end_case_1 = ref_len - (self.ref_left_by_alignment + (self.lib_len - 1) + right_potential_extension)
                right_end_case_2 = ref_len - (self.pos_next + (self.seq_len - 1) + right_potential_extension)
            self.left_to_end = max(min(left_end_case_1, left_end_case_2), 0)
            self.right_to_end = max(min(right_end_case_1, right_end_case_2), 0)
            self.min_dist_to_end = min(self.left_to_end, self.right_to_end)
        else:
            self.left_to_end = max((self.ref_left_by_alignment - left_potential_extension) - 1, 0)
            self.right_to_end = max(ref_len - (self.ref_left_by_alignment + (self.seq_len - 1) + right_potential_extension), 0)
            if self.left_to_end == 0 or self.right_to_end == 0:
                self.min_dist_to_end = 0
            else:
                if self.is_reverse_complementary:
                    self.min_dist_to_end = self.left_to_end
                else:
                    self.min_dist_to_end = self.right_to_end

    def get_query_end_info(self):
        if self.is_reverse_complementary:
            for op_len, operation in self.cigar[::-1]:
                if operation in CIGAR_EQUAL:
                    break
                elif operation in CIGAR_ADD or operation in CIGAR_CONSTANT:
                    self.query_start_align += op_len
            for op_len, operation in self.cigar:
                if operation in CIGAR_EQUAL:
                    break
                elif operation in CIGAR_ADD or operation in CIGAR_CONSTANT:
                    self.query_end_align -= op_len
        else:
            for op_len, operation in self.cigar:
                if operation in CIGAR_EQUAL:
                    break
                elif operation in CIGAR_ADD or operation in CIGAR_CONSTANT:
                    self.query_start_align += op_len
            for op_len, operation in self.cigar[::-1]:
                if operation in CIGAR_EQUAL:
                    break
                elif operation in CIGAR_ADD or operation in CIGAR_CONSTANT:
                    self.query_end_align -= op_len


class ReadLibrary:
    def __init__(self, sam_file):
        self.sam_file = sam_file
        self.references = {}
        self.queries = {}
        self.library_len = 0
        self.lib_mean = 0
        self.lib_std = 0
        self.remains = 0
        self.dropped = 0
        self.parse_sam(sam_file)

    def parse_sam(self, sam_file):

        candidate_library_lengths = []
        with open(sam_file, "r") as open_sam:
            for aligned_line in open_sam:
                if aligned_line.startswith("@SQ\tSN"):
                    ref_name, ref_len = aligned_line.strip().split("\t")[1:3]
                    self.references[ref_name.split(":")[-1]] = {"len": int(ref_len.split(":")[-1])}
                elif not aligned_line.startswith("@"):
                    aligned_split = aligned_line.strip().split("\t")
                    q_name, flag, ref_name, pos, map_quality, cigar, ref_next, pos_next, lib_len, seq, seq_quality = \
                        aligned_split[:11]
                    optional_fields = {}
                    for flag_type_val in aligned_split[11:]:
                        op_flag, op_type, op_val = flag_type_val.split(":")
                        if op_type == "i":
                            optional_fields[op_flag] = int(op_val)
                    flag_fields = [bool(int(bit)) for bit in '{:012b}'.format(int(flag))][::-1]
                    is_paired, is_all_mp, not_mp, not_next_mp, is_reversed, is_next_reversed, is_first, is_last = \
                        flag_fields[:8]
                    if (q_name, is_first) not in self.queries:
                        self.queries[(q_name, is_first)] = []
                    # if lib_len == "40":
                    #     print(q_name, is_first)
                    #     print(aligned_line)
                    #     print(aligned_split)
                    self.queries[(q_name, is_first)].append(
                        MapRecord(ref=ref_name,
                                  position=int(pos),
                                  cigar=cigar,
                                  ref_next=ref_next,
                                  pos_next=int(pos_next),
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
                                  ref_len=self.references[ref_name]["len"]))

                    # all segments are mapped to the same reference
                    if int(lib_len) > 0:
                        candidate_library_lengths.append(int(lib_len))

        self.lib_mean, self.lib_std = mean_and_std_pauta_criterion(candidate_library_lengths, SIGMA_MULTIPLE)

    def remove_records_in_the_middle(self):
        # libs outside this scope would not be took into consideration during scaffolding
        max_dev_with_lib_across_end = 2 * SIGMA_MULTIPLE * self.lib_std - 2
        for q_name_with_direction in list(self.queries):
            go_to_record = 0
            while go_to_record < len(self.queries[q_name_with_direction]):
                map_record = self.queries[q_name_with_direction][go_to_record]
                # map_record = MapRecord()
                if map_record.min_dist_to_end < max_dev_with_lib_across_end:
                    self.remains += 1
                    go_to_record += 1
                else:
                    del self.queries[q_name_with_direction][go_to_record]
                    self.dropped += 1
            if not self.queries[q_name_with_direction]:
                del self.queries[q_name_with_direction]


# assembly = Assembly("/Users/Kinggerm/Documents/ResearchLife3/Teamworks/Zhouzhuo/2018-10-18-cp-assembly"
#                     "/plastome-300/K125/target.graph1.selected_graph.gfa")
# a = Scaffolding("/Users/Kinggerm/Documents/ResearchLife3/Teamworks/Zhouzhuo/2018-10-18-cp-assembly/800-a-n.sam")
# print(len(cigars))
# print(sorted(cigars))
# print(set(str(sorted(cigars))))




# read_len = 150
# kmer = 125
# for vertex in assembly.vertex_info:
#     for this_end in (False, True):

