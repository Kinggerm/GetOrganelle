import os
import sys
import math
import re
major_version, minor_version = sys.version_info[:2]
if major_version == 2 and minor_version >= 7:
    python_version = "2.7+"
elif major_version == 3 and minor_version >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)


if python_version == "2.7+":
    # python2
    import string
    translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


    def complementary_seq(input_seq):
        return string.translate(input_seq, translator)[::-1]

else:
    # python3
    translator = str.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


    def complementary_seq(input_seq):
        return str.translate(input_seq, translator)[::-1]


def complementary_seqs(input_seq_iter):
    return tuple([complementary_seq(seq) for seq in input_seq_iter])


class Sequence:
    def __init__(self, label, seq):
        self.label = label
        self.seq = seq
        
    def __len__(self):
        return len(self.seq)

    def fasta_str(self, interleaved=False):
        out_str = []
        if interleaved:
            out_str.extend(['>', self.label, '\n'])
            j = interleaved
            while j < len(self):
                out_str.append(self.seq[(j - interleaved):j])
                out_str.append('\n')
                j += interleaved
            out_str.append(self.seq[(j - interleaved):j])
            out_str.append('\n')
        else:
            out_str = ['>', self.label, '\n', self.seq, "\n"]
        return "".join(out_str)


class SequenceList:
    def __init__(self, input_fasta_file=None):
        self.sequences = []
        self.interleaved = False
        if input_fasta_file:
            self.read_fasta(input_fasta_file)
        
    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        for seq in self.sequences:
            yield seq
    
    def append(self, sequence):
        self.sequences.append(sequence)
        
    def read_fasta(self, fasta_file):
        fasta_file = open(fasta_file, 'rU')
        this_line = fasta_file.readline()
        interleaved = 0
        while this_line:
            if this_line.startswith('>'):
                this_name = this_line[1:].strip()
                this_seq = ''
                this_line = fasta_file.readline()
                seq_line_count = 0
                while this_line and not this_line.startswith('>'):
                    if seq_line_count == 1:
                        interleaved = len(this_seq)
                    this_seq += this_line.strip()
                    this_line = fasta_file.readline()
                    seq_line_count += 1
                self.append(Sequence(this_name, this_seq))
            else:
                this_line = fasta_file.readline()
        fasta_file.close()
        self.interleaved = interleaved

    def write_fasta(self, fasta_file, overwrite=True):
        if not overwrite:
            while os.path.exists(fasta_file):
                fasta_file = '.'.join(fasta_file.split('.')[:-1]) + '_.' + fasta_file.split('.')[-1]
        fasta_file_handler = open(fasta_file, 'w')
        for seq in self:
            fasta_file_handler.write(seq.fasta_str(self.interleaved))
        fasta_file_handler.close()


def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    interleaved = 0
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = fasta_file.readline()
            seq_line_count = 0
            while this_line and not this_line.startswith('>'):
                if seq_line_count == 1:
                    interleaved = len(this_seq)
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
                seq_line_count += 1
            seqs.append(this_seq)
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs, interleaved]


def write_fasta(out_dir, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_dir):
            out_dir = '.'.join(out_dir.split('.')[:-1])+'_.'+out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
                j += matrix[2]
            fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
    else:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            fasta_file.write(matrix[1][i]+'\n')
    fasta_file.close()


def read_fasta_as_list(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    interleaved = 0
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = fasta_file.readline()
            seq_line_count = 0
            while this_line and not this_line.startswith('>'):
                if seq_line_count == 1:
                    interleaved = len(this_seq)
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
                seq_line_count += 1
            seqs.append(list(this_seq))
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs, interleaved]


def write_fasta_with_list(out_dir, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_dir):
            out_dir = '.'.join(out_dir.split('.')[:-1])+'_.'+out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(''.join(matrix[1][i][(j-matrix[2]):j])+'\n')
                j += matrix[2]
            fasta_file.write(''.join(matrix[1][i][(j-matrix[2]):j])+'\n')
    else:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            fasta_file.write(''.join(matrix[1][i])+'\n')
    fasta_file.close()


matching_char = {}
if python_version == "2.7+":
    for char in string.ascii_lowercase:
        matching_char[(char, char)] = 0
        matching_char[(char, char.upper())] = 0
        matching_char[(char.upper(), char)] = 0
        matching_char[(char.upper(), char.upper())] = 0
else:
    for char in "atgcrmykhbdvn-":
        matching_char[(char, char)] = 0
        matching_char[(char, char.upper())] = 0
        matching_char[(char.upper(), char)] = 0
        matching_char[(char.upper(), char.upper())] = 0


def find_string_difference(this_string, this_reference, dynamic_span=2.0):
    len_str = len(this_string)
    len_ref = len(this_reference)
    if dynamic_span == 0:
        difference = sum([not (this_string[i], this_reference[i]) in matching_char for i in range(min(len_ref, len_str))]) + abs(len_ref - len_str)
        proper_end = this_string[-1] == this_reference[-1]
        return difference, proper_end
    else:
        dynamic_span = max(abs(len(this_string)-len(this_reference))+1, dynamic_span)
        this_match = int(not (this_string[0], this_reference[0]) in matching_char)
        this_matrix = {(0, 0): {'state': this_match}}
        # calculate the first column
        for i in range(1, min(int(math.ceil(dynamic_span))+1, len_str)):
            this_matrix[(i, 0)] = {'right_out': this_match+i, 'state': this_match+i}
        # calculate the first line
        for j in range(1, min(int(math.ceil(dynamic_span))+1, len_ref)):
            this_matrix[(0, j)] = {'right_out': this_match+j, 'state': this_match+j}
        # calculate iteratively
        start = 0
        for i in range(1, len_str):
            start = max(1, int(i-dynamic_span))
            end = min(len_ref, int(math.ceil(i+dynamic_span)))
            # start: no right_in
            this_match = int(not (this_string[i], this_reference[start]) in matching_char)
            this_matrix[(i, start)] = {'diagonal_out': this_matrix[(i-1, start-1)]['state'] + this_match,
                                       'down_out': this_matrix[(i-1, start)]['state'] + 1}
            this_matrix[(i, start)]['state'] = min(this_matrix[(i, start)].values())
            # middle
            for j in range(start+1, end-1):
                this_match = not (this_string[i], this_reference[j]) in matching_char
                this_matrix[(i, j)] = {'diagonal_out': this_matrix[(i-1, j-1)]['state'] + this_match,
                                       'down_out': this_matrix[(i-1, j)]['state'] + 1,
                                       'right_out': this_matrix[(i, j-1)]['state'] + 1}
                this_matrix[(i, j)]['state'] = min(this_matrix[(i, j)].values())
            # end
            this_match = not (this_string[i], this_reference[end - 1]) in matching_char
            this_matrix[(i, end-1)] = {'diagonal_out': this_matrix[(i-1, end-2)]['state'] + this_match}
            if (i, end-2) in this_matrix:
                this_matrix[(i, end-1)]['right_out'] = this_matrix[(i, end-2)]['state'] + 1
            this_matrix[(i, end-1)]['state'] = min(this_matrix[(i, end-1)].values())
        # print time.time()-time0
        difference = this_matrix[(len_str-1, len_ref-1)]['state']
        proper_end = True
        for j in range(start, len_ref):
            try:
                if this_matrix[(len_str-1, j)]['state'] < difference:
                    proper_end = False
                    break
            except KeyError:
                pass
        for i in range(max(0, len_str-len_ref+start), len_str):
            try:
                if this_matrix[(i, len_ref-1)]['state'] < difference:
                    proper_end = False
                    break
            except KeyError:
                pass
        return difference, proper_end


DEGENERATE_BASES = {"N", "V", "H", "D", "B", "Y", "R", "K", "M", "S", "W"}


DEGENERATE_DICT = {  # degenerate
    "N": ["A", "C", "G", "T"],
    "V": ["A", "C", "G"], "H": ["A", "C", "T"], "D": ["A", "G", "T"], "B": ["C", "G", "T"],
    "Y": ["C", "T"], "R": ["A", "G"], "K": ["G", "T"], "M": ["A", "C"],
    "S": ["C", "G"], "W": ["A", "T"],
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
    #  consensus
    ('A', 'C', 'G', 'T'): 'N',
    ('A', 'C', 'G'): 'V', ('A', 'C', 'T'): 'H', ('A', 'G', 'T'): 'D', ('C', 'G', 'T'): 'B',
    ('C', 'T'): 'Y', ('A', 'G'): 'R', ('G', 'T'): 'K', ('A', 'C'): 'M',
    ('C', 'G'): 'S', ('A', 'T'): 'W',
    ('A',): 'A', ('C',): 'C', ('G',): 'G', ('T',): 'T'}

DEGENERATE_DICT_DIGIT = {  # degenerate
    "N": [1, 2, 4, 8],
    "V": [1, 2, 4], "H": [1, 2, 8], "D": [1, 4, 8], "B": [2, 4, 8],
    "Y": [2, 8], "R": [1, 4], "K": [4, 8], "M": [1, 2],
    "S": [2, 4], "W": [1, 8],
    "A": [1], "C": [2], "G": [4], "T": [8],
    #  consensus
    15: 'N',
    7: 'V', 11: 'H', 13: 'D', 14: 'B',
    10: 'Y', 5: 'R', 12: 'K', 3: 'M',
    6: 'S', 9: 'W',
    1: "A", 2: "C", 4: "G", 8: "T"}


def generate_consensus(*seq_str):
    consensus_res = []
    seq_num = len(seq_str)
    for go_to_base in range(len(seq_str[0])):
        this_base_set = set()
        for go_to_seq in range(seq_num):
            this_base_set.update(DEGENERATE_DICT_DIGIT.get(seq_str[go_to_seq][go_to_base], []))
        consensus_res.append(DEGENERATE_DICT_DIGIT[sum(this_base_set)])
    return "".join(consensus_res)


def split_seq_by_quality_pattern(sequence, quality_str, low_quality_pattern, min_sub_seq=0):
    seq_list = []
    start = 0
    for splicer in re.finditer(low_quality_pattern, quality_str):
        end, new_start = splicer.span()
        seq_part = sequence[start:end]
        if len(seq_part) > min_sub_seq:
            seq_list.append(seq_part)
        start = new_start
    seq_part = sequence[start:]
    if len(seq_part) > min_sub_seq:
        seq_list.append(seq_part)
    return tuple(seq_list)


def fq_seq_simple_generator(fq_dir_list, go_to_line=1, split_pattern=None, min_sub_seq=0, max_n_reads=float("inf")):
    if not ((type(fq_dir_list) is list) or (type(fq_dir_list) is tuple)):
        fq_dir_list = [fq_dir_list]
    max_n_lines = 4 * max_n_reads
    if split_pattern and len(split_pattern) > 2:
        for fq_dir in fq_dir_list:
            count = 0
            with open(fq_dir, 'rU') as fq_handler:
                seq_line = fq_handler.readline()
                while seq_line:
                    if count % 4 == go_to_line:
                        fq_handler.readline()
                        quality_str = fq_handler.readline()[:-1]
                        count += 2
                        yield split_seq_by_quality_pattern(seq_line[:-1], quality_str, split_pattern, min_sub_seq)
                    count += 1
                    if count >= max_n_lines:
                        break
                    seq_line = fq_handler.readline()
    else:
        for fq_dir in fq_dir_list:
            count = 0
            with open(fq_dir, 'rU') as fq_handler:
                for fq_line in fq_handler:
                    if count % 4 == go_to_line:
                        yield fq_line[:-1]
                    if count >= max_n_lines:
                        break
                    count += 1


def chop_seqs(seq_generator_or_list, word_size):
    return_words = set()
    for seed in seq_generator_or_list:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            for i in range(0, this_seq_len - word_size + 1):
                forward = seed[i:i + word_size]
                return_words.add(forward)
                reverse = cpt_seed[i:i + word_size]
                return_words.add(reverse)
    return return_words


def chop_seqs_as_empty_dict(seq_generator_or_list, word_size):
    return_words = dict()
    for seed in seq_generator_or_list:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            for i in range(0, this_seq_len - word_size + 1):
                forward = seed[i:i + word_size]
                return_words[forward] = 0
                reverse = cpt_seed[i:i + word_size]
                return_words[reverse] = 0
    return return_words


def chop_seq_list(seq_generator_or_list, word_size):
    return_words = set()
    for seed in seq_generator_or_list:
        for seq_part in seed:
            this_seq_len = len(seq_part)
            if this_seq_len >= word_size:
                cpt_seed = complementary_seq(seq_part)
                for i in range(0, this_seq_len - word_size + 1):
                    forward = seq_part[i:i + word_size]
                    return_words.add(forward)
                    reverse = cpt_seed[i:i + word_size]
                    return_words.add(reverse)
    return return_words


def counting_words(seq_generator, words_initial_dict, word_size):
    for seq in seq_generator:
        for i in range(0, len(seq) - word_size + 1):
            this_word = seq[i: i + word_size]
            if this_word in words_initial_dict:
                words_initial_dict[this_word] += 1
    return words_initial_dict


def check_fasta_seq_names(original_fas, log=None):
    fas_matrix = read_fasta(original_fas)
    short_names = [s_n.split(" ")[0] for s_n in fas_matrix[0]]
    if len(short_names) == len(set(short_names)):
        return original_fas
    else:
        new_fas = original_fas + ".modified"
        if os.path.exists(new_fas):
            return check_fasta_seq_names(new_fas)
        else:
            if log:
                log.info("Setting '-s " + new_fas + "'")
            for go_to, name in enumerate(fas_matrix[0]):
                fas_matrix[0][go_to] = name.replace(" ", "_")
            if len(fas_matrix[0]) != len(set(fas_matrix[0])):
                existed = {}
                for go_to, name in enumerate(fas_matrix[0]):
                    if name in existed:
                        existed[name] += 1
                    else:
                        existed[name] = 1
                    fas_matrix[0][go_to] += "-" + str(existed[name])
            write_fasta(new_fas, fas_matrix, True)
            return new_fas


def get_read_len_mean_max_count(fq_files, maximum_n_reads, sampling_percent=1.):
    read_lengths = []
    all_count = 0
    if sampling_percent == 1:
        for fq_f in fq_files:
            count_r = 0
            for seq in fq_seq_simple_generator(fq_f):
                count_r += 1
                read_lengths.append(len(seq.strip("N")))
                if count_r >= maximum_n_reads:
                    break
            all_count += count_r
    else:
        sampling_percent = int(1 / sampling_percent)
        for fq_f in fq_files:
            count_r = 0
            for seq in fq_seq_simple_generator(fq_f):
                count_r += 1
                if count_r % sampling_percent == 0:
                    read_lengths.append(len(seq.strip("N")))
                if count_r >= maximum_n_reads:
                    break
            all_count += count_r
    return sum(read_lengths)/len(read_lengths), max(read_lengths), all_count


def get_read_quality_info(fq_files, maximum_n_reads, min_quality_score, log,
                          maximum_ignore_percent=0.05, sampling_percent=0.1):
    sampling_percent = int(1 / sampling_percent)
    all_quality_chars_list = []
    record_fq_beyond_read_num_limit = []
    for fq_f in fq_files:
        count_r = 0
        record_fq_beyond_read_num_limit.append(False)
        this_fq_generator = fq_seq_simple_generator(fq_f, go_to_line=3)
        for quality_str in this_fq_generator:
            if count_r % sampling_percent == 0:
                all_quality_chars_list.append(quality_str)
            count_r += 1
            if count_r >= maximum_n_reads:
                break
        for quality_str in this_fq_generator:
            if quality_str:
                log.info("Number of reads exceeded " + str(int(maximum_n_reads)) + " in " + os.path.basename(fq_f)
                         + ", only top " + str(int(maximum_n_reads))
                         + " reads are used in downstream analysis.")
                record_fq_beyond_read_num_limit[-1] = True
            break
    all_quality_chars = "".join(all_quality_chars_list)
    len_quality_chars_total = float(len(all_quality_chars))
    max_quality = max(all_quality_chars)
    min_quality = min(all_quality_chars)
    max_quality = ord(max_quality)
    min_quality = ord(min_quality)
    decision_making = []
    for type_name, char_min, char_max, score_min, score_max in [("Sanger", 33, 73, 0, 40),
                                                                ("Solexa", 59, 104, -5, 40),
                                                                ("Illumina 1.3+", 64, 104, 0, 40),
                                                                ("Illumina 1.5+", 67, 105, 3, 41),
                                                                ("Illumina 1.8+", 33, 74, 0, 41)]:
        decision_making.append((type_name, char_min, char_max, score_min, score_max,
                                (max_quality - char_max) ** 2 + (min_quality - char_min) ** 2))
    the_form, the_c_min, the_c_max, the_s_min, the_s_max, deviation = sorted(decision_making, key=lambda x: x[-1])[0]
    log.info("Identified quality encoding format = " + the_form)
    if max_quality > the_c_max:
        log.warning("Max quality score " + repr(chr(max_quality)) +
                    "(" + str(max_quality) + ":" + str(max_quality - the_c_min + the_s_min) +
                    ") in your fastq file exceeds the usual boundary " + str((the_c_min, the_c_max)))
    if min_quality < the_c_min:
        log.warning("Min quality score " + repr(chr(min_quality)) +
                    "(" + str(min_quality) + ":" + str(min_quality - the_c_min + the_s_min) +
                    ") in your fastq file is under the usual lower boundary " + str((the_c_min, the_c_max)))

    # increase --min-quality-score if ignoring too much data.
    low_quality_score = min(the_c_min, min_quality)
    ignore_percent = 0
    trimmed_quality_chars = []

    while low_quality_score < the_c_min + min_quality_score - the_s_min:
        ignore_this_time = all_quality_chars.count(chr(low_quality_score)) / len_quality_chars_total
        ignore_percent += ignore_this_time
        if ignore_percent > maximum_ignore_percent:
            ignore_percent -= ignore_this_time
            break
        else:
            trimmed_quality_chars.append(chr(low_quality_score))
        low_quality_score += 1
    if low_quality_score < the_c_min + min_quality_score - the_s_min:
        log.info("Resetting '--min-quality-score " + str(low_quality_score + the_s_min - the_c_min) + "'")
    if trimmed_quality_chars:
        log.info("Trimming bases with qualities (" + "%.2f" % (ignore_percent * 100) + "%): " +
                 str(ord(min("".join(trimmed_quality_chars)))) + ".." + str(ord(max("".join(trimmed_quality_chars)))) +
                 "  " + "".join(trimmed_quality_chars))
    trimmed_quality_chars = "".join(trimmed_quality_chars)

    # calculate mean error rate
    all_quality_char_dict = {in_quality_char: all_quality_chars.count(in_quality_char)
                             for in_quality_char in set(all_quality_chars)}
    error_prob_func = chose_error_prob_func[the_form]
    mean_error_rate = sum([error_prob_func(in_quality_char) * all_quality_char_dict[in_quality_char]
                           for in_quality_char in all_quality_char_dict]) / len_quality_chars_total
    log.info("Mean error rate = " + str(round(mean_error_rate, 4)))

    return "[" + trimmed_quality_chars + "]", mean_error_rate, record_fq_beyond_read_num_limit  # , post_trimming_mean


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
    return round(mean - down_dev, 2), round(mean, 2), round(mean + up_dev, 2)


# Q = -10*log10(P)
def score_2_error_prob_sanger_quality(quality_score):
    return 10**(-quality_score / 10.)


# Q = -10*log10(P/(1-P))
def score_2_error_prob_solexa_quality(quality_score):
    this_sanger_prob = score_2_error_prob_sanger_quality(quality_score)
    return this_sanger_prob/(1 + this_sanger_prob)


# solexa_quality_str = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
all_quality_str = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
phred_33_quality_trans = {in_trans_char: ord(in_trans_char) - 33
                          for in_trans_char in all_quality_str}
phred_64_quality_trans = {in_trans_char: ord(in_trans_char) - 64
                          for in_trans_char in all_quality_str}


def sanger_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_33_quality_trans.get(quality_char, 0))


def solexa_error_prob(quality_char):
    return score_2_error_prob_solexa_quality(phred_64_quality_trans.get(quality_char, -31))


def illumina_1_3_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_64_quality_trans.get(quality_char, -31))


def illumina_1_5_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_64_quality_trans.get(quality_char, -31))


def illumina_1_8_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_33_quality_trans.get(quality_char, 0))


chose_error_prob_func = {"Sanger": sanger_error_prob,
                         "Solexa": solexa_error_prob,
                         "Illumina 1.3+": illumina_1_3_error_prob,
                         "Illumina 1.5+": illumina_1_5_error_prob,
                         "Illumina 1.8+": illumina_1_8_error_prob}