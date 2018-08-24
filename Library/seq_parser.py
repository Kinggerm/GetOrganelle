import os
import sys
import math
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
            names.append(this_line[1:].strip('\n').strip('\r'))
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


degenerate_dict = {  # degenerate
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

degenerate_dict_digit = {  # degenerate
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


# def consensus_seq(*seq_str):
#     consensus_res = []
#     seq_num = len(seq_str)
#     for go_to_base in range(len(seq_str[0])):
#         this_base_set = set()
#         for go_to_seq in range(seq_num):
#             this_base_set.update(degenerate_dict.get(seq_str[go_to_seq][go_to_base], []))
#         consensus_res.append(degenerate_dict[tuple(sorted(this_base_set))])
#     return "".join(consensus_res)


def generate_consensus(*seq_str):
    consensus_res = []
    seq_num = len(seq_str)
    for go_to_base in range(len(seq_str[0])):
        this_base_set = set()
        for go_to_seq in range(seq_num):
            this_base_set.update(degenerate_dict_digit.get(seq_str[go_to_seq][go_to_base], []))
        consensus_res.append(degenerate_dict_digit[sum(this_base_set)])
    return "".join(consensus_res)

