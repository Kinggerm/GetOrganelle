import os
import sys
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
