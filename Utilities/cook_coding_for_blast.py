import sys


def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = fasta_file.readline()
            while this_line and not this_line.startswith('>'):
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
            seqs.append(this_seq)
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs]
    
in_fasta = input('fasta:').strip()
f_matrix = read_fasta(in_fasta)
for i in range(len(f_matrix[0])):
    temp = [x.strip() for x in  f_matrix[0][i].split('-')]
    f_matrix[0][i] = temp[1]+' - '+temp[0]
i = 0
del_count = 0
seq_sets = set()
while i < len(f_matrix[0]):
    if f_matrix[1][i] in seq_sets or len(f_matrix[1][i])<20 or f_matrix[0][i].startswith('gene'):
        del f_matrix[0][i]
        del f_matrix[1][i]
        del_count += 1
    else:
        seq_sets.add(f_matrix[1][i])
        i += 1
sys.stdout.write('delete '+str(del_count)+'\n')
out_fasta = open(in_fasta+'.new.fasta', 'w')
for i in range(len(f_matrix[0])):
    out_fasta.write('>'+f_matrix[0][i]+'\n'+f_matrix[1][i]+'\n')
out_fasta.close()