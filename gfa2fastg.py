import string
import os

translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")
direction = {'+': True, '-': False}


def complementary_seq(input_seq):
    return string.translate(input_seq, translator)[::-1]


def write_fasta(out_dir, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_dir):
            out_dir = '.'.join(out_dir.split('.')[:-1])+'_.'+out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'wb')
    # if interleaved
    if matrix[2]:
        for i in xrange(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
                j += matrix[2]
            fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
    else:
        for i in xrange(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            fasta_file.write(matrix[1][i]+'\n')
    fasta_file.close()


def read_gfa_as_fastg(gfa_file):
    edges = {}
    count_edge = 0
    for line in open(gfa_file, 'rU'):
        line_split = line.rstrip().split('\t')
        if line_split[0] == 'S':
            count_edge += 1
            seq_len = int(line_split[3].split(':')[-1])
            coverage = round(float(line_split[4].split(':')[-1])/seq_len, 5)
            edges[line_split[1]] = {'name': 'EDGE_'+str(count_edge)+'_length_'+str(seq_len)+'_cov_'+str(coverage),
                                    ('seq', True): line_split[2],
                                    ('seq', False): complementary_seq(line_split[2]),
                                    True: [],
                                    False: []}
        elif line_split[0] == 'L':
            edges[line_split[1]][direction[line_split[2]]].append((line_split[3], direction[line_split[4]]))
            edges[line_split[3]][not direction[line_split[4]]].append((line_split[1], not direction[line_split[2]]))
    fasta_matrix = [[], [], 70]
    for original_edge_name in edges:
        for this_direction in [True, False]:
            if edges[original_edge_name][this_direction]:
                seq_name = edges[original_edge_name]['name']+(not this_direction)*'\''+':'
                list_next_edge = []
                for next_edge, next_direction in edges[original_edge_name][this_direction]:
                    list_next_edge.append(edges[next_edge]['name']+(not next_direction)*'\'')
                seq_name += ','.join(list_next_edge)+';'
            else:
                seq_name = edges[original_edge_name]['name']+(not this_direction)*'\''+';'
            fasta_matrix[0].append(seq_name)
            fasta_matrix[1].append(edges[original_edge_name][('seq', this_direction)])
    return fasta_matrix


def main():
    gfa_file = raw_input('Please input gfa file:').strip()
    write_fasta(gfa_file+'.fastg', read_gfa_as_fastg(gfa_file), False)


if __name__ == '__main__':
    main()
