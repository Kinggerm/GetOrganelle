#!/usr/bin/env python
import sys
from optparse import OptionParser
import subprocess
try:
    import commands
except:
    pass
import os


def check_db(reference_fa_base):
    in_index = reference_fa_base + '.index'
    if min([os.path.exists(in_index+postfix) for postfix in ('.nhr', '.nin', '.nsq')]):
        pass
    elif reference_fa_base:
        print('Making BLAST db ... ')
        makedb_result = subprocess.Popen('makeblastdb -dbtype nucl -in '+reference_fa_base+' -out '+in_index,
                                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, error = makedb_result.communicate()
        if 'Error' in str(output) or 'error' in str(output):
            print('ERROR: Blast terminated with following info:\n'+output.decode('utf-8'))
            exit()
        print('Making BLAST db finished.\n')
    else:
        print('ERROR: No reference input!')
        exit()
    return in_index


def execute_blast(query, blast_db, output, outfmt):
    print("Executing BLAST ...")
    make_blast = subprocess.Popen(
        'blastn -evalue 1E-30 -num_threads 4 -word_size 10 -query ' + query + ' -db ' + blast_db + ' -out ' + output + ' -outfmt '+str(outfmt),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = make_blast.communicate()
    if "(ERR)" in str(output) or "Error:" in str(output) or 'error' in str(output):
        print('\nERROR:\n', output.decode('utf-8'))
        exit()
    print("Executing BLAST finished.")


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


def require_options():
    try:
        # python3
        blast_in_path = subprocess.getstatusoutput('blastn')
    except AttributeError:
        # python2
        blast_in_path = commands.getstatusoutput('blastn')
    if blast_in_path[0] == 32512:
        sys.stdout.write('\nError: blastn not in the path!')
        exit()
    try:
        # python3
        makeblastdb_in_path = subprocess.getstatusoutput('makeblastdb')
    except AttributeError:
        # python2
        makeblastdb_in_path = commands.getstatusoutput('makeblastdb')
    if makeblastdb_in_path[0] == 32512:
        sys.stdout.write('\nError: makeblastdb not in the path!')
        exit()
    usage = "Usage: rm_low_coverage_duplicated_contigs.py *.fastg"
    parser = OptionParser(usage=usage)
    parser.add_option('--cov-t', dest='coverage_threshold', default=0.12,
                      help='With ratio (coverage of query/coverage of subject) below which, '
                           'the query would be exposed to discarded. Default: 0.12')
    parser.add_option('--len-t', dest='length_threshold', default=0.9,
                      help='With overlap (length of hit of query/ length of query) above which, '
                           'the query would be exposed to discarded. Default: 0.9')
    parser.add_option('--blur', dest='blur_bases', default=False, action='store_true',
                      help='Replace hit low-coverage bases with N.')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true',
                      help='Keep temp blast files.')
    options, args = parser.parse_args()
    if not args:
        parser.print_help()
        sys.stdout.write('\n######################################\nERROR: Insufficient REQUIRED arguments!\n\n')
        exit()
    return options, args


def purify_fastg(fastg_file, cov_threshold, len_threshold, blur, keep_temp):
    index_base = check_db(fastg_file)
    execute_blast(fastg_file, index_base, fastg_file + '.blast', 7)
    blast_result = open(fastg_file + '.blast', 'rU')
    suspicious_nodes = {}
    for line in blast_result:
        if not line.startswith("#") and line.strip():
            records = line.strip().split('\t')
            query_id = records[0]
            subject_id = records[1]
            if query_id != subject_id:
                query_cov = float(query_id.split(':')[0].rstrip('\'').split('cov_')[1])
                subject_cov = float(subject_id.split(':')[0].rstrip('\'').split('cov_')[1])
                node_name = query_id.split(':')[0]
                if query_cov/subject_cov < cov_threshold:
                    if node_name in suspicious_nodes:
                        for base in range(int(records[6]), int(records[7]) + 1):
                            suspicious_nodes[node_name].add(base)
                    else:
                        suspicious_nodes[node_name] = set(range(int(records[6]), int(records[7]) + 1))
    fastg_matrix = read_fasta_as_list(fastg_file)
    i = 0
    while i < len(fastg_matrix[0]):
        node_name = fastg_matrix[0][i].split(':')[0].strip(';').strip('-')
        if node_name in suspicious_nodes:
            if len(suspicious_nodes[node_name])/float(len(fastg_matrix[1][i])) > len_threshold:
                del fastg_matrix[0][i]
                del fastg_matrix[1][i]
                continue
            elif blur:
                for base_to_blur in suspicious_nodes[node_name]:
                    fastg_matrix[1][i][base_to_blur-1] = 'N'
        i += 1
    write_fasta_with_list(fastg_file+'.purified.fastg', fastg_matrix, True)
    blast_result.close()
    if not keep_temp:
        os.remove(fastg_file + '.blast')
        for postfix in ('.nhr', '.nin', '.nsq'):
            os.remove(fastg_file + '.index' + postfix)


def main():
    options, args = require_options()
    for fastg_file in args:
        purify_fastg(fastg_file, options.coverage_threshold, options.length_threshold,
                     options.blur_bases, options.keep_temp)


if __name__ == '__main__':
    main()


"""Copyright 2016 Jianjun Jin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License."""