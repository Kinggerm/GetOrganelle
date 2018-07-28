#!/usr/bin/env python
# coding: utf8
import os
import sys
import string
import math
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *


def require_commands():
    usage = 'python '+str(os.path.basename(__file__))+' input_files'
    parser = OptionParser(usage=usage)
    parser.add_option('-k', dest='minimum_match', help='minimum continuous match above which as overlap. Default: 10', default=10, type=int)
    parser.add_option('-s', dest='similarity', help='similarity threshold above which took as overlap. Default: 0.8', default=0.8, type=float)
    options, args = parser.parse_args()
    if not args:
        raise FileNotFoundError("No input fasta file found!")
    else:
        return options, args


transfer = {}
for char in string.ascii_lowercase:
    transfer[(char, char)] = 0
    transfer[(char, char.upper())] = 0
    transfer[(char.upper(), char)] = 0
    transfer[(char.upper(), char.upper())] = 0


def find_string_difference(this_string, this_reference, dynamic_span=2.0):
    len_str = len(this_string)
    len_ref = len(this_reference)
    if dynamic_span == 0:
        difference = sum([not (this_string[i], this_reference[i]) in transfer for i in range(min(len_ref, len_str))])+abs(len_ref-len_str)
        proper_end = this_string[-1] == this_reference[-1]
        return difference, proper_end
    else:
        dynamic_span = max(abs(len(this_string)-len(this_reference))+1, dynamic_span)
        this_match = int(not (this_string[0], this_reference[0]) in transfer)
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
            this_match = int(not (this_string[i], this_reference[start]) in transfer)
            this_matrix[(i, start)] = {'diagonal_out': this_matrix[(i-1, start-1)]['state'] + this_match,
                                       'down_out': this_matrix[(i-1, start)]['state'] + 1}
            this_matrix[(i, start)]['state'] = min(this_matrix[(i, start)].values())
            # middle
            for j in range(start+1, end-1):
                this_match = not (this_string[i], this_reference[j]) in transfer
                this_matrix[(i, j)] = {'diagonal_out': this_matrix[(i-1, j-1)]['state'] + this_match,
                                       'down_out': this_matrix[(i-1, j)]['state'] + 1,
                                       'right_out': this_matrix[(i, j-1)]['state'] + 1}
                this_matrix[(i, j)]['state'] = min(this_matrix[(i, j)].values())
            # end
            this_match = not (this_string[i], this_reference[end - 1]) in transfer
            this_matrix[(i, end-1)] = {'diagonal_out': this_matrix[(i-1, end-2)]['state'] + this_match}
            if (i, end-2) in this_matrix:
                this_matrix[(i, end-1)]['right_out'] = this_matrix[(i, end-2)]['state'] + 1
            this_matrix[(i, end-1)]['state'] = min(this_matrix[(i, end-1)].values())
        # print time.time()-time0
        difference = this_matrix[(len_str-1, len_ref-1)]['state']
        proper_end = True
        try:
            for j in range(start, len_ref):
                if this_matrix[(len_str-1, j)]['state'] < difference:
                    proper_end = False
                    break
            for i in range(max(0, len_str-len_ref+start), len_str):
                if this_matrix[(i, len_ref-1)]['state'] < difference:
                    proper_end = False
                    break
        except KeyError:
            return difference, False
        else:
            return difference, proper_end


def check_similarity_and_continuity(string1, string2, k_mer, similarity):
    overlap_len = len(string1)
    # continuity
    words_string1 = set([string1[x:x+k_mer] for x in range(overlap_len-k_mer+1)])
    for i in range(overlap_len-k_mer+1):
        if string2[i:i+k_mer] in words_string1:
            return True
    # similarity
    minimum_dif = (1-similarity)*overlap_len
    if find_string_difference(string1, string2, math.ceil(minimum_dif))[0] > minimum_dif:
        return False
    else:
        return True


def main():
    options, argv = require_commands()
    for in_fasta_file in argv:
        if not os.path.exists(in_fasta_file):
            raise FileExistsError(argv+" not found!")
        else:
            fasta_matrix = read_fasta(in_fasta_file)
            first_seq = fasta_matrix[1][0]
            new_len = len(first_seq) - (len(first_seq) - len(first_seq.replace('?', '')))*2
            new_fasta_matrix = [[], [], 70]
            count_contig = 1
            i = 0
            former_start = 0
            if first_seq[0] == '?' or first_seq[-1] == '?':
                sys.stdout.write('\nError: overlap cannot be in the two ends!')
                os._exit(0)
            while i < len(first_seq):
                if first_seq[i] == '?':
                    this_start = i
                    while first_seq[i] == '?':
                        this_end = i
                        i += 1
                    this_overlap_len = this_end-this_start+1
                    if check_similarity_and_continuity(first_seq[this_start-this_overlap_len:this_start],
                                                       first_seq[this_end+1:this_end+this_overlap_len+1],
                                                       k_mer=options.minimum_match,
                                                       similarity=options.similarity):
                        new_fasta_matrix[0].append('contig_'+str(count_contig))
                        new_fasta_matrix[1].append(former_start*'-'+first_seq[:this_start]+(new_len-this_start-former_start)*'-')
                        first_seq = first_seq[this_end+1:]
                        former_start += this_start - this_overlap_len
                        count_contig += 1
                        i = 0
                    else:
                        new_fasta_matrix[0].append('contig_'+str(count_contig))
                        new_fasta_matrix[1].append(former_start*'-'+first_seq[:this_start]+(new_len-this_start-former_start)*'-')
                        first_seq = first_seq[this_end+1:]
                        former_start += this_start - this_overlap_len
                        count_contig += 1
                        i = 0
                        sys.stdout.write('\nUnrecognized overlap at alignment!')  # ('+str(former_start+1)+'-'+str(former_start+this_overlap_len)+')!'
                        pass
                else:
                    i += 1
            new_fasta_matrix[0].append('contig_'+str(count_contig))
            new_fasta_matrix[1].append(former_start*'-'+first_seq)
            write_fasta(in_fasta_file+'.contig_alignments.fasta', new_fasta_matrix, False)


if __name__ == '__main__':
    main()
