#!/usr/bin/env python
"""This script converts a gfa (Graphical Fragment Assembly) file into a fasta file"""
import sys
import os
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
from GetOrganelleLib.seq_parser import *
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]


def read_gfa_as_fasta(gfa_file):
    fasta_matrix = [[], [], 70]
    for line in open(gfa_file):
        line_split = line.rstrip().split('\t')
        if line_split[0] == 'S':
            seq_len = int(line_split[3].split(':')[-1])
            coverage = round(float(line_split[4].split(':')[-1])/seq_len, 5)
            fasta_matrix[0].append(line_split[1] + ' length_' + str(seq_len) + ' cov_' + str(coverage))
            fasta_matrix[1].append(line_split[2])
    return fasta_matrix


def main():
    if len(sys.argv) > 1:
        for i in sys.argv:
            if '-h' in i or 'help' in i:
                print("Usage: gfa2fasta.py *.gfa")
                break
        else:
            for gfa_file in sys.argv[1:]:
                write_fasta(gfa_file + '.fasta', read_gfa_as_fasta(gfa_file), False)
    else:
        if type(2/1) == float:
            gfa_file = input('Please input gfa file:').strip()
        else:
            gfa_file = raw_input('Please input gfa file:').strip() # type: ignore
        if gfa_file.strip():
            write_fasta(gfa_file +'.fasta', read_gfa_as_fasta(gfa_file), False)


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