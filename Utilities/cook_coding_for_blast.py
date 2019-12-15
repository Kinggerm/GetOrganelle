#!/usr/bin/env python
import sys
import os
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.seq_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]

if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
    sys.stdout.write("usage: python cook_coding_for_blast.py fasta_files\n")
else:
    for in_fasta in sys.argv[1:]:
        try:
            f_matrix = read_fasta(in_fasta)
            for i in range(len(f_matrix[0])):
                temp = [x.strip() for x in f_matrix[0][i].split('-')]
                f_matrix[0][i] = temp[1]+' - '+temp[0]
            i = 0
            del_count = 0
            seq_sets = set()
            while i < len(f_matrix[0]):
                if f_matrix[1][i] in seq_sets or len(f_matrix[1][i]) < 20 or f_matrix[0][i].startswith('gene'):
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
        except IndexError:
            sys.stdout.write("Error: " + in_fasta + " is in illegal format!\n"
                             "See samples in GetOrganelleLib/LabelDatabase/*.fasta\n")
