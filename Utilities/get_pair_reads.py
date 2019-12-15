#!/usr/bin/env python
import sys
import time
import os
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.seq_parser import get_paired_and_unpaired_reads
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
time0 = time.time()
if len(sys.argv) == 3:
    try:
        fq_file_1 = sys.argv[1].strip().strip("\"").strip("\'")
        fq_file_2 = sys.argv[2].strip().strip("\"").strip("\'")
        if fq_file_1.endswith(".fastq"):
            postfix = ".fastq"
        elif fq_file_1.endswith(".fq"):
            postfix = ".fq"
        else:
            postfix = ""
        out_paired_1 = fq_file_1[:-len(postfix)] + "_paired.fq"
        out_paired_2 = fq_file_2[:-len(postfix)] + "_paired.fq"
        out_unpaired_1 = fq_file_1[:-len(postfix)] + "_unpaired.fq"
        out_unpaired_2 = fq_file_2[:-len(postfix)] + "_unpaired.fq"
        get_paired_and_unpaired_reads(fq_file_1, fq_file_2, out_paired_1, out_paired_2, out_unpaired_1, out_unpaired_2)
        sys.stdout.write("\nCost"+str(round(time.time()-time0, 3))+"\nThank you!\n\n")
    except IOError:
        sys.stdout.write("\nUsage: get_pair_reads.py *_1.fq *_2.fq\n\n")
        exit()
else:
    sys.stdout.write("\nUsage: get_pair_reads.py *_1.fq *_2.fq\n\n")
    exit()
