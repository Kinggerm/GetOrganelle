#!/usr/bin/env python
import sys
import time
time0 = time.time()
if len(sys.argv) == 3:
    try:
        fq1 = sys.argv[1].strip().strip("\"").strip("\'")
        fq2 = sys.argv[2].strip().strip("\"").strip("\'")
        file1 = open(fq1, 'rU').readlines()
        file2 = open(fq2, 'rU').readlines()
        names = {file1[i].split()[0].split('#')[0]: i for i in range(0, len(file1), 4)}
        outp_1 = open(fq1.replace('.fq', '')+'_paired.fq', 'w')
        outp_2 = open(fq2.replace('.fq', '')+'_paired.fq', 'w')
        outu_1 = open(fq1.replace('.fq', '')+'_unpaired.fq', 'w')
        outu_2 = open(fq2.replace('.fq', '')+'_unpaired.fq', 'w')
        for i in range(0, len(file2), 4):
            this_name = file2[i].split()[0]
            if this_name in names:
                id = names[this_name]
                outp_1.writelines(file1[id:id+4])
                outp_2.writelines(file2[i:i+4])
                del names[this_name]
            else:
                outu_2.writelines(file2[i:i+4])
        for i in range(0, len(file1), 4):
            this_name = file1[i].split()[0]
            if this_name in names:
                outu_1.writelines(file1[i:i+4])
        sys.stdout.write("\nCost"+str(round(time.time()-time0, 3))+"\nThank you!\n\n")
    except IOError:
        sys.stdout.write("\nUsage: get_pair_reads.py *_1.fq *_2.fq\n\n")
        exit()
else:
    sys.stdout.write("\nUsage: get_pair_reads.py *_1.fq *_2.fq\n\n")
    exit()
