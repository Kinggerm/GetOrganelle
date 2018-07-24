#!/usr/bin/env python
import sys
import os
import time
time0 = time.time()
if len(sys.argv) == 3:
    try:
        fq1 = sys.argv[1].strip().strip("\"").strip("\'")
        fq2 = sys.argv[2].strip().strip("\"").strip("\'")
        file1 = open(fq1, 'rU').readlines()
        file2 = open(fq2, 'rU')
        names = {}
        # common_parts = [file1[0].split()[0].split('#')[0].split(".")]
        # len_parts = len(common_parts)
        first_n = file1[0].split()[0].split('#')[0].split(".")
        split_by_dot = len(first_n) > 2
        if split_by_dot:
            for i in range(0, len(file1), 4):
                this_n = ".".join(file1[i].split()[0].split('#')[0].split(".")[:2])
                names[this_n] = i
        else:
            for i in range(0, len(file1), 4):
                this_n = file1[i].split()[0].split('#')[0]
                names[this_n] = i
        outp_1 = open(fq1.rstrip('.fq')+'_paired.temp', 'w')
        outp_2 = open(fq2.rstrip('.fq')+'_paired.temp', 'w')
        outu_1 = open(fq1.rstrip('.fq')+'_unpaired.temp', 'w')
        outu_2 = open(fq2.rstrip('.fq')+'_unpaired.temp', 'w')
        this_line = file2.readline()
        while this_line:
            if split_by_dot:
                this_name = ".".join(this_line.split()[0].split('#')[0].split(".")[:2])
            else:
                this_name = this_line.split()[0].split('#')[0]
            if this_name in names:
                here_id = names[this_name]
                outp_1.writelines(file1[here_id:here_id + 4])
                outp_2.write(this_line)
                for k in range(3):
                    outp_2.write(file2.readline())
                this_line = file2.readline()
                del names[this_name]
            else:
                outu_2.write(this_line)
                for k in range(3):
                    outu_2.write(file2.readline())
                this_line = file2.readline()
        for i in range(0, len(file1), 4):
            this_name = file1[i].split()[0]
            if this_name in names:
                outu_1.writelines(file1[i:i+4])
        sys.stdout.write("\nCost"+str(round(time.time()-time0, 3))+"\nThank you!\n\n")
    except IOError:
        sys.stdout.write("\nUsage: get_pair_reads.py *_1.fq *_2.fq\n\n")
        exit()
    else:
        os.rename(fq1.rstrip('.fq')+'_paired.temp', fq1.rstrip('.fq')+'_paired.fq')
        os.rename(fq2.rstrip('.fq')+'_paired.temp', fq2.rstrip('.fq')+'_paired.fq')
        os.rename(fq1.rstrip('.fq')+'_unpaired.temp', fq1.rstrip('.fq')+'_unpaired.fq')
        os.rename(fq2.rstrip('.fq')+'_unpaired.temp', fq2.rstrip('.fq')+'_unpaired.fq')
else:
    sys.stdout.write("\nUsage: get_pair_reads.py *_1.fq *_2.fq\n\n")
    exit()
