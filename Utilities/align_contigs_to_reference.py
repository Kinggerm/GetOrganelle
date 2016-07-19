#!/usr/bin/env python
from optparse import OptionParser
import commands, subprocess
import string
import logging
import sys
import os


translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


def complementary_seq(input_seq):
    return string.translate(input_seq, translator)[::-1]


def simple_log(log, output_base):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'align_contigs.log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_simple.addHandler(console)
    log_simple.addHandler(logfile)
    return log_simple


def complicated_log(log, output_base):
    log_timed = log
    for handler in list(log_timed.handlers):
        log_timed.removeHandler(handler)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'align_contigs.log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_timed.addHandler(console)
    log_timed.addHandler(logfile)
    return log_timed


def require_options(print_title):
    usage = "python this_script.py -r reference.fasta -i *.contigs.fasta -o output"
    parser = OptionParser(usage=usage)
    parser.add_option('-r', dest='reference_fasta',
                      help='Input reference sequence or matrix in fasta format.'
                           'Only the first sequence would be taken as the reference.')
    parser.add_option('-i', dest='input_contigs_fasta',
                      help='Input contigs in fasta format. Could be bunch of files or *.fasta.')
    parser.add_option('-o', dest='output_directory',
                      help='Output directory.')
    parser.add_option('--mauve-path', default='progressiveMauve', dest='mauve_path',
                      help='If you have executable progressiveMauve but it is not in the path, please add it here')
    options, args = parser.parse_args()
    if not (options.reference_fasta and options.input_contigs_fasta and options.output_directory):
        parser.print_help()
        print '\n######################################\nERROR: Insufficient REQUIRED arguments!\n'
        exit()
    if commands.getstatusoutput(options.mauve_path+' -h')[0] in {32512, 32256}:
        print '\nERROR: cannot find executable progressiveMauve!'
        exit()
    if not os.path.exists(options.output_directory):
        os.mkdir(options.output_directory)
    log = simple_log(logging.getLogger(), options.output_directory)
    log.info(print_title)
    log.info(' '.join(sys.argv) + '\n')
    log = complicated_log(log, options.output_directory)
    return options, args, log


def execute_progressivemauve(mauve_path, input_files, output_file, log):
    log.info("Executing progressiveMauve ...")
    make_alignment = subprocess.Popen(mauve_path+" "+" ".join(input_files)+" --output "+output_file, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = make_alignment.communicate()
    log.info("\n"+output.strip())
    log.info("Executing progressiveMauve finished.")


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


def read_mauve_result(mauve_result):
    result_file = open(mauve_result, 'rU')
    names = []
    seqs = []
    this_line = result_file.readline()
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = result_file.readline()
            while this_line and not this_line.startswith('>'):
                this_seq += this_line.strip()
                this_line = result_file.readline()
            seqs.append(this_seq)
        else:
            this_line = result_file.readline()
    result_file.close()
    return [names, seqs]


def main():
    print_title = "\nAlign contigs to reference by Jianjun Jin on Jul 2 2016" \
                  "\n" \
                  "\nThis scripts would produce both mapping-like alignmnets and consensus alignments for phylogenetic analysis." \
                  "\n"
    options, args, log = require_options(print_title)

    log.info("Mapping finished.")


if __name__ == "__main__":
    main()