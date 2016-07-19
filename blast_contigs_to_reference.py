#!/usr/bin/env python
from optparse import OptionParser
import commands, subprocess
import logging
import sys
import os


def simple_log(log, output_base):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'blast_contigs.log.txt'), mode='a')
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
    logfile = logging.FileHandler(os.path.join(output_base, 'blast_contigs.log.txt'), mode='a')
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
                      help='Input reference sequence or matrix in fasta format.'
                           'Only the first sequence would be taken as the reference.')
    options, args = parser.parse_args()
    if not (options.reference_fasta and options.input_contigs_fasta and options.output_directory):
        parser.print_help()
        print '\n######################################\nERROR: Insufficient REQUIRED arguments!\n'
        exit()
    if commands.getstatusoutput('blastn')[0] == 32512:
        print '\nERROR: blastn not in the PATH!'
        exit()
    if commands.getstatusoutput('makeblastdb')[0] == 32512:
        print '\nERROR: makeblastdb not in the PATH!'
        exit()
    log = simple_log(logging.getLogger(), options.output_directory)
    log.info(print_title)
    log.info(' '.join(sys.argv) + '\n')
    log = complicated_log(log, options.output_directory)
    return options, args, log


def makeblastdb(in_fasta, out_base, log):
    log.info("Making blast database ...")
    making_db = subprocess.Popen("makeblastdb -dbtype nucl -in "+in_fasta+" -out "+out_base, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = making_db.communicate()
    log.info("\n"+output.strip())
    log.info("Making blast database finished.")


def blastn(query_fasta, database, log):
    log.info("Making blast for "+query_fasta+" ...")
    making_blast = subprocess.Popen('blastn -num_threads 4 -query '+query_fasta+' -db '+database+' -out '+query_fasta+'.blast_result -outfmt 1', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = making_blast.communicate()
    log.info("\n" + output.strip())


def main():
    print_title = "\nBlast contigs to reference by Jianjun Jin on Jul 2 2016" \
                  "\n" \
                  "\nThis scripts would produce both mapping-like alignmnets and consensus alignments." \
                  "\n"
    options, args, log = require_options(print_title)

    log.info("Mapping finished.")


if __name__ == "__main__":
    main()