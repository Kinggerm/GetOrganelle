#!/usr/bin/env python
import sys
import os
from argparse import ArgumentParser
import time
# from https://github.com/Kinggerm/PersonalUtilities
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
from GetOrganelleLib.seq_parser import read_fasta, detect_plastome_architecture
from GetOrganelleLib.versions import get_versions
# from GetOrganelleLib.sam_parser import *
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]


def get_options(description):
    parser = ArgumentParser(description=description, usage="plastome_arch_info.py fasta_format_sequence_file(s)")
    parser.add_argument("sequences", metavar="sequences", type=str, nargs="+",
                        help="Input fasta format sequences (split the files by spaces).")
    parser.add_argument("-o", dest="output",
                        help="output file.")
    parser.add_argument("-r", dest="min_ir_length", default=5000, type=int,
                        help="The minimum repeat length treated as the IR region of plastome. Default: [%(default)s]")
    parser.add_argument("-v", dest="valid_bases", default="ATGCRMYKHBDVatgcrmykhbdv",
                        help="Valid bases. Default: ATGCRMYKHBDVatgcrmykhbdv")
    parser.add_argument("--version", action="version",
                        version="GetOrganelle v{version}".format(version=get_versions()))
    options = parser.parse_args()
    if not len(options.sequences):
        parser.print_help()
        sys.exit()
    else:
        for f in options.sequences:
            if not os.path.isfile(f):
                raise IOError(f + " not found/valid!")
        options.valid_bases = set(list(options.valid_bases))
        return options, options.sequences


def main():
    time0 = time.time()
    description = """\n## This script helps you count the LSC/SSC/IR-DR lengths from a batch of plastome sequences.\n
                  ## Jianjun\n\n"""
    options, sequence_files = get_options(description=description)
    sys.stdout.write(
        "file_name\tsequence_name\ttotal_length\tLSC_length\tSSC_length\tIR/DR_length\tarch_Notes\tGC_content\n")
    if options.output:
        out_handler = open(options.output, "w")
        out_handler.close()
    for this_f in sequence_files:
        this_matrix = read_fasta(this_f)
        for i in range(len(this_matrix[0])):
            arch = detect_plastome_architecture(this_matrix[1][i], options.min_ir_length, options.valid_bases)
            this_upper = this_matrix[1][i].upper()
            this_gc = this_upper.count("G") + this_upper.count("C")
            this_at = this_upper.count("A") + this_upper.count("T")
            if this_gc + this_at == len(this_upper):
                gc_content = "%.4f" % (this_gc / float(len(this_upper)))
            else:
                gc_content = "ca. " + "%.8f" % (this_gc / float(len(this_upper)))
            out_line = "\t".join([this_f, this_matrix[0][i], str(len(this_matrix[1][i]))] +
                                 [str(x) for x in arch] + [gc_content]) + "\n"
            sys.stdout.write(out_line)
            if options.output:
                open(options.output, "a").write(out_line)
    sys.stdout.write("\n## Cost: " + str(round(time.time() - time0, 2)) + "s\n\n")


if __name__ == '__main__':
    main()
