#!/usr/bin/env python
"""This script uses an naive De Bruijn approach to convert sequence back into an assembly graph file,
   such as a gfa (Graphical Fragment Assembly) or a fastg file"""
import sys
import os
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.assembly_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


def get_options():
    parser = OptionParser(
        description="This script uses an naive De Bruijn approach to convert sequence back into an "
                    "assembly graph file, such as a gfa (Graphical Fragment Assembly) or a fastg file.",
        usage="reconstruct_graph_from_fasta.py -i fasta_file -o out.gfa")
    parser.add_option("-i", dest="input",
                      help="Input fasta file.")
    parser.add_option("-o", dest="output",
                      help="Output graph file. The output format is GFA by default, but FASTG only when "
                           "indicated with postfix '.fastg'.")
    parser.add_option("-k", dest="kmer", default=55,
                      help="kmer for reconstructing De Bruijn graph. Default:%default")
    parser.add_option("-c", "--circular", dest="circular", default="auto",
                      help="Sequences in input fasta file are all circular (yes/no/auto). "
                           "The auto mode enables detection by checking the existence of '(circular)' in "
                           "the end of the header of each sequence. Default:%default")
    options, argv = parser.parse_args()
    if not (options.output and options.input):
        parser.print_help()
        sys.stdout.write("Insufficient arguments!\n")
        sys.exit()
    elif options.circular not in ("auto", "yes", "no"):
        # parser.print_help()
        sys.stdout.write("Illegal -c input! circular mode must be one of yes/no/auto!\n")
        sys.exit()
    elif options.kmer % 2 == 0:
        # parser.print_help()
        sys.stdout.write("Illegal -k input! kmer must be an odd number!\n")
    return options, argv


def main():
    options, argv = get_options()
    # detect postfix
    de_burijn_graph = NaiveDeBruijnGraph(options.input, kmer_len=options.kmer, circular=options.circular)
    if options.output.endswith(".fastg"):
        sys.stdout.warning("Fastg is not recommended!\n")
        de_burijn_graph.write_to_fastg(options.output)
    else:
        de_burijn_graph.write_to_gfa(options.output)


if __name__ == '__main__':
    main()
