#!/usr/bin/env python
"""This script uses an naive De Bruijn approach to convert sequence back into an assembly graph file,
   such as a gfa (Graphical Fragment Assembly) or a fastg file"""
import sys
import os
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.assembly_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


def get_options():
    parser = OptionParser(
        description="This script uses an naive De Bruijn approach to convert sequence back into an "
                    "assembly graph file, such as a gfa (Graphical Fragment Assembly) or a fastg file.",
        usage="reconstruct_graph_from_fasta.py -i fasta_file -o out.gfa")
    parser.add_option("-i", dest="input",
                      help="Input fasta file.")
    parser.add_option("-o", dest="output", default="",
                      help="Output graph file. The output format is GFA by default, but FASTG only when "
                           "indicated with postfix '.fastg'.")
    parser.add_option("-L", "--overlap", dest="overlap", default=55, type=int,
                      help="overlap for reconstructing De Bruijn graph. Default:%default")
    parser.add_option("-c", "--circular", dest="circular", default="auto",
                      help="Sequences in input fasta file are all circular (yes/no/auto). "
                           "The auto mode enables detection by checking the existence of '(circular)' in "
                           "the end of the header of each sequence. Default:%default")
    parser.add_option("--single-chain", dest="single_chain", default=False, action="store_true",
                      help="The input sequence(s) was by default treated as DNA double-chain with its complementary "
                           "sequence. Choose this flag to turn off.")
    parser.add_option("--out-kg", dest="out_kg", default="",
                      help="Output kmer node graph.")
    options, argv = parser.parse_args()
    if not ((options.output or options.out_kg) and options.input):
        parser.print_help()
        sys.stdout.write("Insufficient arguments!\n")
        sys.exit()
    elif options.circular not in ("auto", "yes", "no"):
        # parser.print_help()
        sys.stdout.write("Illegal -c input! circular mode must be one of yes/no/auto!\n")
        sys.exit()
    elif options.overlap % 2 == 0:
        # parser.print_help()
        sys.stdout.write("Illegal -k input! kmer must be an odd number!\n")
    return options, argv


def main():
    import time
    time_0 = time.time()
    options, argv = get_options()
    # detect postfix
    kmer_node_graph = NaiveKmerNodeGraph(options.input, kmer_len=options.overlap,
                                         circular=options.circular, single_chain=options.single_chain)
    if options.output:
        assembly_graph = kmer_node_graph.generate_assembly_graph()
        if options.output.endswith(".fastg"):
            sys.stdout.warning("Fastg is not recommended!\n")
            assembly_graph.write_to_fastg(options.output)
        else:
            assembly_graph.write_to_gfa(options.output)
    if options.out_kg:
        kmer_node_graph.write_to_gfa(options.out_kg)
    sys.stdout.write("Took " + "%.4f" % (time.time() - time_0) + "s in generating " +
                     options.output * int(bool(options.output)) + ", " * int(bool(options.output and options.out_kg)) +
                     options.out_kg * int(bool(options.out_kg)) + "\n")


if __name__ == '__main__':
    main()
