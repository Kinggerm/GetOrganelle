#!/usr/bin/env python
# coding: utf8
import os
import sys
import csv
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *
from Library.pipe_control_func import LogInfo
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


def get_options():
    parser = OptionParser(usage="summary_get_organelle_output.py list_of_sample_folders -o tab_file")
    parser.add_option("-o", dest="output",
                      help="Output csv file.")
    # parser.add_option("--verbose", dest="verbose", default=False, action="store_true",
    #                   help="Verbose style.")
    options, argv = parser.parse_args()
    if not options.output or not len(argv):
        parser.print_help()
        sys.stdout.write("Insufficient arguments!\n")
        sys.exit()
    return options, argv


def main():
    options, argv = get_options()
    first_sample = LogInfo(argv[0])
    header = first_sample.header
    # if options.verbose:
    #     header.remove()
    with open(options.output, "w") as output_handler:
        res_out = csv.DictWriter(output_handler, fieldnames=header, delimiter="\t", extrasaction="ignore")
        res_out.writeheader()
        for sample_d in argv:
            res_out.writerow(LogInfo(sample_d).__dict__)


if __name__ == '__main__':
    main()

