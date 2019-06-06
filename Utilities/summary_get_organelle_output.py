#!/usr/bin/env python
# coding: utf8
import os
import sys
import csv
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.pipe_control_func import LogInfo
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
    # detect prefix
    dir_prefix_pairs = []
    for folder in argv:
        if os.path.isdir(folder):
            prefix_list = [log_f[:-len("get_org.log.txt")]
                           for log_f in os.listdir(folder) if log_f.endswith("get_org.log.txt")]
            for this_prefix in prefix_list:
                dir_prefix_pairs.append((folder, this_prefix))
    first_sample = LogInfo(sample_out_dir=dir_prefix_pairs[0][0], prefix=dir_prefix_pairs[0][1])
    header = first_sample.header
    # if options.verbose:
    #     header.remove()
    with open(options.output, "w") as output_handler:
        res_out = csv.DictWriter(output_handler, fieldnames=header, delimiter="\t", extrasaction="ignore")
        res_out.writeheader()
        for sample_d, sample_p in dir_prefix_pairs:
            res_out.writerow(LogInfo(sample_out_dir=sample_d, prefix=sample_p).__dict__)


if __name__ == '__main__':
    main()

