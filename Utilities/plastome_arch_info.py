#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import time
# from https://github.com/Kinggerm/PersonalUtilities
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.sam_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


def get_options():
    parser = OptionParser(usage="plastome_arch_info.py fasta_format_sequence_file(s)")
    parser.add_option("-o", dest="output",
                      help="output file.")
    parser.add_option("-r", dest="min_ir_length", default=5000, type=int,
                      help="The minimum repeat length treated as the IR region of plastome. Default: [%default]")
    parser.add_option("-v", dest="valid_bases", default="ATGCRMYKHBDVatgcrmykhbdv",
                      help="Valid bases. Default: ATGCRMYKHBDVatgcrmykhbdv")
    options, argv = parser.parse_args()
    if not len(argv):
        parser.print_help()
        sys.exit()
    else:
        for f in argv:
            if not os.path.isfile(f):
                raise IOError(f + " not found/valid!")
        options.valid_bases = set(list(options.valid_bases))
        return options, argv


def detect_architecture(sequence, min_repeat_length, accepted_char):
    # assume the longest is ir
    all_repeats = find_exact_repeats(sequence, min_repeat_length, True, accepted_char)
    if all_repeats:
        # Sorting makes:
        # direct1==1 and direct2==-1
        # start1 be the smallest forward start
        ir_locations_1 = sorted(all_repeats[0], key=lambda x: (-x["direction"], x["start"]))
        ir_locations_2 = sorted(reverse_repeats_info(ir_locations_1), key=lambda x: (-x["direction"], x["start"]))
        ir_locations = sorted([ir_locations_1, ir_locations_2],
                              key=lambda x: (-max([y["direction"] for y in x]), x[0]["start"]))[0]
        if len(ir_locations) != 2:
            return "-", "-", "-", "not canonical IR"
        elif ir_locations[0]["direction"] == ir_locations[1]["direction"]:
            start1, end1, direct1 = ir_locations[0]["start"], ir_locations[0]["end"], ir_locations[0]["direction"]
            start2, end2, direct2 = ir_locations[1]["start"], ir_locations[1]["end"], ir_locations[1]["direction"]
            # cross the end, meaning site:seq_len in (DR1)
            if end1 < start1:
                if end2 >= start1:
                    return 0, 0, ir_locations[0]["length"], "DR detected and overlaps"
                else:
                    return start1 - end2 - 1, 0, ir_locations[0]["length"], "DR detected and overlaps"
            elif end2 < start2:
                if end2 >= start1:
                    if end1 >= start2:
                        return 0, 0, ir_locations[0]["length"], "DR detected and overlaps"
                    else:
                        return start2 - end1 - 1, 0, ir_locations[0]["length"], "DR detected and overlaps"
                elif end1 >= start2:
                    return start1 - end2 - 1, 0, ir_locations[0]["length"], "DR detected and overlaps"
                else:
                    ssc, lsc = sorted([start1 - end2 - 1, start2 - end1 - 1])
                    return lsc, ssc, ir_locations[0]["length"], "DR detected"
            else:
                ssc, lsc = sorted([start2 - end1 - 1, len(sequence) + start1 - end2 - 1])
                return lsc, ssc, ir_locations[0]["length"], "DR detected"
        else:
            start1, end1, direct1 = ir_locations[0]["start"], ir_locations[0]["end"], ir_locations[0]["direction"]
            start2, end2, direct2 = ir_locations[1]["start"], ir_locations[1]["end"], ir_locations[1]["direction"]
            # cross the end, meaning site:seq_len in (IR1)
            if end1 < start1:
                # seq_len >= start2 >= start1
                if start2 >= start1:
                    return 0, 0, ir_locations[0]["length"], "IR overlaps"
                else:
                    return start1 - start2 - 1, 0, ir_locations[0]["length"], "IR overlaps"
            elif start2 < end2:
                if start2 >= start1:
                    if end1 >= end2:
                        return 0, 0, ir_locations[0]["length"], "IR overlaps"
                    else:
                        return end2 - end1 - 1, 0, ir_locations[0]["length"], "IR overlaps"
                else:
                    ssc, lsc = sorted([end2 - end1 - 1, start1 - start2 - 1])
                    return lsc, ssc, ir_locations[0]["length"], "IR detected"
            else:
                ssc, lsc = sorted([end2 - end1 - 1, len(sequence) + start1 - start2 - 1])
                return lsc, ssc, ir_locations[0]["length"], "IR detected"
    else:
        return "-", "-", "-", "no IR found"


def main():
    time0 = time.time()
    sys.stdout.write("\n"
                     "## This script helps you count the LSC/SSC/IR-DR lengths from a batch of plastome sequences.\n"
                     "## by jinjianjun@mail.kib.ac.cn\n\n")
    options, argv = get_options()
    sys.stdout.write("file_name\tsequence_name\ttotal_length\tLSC_length\tSSC_length\tIR/DR_length\tNotes\n")
    if options.output:
        out_handler = open(options.output, "w")
        out_handler.close()
    for this_f in argv:
        this_matrix = read_fasta(this_f)
        for i in range(len(this_matrix[0])):
            arch = detect_architecture(this_matrix[1][i], options.min_ir_length, options.valid_bases)
            out_line = "\t".join([this_f, this_matrix[0][i], str(len(this_matrix[1][i]))] +
                                 [str(x) for x in arch]) + "\n"
            sys.stdout.write(out_line)
            if options.output:
                open(options.output, "a").write(out_line)
    sys.stdout.write("\n## Cost: " + str(round(time.time() - time0, 2)) + "s\n\n")


if __name__ == '__main__':
    main()
