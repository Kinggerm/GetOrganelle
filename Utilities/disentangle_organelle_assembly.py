#!/usr/bin/env python
# coding: utf8
import time
import os
import sys
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.assembly_parser import Assembly
from Library.seq_parser import *
from Library.pipe_control_func import logging, timed_log, simple_log
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


def disentangle_circular_assembly(fastg_file, tab_file, prefix, weight_factor, mode="cp",
                                  log_hard_cov_threshold=10.,
                                  contamination_depth=5., contamination_similarity=5.,
                                  degenerate=True, degenerate_depth=1.5, degenerate_similarity=1.5,
                                  min_sigma_factor=0.1, keep_temp=False,
                                  verbose=False, log=None, debug=False):
    time_a = time.time()
    if log:
        log.info(">>> Parsing " + fastg_file + " ..")
    else:
        sys.stdout.write("Parsing " + fastg_file + " ..\n")
    input_graph = Assembly(fastg_file)
    time_b = time.time()
    if log:
        log.info(">>> Parsing input fastg file finished: " + str(round(time_b - time_a, 4)) + "s")
    else:
        sys.stdout.write("\n>>> Parsing input fastg file finished: " + str(round(time_b - time_a, 4)) + "s\n")
    temp_graph = prefix + ".temp.fastg" if keep_temp else None

    copy_results = input_graph.find_target_graph(tab_file, mode=mode, weight_factor=weight_factor,
                                                 log_hard_cov_threshold=log_hard_cov_threshold,
                                                 contamination_depth=contamination_depth,
                                                 contamination_similarity=contamination_similarity,
                                                 degenerate=degenerate, degenerate_depth=degenerate_depth,
                                                 degenerate_similarity=degenerate_similarity,
                                                 min_sigma_factor=min_sigma_factor,
                                                 temp_graph=temp_graph, verbose=verbose, log_handler=log, debug=debug)
    time_c = time.time()
    if log:
        log.info(">>> Detecting target graph finished: " + str(round(time_c - time_b, 4)) + "s")
        log.info(str(len(copy_results)) + " set(s) of graph detected.")
    else:
        sys.stdout.write("\n\n>>> Detecting target graph finished: " + str(round(time_c - time_b, 4)) + "s\n")
        sys.stdout.write(str(len(copy_results)) + " set(s) of graph detected.\n")

    for go_res, copy_res in enumerate(copy_results):
        idealized_graph = copy_res["graph"]
        # should add making one-step-inversion pairs for paths,
        # which would be used to identify existence of a certain isomer using mapping information
        count_path = 0
        for this_path, other_tag in idealized_graph.get_all_circular_paths(mode=mode, log_handler=log):
            count_path += 1
            open(prefix + ".graph" + str(go_res + 1) + other_tag + "." + str(count_path) + ".path_sequence.fasta", "w"). \
                write(idealized_graph.export_path(this_path).fasta_str())
        idealized_graph.write_to_gfa(prefix + ".graph" + str(go_res + 1) + ".selected_graph.gfa")
    time_d = time.time()
    if log:
        log.info(">>> Solving and unfolding graph finished: " + str(round(time_d - time_c, 4)) + "s")
    else:
        sys.stdout.write("\n\n>>> Solving and unfolding graph finished: " + str(round(time_d - time_c, 4)) + "s\n")


def get_options(print_title):
    parser = OptionParser("disentangle_organelle_assembly.py -g input.fastg -t input.tab -o output_dir")
    parser.add_option("-g", dest="fastg_file",
                      help="input fastg format file.")
    parser.add_option("-t", dest="tab_file",
                      help="input tab format file produced by slim_fastg.py.")
    parser.add_option("-o", dest="output_directory",
                      help="output directory.")
    parser.add_option("-m", dest="mode", default="cp",
                      help="organelle mode: cp, mt, or nt. Default:%default")
    parser.add_option("--weight-f", dest="weight_factor", type=float, default=100.0,
                      help="weight factor for excluding non-target contigs. Default:%default")
    parser.add_option("--depth-f", dest="depth_factor", type=float, default=10.,
                      help="Depth factor for excluding non-target contigs. Default:%default")
    parser.add_option("--contamination-depth", dest="contamination_depth", default=5., type=float,
                      help="Depth factor for confirming contaminating contigs. Default:%default")
    parser.add_option("--contamination-similarity", dest="contamination_similarity", default=0.9, type=float,
                      help="Similarity threshold for confirming contaminating contigs. Default:%default")
    parser.add_option("--no-degenerate", dest="degenerate", default=True, action="store_false",
                      help="Disable making consensus from parallel contig based on nucleotide degenerate table.")
    parser.add_option("--degenerate-depth", dest="degenerate_depth", default=1.5, type=float,
                      help="Depth factor for confirming parallel contigs. Default:%default")
    parser.add_option("--degenerate-similarity", dest="degenerate_similarity", default=0.95, type=float,
                      help="Similarity threshold for confirming parallel contigs. Default:%default")
    parser.add_option("--min-sigma", dest="min_sigma_factor", type=float, default=0.1,
                      help="Minimum deviation factor for excluding non-target contigs. Default:%default")
    parser.add_option("--keep-temp", dest="keep_temp_graph", default=False, action="store_true",
                      help="export intermediate graph file.")
    parser.add_option("--verbose", dest="verbose", default=False, action="store_true",
                      help="verbose logging.")
    parser.add_option("--debug", dest="debug", default=False, action="store_true",
                      help="for debug.")
    options, argv = parser.parse_args()
    if (options.fastg_file is None) or (options.tab_file is None) or (options.output_directory is None):
        parser.print_help()
        sys.stdout.write("Insufficient arguments!\n")
        sys.exit()
    else:
        if options.output_directory and not os.path.exists(options.output_directory):
            os.mkdir(options.output_directory)
        log = simple_log(logging.getLogger(), options.output_directory, "disentangle.")
        log.info(print_title)
        log.info(' '.join(sys.argv) + '\n')
        log = timed_log(log, options.output_directory, "disentangle.")
        return options, log


def main():
    time0 = time.time()
    print_title = "\nThis is a script for extracting circular organelle genome from assembly result (fastg). " + \
                  "\nBy jinjianjun@mail.kib.ac.cn\n\n"
    options, log = get_options(print_title)
    try:
        disentangle_circular_assembly(options.fastg_file, options.tab_file,
                                      os.path.join(options.output_directory, "target"),
                                      mode=options.mode,
                                      weight_factor=options.weight_factor,
                                      log_hard_cov_threshold=options.depth_factor,
                                      contamination_depth=options.contamination_depth,
                                      contamination_similarity=options.contamination_similarity,
                                      degenerate=options.degenerate, degenerate_depth=options.degenerate_depth,
                                      degenerate_similarity=options.degenerate_similarity,
                                      min_sigma_factor=options.min_sigma_factor,
                                      keep_temp=options.keep_temp_graph,
                                      log=log, verbose=options.verbose, debug=options.debug)
        log = simple_log(logging.getLogger(), options.output_directory, "disentangle.")

        log.info('\nTotal cost: ' + str(round(time.time() - time0, 4)) + 's\n')
    except:
        log.exception("")
        log = simple_log(log, options.output_directory, "disentangle.")
        log.info("\nTotal cost " + str(time.time() - time0))
        log.info("Please email jinjianjun@mail.kib.ac.cn if you find bugs!\n")
    logging.shutdown()


if __name__ == '__main__':
    main()


"""Copyright 2018 Jianjun Jin"""
