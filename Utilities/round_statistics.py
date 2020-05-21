#! /usr/bin/env python

from optparse import OptionParser
import os
import sys
from math import ceil
from shutil import copyfile
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(path_of_this_script, ".."))
import GetOrganelleLib
from GetOrganelleLib.pipe_control_func import *
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.sam_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
import platform
SYSTEM_NAME = ""
if platform.system() == "Linux":
    SYSTEM_NAME = "linux"
elif platform.system() == "Darwin":
    SYSTEM_NAME = "macOS"
else:
    sys.stdout.write("Error: currently GetOrganelle is not supported for " + platform.system() + "! ")
    exit()
GO_LIB_PATH = os.path.split(GetOrganelleLib.__file__)[0]
GO_DEP_PATH = os.path.realpath(os.path.join(GO_LIB_PATH, "..", "GetOrganelleDep", SYSTEM_NAME))


def get_options():
    parser = OptionParser("round_statistics.py -f fasta_file -d output_per_round_folder -i Initial_mapped.fq -o output")
    parser.add_option("-f", dest="fasta",
                      help="input fasta file.")
    parser.add_option("-d", dest="output_per_round_dir",
                      help="output per round directory.")
    parser.add_option("-i", dest="initial_mapped",
                      help="seed fastq.")
    parser.add_option("-o", dest="output_base",
                      help="output folder.")
    parser.add_option("-R", dest="round", type=int,
                      help="rounds to check. default:automatic stop!")
    parser.add_option("-t", dest="threads", type=int, default=2,
                      help="threads.")
    parser.add_option("--which-bowtie2", dest="which_bowtie2", default="",
                      help="Assign the path to Bowtie2 binary files if not added to the path. "
                           "Default: try GetOrganelleDep/" + SYSTEM_NAME + "/bowtie2 first, then $PATH")
    parser.add_option('--random-seed', dest="random_seed", type=int, default=12345,
                      help="seed for random generator for bowtie2. Default: %default")
    parser.add_option("--threshold", dest="threshold", default="0,10",
                      help="sites with coverage above the threshold would be marked as covered. default: %default")
    parser.add_option("--continue", dest="resume", default=False, action="store_true")
    parser.add_option("--keep-temp", dest="keep_temp", default=False, action="store_true")
    parser.add_option("--draw", dest="draw_plot", default=False, action="store_true",
                      help="Draw density plot using matplotlib, which should be installed.")
    parser.add_option("--max-coverage-tick", dest="max_cov_tick")
    # parser.add_option("--average", default=False, action="store_true",
    #                   help="output average coverage.")
    parser.add_option("--debug", dest="debug",
                      help="Debug mode.")
    options, argv = parser.parse_args()
    if not (options.fasta and options.initial_mapped and options.output_base and options.output_per_round_dir):
        sys.stderr.write("Insufficient arguments!\n")
        sys.exit()
    if not os.path.isdir(options.output_base):
        os.mkdir(options.output_base)
    if options.debug:
        log_level = "DEBUG"
    else:
        log_level = "INFO"
    log_handler = simple_log(logging.getLogger(), options.output_base, "", log_level=log_level)
    log_handler.info("")
    log_handler.info(" ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n")
    if not options.which_bowtie2:
        try_this_bin = os.path.join(GO_DEP_PATH, "bowtie2", "bowtie2")
        if os.path.isfile(try_this_bin) and executable(try_this_bin):
            options.which_bowtie2 = os.path.split(try_this_bin)[0]
    if not executable(os.path.join(options.which_bowtie2, "bowtie2")):
        log_handler.error(os.path.join(options.which_bowtie2, "bowtie2") + " not accessible!")
        exit()
    if not executable(os.path.join(options.which_bowtie2, "bowtie2-build") + " --large-index"):
        log_handler.error(os.path.join(options.which_bowtie2, "bowtie2-build") + " not accessible!")
        exit()
    # if not executable(os.path.join(options.which_bowtie2, "bowtie2-build-l")):
    #     log_handler.error(os.path.join(options.which_bowtie2, "bowtie2-build-l") + " not accessible!")
    #     exit()
    log_handler = timed_log(log_handler, options.output_base, "", log_level=log_level)
    return options, log_handler


def count_fq_reads(fq_files):
    if type(fq_files) == str:
        fq_files = [fq_files]
    count = 0
    for fq_f in fq_files:
        for line in open(fq_f):
            count += 1
    return int(count/4.)


def main():
    options, log_handler = get_options()
    out_base = options.output_base
    all_fq_files = []
    for fq_f in os.listdir(options.output_per_round_dir):
        if fq_f.endswith("_1.fq") \
                and os.path.exists(os.path.join(options.output_per_round_dir, fq_f.replace("_1.fq", "_2.fq"))):
            all_fq_files.append((fq_f, fq_f.replace("_1.fq", "_2.fq")))
    all_fq_files.sort(key=lambda x: int(x[0].split(".")[1].split("_")[0]))
    len_ref_seq = len(read_fasta(options.fasta)[1][0])
    all_coverages = {}
    if options.round:
        log_handler.info(str(len(all_fq_files)) + " paired fastq files found. " + str(options.round) + " rounds required.")
    else:
        log_handler.info(str(len(all_fq_files)) + " paired fastq files found.")
    thresholds = [int(thr) for thr in options.threshold.split(",")]
    if all_fq_files:
        log_handler.info("\t".join(
            ["#Round"] +
            ["Round_Covered_T" + str(x) for x in thresholds] +
            ["Cumulative_Covered_T" + str(x) for x in thresholds] +
            ["Other_reads", "Other_percentage"]
        ))
    results_to_draw = []
    real_fq_files = [[os.path.join(options.output_per_round_dir, in_fq_f) for in_fq_f in fq_p] for fq_p in all_fq_files]
    all_fq_files = [(os.path.basename(options.initial_mapped), )] + all_fq_files
    real_fq_files = [(options.initial_mapped, )] + real_fq_files
    for go_to, fq_pairs in enumerate(all_fq_files):
        if options.round and go_to > options.round:
            log_handler.info("Hit required rounds! Exiting ..")
            break
        real_fq = real_fq_files[go_to]
        bowtie_base = os.path.join(out_base, fq_pairs[0].replace("_1.fq", ""))
        seed_file = os.path.join(out_base, os.path.basename(options.fasta))
        copyfile(options.fasta, seed_file)
        map_with_bowtie2(seed_file=seed_file, original_fq_files=real_fq, bowtie_out=bowtie_base,
                         resume=options.resume, threads=options.threads, random_seed=options.random_seed,
                         which_bowtie2=options.which_bowtie2,
                         silent=True, log_handler=log_handler, generate_fq=True, bowtie2_mode="--very-fast-local")
        this_result = [fq_pairs[0].replace("_1.fq", "")]
        # this round coverage
        this_records = MapRecords(bowtie_base + ".sam")
        this_records.update_coverages()
        this_coverages = this_records.coverages
        if this_coverages:
            ref_bowtie = sorted(this_coverages.keys())[0]
            for threshold in thresholds:
                count_site = len([site_cov for site_cov in this_coverages[ref_bowtie] if site_cov > threshold])
                this_result.append(float(count_site) / len_ref_seq)
            results_to_draw.append(this_coverages)
            # merge
            for ref in this_coverages:
                if ref not in all_coverages:
                    all_coverages[ref] = list(this_coverages[ref])
                else:
                    for go_s, site_cov in enumerate(this_coverages[ref]):
                        all_coverages[ref][go_s] += site_cov
            ref_bowtie = sorted(all_coverages.keys())[0]
            for threshold in thresholds:
                count_site = len([site_cov for site_cov in all_coverages[ref_bowtie] if site_cov > threshold])
                this_result.append(float(count_site) / len_ref_seq)
        elif not this_coverages and not options.round:
            log_handler.info("No more target found! Exiting ..")
            if not options.keep_temp:
                os.remove(bowtie_base + ".fq")
                os.remove(bowtie_base + ".sam")
            break
        else:
            for threshold in thresholds:
                this_result.append("-")
            for threshold in thresholds:
                this_result.append("-")
        # all reads
        all_reads_num = count_fq_reads(real_fq)
        other_reads_num = all_reads_num - count_fq_reads(bowtie_base + ".fq")
        this_result.append(other_reads_num)
        if all_reads_num:
            this_result.append(other_reads_num/float(all_reads_num))
        else:
            this_result.append("-")
        if not options.keep_temp:
            os.remove(bowtie_base + ".fq")
            os.remove(bowtie_base + ".sam")
        log_handler.info("\t".join([str(val) for val in this_result]))
    # draw
    if options.draw_plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        len_round = len(results_to_draw)
        figure = plt.figure(figsize=(50, len_round*3))
        width = 50
        x_values = list(range(1, len_ref_seq + 1))
        x_trans_values = [x_values[n] for n in range(0, len(x_values), width)]
        y_values_list = []
        for i in range(len_round):
            ref_bowtie = sorted(results_to_draw[i].keys())[0]
            y_values_list.append([results_to_draw[i][ref_bowtie][site_here-1] for site_here in x_values])
            y_values_list[-1] = [sum(y_values_list[-1][n:n + width]) / float(len(y_values_list[-1][n:n + width]))
                                 for n in range(0, len(y_values_list[-1]), width)]
        if options.max_cov_tick:
            y_tick_max = int(options.max_cov_tick)
            # set larger value to max_cov_tick
            for go_r in range(len_round):
                for check_go in range(len(y_values_list[go_r])):
                    y_values_list[go_r][check_go] = min(y_tick_max, y_values_list[go_r][check_go])
        else:
            y_tick_max = ceil(max([max(y_v) for y_v in y_values_list])/100.) * 100
        y_ticks = [int(y_tick_max/3), int(y_tick_max*2/3), y_tick_max]
        for i in range(len_round):
            plt.subplot(len_round, 1, i + 1)
            plt.bar(x_trans_values, y_values_list[i], width=width, color="darkgreen")
            # plt.hist([1, 2, 3, 5, 2, 1], color="darkgreen")
            plt.yticks(y_ticks, [str(y_t) + "x" for y_t in y_ticks], fontsize=40)
            if i < len_round - 1:
                plt.xticks([])
            else:
                plt.xticks(fontsize=40)
        figure.savefig(os.path.join(out_base, "coverage_per_round.pdf"), bbox_inches="tight")


if __name__ == '__main__':
    main()
