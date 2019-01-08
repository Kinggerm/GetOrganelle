#! /usr/bin/env python

from optparse import OptionParser
import os
import sys
from math import ceil
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.pipe_control_func import *
from Library.seq_parser import *
from Library.sam_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


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
    parser.add_option("--threshold", dest="threshold", default="0,10",
                      help="sites with coverage above the threshold would be marked as covered. default:[%default]")
    parser.add_option("--continue", dest="resume", default=False, action="store_true")
    parser.add_option("--keep-temp", dest="keep_temp", default=False, action="store_true")
    parser.add_option("--draw", dest="draw_plot", default=False, action="store_true",
                      help="Draw density plot using matplotlib, which should be installed.")
    parser.add_option("--max-coverage-tick", dest="max_cov_tick")
    # parser.add_option("--average", default=False, action="store_true",
    #                   help="output average coverage.")
    options, argv = parser.parse_args()
    if not (options.fasta and options.initial_mapped and options.output_base and options.output_per_round_dir):
        sys.stderr.write("Insufficient arguments!\n")
        sys.exit()
    if not os.path.isdir(options.output_base):
        os.mkdir(options.output_base)
    log = simple_log(logging.getLogger(), options.output_base, "")
    log.info("")
    log.info(' '.join(sys.argv) + '\n')
    log = timed_log(log, options.output_base, "")
    return options, log


def mapping_with_bowtie2(seed_file, original_fq_files, bowtie_out, resume, threads, log):
    if not (os.path.exists(seed_file + '.index.1.bt2l')):
        # log.info("Making seed bowtie2 index ...")
        build_seed_index = subprocess.Popen("bowtie2-build --large-index " + seed_file + " " + seed_file + '.index',
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = build_seed_index.communicate()
        if "unrecognized option" in str(output):
            build_seed_index = subprocess.Popen("bowtie2-build " + seed_file + " " + seed_file + '.index',
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = build_seed_index.communicate()
        if "(ERR)" in str(output) or "Error:" in str(output):
            log.error('\n' + str(output))
            exit()
        # if verbose_log:
        #     log.info("\n" + str(output).strip())
        # log.info("Making seed bowtie2 index finished.")
    seed_index_base = seed_file + '.index'
    res_path_name, res_base_name = os.path.split(bowtie_out)
    total_seed_file = [os.path.join(res_path_name, x + res_base_name + ".fq") for x in ("temp.", "")]
    total_seed_sam = [os.path.join(res_path_name, x + res_base_name + ".sam") for x in ("temp.", "")]
    if not (resume and os.path.exists(total_seed_file[1])):
        # log.info("Mapping reads to seed bowtie2 index ...")
        make_seed_bowtie2 = subprocess.Popen(
            "bowtie2 --mm -p " + str(threads) + " --very-fast-local --al " + total_seed_file[
                0] + " -x " + seed_index_base + " -U " +
            ",".join(original_fq_files) + " -S " + total_seed_sam[0] + " --no-unal --no-hd --no-sq --omit-sec-seq -t",
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in str(output) or "Error:" in str(output):
            log.error('\n' + str(output))
            exit()
        # if verbose_log:
        #     log.info("\n" + str(output).strip())
        if os.path.exists(total_seed_sam[0]):
            os.rename(total_seed_sam[0], total_seed_sam[1])
            os.rename(total_seed_file[0], total_seed_file[1])
            # log.info("Mapping finished.")
        else:
            log.error("Cannot find bowtie2 result!")
            exit()


def count_fq_reads(fq_files):
    if type(fq_files) == str:
        fq_files = [fq_files]
    count = 0
    for fq_f in fq_files:
        for line in open(fq_f):
            count += 1
    return int(count/4.)


def main():
    options, log = get_options()
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
        log.info(str(len(all_fq_files)) + " paired fastq files found. " + str(options.round) + " rounds required.")
    else:
        log.info(str(len(all_fq_files)) + " paired fastq files found.")
    thresholds = [int(thr) for thr in options.threshold.split(",")]
    if all_fq_files:
        log.info("\t".join(
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
            log.info("Hit required rounds! Exiting ..")
            break
        real_fq = real_fq_files[go_to]
        bowtie_base = os.path.join(out_base, fq_pairs[0].replace("_1.fq", ""))
        mapping_with_bowtie2(options.fasta, real_fq, bowtie_base, options.resume, options.threads, log)
        this_result = [fq_pairs[0].replace("_1.fq", "")]
        # this round coverage
        this_coverage = get_coverage_from_sam_fast(bowtie_base + ".sam")
        if this_coverage:
            if not this_coverage and not options.round:
                log.info("No more target found! Exiting ..")
                if not options.keep_temp:
                    os.remove(bowtie_base + ".fq")
                    os.remove(bowtie_base + ".sam")
                break
            ref_bowtie = sorted(this_coverage.keys())[0]
            for threshold in thresholds:
                count_site = 0
                for site in range(1, len_ref_seq + 1):
                    if site in this_coverage[ref_bowtie] and this_coverage[ref_bowtie][site] > threshold:
                        count_site += 1
                this_result.append(round(float(count_site) / len_ref_seq, 3))
            results_to_draw.append(this_coverage)
            # merge
            for ref in this_coverage:
                if ref not in all_coverages:
                    all_coverages[ref] = {}
                    for site in this_coverage[ref]:
                        all_coverages[ref][site] = this_coverage[ref][site]
                else:
                    for site in this_coverage[ref]:
                        if site in all_coverages[ref]:
                            all_coverages[ref][site] += this_coverage[ref][site]
                        else:
                            all_coverages[ref][site] = this_coverage[ref][site]
            ref_bowtie = sorted(all_coverages.keys())[0]
            for threshold in thresholds:
                count_site = 0
                for site in range(1, len_ref_seq + 1):
                    if site in all_coverages[ref_bowtie] and all_coverages[ref_bowtie][site] > threshold:
                        count_site += 1
                this_result.append(round(float(count_site) / len_ref_seq, 3))
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
        log.info("\t".join([str(val) for val in this_result]))
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
            y_values_list.append([results_to_draw[i][ref_bowtie].get(site_here, 0) for site_here in x_values])
            y_values_list[-1] = [sum(y_values_list[-1][n:n + width]) / float(len(y_values_list[-1][n:n + width]))
                                 for n in range(0, len(y_values_list[-1]), width)]
        if options.max_cov_tick:
            y_tick_max = int(options.max_cov_tick)
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
