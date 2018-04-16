#! /usr/bin/env python

from optparse import OptionParser
import os
import sys
import subprocess
import logging


def get_options():
    parser = OptionParser("this_script.py -f fasta_file -s sam_file")
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
    parser.add_option("--draw", dest="draw_plot", default=False, action="store_true")
    # parser.add_option("--average", default=False, action="store_true",
    #                   help="output average coverage.")
    options, argv = parser.parse_args()
    if not (options.fasta and options.initial_mapped and options.output_base and options.output_per_round_dir):
        sys.stderr.write("Insufficient arguments!\n")
        sys.exit()
    if not os.path.isdir(options.output_base):
        os.mkdir(options.output_base)
    log = simple_log(logging.getLogger(), options.output_base)
    log.info("")
    log.info(' '.join(sys.argv) + '\n')
    log = timed_log(log, options.output_base)
    return options, log


def simple_log(log, output_base):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_simple.addHandler(console)
    log_simple.addHandler(logfile)
    return log_simple


def timed_log(log, output_base):
    log_timed = log
    for handler in list(log_timed.handlers):
        log_timed.removeHandler(handler)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_timed.addHandler(console)
    log_timed.addHandler(logfile)
    return log_timed


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
            "bowtie2 -p " + str(threads) + " --very-fast-local --al " + total_seed_file[
                0] + " -x " + seed_index_base + " -U " +
            ",".join(original_fq_files) + " -S " + total_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
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


# does not deal with multiple hits!
def get_coverage(bowtie_sam_file):
    coverage = {}
    for line in open(bowtie_sam_file):
        if line.strip() and not line.startswith('@'):
            line_split = line.strip().split('\t')
            start_position = int(line_split[3])
            if start_position:
                flag = int(line_split[1])
                reference = line_split[2]
                direction = 1 if flag % 32 < 16 else -1
                read_len = len(line_split[9])
                if reference in coverage:
                    for position in range(max(1, start_position), max(1, start_position + read_len*direction), direction):
                        if position in coverage[reference]:
                            coverage[reference][position] += 1
                        else:
                            coverage[reference][position] = 1
                else:
                    coverage[reference] = {}
                    for position in range(max(1, start_position), max(1, start_position + read_len * direction), direction):
                        coverage[reference][position] = 1
    return coverage


def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    interleaved = 0
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip('\n').strip('\r'))
            this_seq = ''
            this_line = fasta_file.readline()
            seq_line_count = 0
            while this_line and not this_line.startswith('>'):
                if seq_line_count == 1:
                    interleaved = len(this_seq)
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
                seq_line_count += 1
            seqs.append(this_seq)
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs, interleaved]


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
    count_round = 0
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
        if options.round and count_round > options.round:
            log.info("Hit required rounds! Exiting ..")
            break
        real_fq = real_fq_files[go_to]
        bowtie_base = os.path.join(out_base, fq_pairs[0].replace("_1.fq", ""))
        mapping_with_bowtie2(options.fasta, real_fq, bowtie_base, options.resume, options.threads, log)
        this_result = [fq_pairs[0].replace("_1.fq", "")]
        # this round coverage
        this_coverage = get_coverage(bowtie_base + ".sam")
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
        this_result.append(other_reads_num/float(all_reads_num))
        if not options.keep_temp:
            os.remove(bowtie_base + ".fq")
            os.remove(bowtie_base + ".sam")
        log.info("\t".join([str(val) for val in this_result]))
    # draw
    if options.draw_plot:
        import matplotlib.pyplot as plt
        len_round = len(results_to_draw)
        figure = plt.figure(figsize=(50, len_round*3))
        width = 50
        for i in range(len_round):
            plt.subplot(len_round, 1, i + 1)
            ref_bowtie = sorted(results_to_draw[i].keys())[0]
            X = list(range(1, len_ref_seq + 1))
            Y = [results_to_draw[i][ref_bowtie].get(site_here, 0) for site_here in X]
            X = [X[n] for n in range(0, len(X), width)]
            Y = [sum(Y[n:n+width])/float(len(Y[n:n+width])) for n in range(0, len(Y), width)]
            plt.bar(X, Y, width=width, color="darkgreen")
            # plt.hist([1, 2, 3, 5, 2, 1], color="darkgreen")
            plt.yticks([100, 200, 300, 400, 500])
        figure.savefig(os.path.join(out_base, "coverage_per_round.pdf"), bbox_inches="tight")


if __name__ == '__main__':
    main()
