#! /usr/bin/env python

from optparse import OptionParser
import os
import sys
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.pipe_control_func import *
from Library.seq_parser import *
from Library.sam_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
from math import ceil


def get_options():
    parser = OptionParser("evaluate_assembly_using_mapping.py -f fasta_file -1 RAW_1.fq -2 RAW_2.fq -o output")
    parser.add_option("-f", dest="fasta",
                      help="input assembly fasta file.")
    parser.add_option("-1", dest="original_fq_1")
    parser.add_option("-2", dest="original_fq_2")
    parser.add_option("--max-lib-len", dest="max_lib_len", type=int, default=1200,
                      help="default: %default.")
    parser.add_option("-c", dest="is_circular", default=False, action="store_true",
                      help="input assembly result is circular.")  # detect by endswith "circular"
    parser.add_option("-o", dest="output_base",
                      help="output folder.")
    parser.add_option("-t", dest="threads", type=int, default=2,
                      help="threads.")
    parser.add_option("--continue", dest="resume", default=False, action="store_true")
    parser.add_option("--draw", dest="draw_plot", default=False, action="store_true",
                      help="Draw density plot using matplotlib, which should be installed.")
    parser.add_option("--plot-format", dest="plot_format", default="pdf,png",
                      help='Default: pdf,png')
    parser.add_option("--plot-title", dest="plot_title",
                      help="Default: `the file name of the input fasta`")
    parser.add_option("--plot-transparent", dest="plot_transparent", default=False,
                      help="Default: False")
    parser.add_option("--max-coverage-tick", dest="max_cov_tick", type=int,
                      help="default: auto-estimated.")
    options, argv = parser.parse_args()
    if not (options.fasta and options.original_fq_1 and options.original_fq_2 and options.output_base):
        sys.stderr.write("Insufficient arguments!\n")
        sys.exit()
    if not os.path.isdir(options.output_base):
        os.mkdir(options.output_base)
    log = simple_log(logging.getLogger(), options.output_base, "")
    log.info("")
    log.info(' '.join(sys.argv) + '\n')
    log = timed_log(log, options.output_base, "")
    return options, log


def mapping_with_bowtie2(seed_file, original_fq_1, original_fq_2, bowtie_out, max_lib_len, resume, threads, log):
    if not (os.path.exists(seed_file + '.index.1.bt2l')):
        build_seed_index = subprocess.Popen("bowtie2-build --large-index " + seed_file + " " + seed_file + '.index',
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = build_seed_index.communicate()
        if "unrecognized option" in str(output):
            build_seed_index = subprocess.Popen("bowtie2-build " + seed_file + " " + seed_file + '.index',
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = build_seed_index.communicate()
        if "(ERR)" in str(output) or "Error:" in str(output) or "error" in str(output):
            log.error('\n' + str(output))
            exit()
    seed_index_base = seed_file + '.index'
    res_path_name, res_base_name = os.path.split(bowtie_out)
    total_seed_sam = [os.path.join(res_path_name, x + res_base_name + ".sam") for x in ("temp.", "")]
    if not (resume and os.path.exists(total_seed_sam[1])):
        make_seed_bowtie2 = subprocess.Popen(
            "bowtie2 --mm -p " + str(threads) + " -X " + str(max_lib_len) + " --no-discordant --dovetail" +
            " --sensitive -x " + seed_index_base + " -1 " + original_fq_1 + " -2 " + original_fq_2 +
            " -S " + total_seed_sam[0] + " --no-unal --omit-sec-seq -t",
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in str(output) or "Error:" in str(output):
            log.error('\n' + str(output))
            exit()
        if os.path.exists(total_seed_sam[0]):
            os.rename(total_seed_sam[0], total_seed_sam[1])
        else:
            log.error("Cannot find bowtie2 result!")
            exit()


def modify_fasta(original_fasta, new_fasta, is_circular, max_lib_len):
    if is_circular:
        fasta_ob = SequenceList(original_fasta)
        for record in fasta_ob:
            if len(record.seq):
                to_add = record.seq[:max_lib_len]
                added_len = len(to_add)
                record.seq += to_add
                # in case ref is extremely short
                while added_len < max_lib_len:
                    to_add = record.seq[:(max_lib_len - added_len)]
                    added_len += len(to_add)
                    record.seq += to_add
        fasta_ob.write_fasta(new_fasta)
    else:
        os.system("cp " + original_fasta + " " + new_fasta)


def adjust_horizontally(max_val, min_val=0, soft_min_gap=20, *y_positions):
    record_original_order = [[this_y_pos, raw_order] for raw_order, this_y_pos in enumerate(y_positions)]
    record_original_order.sort()
    len_val = len(record_original_order)
    go_val = 0
    while go_val + 1 < len_val:
        if record_original_order[go_val][0] + soft_min_gap > record_original_order[go_val + 1][0]:
            record_original_order[go_val + 1][0] = record_original_order[go_val][0] + soft_min_gap
            go_change = go_val + 2
            while go_change < len_val:
                if record_original_order[go_change - 1][0] + soft_min_gap > record_original_order[go_change][0]:
                    record_original_order[go_change][0] = record_original_order[go_change - 1][0] + soft_min_gap
                    go_change += 1
                else:
                    break
        go_val += 1
    # push out
    if record_original_order[-1][0] > max_val:
        record_original_order[-1][0] = max_val
        record_original_order.sort(reverse=True)
        go_val = 0
        while go_val + 1 < len_val:
            if record_original_order[go_val][0] - soft_min_gap < record_original_order[go_val + 1][0]:
                record_original_order[go_val + 1][0] = record_original_order[go_val][0] - soft_min_gap
                go_change = go_val + 2
                while go_change < len_val:
                    if record_original_order[go_change - 1][0] - soft_min_gap < record_original_order[go_change][0]:
                        record_original_order[go_change][0] = record_original_order[go_change - 1][0] - soft_min_gap
                        go_change += 1
                    else:
                        break
            go_val += 1
        # push back, mean
        if record_original_order[0][0] < min_val:
            mean_dist = float(max_val - min_val) / (len_val - 1)
            record_original_order[0][0] = max_val
            record_original_order[-1][0] = min_val
            for go_val in range(1, len_val - 1):
                record_original_order[go_val][0] = max_val - go_val * mean_dist
    # sort by original order
    record_original_order.sort(key=lambda x: x[1])
    return [new_val[0] for new_val in record_original_order]


def main():
    options, log_handler = get_options()
    new_fasta = os.path.join(options.output_base, "modified.fasta")
    if not (options.resume and os.path.exists(new_fasta)):
        modify_fasta(options.fasta, new_fasta, options.is_circular, max_lib_len=options.max_lib_len)
    mapping_with_bowtie2(seed_file=new_fasta, original_fq_1=options.original_fq_1, original_fq_2=options.original_fq_2,
                         bowtie_out=os.path.join(options.output_base, "check"), max_lib_len=options.max_lib_len,
                         resume=options.resume, threads=options.threads, log=log_handler)
    ref_lengths = {record.label.split()[0]: len(record.seq) for record in SequenceList(options.fasta)}
    mapping_records = MapRecords(sam_file=os.path.join(options.output_base, "check.sam"), ref_real_len_dict=ref_lengths)
    sequence_statistics = mapping_records.get_customized_mapping_characteristics()

    if options.draw_plot:
        import matplotlib.pyplot as plt
        # make data and default settings
        gap_len = 1000
        slide_window_size = 1
        x_data_len = gap_len * (len(mapping_records.references) - 1) \
                     + sum([ceil(mapping_records.references[ref]["real_len"]/float(slide_window_size))
                            for ref in mapping_records.references])
        fig_width, fig_height = x_data_len * slide_window_size / 10000, 5
        plot_area_l, plot_area_r, plot_area_b, plot_area_t = 0.06, 0.88, 0.09, 0.91
        cigar_chars = ["M", "X", "I", "D"]
        cigar_char_dict = {"M": "Aligned", "X": "Mismatched", "I": "Inserted", "D": "Deleted"}
        color_used = {"M": [(0.133, 0.616, 0.361), 0.5],
                      "X": [(0.145, 0.651, 0.961), 0.3],
                      "I": [(0.996, 0.804, 0.322), 0.8],
                      "D": [(0.831, 0.310, 0.275), 0.5]}

        x_data = np.array(range(x_data_len))
        y_data = {}
        # log mean and std
        y_stat = {}
        for cigar_char in cigar_chars:
            y_stat[cigar_char] = {"subset": [], "all": []}
            this_cover = []
            this_cover_for_stat = []
            ref_list = sorted(list(sequence_statistics[cigar_char]))
            while ref_list:
                ref = ref_list.pop(0)
                averaged_cover = sequence_statistics[cigar_char][ref]
                this_mean, this_std = round(np.average(averaged_cover), 2), round(float(np.std(averaged_cover)), 2)
                y_stat[cigar_char]["subset"].append((this_mean, this_std, len(averaged_cover)))
                if slide_window_size != 1:
                    averaged_cover = np.array([np.average(averaged_cover[j: j + slide_window_size])
                                               for j in range(0, len(averaged_cover), slide_window_size)])
                this_cover.extend(averaged_cover)
                this_cover_for_stat.extend(averaged_cover)
                if ref_list:
                    this_cover.extend([0] * gap_len)
            y_data[cigar_char] = np.ma.masked_where(np.array(this_cover) <= 0, this_cover)
            y_stat[cigar_char]["all"] = \
                round(np.average(this_cover_for_stat), 2), round(float(np.std(this_cover_for_stat)), 2)
        max_y_dat = max([max(y_data[temp_sub]) for temp_sub in y_data])

        # create figure
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        fig.subplots_adjust(left=plot_area_l, right=plot_area_r, bottom=plot_area_b, top=plot_area_t)

        lines = plt.plot(x_data, y_data["M"], 'o',
                         x_data, y_data["X"], 'o',
                         x_data, y_data["I"], 'o',
                         x_data, y_data["D"], 'o')
        for go_to, this_char in enumerate(cigar_chars):
            plt.setp(lines[go_to], color=color_used[this_char][0], alpha=color_used[this_char][1], markersize=0.18)  # linewidth=0.3
        plt.title(options.plot_title if options.plot_title else options.fasta, fontsize=18)

        # write to log
        line_labels = {c_: full_name + ": " + "%.2f" % y_stat[c_]["all"][0] + "±" + "%.2f" % y_stat[c_]["all"][1]
                       for c_, full_name in cigar_char_dict.items()}
        echo_statistics = [
            line_labels[cigar_char] +
            " (" + ", ".join(["%.2f" % here_mean + "±" + "%.2f" % here_std
                              for here_mean, here_std, here_len in y_stat[cigar_char]["subset"]]) + ")"
            for cigar_char in cigar_chars]
        for log_line in echo_statistics:
            log_handler.info(log_line)

        # plot txt
        new_y_pos = adjust_horizontally(max_y_dat, 0, 20,
                                        y_stat["M"]["subset"][-1][0],
                                        y_stat["X"]["subset"][-1][0],
                                        y_stat["I"]["subset"][-1][0],
                                        y_stat["D"]["subset"][-1][0])
        for go_to, this_char in enumerate(cigar_chars):
            plt.text(x_data_len * 1.01, new_y_pos[go_to], line_labels[this_char],
                     color=color_used[this_char][0], fontsize=12)
        middle_pos = (y_stat["M"]["subset"][-1][0] - y_stat["M"]["subset"][-1][1] +
                      y_stat["X"]["subset"][-1][0] + y_stat["X"]["subset"][-1][1]) / 2
        middle_set = middle_pos + 27, middle_pos + 9, middle_pos - 9, middle_pos - 27
        for go_to, this_char in enumerate(cigar_chars):
            accumulated_pos = 0
            for this_mean_cov, this_std_cov, this_len in y_stat[this_char]["subset"]:
                plt.text(accumulated_pos + this_len / 2, middle_set[go_to],
                         "%.2f" % this_mean_cov + "±" + "%.2f" % this_std_cov,
                         color=color_used[this_char][0], fontsize=9)
                accumulated_pos += this_len + gap_len

        # density plot could also be added follow this
        # a = plt.axes([.65, .6, .2, .2], facecolor='k')
        # n, bins, patches = plt.hist(s, 400, density=True)
        # plt.title('Probability')
        # plt.xticks([])
        # plt.yticks([])

        for plot_format in options.plot_format.split(","):
            plt.savefig(os.path.join(options.output_base, "mapping." + plot_format),
                        transparent=options.plot_transparent, dpi=300)

    log_handler = simple_log(log_handler, options.output_base, "")
    log_handler.info("")


if __name__ == '__main__':
    main()
