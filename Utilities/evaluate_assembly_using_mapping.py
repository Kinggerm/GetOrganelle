#! /usr/bin/env python
# coding:utf8

from optparse import OptionParser
import os
import sys
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
import GetOrganelleLib
from GetOrganelleLib.pipe_control_func import *
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.sam_parser import *
from GetOrganelleLib.statistical_func import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
from sympy import Interval
import sys
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

try:
    # python2 UnicodeDecodeError ±
    reload(sys)
    sys.setdefaultencoding('utf8')
except NameError:
    pass


def get_options():
    parser = OptionParser("evaluate_assembly_using_mapping.py -f fasta_file -1 RAW_1.fq -2 RAW_2.fq -o output")
    parser.add_option("-f", dest="fasta",
                      help="input assembly fasta file.")
    parser.add_option("-1", dest="original_fq_1")
    parser.add_option("-2", dest="original_fq_2")
    parser.add_option("-u", dest="unpaired_fq_files", default="",
                      help="Input file(s) with unpaired (single-end) reads to be added to the pool. "
                           "files could be comma-separated lists such as 'seq1,seq2'.")
    parser.add_option("-X", "--max-lib-len", dest="max_lib_len", type=int, default=1200,
                      help="Corresponding to '-X' option in Bowtie2. Default: %default.")
    parser.add_option("-c", dest="is_circular", default="auto",
                      help="(yes/no/auto) input fasta is circular. "
                           "If auto was chosen, the input fasta would be treated as circular when the sequence name "
                           "ends with '(circular)'. "
                           "Default: auto")
    parser.add_option("-o", dest="output_base",
                      help="output folder.")
    parser.add_option("-t", dest="threads", type=int, default=2,
                      help="threads.")
    parser.add_option("--continue", dest="resume", default=False, action="store_true")
    parser.add_option("--seed", dest="random_seed", default=12345, type=int,
                      help="Seed for random number generator. Default: %default")
    parser.add_option("--draw", dest="draw_plot", default=False, action="store_true",
                      help="Draw density plot using matplotlib, which should be installed.")
    parser.add_option("--plot-format", dest="plot_format", default="pdf,png",
                      help='Default: pdf,png')
    parser.add_option("--plot-title", dest="plot_title",
                      help="Default: `the file name of the input fasta`")
    parser.add_option("--plot-subtitle", dest="plot_subtitle", default="",
                      help="A 4-space indicates a line break. Default: None")
    parser.add_option("--plot-transparent", dest="plot_transparent", default=False, action="store_true",
                      help="Default: False")
    parser.add_option("--plot-x-density", dest="plot_x_density", default=12000., type=float,
                      help="Default: %default")
    # parser.add_option("--plot-x-sliding-window", dest="sliding_window_size", default=1, type=int,
    #                   help="Default: %default")
    parser.add_option("--plot-x-gap-dots", dest="gap_len", default=3000, type=int,
                      help="Number of sites added in-between isolated contigs. Default: %default")
    parser.add_option("--plot-figure-height", dest="figure_height", default=5., type=float,
                      help="Default: %default")
    parser.add_option("--plot-y-lim", dest="y_lim", type=float,
                      help="Y axis value limit. ")
    # parser.add_option("--plot-figure-extra-width", dest="extra_width", default=3., type=float,
    #                   help="Default: %default")
    parser.add_option("--plot-font", dest="plot_font", default=None,
                      help="For plot of unicode characters in some environments. Use 'Times New Roman','Arial' etc. "
                           "Default: %default.")
    parser.add_option("--disable-customized-error-rate", dest="customized_error_rate", default=True,
                      action="store_true")
    parser.add_option("--which-bowtie2", dest="which_bowtie2", default="",
                      help="Assign the path to Bowtie2 binary files if not added to the path. "
                      "Default: try GetOrganelleDep/" + SYSTEM_NAME + "/bowtie2 first, then $PATH")
    parser.add_option("--bowtie2-mode", dest="bowtie2_mode", default="--sensitive",
                      help="Default: %default")
    parser.add_option("--bowtie2-options", dest="other_bowtie2_options", default="--no-discordant --dovetail",
                      help="Default: %default")
    parser.add_option("--stat-mode", dest="stat_mode", default="best",
                      help="Statistical mode for counting multiple hits of a single read: best/all. "
                           "The all mode is meaningful only when '-k <INT>' was included in '--bowtie2-options'. "
                           "Default: %default")
    parser.add_option("--debug", dest="debug_mode", default=False, action="store_true",
                      help="Turn on debug mode.")
    options, argv = parser.parse_args()
    if not (options.fasta and
            ((options.original_fq_1 and options.original_fq_2) or options.unpaired_fq_files)
            and options.output_base):
        sys.stderr.write("Insufficient arguments!\n")
        sys.exit()
    if not os.path.isdir(options.output_base):
        os.mkdir(options.output_base)
    if options.debug_mode:
        log_level = "DEBUG"
    else:
        log_level = "INFO"
    assert options.stat_mode in ("best", "all")
    log_handler = simple_log(logging.getLogger(), options.output_base, "", log_level=log_level)
    log_handler.info("")
    log_handler.info("Python " + str(sys.version).replace("\n", " "))
    # log versions of python libs
    lib_versions_info = []
    if options.draw_plot:
        try:
            import matplotlib
        except ImportError:
            pass
        else:
            lib_versions_info.append("matplotlib " + matplotlib.__version__)
    lib_versions_info.append("GetOrganelleLib " + GetOrganelleLib.__version__)
    log_handler.info("PYTHON LIBS: " + "; ".join(lib_versions_info))
    # log versions of dependencies
    dep_versions_info = []
    if not options.which_bowtie2:
        try_this_bin = os.path.join(GO_DEP_PATH, "bowtie2", "bowtie2")
        if os.path.isfile(try_this_bin) and executable(try_this_bin):
            options.which_bowtie2 = os.path.split(try_this_bin)[0]
    if not executable(os.path.join(options.which_bowtie2, "bowtie2")):
        log_handler.error(os.path.join(options.which_bowtie2, "bowtie2") + " not accessible!")
        exit()
    else:
        output, err = subprocess.Popen(
            os.path.join(options.which_bowtie2, "bowtie2") + " --version", stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True).communicate()
        this_lines = output.decode("utf8").split("\n")[:3]
        dep_versions_info.append("Bowtie2 " + this_lines[0].split()[-1].strip())
    if not executable(os.path.join(options.which_bowtie2, "bowtie2-build") + " --large-index"):
        log_handler.error(os.path.join(options.which_bowtie2, "bowtie2-build") + " not accessible!")
        exit()
    log_handler.info("DEPENDENCIES: " + "; ".join(dep_versions_info))
    log_handler.info("WORKING DIR: " + os.getcwd())
    # if not executable(os.path.join(options.which_bowtie2, "bowtie2-build-l")):
    #     log_handler.error(os.path.join(options.which_bowtie2, "bowtie2-build-l") + " not accessible!")
    #     exit()
    log_handler.info(" ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n")
    log_handler = timed_log(log_handler, options.output_base, "", log_level=log_level)
    return options, log_handler


def modify_fasta(original_fasta, new_fasta, is_circular, max_lib_len):
    # duplicated seq names would cause error in downstream analysis
    count_name_freq = {}
    fasta_ob = SequenceList(original_fasta)
    for record in fasta_ob:
        if record.label in count_name_freq:
            count_name_freq[record.label] += 1
        else:
            count_name_freq[record.label] = 1
        record.seq = record.seq.replace("*", "N").replace("-", "N")
    duplicated_name_go = {seq_name: 0 for seq_name in count_name_freq if count_name_freq[seq_name] > 1}

    #
    if is_circular == "yes":
        for record in fasta_ob:
            if len(record.seq):
                record.seq = re_linear_circular_seqs(record.seq)
                to_add = record.seq[:max_lib_len]
                added_len = len(to_add)
                record.seq += to_add
                # in case ref is extremely short
                while added_len < max_lib_len:
                    to_add = record.seq[:(max_lib_len - added_len)]
                    added_len += len(to_add)
                    record.seq += to_add
            if record.label in duplicated_name_go:
                duplicated_name_go[record.label] += 1
                record.label = record.label.split(" ")[0] + "--" + str(duplicated_name_go[record.label])  # \
                # + " ".join(record.label.split(" ")[1:])
        fasta_ob.write_fasta(new_fasta)
    elif is_circular == "auto":
        for record in fasta_ob:
            if len(record.seq) and record.label.endswith("(circular)"):
                record.seq = re_linear_circular_seqs(record.seq)
                to_add = record.seq[:max_lib_len]
                added_len = len(to_add)
                record.seq += to_add
                # in case ref is extremely short
                while added_len < max_lib_len:
                    to_add = record.seq[:(max_lib_len - added_len)]
                    added_len += len(to_add)
                    record.seq += to_add
            if record.label in duplicated_name_go:
                duplicated_name_go[record.label] += 1
                record.label = record.label.split(" ")[0] + "--" + str(duplicated_name_go[record.label])  # \
                # + " ".join(record.label.split(" ")[1:])
        fasta_ob.write_fasta(new_fasta)
    else:
        for record in fasta_ob:
            if record.label in duplicated_name_go:
                duplicated_name_go[record.label] += 1
                record.label = record.label.split(" ")[0] + "--" + str(duplicated_name_go[record.label])  # \
                # + " ".join(record.label.split(" ")[1:])
        fasta_ob.write_fasta(new_fasta)


def get_lengths_with_seq_names_modified(raw_fasta_file, log_handler=None):
    # duplicated seq names would cause error in downstream analysis
    count_name_freq = {}
    fasta_ob = SequenceList(raw_fasta_file)
    for record in fasta_ob:
        if record.label in count_name_freq:
            count_name_freq[record.label] += 1
        else:
            count_name_freq[record.label] = 1
    duplicated_name_go = {seq_name: 0 for seq_name in count_name_freq if count_name_freq[seq_name] > 1}
    for record in fasta_ob:
        if record.label in duplicated_name_go:
            duplicated_name_go[record.label] += 1
            record.label = record.label.split(" ")[0] + "--" + str(duplicated_name_go[record.label])
        else:
            record.label = record.label.split(" ")[0]
    if log_handler:
        lengths = [len(rc.seq) for rc in fasta_ob]
        log_handler.info("Reference length: " + str(sum(lengths)) + " (" + ", ".join([str(l) for l in lengths]) + ")")
    return {record.label: len(record.seq) for record in fasta_ob}


def adjust_vertically_in_one_line(max_val, min_val=0, soft_min_gap=20., *y_positions):
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


def adjust_vertically_in_different_lines(middle_y, min_graph_dist, x_factor=1., y_factor=1., *sorted_x_y_positions):
    x_y_positions = deepcopy(sorted_x_y_positions)
    for go_p, (x_pos, y_pos) in enumerate(x_y_positions):
        if go_p > 0 and (x_pos - x_y_positions[go_p - 1][0]) * x_factor < min_graph_dist:
            go_back = go_p - 1
            constraints = Interval(-inf, inf)
            while go_back >= 0 and (x_pos - x_y_positions[go_back][0]) * x_factor < min_graph_dist:
                prev_x = x_y_positions[go_back][0]
                if len(x_y_positions[go_back][1]):
                    prev_y = sum(x_y_positions[go_back][1]) / len(x_y_positions[go_back][1])
                else:
                    prev_y = 0.
                y_move = (min_graph_dist ** 2 - ((x_pos - prev_x) * x_factor) ** 2) ** 0.5 / y_factor
                constraints &= (Interval(-inf, prev_y - y_move) | Interval(prev_y + y_move, inf))
                go_back -= 1
            if middle_y not in constraints:
                new_average_y = sorted(constraints.boundary, key=lambda cons_y: abs(cons_y-middle_y))[0]
                if len(y_pos):
                    old_y = sum(y_pos) / len(y_pos)
                else:
                    old_y = 0.
                x_y_positions[go_p][1] = x_y_positions[go_p][1] - (old_y - new_average_y)
    return x_y_positions


def main():
    options, log_handler = get_options()
    try:
        new_fasta = os.path.join(options.output_base, "modified.fasta")
        if not (options.resume and os.path.exists(new_fasta)):
            modify_fasta(options.fasta, new_fasta, options.is_circular, max_lib_len=options.max_lib_len)
        unpaired_fq_files = []
        if options.unpaired_fq_files:
            unpaired_fq_files = options.unpaired_fq_files.split(",")
        other_bowtie2_options = " -X " + str(options.max_lib_len) + " " + options.other_bowtie2_options + " "
        bowtie2_mode = " " + options.bowtie2_mode + " "
        map_with_bowtie2(seed_file=new_fasta,
                         original_fq_files=unpaired_fq_files,
                         bowtie_out=os.path.join(options.output_base, "check"),
                         resume=options.resume, threads=options.threads, random_seed=options.random_seed,
                         silent=False or options.debug_mode, log_handler=log_handler, verbose_log=options.debug_mode,
                         which_bowtie2=options.which_bowtie2, bowtie2_other_options=other_bowtie2_options,
                         fq_1=options.original_fq_1, fq_2=options.original_fq_2, bowtie2_mode=bowtie2_mode)
        ref_lengths = get_lengths_with_seq_names_modified(options.fasta, log_handler)
        mapping_records = MapRecords(sam_file=os.path.join(options.output_base, "check.sam"),
                                     ref_real_len_dict=ref_lengths)
        sequence_statistics = \
            mapping_records.get_customized_mapping_characteristics(multiple_hits_mode=options.stat_mode)
        num_mapped_reads = mapping_records.get_number_of_mapped_reads()
        num_paired, num_single = num_mapped_reads["paired"], num_mapped_reads["single"]

        if options.draw_plot:
            import matplotlib
            matplotlib.use('Agg')
            if options.plot_font:
                matplotlib.rc('font', family=options.plot_font)
            import matplotlib.pyplot as plt
        else:
            plt = None

        # make data and default settings
        gap_len = options.gap_len
        extra_width = 3.  # options.extra_width
        sliding_w_size = 1   # options.sliding_window_size
        x_data_len = gap_len * (len(mapping_records.references) - 1) \
                     + sum([mapping_records.references[ref]["real_len"] for ref in mapping_records.references])
        fig_width = extra_width + x_data_len / options.plot_x_density
        fig_height = options.figure_height
        extra_percent = extra_width / fig_width
        title_height_percent = 0.09
        add_extra_to_left = 0.27   # for extra_width==3
        plot_area_l, plot_area_r = extra_percent * add_extra_to_left, 1 - extra_percent * (1 - add_extra_to_left)
        plot_area_b, plot_area_t = title_height_percent, 1 - title_height_percent
        cigar_chars = ["M", "X", "I", "D"]
        cigar_char_dict = {"M": "Aligned", "X": "Mismatched", "I": "Inserted", "D": "Deleted"}
        color_used = {"M": [(0.133, 0.616, 0.361), 0.5],
                      "X": [(0.145, 0.651, 0.961), 0.3],
                      "I": [(0.996, 0.804, 0.322), 0.8],
                      "D": [(0.831, 0.310, 0.275), 0.5]}

        x_data = np.array(range(x_data_len))
        y_data = {}
        # log mean and std
        y_stat = dict()
        max_y_dat = 0
        # start @20190304: add extra error rate @20190304
        err_all_cover = np.array([])
        err_subset_cover = [np.array([]) for foo in range(len(sequence_statistics["M"]))]
        # end @20190304
        for cigar_char in cigar_chars:
            y_stat[cigar_char] = {"subset": [], "all": []}
            this_cover = []
            this_cover_for_stat = []
            ref_list = sorted(list(sequence_statistics[cigar_char]))
            go_to_subset = 0
            while ref_list:
                ref = ref_list.pop(0)
                smoothed_cover_per_site = sequence_statistics[cigar_char][ref]
                this_mean, this_std = float(np.average(smoothed_cover_per_site)), float(np.std(smoothed_cover_per_site))

                y_stat[cigar_char]["subset"].append((this_mean, this_std, len(smoothed_cover_per_site)))
                this_cover_for_stat.extend(smoothed_cover_per_site)
                # start @20190304
                if options.customized_error_rate:
                    if cigar_char in {"X", "I", "D"}:
                        if len(err_subset_cover[go_to_subset]):
                            err_subset_cover[go_to_subset] += np.array(smoothed_cover_per_site)
                        else:
                            err_subset_cover[go_to_subset] = np.array(smoothed_cover_per_site)
                # end @20190304
                if sliding_w_size != 1:
                    new_averaged_cover = []
                    for j in range(len(smoothed_cover_per_site)):
                        if j % sliding_w_size:
                            new_averaged_cover.append(0)
                        else:
                            new_averaged_cover.append(np.average(smoothed_cover_per_site[j: j + sliding_w_size]))
                    smoothed_cover_per_site = np.array(new_averaged_cover)
                this_cover.extend(smoothed_cover_per_site)
                max_y_dat = max(max(smoothed_cover_per_site), max_y_dat)
                if ref_list:
                    this_cover.extend([0] * gap_len)
                go_to_subset += 1
            y_data[cigar_char] = np.ma.masked_where(np.array(this_cover) <= 0, this_cover)
            y_stat[cigar_char]["all"] = float(np.average(this_cover_for_stat)), float(np.std(this_cover_for_stat))
            # start @20190304
            if options.customized_error_rate:
                if cigar_char in {"X", "I", "D"}:
                    if len(err_all_cover):
                        err_all_cover += np.array(this_cover_for_stat)
                    else:
                        err_all_cover = np.array(this_cover_for_stat)
            # end @20190304

        if not max_y_dat:
            raise ValueError("No mapped reads found!")
        if options.y_lim:
            max_y_dat = options.y_lim

        # start @20190304
        if options.customized_error_rate:
            y_stat["error"] = {"all": [np.average(err_all_cover) / y_stat["M"]["all"][0],
                                       np.std(err_all_cover) / y_stat["M"]["all"][0]],
                               "subset": [[np.average(err_subset_cover[go_to_sb]) / y_stat["M"]["subset"][go_to_sb][0],
                                           np.std(err_subset_cover[go_to_sb]) / y_stat["M"]["subset"][go_to_sb][0],
                                           y_stat["M"]["subset"][go_to_sb][2]]
                                          for go_to_sb in range(len(sequence_statistics["M"]))]}
        # end @20190304

        # create figure
        if options.draw_plot:
            fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
            ax.spines['top'].set_visible(False)
            # ax.spines['bottom'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # ax.spines['left'].set_visible(False)
            plt.setp(ax.get_xticklabels(), fontsize=12)
            plt.setp(ax.get_yticklabels(), fontsize=12)
            fig.subplots_adjust(left=plot_area_l, right=plot_area_r, bottom=plot_area_b, top=plot_area_t)

            lines = plt.plot(x_data, y_data["M"], 'o',
                             x_data, y_data["X"], 'o',
                             x_data, y_data["I"], 'o',
                             x_data, y_data["D"], 'o')
            for go_to, this_char in enumerate(cigar_chars):
                plt.setp(lines[go_to], color=color_used[this_char][0], alpha=color_used[this_char][1], markersize=0.2)
            plt.title("  " + (options.plot_title if options.plot_title else options.fasta), fontsize=18, loc='left')
            subtitle_x_pos = x_data_len + (1 - plot_area_r) * fig_width * options.plot_x_density
            subtitle_y_pos = max_y_dat * (1 + (1. / (1 - 2 * title_height_percent) - 1) / 8)

            mapped_str = []
            if num_paired:
                mapped_str.append("# mapped pairs: " + str(int(num_paired / 2)))
            if num_single:
                mapped_str.append("# mapped reads: " + str(int(num_paired / 2)) + "×2+" + str(num_single))
            for subtitle_str in [sub_str.strip() for sub_str in options.plot_subtitle.split("    ")] + mapped_str:
                plt.text(subtitle_x_pos, subtitle_y_pos, subtitle_str, fontsize=12, alpha=0.7,
                         horizontalalignment='right', verticalalignment='center', multialignment='center')
                subtitle_y_pos -= max_y_dat / options.figure_height / 4.  # for fontsize==2

        # write to log
        line_labels = {c_: full_name + ": " + "%.2f" % y_stat[c_]["all"][0] + "±" + "%.2f" % y_stat[c_]["all"][1]
                       for c_, full_name in cigar_char_dict.items()}
        echo_statistics = [
            line_labels[cigar_char] +
            " (" + ", ".join(["%.2f" % here_mean + "±" + "%.2f" % here_std
                              for here_mean, here_std, here_len in y_stat[cigar_char]["subset"]]) + ")"
            for cigar_char in cigar_chars]
        if options.customized_error_rate:
            echo_statistics.append("Customized error rate: " +
                                   "%.4f" % y_stat["error"]["all"][0] + "±" + "%.4f" % y_stat["error"]["all"][1] +
                                   " (" +
                                   ", ".join(["%.4f" % here_mean + "±" + "%.4f" % here_std
                                              for here_mean, here_std, here_len in y_stat["error"]["subset"]]) +
                                   ")")
        echo_statistics.insert(0, "# mapped pairs: " + str(int(num_paired / 2)))
        echo_statistics.insert(0, "# mapped reads: " + str(int(num_paired / 2)) + "×2+" + str(num_single))
        for log_line in echo_statistics:
            log_handler.info(log_line)

        # plot txt
        if options.draw_plot:
            new_y_pos = adjust_vertically_in_one_line(max_y_dat, 0, max_y_dat / 20,
                                                      y_stat["M"]["subset"][-1][0],
                                                      y_stat["X"]["subset"][-1][0],
                                                      y_stat["I"]["subset"][-1][0],
                                                      y_stat["D"]["subset"][-1][0])
            for go_to, this_char in enumerate(cigar_chars):
                plt.text(x_data_len + 0.05 * options.plot_x_density, new_y_pos[go_to], line_labels[this_char],
                         color=color_used[this_char][0], fontsize=12)
            if y_stat["M"]["subset"][-1][0] and max_y_dat / y_stat["M"]["subset"][-1][0] > 2.:
                middle_pos = (y_stat["M"]["subset"][-1][0] + y_stat["M"]["subset"][-1][1] + max_y_dat) / 2
            else:
                middle_pos = (y_stat["M"]["subset"][-1][0] - y_stat["M"]["subset"][-1][1] +
                              y_stat["X"]["subset"][-1][0] + y_stat["X"]["subset"][-1][1]) / 2
            middle_set = np.array([middle_pos + max_y_dat/15, middle_pos + max_y_dat/45,
                                   middle_pos - max_y_dat/45, middle_pos - max_y_dat/15])
            x_y_pos = []
            this_accumulated_pos = 0
            for this_mean_cov, this_std_cov, this_len in y_stat["M"]["subset"]:
                x_y_pos.append([this_accumulated_pos + this_len / 2, deepcopy(middle_set)])
                this_accumulated_pos += this_len + gap_len
            width_over_x = (fig_width - extra_width) / x_data_len
            height_over_y = fig_height / max_y_dat
            new_x_y_pos = adjust_vertically_in_different_lines(np.average(middle_set),
                                                               (fig_height**2 + (fig_width - extra_width)**2) ** 0.5 / 10,
                                                               width_over_x, height_over_y,
                                                               *x_y_pos)
            label_x_offset = -0.3 / width_over_x
            for go_to, this_char in enumerate(cigar_chars):
                # accumulated_pos = 0
                for go_subset, (this_mean_cov, this_std_cov, this_len) in enumerate(y_stat[this_char]["subset"]):
                    this_x, this_y = new_x_y_pos[go_subset]
                    plt.text(this_x + label_x_offset, this_y[go_to],
                             "%.2f" % this_mean_cov + "±" + "%.2f" % this_std_cov,
                             color=color_used[this_char][0], fontsize=9)
                    # accumulated_pos += this_len + gap_len
            for go_subset, (this_mean_cov, this_std_cov, this_len) in enumerate(y_stat["M"]["subset"]):
                this_x, this_y = new_x_y_pos[go_subset]
                plt.text(this_x + label_x_offset, this_y[0] + max_y_dat * 2/45,
                         str(this_len) + " bp",
                         color="black", fontsize=9)

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
    except Exception as e:
        log_handler.exception(str(e))
    logging.shutdown()


if __name__ == '__main__':
    main()
