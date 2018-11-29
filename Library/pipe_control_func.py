import signal
import time
import logging
import sys
import os
from multiprocessing import Pool

major_version, minor_version = sys.version_info[:2]
if major_version == 2 and minor_version >= 7:
    python_version = "2.7+"
elif major_version == 3 and minor_version >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)

if python_version == "2.7+":
    from commands import getstatusoutput
else:
    from subprocess import getstatusoutput
import subprocess
dead_code = {"2.7+": 32512, "3.5+": 127}[python_version]

path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]


def simple_log(log, output_base, prefix):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, prefix + 'log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_simple.addHandler(console)
    log_simple.addHandler(logfile)
    return log_simple


def timed_log(log, output_base, prefix):
    log_timed = log
    for handler in list(log_timed.handlers):
        log_timed.removeHandler(handler)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, prefix + 'log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_timed.addHandler(console)
    log_timed.addHandler(logfile)
    return log_timed


def set_time_limit(num):
    def wrap(func):
        def handle(sig_num, interrupted_stack_frame):
            raise RuntimeError("\n\nIncrease '--time-limit' to meet the specific need of data!")

        def func_modified(*args, **kwargs):
            signal.signal(signal.SIGALRM, handle)
            signal.alarm(num)
            r = func(*args, **kwargs)
            signal.alarm(0)
            return r

        return func_modified

    return wrap


def pool_multiprocessing(target, iter_args, constant_args, num_process):
    # parse args
    if iter_args:
        if type(iter_args) in {list, tuple}:
            if type(iter_args[0]) not in {list, tuple}:
                iter_args = [[each_arg] for each_arg in iter_args]
        else:
            sys.stderr.write("iter_args must be list/tuple!\n")
    else:
        return
    if constant_args:
        if type(constant_args) not in {list, tuple}:
            constant_args = [constant_args]
    else:
        constant_args = []
    pool = Pool(processes=num_process)
    for this_arg in iter_args:
        pool.apply_async(target, tuple(list(this_arg) + list(constant_args)))
    pool.close()
    try:
        pool.join()
    except KeyboardInterrupt:
        pool.terminate()
        raise KeyboardInterrupt


# test whether an external binary is executable
def executable(test_this):
    return True if os.access(test_this, os.X_OK) or getstatusoutput(test_this)[0] != dead_code else False


def mapping_with_bowtie2_for_library(graph_fasta, fq_1, fq_2, out_base, threads, resume, log, verbose_log):

    log.info("Making seed - bowtie2 index ...")
    build_seed_index = subprocess.Popen("bowtie2-build --large-index " + graph_fasta + " " + graph_fasta + '.index',
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log:
        log.info("bowtie2-build --large-index " + graph_fasta + " " + graph_fasta + '.index')
    output, err = build_seed_index.communicate()
    if "unrecognized option" in output.decode("utf8"):
        build_seed_index = subprocess.Popen("bowtie2-build " + graph_fasta + " " + graph_fasta + '.index',
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = build_seed_index.communicate()
    if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
        log.error('\n' + output.decode("utf8"))
        exit()
    if verbose_log:
        log.info("\n" + output.decode("utf8").strip())
    log.info("Making seed - bowtie2 index finished.")

    mapping_sam = [os.path.join(os.path.split(out_base)[0], x + os.path.split(out_base)[1] + "scaffold_bowtie.sam")
                   for x in ("temp.", "")]
    if resume and os.path.exists(mapping_sam[1]):
        log.info("Mapping result existed!")
        return mapping_sam[1]
    else:
        log.info("Mapping reads to graph - bowtie2 index ...")
        "-S 800-a.sam -1 800bp_23_1.fq -2 800bp_23_2.fq "

        this_command = "bowtie2 -p " + str(threads) + " --local -a --no-discordant --no-contain -X 1000000 " + \
                       " -x " + graph_fasta + ".index -1 " + fq_1 + " -2 " + fq_2 + " -S " + mapping_sam[0] + \
                       " --no-unal -t"
        make_seed_bowtie2 = subprocess.Popen(this_command,
                                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log:
            log.info(this_command)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            log.error('\n' + output.decode("utf8"))
            return False
        if verbose_log:
            log.info("\n" + output.decode("utf8").strip())

        if os.path.exists(mapping_sam[0]):
            os.rename(mapping_sam[0], mapping_sam[1])
            log.info("Mapping finished.")
            return mapping_sam[1]
        else:
            log.error("Cannot find bowtie2 result!")
            return False


MEM_TRANS = {"K": 1000, "M": 1000000, "G":1000000000}


class LogInfo:
    def __init__(self, sample_out_dir):
        self.header = ["sample", "fastq_format", "mean_error_rate", "trim_percent", "trim_int", "trim_chars",
                       "mean_read_len", "max_read_len", "estimated_cp_base_coverage",
                       "w", "k", "pre_w", "mem_dup", "mem_used", "n_candidates", "n_reads", "n_bases", "dup_used",
                       "mem_group", "rounds", "accepted_lines", "accepted_words", "mem_extending", "circular",
                       "degenerate_base_used", "library_size", "library_deviation", "library_left", "library_right",
                       "res_kmer", "res_kmer_cov", "res_base_cov", "res_path_count",
                       "res_len", "time", "mem_spades", "mem_max"]
        this_record = {"sample": sample_out_dir, "circular": "no", "mem_max": 0, "degenerate_base_used": "no",
                       "time": 0.}
        if os.path.exists(os.path.join(sample_out_dir, "get_org.log.txt")):
            log_contents = open(os.path.join(sample_out_dir, "get_org.log.txt"), "r").read().split("\n\n")
            for log_part in log_contents:
                if "get_organelle_reads" in log_part:
                    command_line = log_part.split()
                    if "-w" in command_line:
                        this_record["w"] = command_line[command_line.index("-w") + 1]
                    if "-k" in command_line:
                        this_record["k"] = command_line[command_line.index("-k") + 1]
                if "Pre-reading fastq" in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Identified quality encoding format = "):
                                this_record["fastq_format"] = detail_record.split(" = ")[-1]
                            elif detail_record.startswith("Trimming bases with qualities"):
                                context, trimmed = detail_record.split(":")
                                this_record["trim_percent"] = context.split("(")[-1].split(")")[0]
                                this_record["trim_int"], this_record["trim_chars"] = trimmed.strip().split("  ")
                            elif detail_record.startswith("Mean error rate = "):
                                this_record["mean_error_rate"] = detail_record.split(" = ")[-1]
                            elif detail_record.startswith("Mean = "):
                                this_record["mean_read_len"] = detail_record.split(" = ")[1].split(" bp")[0]
                                this_record["max_read_len"] = detail_record.split(" = ")[-1].split(" bp")[0]
                elif "Checking seed reads and parameters" in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Estimated cp base-coverage = "):
                                this_record["estimated_cp_base_coverage"] = detail_record.split("coverage = ")[-1]
                            elif detail_record.startswith("Setting '-w "):
                                this_record["w"] = detail_record.split("-w ")[-1][:-1]
                            elif detail_record.startswith("Setting '-k "):
                                this_record["k"] = detail_record.split("-k ")[-1][:-1]
                elif "Making read index" in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if "candidates in all" in detail_record:
                                if "Mem" in detail_record:
                                    mem_part, data_part = detail_record.split(", ")
                                    mem = mem_part.split("Mem ")[-1]
                                    this_record["mem_dup"] = int(float(mem.split()[0]) * MEM_TRANS[mem.split()[1]])
                                    this_record["mem_max"] = max(this_record["mem_max"], this_record["mem_dup"])
                                    data_part_list = data_part.split(" ")
                                else:
                                    data_part_list = detail_record.split(" ")
                                this_record["n_candidates"] = data_part_list[0]
                                this_record["n_reads"] = data_part_list[4]
                                if "mean_read_len" in this_record:
                                    this_record["n_bases"] = str(int(float(this_record["mean_read_len"]) *
                                                                     int(this_record["n_reads"])))
                            elif detail_record.startswith("Setting '--pre-w "):
                                this_record["pre_w"] = detail_record.split("--pre-w ")[-1][:-1]
                            elif detail_record.endswith("used/duplicated"):
                                if "Mem" in detail_record:
                                    mem = detail_record.split("Mem ")[-1].split(",")[0]
                                    this_record["mem_used"] = int(float(mem.split()[0]) * MEM_TRANS[mem.split()[1]])
                                    this_record["mem_max"] = max(this_record["mem_max"], this_record["mem_used"])
                                this_record["dup_used"] = detail_record.split(", ")[-1].split(" ")[0]
                            elif "groups made" in detail_record and "Mem" in detail_record:
                                mem = detail_record.split("Mem ")[-1].split(",")[0]
                                this_record["mem_group"] = int(float(mem.split()[0]) * MEM_TRANS[mem.split()[1]])
                                this_record["mem_max"] = max(this_record["mem_max"], this_record["mem_group"])
                elif "Extending finished" in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Setting '-w "):
                                this_record["w"] = detail_record.split("-w ")[-1][:-1]
                            elif detail_record.startswith("Round "):
                                if "Mem" in detail_record:
                                    flag, r_num, l_count, ai_flag, acc_i, aw_flag, acc_w, mem_flag, mem = \
                                        detail_record.split(" ")
                                    this_record["mem_extending"] = int(float(mem) * MEM_TRANS["G"])
                                    this_record["mem_max"] = max(this_record["mem_max"], this_record["mem_extending"])
                                else:
                                    flag, r_num, l_count, ai_flag, acc_i, aw_flag, acc_w = detail_record.split(" ")
                                this_record["rounds"] = r_num[:-1]
                                this_record["accepted_lines"] = acc_i
                                this_record["accepted_words"] = acc_w
                elif "Assembling using SPAdes ..." in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Setting '-k "):
                                this_record["k"] = detail_record.split("-k ")[-1][:-1]
                            elif detail_record.startswith("Insert size = "):
                                lib_parts = detail_record.split(", ")
                                # try:
                                this_record["library_size"] = float(lib_parts[0].split(" = ")[-1])
                                this_record["library_deviation"] = float(lib_parts[1].split(" = ")[-1])
                                this_record["library_left"] = float(lib_parts[2].split(" = ")[-1])
                                this_record["library_right"] = float(lib_parts[3].split(" = ")[-1])
                                # except IndexError:
                                #     pass
                elif "Slimming and disentangling graph ..." in log_part:
                    these_lines = log_part.split("\n")
                    for go_l, line in enumerate(these_lines):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Average ") and "kmer-coverage" in detail_record:
                                this_record["res_kmer_cov"] = detail_record.split(" = ")[-1]
                            elif detail_record.startswith("Average ") and "base-coverage" in detail_record:
                                this_record["res_base_cov"] = detail_record.split(" = ")[-1]
                elif "Writing output ..." in log_part:
                    these_lines = log_part.split("\n")
                    res_lengths = set()
                    for go_l, line in enumerate(these_lines):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Writing GRAPH to "):
                                if " - INFO: " in these_lines[go_l - 1] and these_lines[go_l - 1][:4].isdigit():
                                    time_point, detail_record = these_lines[go_l - 1].split(" - INFO:")
                                    detail_record = detail_record.strip()
                                    for out_name in os.path.split(detail_record.split()[-1])[-1].split(".")[::-1]:
                                        if out_name.startswith("K") and out_name[1:].isdigit():
                                            this_record["res_kmer"] = int(out_name[1:])
                                    this_record["res_path_count"] = detail_record.split("Writing PATH")[-1].split(" ")[0]
                                    this_p_file = detail_record.split("Writing PATH")[-1].split(" to ")[-1]
                                    this_p_file = os.path.join(os.path.split(sample_out_dir)[0], this_p_file)
                                    if os.path.exists(this_p_file):
                                        res_lengths.add(tuple(sorted([len(this_seq.seq)
                                                                      for this_seq in SequenceList(this_p_file)])))
                            elif detail_record.startswith("Result status"):
                                if "circular genome" in detail_record:
                                    this_record["circular"] = "yes"
                        elif " - WARNING: " in line and line[:4].isdigit():
                            if "Degenerate base(s) used!" in line:
                                this_record["degenerate_base_used"] = "yes"
                    this_record["res_len"] = ";".join([",".join([str(l) for l in list(l_l)]) for l_l in sorted(res_lengths)])
                elif "Total cost " in log_part:
                    for line in log_part.split("\n"):
                        if "Total cost " in line:
                            this_record["time"] += float(line.split("Total cost ")[-1].split(" ")[0])
        else:
            raise FileNotFoundError(os.path.join(sample_out_dir, "get_org.log.txt"))
        if os.path.exists(os.path.join(sample_out_dir, "filtered_spades", "spades.log")):
            with open(os.path.join(sample_out_dir, "filtered_spades", "spades.log"), "r") as spades_log:
                for line in spades_log:
                    line = line.strip()
                    if line.count(":") > 2:
                        line = line.split()
                        if line[0].replace(":", "").replace(".", "").isdigit():
                            this_record["mem_spades"] = max(this_record.get("mem_spades", 0),
                                                            int(line[1][:-1]) * MEM_TRANS[line[1][-1]])
                if "mem_spades" in this_record:
                    this_record["mem_max"] = max(this_record["mem_max"], this_record["mem_spades"])
        self.__dict__.update(this_record)

    def get_all_values_as_tab_line(self):
        return "\t".join([self.__dict__.get(this_key, "") for this_key in self.header])


