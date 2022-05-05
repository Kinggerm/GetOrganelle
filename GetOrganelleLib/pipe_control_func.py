import signal
import time
import logging
import sys
import os
import hashlib
from multiprocessing import Pool


MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    python_version = "2.7+"
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)

if python_version == "2.7+":
    from commands import getstatusoutput
    class ConnectionRefusedError(Exception):
        def __init__(self, value=""):
            self.value = value

        def __str__(self):
            return repr(self.value)

    class FileNotFoundError:
        def __init__(self, value=""):
            self.value = value

        def __str__(self):
            return repr(self.value)
else:
    from subprocess import getstatusoutput
import subprocess
# aborted
# directory
# command not found
DEAD_CODES = {"2.7+": (512, 32256, 32512),
              "3.5+": (2, 126, 127)}[python_version]

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
from GetOrganelleLib.seq_parser import phred_offset_table

ORGANELLE_TYPE_SET = {"embplant_pt", "embplant_mt", "embplant_nr", "fungus_mt", "fungus_nr", "animal_mt", "other_pt"}
ORGANELLE_TYPE_LIST = ["embplant_pt", "embplant_mt", "embplant_nr", "fungus_mt", "fungus_nr", "animal_mt", "other_pt"]
MAX_SLIM_EXTENDING_LENS = {"embplant_pt": 15000,
                           "embplant_mt": 50000,
                           "embplant_nr": 12500,
                           "fungus_mt": 50000,
                           "fungus_nr": 12500,
                           "animal_mt": 12500,
                           "other_pt": 50000,
                           "anonym": float("inf")}


if "GETORG_PATH" in os.environ:
    GO_PATH = os.path.expanduser(os.environ["GETORG_PATH"])
else:
    GO_PATH = os.path.expanduser("~/.GetOrganelle")
LBL_NAME = "LabelDatabase"
SEQ_NAME = "SeedDatabase"
LBL_DB_PATH = os.path.join(GO_PATH, LBL_NAME)
SEQ_DB_PATH = os.path.join(GO_PATH, SEQ_NAME)

HEAD_MAXIMUM_LINES = 2147483647

INVALID_PATH_CHAR_RANGES = [
    [u'\u4e00', u'\u9fff'],  # chinese characters are not accepted by SPAdes
    [" ", " "]
]


def simple_log(log, output_base, prefix, log_level="NOTSET"):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    if log_level == "INFO":
        log_simple.setLevel(logging.INFO)
    else:
        log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    if log_level == "INFO":
        console.setLevel(logging.INFO)
    else:
        console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, prefix + 'log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(message)s'))
    if log_level == "INFO":
        logfile.setLevel(logging.INFO)
    else:
        logfile.setLevel(logging.NOTSET)
    log_simple.addHandler(console)
    log_simple.addHandler(logfile)
    return log_simple


def timed_log(log, output_base, prefix, log_level="NOTSET"):
    log_timed = log
    for handler in list(log_timed.handlers):
        log_timed.removeHandler(handler)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    if log_level == "INFO":
        console.setLevel(logging.INFO)
    else:
        console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, prefix + 'log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    if log_level == "INFO":
        logfile.setLevel(logging.INFO)
    else:
        logfile.setLevel(logging.NOTSET)
    log_timed.addHandler(console)
    log_timed.addHandler(logfile)
    return log_timed


if MAJOR_VERSION == 2:
    class TimeoutError(Exception):
        def __init__(self, value=""):
            self.value = value

        def __str__(self):
            return repr(self.value)


def set_time_limit(num, flag_str="'--time-limit'"):
    def wrap(func):
        def handle(sig_num, interrupted_stack_frame):
            raise TimeoutError("\n\nIncrease " + flag_str + " to meet the specific need of data!")

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


# test whether an full-path-to binary (os.access & os.path.isfile) is executable:
#      or (os.access(test_this, os.X_OK) and os.path.isfile(test_this.split(" ")[0]))
# test whether a variable with argv (getstatusoutput) is executable
def executable(test_this):
    return True if getstatusoutput(test_this)[0] not in DEAD_CODES else False


def is_valid_path(path_str):
    for char in path_str:
        for down_str, up_str in INVALID_PATH_CHAR_RANGES:
            if down_str <= char <= up_str:
                return False
    return True


def draw_assembly_graph_using_bandage(input_graph_file, output_image_file, assembly_graph_ob,
                                      resume=False, log_handler=None, verbose_log=False, which_bandage=""):
    if resume and os.path.exists(output_image_file):
        return True
    elif executable(os.path.join(which_bandage, "Bandage -v")):
        assert output_image_file.split(".")[-1] in {"png", "jpg", "svg"}, "Output image format must be png/jpg/svg!"
        temp_file = os.path.join(os.path.split(output_image_file)[0], "temp." + os.path.split(output_image_file)[1])
        # preparing for color
        if not assembly_graph_ob.vertex_to_copy:
            assembly_graph_ob.estimate_copy_and_depth_by_cov(mode="all", verbose=verbose_log, log_handler=log_handler)
        all_coverages = [assembly_graph_ob.vertex_info[g_v].cov for g_v in assembly_graph_ob.vertex_info]
        low_cov_threshold = min(all_coverages)
        if set(assembly_graph_ob.copy_to_vertex) == {1}:
            high_cov_threshold = 2 * low_cov_threshold
        else:
            high_cov_threshold = max(all_coverages)
        # preparing for commands
        this_command = os.path.join(which_bandage, "Bandage") + " image " + \
                       input_graph_file + " " + temp_file + " " \
                       "--names --depth --lengths --outline 0 --colour depth " \
                       "--depcollow \"#D8D8D8\" --depcolhi \"#C60E29\" " \
                       "--depvallow " + str(low_cov_threshold) + " --depvalhi " + str(high_cov_threshold) + " " \
                       "--fontsize 6 --iter 4 && mv " + temp_file + " " + output_image_file
        if verbose_log:
            if log_handler:
                log_handler.info(this_command)
            else:
                sys.stdout.write(this_command + "\n")
        draw_image = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = draw_image.communicate()
        if "error" in output.decode("utf8").lower() or "Abort trap" in output.decode("utf8"):
            if verbose_log:
                if log_handler:
                    log_handler.error("\n" + output.decode("utf8").strip())
                else:
                    sys.stdout.write("\n" + output.decode("utf8").strip() + "\n")
            else:
                if log_handler:
                    log_handler.warning("Visualizing the assembly graph failed with Bandage!")
                else:
                    sys.stdout.write("Warning: Visualizing the assembly graph failed with Bandage!\n")
        if os.path.exists(output_image_file):
            if log_handler:
                log_handler.info("Writing GRAPH image to " + output_image_file)
            else:
                sys.stdout.write("Writing GRAPH image to " + output_image_file + "\n")
            return True
        else:
            if log_handler:
                log_handler.warning("Visualizing the assembly graph failed with Bandage!")
            else:
                sys.stdout.write("Warning: Visualizing the assembly graph failed with Bandage!\n")
    else:
        return False


def make_blast_db(input_file, output_base, which_blast="", log_handler=None, verbose_log=False):
    this_command = os.path.join(which_blast, "makeblastdb") + \
                   " -in " + input_file + " -out " + output_base + " -dbtype nucl"
    if verbose_log:
        if log_handler:
            log_handler.info(this_command)
        else:
            sys.stdout.write(this_command + "\n")
    mk_blast_db_run = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = mk_blast_db_run.communicate()
    if "error" in output.decode("utf8").lower() or "(ERR)" in output.decode("utf8"):
        if log_handler:
            log_handler.error(output.decode("utf8"))
        else:
            sys.stdout.write('\n' + output.decode("utf8"))
        exit()


def execute_blast(query, blast_db, output, threads, outfmt=6, e_value="1e-25", word_size=11, other_options="",
                  which_blast="", silent=False, log_handler=None):
    if not silent:
        if log_handler:
            log_handler.info("Executing BLAST ...")
        else:
            sys.stdout.write("Executing BLAST ...\n")
    make_blast = subprocess.Popen(
        os.path.join(which_blast, "blastn") + " -evalue " + str(e_value) + " -num_threads " + str(threads) +
        " -word_size " + str(word_size) + other_options +
        " -query " + query + " -db " + blast_db + " -out " + output + " -outfmt " + str(outfmt),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = make_blast.communicate()
    if "(ERR)" in output.decode("utf-8") or "Error:" in output.decode("utf-8") or "error" in output.decode("utf-8"):
        if log_handler:
            log_handler.error(output.decode("utf-8"))
        else:
            sys.stdout.write("\nERROR:\n" + output.decode("utf-8"))
        exit()
    if not silent:
        if log_handler:
            log_handler.info("Executing BLAST finished.")
        else:
            sys.stdout.write("Executing BLAST finished.\n")


def build_bowtie2_db(seed_file, seed_index_base, which_bowtie2, target_echo_name="", overwrite=False, random_seed=12345,
                     silent=False, log_handler=None, verbose_log=False):
    if not overwrite and (sum([os.path.exists(seed_index_base + postfix)
                               for postfix in
                               (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")]) == 6 or
                          sum([os.path.exists(remove_db_postfix(seed_file) + ".index" + postfix)
                               for postfix in
                               (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")]) == 6):
        if not silent:
            if log_handler:
                log_handler.info((target_echo_name + " ").capitalize() * int(bool(target_echo_name)) +
                                 "bowtie2 index existed!")
            else:
                sys.stdout.write((target_echo_name + " ").capitalize() * int(bool(target_echo_name)) +
                                 "bowtie2 index existed!\n")
    else:
        if not silent:
            if log_handler:
                log_handler.info("Making " + target_echo_name + " - bowtie2 index ...")
            else:
                sys.stdout.write("Making " + target_echo_name + " - bowtie2 index ...\n")
        this_command = os.path.join(which_bowtie2, "bowtie2-build") + " --seed " + str(random_seed) + \
                       " --large-index " + seed_file + " " + seed_index_base
        building = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if not silent and verbose_log:
            if log_handler:
                log_handler.info(this_command)
            else:
                sys.stdout.write(this_command + "\n")
        output, err = building.communicate()
        # if "unrecognized option" in output.decode("utf8"):
        #     this_command = os.path.join(which_bowtie2, "bowtie2-build") + " --seed " + str(random_seed) + \
        #                    " " + seed_file + " " + seed_index_base
        #     building = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        #     if not silent and verbose_log:
        #         if log_handler:
        #             log_handler.info(this_command)
        #         else:
        #             sys.stdout.write(this_command + "\n")
        #     output, err = building.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            if log_handler:
                log_handler.error('\n' + output.decode("utf8"))
            else:
                sys.stdout.write("\nError: \n" + output.decode("utf8"))
            exit()
        if not silent and verbose_log:
            if log_handler:
                log_handler.info("\n" + output.decode("utf8").strip())
            else:
                sys.stdout.write("\n" + output.decode("utf8").strip() + "\n")
        if not silent:
            if log_handler:
                log_handler.info("Making " + target_echo_name + " - bowtie2 index finished.")
            else:
                sys.stdout.write("Making " + target_echo_name + " - bowtie2 index finished.\n")


def map_with_bowtie2(seed_file, original_fq_files, bowtie_out, resume, threads, random_seed, log_handler,
                     target_echo_name="", generate_fq=True, silent=False, verbose_log=False, which_bowtie2="",
                     bowtie2_mode="", bowtie2_other_options="", fq_1=None, fq_2=None):
    res_path_name, res_base_name = os.path.split(bowtie_out)
    total_seed_fq = [os.path.join(res_path_name, x + res_base_name + ".fq") for x in ("temp.", "")]
    total_seed_sam = [os.path.join(res_path_name, x + res_base_name + ".sam") for x in ("temp.", "")]
    if resume and (os.path.exists(total_seed_fq[1]) or not generate_fq) and os.path.exists(total_seed_sam[1]):
        if not silent:
            if log_handler:
                log_handler.info((target_echo_name + " ") * int(bool(target_echo_name)) + "reads existed!")
            else:
                sys.stdout.write((target_echo_name + " ") * int(bool(target_echo_name)) + "reads existed!\n")
    else:
        seed_index_base = remove_db_postfix(seed_file) + '.index'
        build_bowtie2_db(seed_file=seed_file, seed_index_base=seed_index_base, which_bowtie2=which_bowtie2,
                         target_echo_name=target_echo_name, overwrite=False, random_seed=random_seed,
                         log_handler=log_handler, silent=silent)
        if not silent:
            if log_handler:
                log_handler.info("Mapping reads to " + target_echo_name + " bowtie2 index ...")
            else:
                sys.stdout.write("Mapping reads to " + target_echo_name + " bowtie2 index ...\n")
        this_input_fq_data = " "
        if fq_1 and fq_2:
            this_input_fq_data += " -1 " + fq_1 + " -2 " + fq_2
        if original_fq_files:
            this_input_fq_data += " -U " + ",".join(original_fq_files)
        if not this_input_fq_data.strip():
            raise ValueError("Parameters: original_fq_files or (fq_1 and fq_2) is/are needed!")
        this_command = os.path.join(which_bowtie2, "bowtie2") + " --large-index --seed " + str(random_seed) + \
                       " --mm -p " + str(threads) + " " + bowtie2_mode + " " + bowtie2_other_options + \
                       (" --al " + total_seed_fq[0]) * int(bool(generate_fq)) + " -x " + seed_index_base + \
                       this_input_fq_data + " -S " + total_seed_sam[0] + \
                       " --no-unal --omit-sec-seq"
        if not silent and verbose_log:
            if log_handler:
                log_handler.info(this_command)
            else:
                sys.stdout.write(this_command + "\n")
        make_seed_bowtie2 = subprocess.Popen(this_command,
                                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            # if not silent:
            if log_handler:
                log_handler.error("\n" + output.decode("utf8"))
            else:
                sys.stdout.write("\nError: " + output.decode("utf8") + "\n")
            if "has more quality values than read characters" in output.decode("utf8"):
                raise Exception("Please check the integrity of your input reads file (e.g. using md5)!!")
            else:
                raise Exception("")
        elif "No such file" in output.decode("utf8") or "not found" in output.decode("utf8"):
            # if not silent:
            if log_handler:
                log_handler.error("\n" + output.decode("utf8"))
            else:
                sys.stdout.write("\nError: " + output.decode("utf8") + "\n")
            if "perl" in output.decode("utf8"):
                raise Exception("perl is required for the wrapper of bowtie2!")
        if not silent and verbose_log:
            if log_handler:
                log_handler.info("\n" + output.decode("utf8").strip())
            else:
                sys.stdout.write("\n" + output.decode("utf8").strip() + "\n")
        if os.path.exists(total_seed_sam[0]):
            os.rename(total_seed_sam[0], total_seed_sam[1])
            if generate_fq:
                os.rename(total_seed_fq[0], total_seed_fq[1])
            if not silent:
                if log_handler:
                    log_handler.info("Mapping finished.")
                else:
                    sys.stdout.write("Mapping finished.\n")
        else:
            raise Exception("Cannot find bowtie2 result!")


def mapping_with_bowtie2_for_library(graph_fasta, fq_1, fq_2, out_base, threads, resume, rand_seed,
                                     log_handler, verbose_log, which_bowtie2=""):

    log_handler.info("Making seed - bowtie2 index ...")
    this_command = os.path.join(which_bowtie2, "bowtie2-build") + " --seed " + str(rand_seed) + " --large-index " + \
                   graph_fasta + " " + graph_fasta + '.index'
    build_seed_index = subprocess.Popen(this_command,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log:
        log_handler.info(this_command)
    output, err = build_seed_index.communicate()
    if "unrecognized option" in output.decode("utf8"):
        this_command = os.path.join(which_bowtie2, "bowtie2-build") + " --seed " + str(rand_seed) + " " + \
                       graph_fasta + " " + graph_fasta + ".index"
        build_seed_index = subprocess.Popen(this_command,
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log:
            log_handler.info(this_command)
        output, err = build_seed_index.communicate()
    if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
        log_handler.error('\n' + output.decode("utf8"))
        exit()
    if verbose_log:
        log_handler.info("\n" + output.decode("utf8").strip())
    log_handler.info("Making seed - bowtie2 index finished.")

    mapping_sam = [os.path.join(os.path.split(out_base)[0], x + os.path.split(out_base)[1] + "scaffold_bowtie.sam")
                   for x in ("temp.", "")]
    if resume and os.path.exists(mapping_sam[1]):
        log_handler.info("Mapping result existed!")
        return mapping_sam[1]
    else:
        log_handler.info("Mapping reads to graph - bowtie2 index ...")
        "-S 800-a.sam -1 800bp_23_1.fq -2 800bp_23_2.fq "

        this_command = os.path.join(which_bowtie2, "bowtie2") + " --large-index --seed " + str(rand_seed) + \
                       " -p " + str(threads) + " --local -a --no-discordant --no-contain -X 1000000 " + \
                       " -x " + graph_fasta + ".index -1 " + fq_1 + " -2 " + fq_2 + " -S " + mapping_sam[0] + \
                       " --no-unal -t"
        make_seed_bowtie2 = subprocess.Popen(this_command,
                                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log:
            log_handler.info(this_command)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            log_handler.error('\n' + output.decode("utf8"))
            return False
        if verbose_log:
            log_handler.info("\n" + output.decode("utf8").strip())

        if os.path.exists(mapping_sam[0]):
            os.rename(mapping_sam[0], mapping_sam[1])
            log_handler.info("Mapping finished.")
            return mapping_sam[1]
        else:
            log_handler.error("Cannot find bowtie2 result!")
            return False


MEM_TRANS = {"K": 1024, "M": 1048576, "G": 1073741824}


class LogInfo:
    def __init__(self, sample_out_dir, prefix=""):
        from GetOrganelleLib.seq_parser import SequenceList
        self.header = ["sample", "organelle_type", "fastq_format", "mean_error_rate",
                       "trim_percent", "trim_int", "trim_chars",
                       "mean_read_len", "max_read_len", "num_reads_1", "seed_read_filesize",
                       "estimated_base_coverage", "closest_seed", "unmapped_percentage", "unmapped_lengths",
                       "num_reads_2", "w", "k", "max_extending_len", "pre_w",
                       "mem_dup", "mem_used", "n_unique_reads", "n_reads", "n_bases", "dup_used",
                       "mem_group", "rounds", "accepted_lines", "accepted_words", "mem_extending", "circular",
                       "degenerate_base_used", "library_size", "library_deviation", "library_left", "library_right",
                       "res_kmer", "res_kmer_cov", "res_base_cov", "res_path_count",
                       "res_len", "time", "mem_spades", "mem_max"]
        this_record = {"sample": os.path.join(sample_out_dir, prefix).rstrip("/"), "mem_max": 0, "time": 0.}
        command_recorded = False
        if os.path.exists(os.path.join(sample_out_dir, prefix + "get_org.log.txt")):
            log_contents = open(os.path.join(sample_out_dir, prefix + "get_org.log.txt"), "r").read().split("\n\n")
            for log_part in log_contents:
                if "Python" in log_part and "get_organelle_from_reads" in log_part and not command_recorded:
                    command_line = log_part.split()
                    if "-w" in command_line:
                        this_record["w"] = command_line[command_line.index("-w") + 1]
                    if "-k" in command_line:
                        this_record["k"] = command_line[command_line.index("-k") + 1]
                    if "-F" in command_line:
                        this_record["organelle_type"] = command_line[command_line.index("-F") + 1]
                    command_recorded = True
                if "Pre-reading fastq" in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Identified quality encoding format = "):
                                this_record["fastq_format"] = detail_record.split(" = ")[-1]
                                if "record_fq_beyond_limit" not in this_record:
                                    this_record["record_fq_beyond_limit"] = False
                                this_record["phred_offset"] = phred_offset_table[this_record["fastq_format"]]
                            if detail_record.startswith("Phred offset = "):
                                this_record["phred_offset"] = int(detail_record.split(" = ")[-1].strip())
                            elif detail_record.startswith("Trimming bases with qualities"):
                                context, trimmed = detail_record.split(":")
                                this_record["trim_percent"] = context.split("(")[-1].split(")")[0]
                                this_record["trim_int"], this_record["trim_chars"] = trimmed.strip().split("  ")
                            elif detail_record.startswith("Mean error rate = "):
                                this_record["mean_error_rate"] = detail_record.split(" = ")[-1]
                            elif detail_record.startswith("Mean = "):
                                this_record["mean_read_len"] = detail_record.split(" = ")[1].split(" bp")[0]
                                this_record["max_read_len"] = detail_record.split(" = ")[-1].split(" bp")[0]
                            elif detail_record.startswith("Reads used = "):
                                this_record["num_reads_1"] = detail_record.split(" = ")[-1]
                            # if " - WARNING: " in line and line[:4].isdigit():
                            #     time_point, detail_record = line.strip().split(" - WARNING:")
                            #     detail_record = detail_record.strip()
                            elif detail_record.endswith("are used in downstream analysis."):
                                this_record["record_fq_beyond_limit"] = True
                elif "Making seed reads .." in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Seed reads made: "):
                                seed_read_filesize = detail_record.split("(")[-1].split(" ")[0]
                                if "seed_read_filesize" in this_record:
                                    this_record["seed_read_filesize"] += " & " + seed_read_filesize
                                else:
                                    this_record["seed_read_filesize"] = seed_read_filesize
                elif "Checking seed reads and parameters" in log_part:
                    for line in log_part.split("\n"):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Estimated ") and "base-coverage = " in detail_record:
                                this_cov = detail_record.split("coverage = ")[-1]
                                if "estimated_base_coverage" in this_record:
                                    this_record["estimated_base_coverage"] += " & " + this_cov
                                else:
                                    this_record["estimated_base_coverage"] = this_cov
                            elif detail_record.startswith("Closest") and "seed sequence: " in detail_record:
                                closest_seed = detail_record.split(": ")[-1]
                                if "closest_seed" in this_record:
                                    this_record["closest_seed"] += " & " + closest_seed
                                else:
                                    this_record["closest_seed"] = closest_seed
                            elif detail_record.startswith("Unmapped percentage"):
                                unmapped_percentage = detail_record.split("Unmapped percentage ")[-1].split()[0]
                                if "unmapped_percentage" in this_record:
                                    this_record["unmapped_percentage"] += " & " + unmapped_percentage
                                else:
                                    this_record["unmapped_percentage"] = unmapped_percentage
                                unmapped_lens = detail_record.split("unmapped lengths ")[-1].split(" ..")[0]
                                if "unmapped_lens" in this_record:
                                    this_record["unmapped_lengths"] += " & " + unmapped_lens
                                else:
                                    this_record["unmapped_lengths"] = unmapped_lens
                            elif detail_record.startswith("Reads reduced to = "):
                                this_record["num_reads_2"] = detail_record.split(" = ")[-1]
                            elif detail_record.startswith("Setting '-w "):
                                this_record["w"] = detail_record.split("-w ")[-1][:-1]
                            elif detail_record.startswith("Setting '-k "):
                                this_record["k"] = detail_record.split("-k ")[-1][:-1]
                            elif detail_record.startswith("Setting '--max-extending-len"):
                                this_record["max_extending_len"] = detail_record.split("--max-extending-len ")[-1][:-1]
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
                                this_record["n_unique_reads"] = data_part_list[0]
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
                elif "Disentangling " in log_part:  # Writing output ...
                    res_lengths = set()
                    this_degenerate = "no"
                    this_circular = "no"
                    these_lines = log_part.split("\n")
                    for go_l, line in enumerate(these_lines):
                        if " - INFO: " in line and line[:4].isdigit():
                            time_point, detail_record = line.strip().split(" - INFO:")
                            detail_record = detail_record.strip()
                            if detail_record.startswith("Average ") and "kmer-coverage" in detail_record:
                                res_kmer_cov = detail_record.split(" = ")[-1]
                                if "res_kmer_cov" in this_record:
                                    this_record["res_kmer_cov"] += " & " + res_kmer_cov
                                else:
                                    this_record["res_kmer_cov"] = res_kmer_cov
                            elif detail_record.startswith("Average ") and "base-coverage" in detail_record:
                                res_base_cov = detail_record.split(" = ")[-1]
                                if "res_base_cov" in this_record:
                                    this_record["res_base_cov"] += " & " + res_base_cov
                                else:
                                    this_record["res_base_cov"] = res_base_cov
                            elif detail_record.startswith("Writing GRAPH to "):
                                if " - INFO: " in these_lines[go_l - 1] and these_lines[go_l - 1][:4].isdigit():
                                    time_point, detail_record = these_lines[go_l - 1].split(" - INFO:")
                                    detail_record = detail_record.strip()
                                    for out_name in os.path.split(detail_record.split()[-1])[-1].split(".")[::-1]:
                                        if out_name.startswith("K") and out_name[1:].isdigit():
                                            if "res_kmer" in this_record:
                                                this_record["res_kmer"] += " & " + str(int(out_name[1:]))
                                            else:
                                                this_record["res_kmer"] = str(int(out_name[1:]))
                                    res_path_count = detail_record.split("Writing PATH")[-1].split(" ")[0]
                                    if "res_path_count" in this_record:
                                        this_record["res_path_count"] += " & " + res_path_count
                                    else:
                                        this_record["res_path_count"] = res_path_count
                                    this_p_file = detail_record.split("Writing PATH")[-1].split(" to ")[-1]
                                    this_p_file = os.path.join(os.path.split(sample_out_dir)[0], this_p_file)
                                    if os.path.exists(this_p_file):
                                        res_lengths.add(tuple(sorted([len(this_seq.seq)
                                                                      for this_seq in SequenceList(this_p_file)])))
                            elif detail_record.startswith("Result status"):
                                if "circular genome" in detail_record:
                                    this_circular = "yes"
                        elif " - WARNING: " in line and line[:4].isdigit():
                            if "Degenerate base(s) used!" in line:
                                this_degenerate = "yes"
                    if "circular" in this_record:
                        this_record["circular"] += " & " + this_circular
                    else:
                        this_record["circular"] = this_circular
                    if "degenerate_base_used" in this_record:
                        this_record["degenerate_base_used"] += " & " + this_degenerate
                    else:
                        this_record["degenerate_base_used"] = this_degenerate
                    res_len_str = ";".join(["+".join([str(l) for l in list(l_l)]) for l_l in sorted(res_lengths)])
                    if "res_len" in this_record:
                        this_record["res_len"] += " & " + res_len_str
                    else:
                        this_record["res_len"] = res_len_str
                elif "Total cost " in log_part:
                    for line in log_part.split("\n"):
                        if "Total cost " in line:
                            this_record["time"] += float(line.split("Total cost ")[-1].split(" ")[0])
        else:
            sys.stdout.write(os.path.join(sample_out_dir, prefix + "get_org.log.txt") + " not found!\n")
        if os.path.exists(os.path.join(sample_out_dir, prefix + "extended_spades", "spades.log")):
            with open(os.path.join(sample_out_dir, prefix + "extended_spades", "spades.log"), "r") as spades_log:
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


def remove_db_postfix(blast_db_basename):
    if blast_db_basename.endswith(".fasta"):
        blast_db_basename = blast_db_basename[:-6]
    elif blast_db_basename.endswith(".fas"):
        blast_db_basename = blast_db_basename[:-4]
    elif blast_db_basename.endswith(".fa"):
        blast_db_basename = blast_db_basename[:-3]
    return blast_db_basename


def generate_in_ex_info_name(include_indices, exclude_indices, exclude_no_con=True, exclude_no_hit=False):
    is_including = int(bool(include_indices))
    is_excluding = int(bool(exclude_indices))
    joint_char = "." * int(is_including + is_excluding == 2)
    in_ex_info = 'only' * int(exclude_no_hit) + \
                 'extend' * int(exclude_no_con) * int(bool(include_indices)) + \
                 ('-' + "-".join([remove_db_postfix(os.path.split(sub_i)[-1])
                                  for sub_i in include_indices])) * is_including + \
                 joint_char + \
                 ('del-' + "-".join([remove_db_postfix(os.path.split(sub_e)[-1])
                                     for sub_e in exclude_indices])) * is_excluding
    return ('.' + in_ex_info) * int(bool(in_ex_info))


def zip_file(source, target, verbose_log=False, log_handler=None, remove_source=True):
    success = False
    target_temp = target + ".Temp"
    run_command = "tar -C " + os.path.split(source)[0] + " -czvf " + target_temp + " " + os.path.split(source)[1]
    zipping = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log and log_handler:
        log_handler.info(run_command)
    output, err = zipping.communicate()
    if "Unrecognized" not in output.decode("utf8") and \
            "error" not in output.decode("utf8").lower():
        success = True
    if success:
        if verbose_log and log_handler:
            log_handler.info(output.decode("utf8"))
        os.rename(target_temp, target)
        if remove_source:
            os.remove(source)
    else:
        if verbose_log and log_handler:
            log_handler.error("\n" + output.decode("utf8"))
        try:
            os.remove(target_temp)
        except:
            pass


def unzip(source, target, line_limit=HEAD_MAXIMUM_LINES, verbose_log=False, log_handler=None):
    target_temp = target + ".Temp"
    if line_limit == float("inf"):
        try_commands = [
            "tar -x -f " + source + " -O > " + target_temp,
            "gunzip -c " + source + " > " + target_temp]
    else:
        try_commands = [
            "tar -x -f " + source + " -O | head -n " + str(int(line_limit)) + " > " + target_temp,
            "gunzip -c " + source + " | head -n " + str(int(line_limit)) + " > " + target_temp]
    # re-order try commands
    if "tar." not in source:
        try_commands = try_commands[1], try_commands[0]
    if log_handler:
        log_handler.info("Unzipping reads file: " + source + " (" + str(os.path.getsize(source)) + " bytes)")
    success = False
    output = b""
    for run_command in try_commands:
        unzipping = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log and log_handler:
            log_handler.info(run_command)
        output, err = unzipping.communicate()
        if "Unrecognized" not in output.decode("utf8") and \
                "error" not in output.decode("utf8").lower():
            success = True
            break
    if success:
        if verbose_log and log_handler:
            log_handler.info(output.decode("utf8"))
        os.rename(target_temp, target)
    else:
        if verbose_log and log_handler:
            log_handler.error("\n" + output.decode("utf8"))
        try:
            os.remove(target_temp)
        except:
            pass
        raise NotImplementedError("unzipping failed!")


def detect_bowtie2_path(which_bowtie2, dependency_path):
    if not which_bowtie2:
        try_this_bin = os.path.join(dependency_path, "bowtie2", "bowtie2")
        if os.path.isfile(try_this_bin) and executable(try_this_bin):
            which_bowtie2 = os.path.split(try_this_bin)[0]
    return which_bowtie2


def detect_spades_path(which_spades, dependency_path):
    if not which_spades:
        try_this_bin = os.path.join(dependency_path, "SPAdes", "bin", "spades.py")
        if os.path.isfile(try_this_bin) and executable(try_this_bin):
            which_spades = os.path.split(try_this_bin)[0]
    return which_spades


def detect_blast_path(which_blast, dependency_path):
    if not which_blast:
        try_this_bin = os.path.join(dependency_path, "ncbi-blast", "blastn")
        if os.path.isfile(try_this_bin) and executable(try_this_bin):
            output, err = subprocess.Popen(
                try_this_bin + " -version", stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, shell=True).communicate()
            if "not found" in output.decode("utf8"):
                sys.stderr.write(output.decode("utf8") + "\n")
                sys.exit()
            else:
                which_blast = os.path.split(try_this_bin)[0]
    return which_blast


def detect_bowtie2_version(which_bowtie2):
    if executable(os.path.join(which_bowtie2, "bowtie2")):
        output, err = subprocess.Popen(
            os.path.join(which_bowtie2, "bowtie2") + " --version", stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True).communicate()
        output = output.decode("utf8")
        if "(ERR)" in output:
            return "Bowtie2 ERROR"
        this_lines = output.split("\n")[:3]
        return "Bowtie2 " + this_lines[0].split()[-1].strip()
    else:
        return "Bowtie2 N/A"


def detect_spades_version(which_spades):
    if executable(os.path.join(which_spades, "spades.py")):
        output, err = subprocess.Popen(
            os.path.join(which_spades, "spades.py") + " -v", stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True).communicate()
        return output.decode("utf8").replace("v", "").replace("genome assembler ", "").strip()
    else:
        return "SPAdes N/A"


def detect_blast_version(which_blast):
    if executable(os.path.join(which_blast, "blastn")):
        try:
            output, err = subprocess.Popen(
                os.path.join(which_blast, "blastn") + " -version", stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, shell=True).communicate()
            this_lines = output.decode("utf8").split("\n")[:2]
            return "Blast " + this_lines[1].strip().split()[2].replace(",", "").strip()
        except IndexError:
            return "Blast N/A"
    else:
        return "Blast N/A"


def detect_bandage_version(which_bandage):
    if executable(os.path.join(which_bandage, "Bandage -v")):
        output, err = subprocess.Popen(
            os.path.join(which_bandage, "Bandage") + " -v", stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True).communicate()
        return "Bandage " + output.decode("utf8").strip().split()[-1]
    else:
        return "Bandage N/A"


def monitor_spades_log(spades_running, log_handler, sensitively_stop=False, silent=False):
    output = []
    for log_block in spades_running.stdout:
        log_block = log_block.decode("utf8")
        output.append(log_block)
        if "not recognized" in log_block:
            if not silent:
                if log_handler:
                    log_handler.error("Error with running SPAdes: " + log_block.strip("\n"))
                else:
                    sys.stdout.write("Error with running SPAdes: " + log_block.strip("\n") + "\n")
            spades_running.terminate()
            break
        elif "terminated by segmentation fault" in log_block:
            if not silent:
                if log_handler:
                    log_handler.error("Error with running SPAdes: " + log_block.strip("\n"))
                else:
                    sys.stdout.write("Error with running SPAdes: " + log_block.strip("\n") + "\n")
            spades_running.terminate()
            break
        elif "== Error ==" in log_block:
            if not silent:
                if log_handler:
                    log_handler.error("Error with running SPAdes: " + log_block.strip("\n"))
                else:
                    sys.stdout.write("Error with running SPAdes: " + log_block.strip("\n") + "\n")
            if sensitively_stop:
                spades_running.terminate()
                break
    return "".join(output)


def get_static_html_context(remote_url, try_times=5, timeout=10, verbose=False, log_handler=None,
                            alternative_url_list=None):
    import requests
    remote_urls = [remote_url]
    if alternative_url_list:
        remote_urls.extend(alternative_url_list)
    response = False
    count = 0
    while not response and count < try_times:
        try:
            if verbose:
                if log_handler:
                    log_handler.info("Connecting to " + remote_urls[count % len(remote_urls)])
                else:
                    sys.stdout.write("Connecting to " + remote_urls[count % len(remote_urls)] + "\n")
            response = requests.get(remote_urls[count % len(remote_urls)], timeout=timeout)
        except (requests.exceptions.ConnectTimeout,
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                ConnectionRefusedError) as e:
            if count + 1 == try_times:
                return {"status": False, "info": "timeout", "content": ""}
        if response and response.status_code == requests.codes.ok:
            # here_data_str = base64.decodebytes(remote_data_api.json()["content"].encode("utf-8")).decode("utf-8")
            return {"status": True, "info": "", "content": response.content.decode("utf-8")}
        time.sleep(1)
        count += 1
    return {"status": False, "info": "unknown", "content": ""}


def download_file_with_progress(remote_url, output_file, log_handler=None, allow_empty=False,
                                sha256_v=None, try_times=5, timeout=100000, alternative_url_list=None, verbose=False):
    time_0 = time.time()
    import requests
    remote_urls = [remote_url]
    if alternative_url_list:
        remote_urls.extend(alternative_url_list)
    temp_file = output_file + ".Temp"
    count_try = 0
    info_list = []
    sys.stdout.write("Downloading %s \n" % os.path.basename(output_file))
    while count_try < try_times:
        count_try += 1
        with open(temp_file, "wb") as file_h:
            try:
                if verbose:
                    sys.stdout.write("Connecting to " + remote_urls[(count_try - 1) % len(remote_urls)] + "\n")
                response = requests.get(remote_urls[(count_try - 1) % len(remote_urls)], stream=True, timeout=timeout)
                if response.status_code == requests.codes.ok:
                    total_length = response.headers.get("content-length")
                    if total_length is None:
                        file_h.write(response.content)
                    else:
                        total_length = int(total_length)
                        chunk_num = 50
                        chunk_size = int(total_length / chunk_num) + 1
                        go_chunk = 0
                        for data_chunk in response.iter_content(chunk_size=chunk_size):
                            file_h.write(data_chunk)
                            go_chunk += 1
                            sys.stdout.write("\rDownloading %s [%s%s] %i%%" %
                                             (os.path.basename(output_file),
                                              "=" * go_chunk,
                                              " " * (chunk_num - go_chunk),
                                              go_chunk * 2))
                            sys.stdout.flush()
                        sys.stdout.write("\n")
                else:
                    info_list.append("404")
                    continue
            except (requests.exceptions.ConnectionError,
                    requests.exceptions.ConnectTimeout,
                    requests.exceptions.ReadTimeout,
                    ConnectionRefusedError) as e:
                info_list.append("timeout")
                os.remove(temp_file)
                if count_try == try_times:
                    sys.stdout.write("Downloading %s failed, cost %.2f s\n" %
                                     (os.path.basename(output_file), time.time() - time_0))
                    return {"status": False, "info": ",".join(info_list), "content": ""}
                else:
                    time.sleep(1)
                    continue
        # check sha256
        if sha256_v:
            this_sha256 = cal_f_sha256(temp_file)
            if this_sha256 != sha256_v:
                os.remove(temp_file)
                info_list.append("sha256_unmatch")
                if count_try == try_times:
                    sys.stdout.write("Downloading %s failed, cost %.2f s\n" %
                                     (os.path.basename(output_file), time.time() - time_0))
                    return {"status": False, "info": ",".join(info_list), "content": ""}
            else:
                break
        else:
            break
        time.sleep(1)

    # check size
    if os.path.exists(temp_file) and (os.path.getsize(temp_file) or allow_empty):
        os.rename(temp_file, output_file)
        if log_handler:
            log_handler.info("Downloaded %s (%i bytes), cost %.2f s" %
                             (os.path.basename(output_file), os.path.getsize(output_file), time.time() - time_0))
        else:
            sys.stdout.write("Downloaded %s (%i bytes), cost %.2f s\n" %
                             (os.path.basename(output_file), os.path.getsize(output_file), time.time() - time_0))
        return {"status": True, "info": "", "content": ""}
    else:
        if os.path.exists(temp_file):
            os.remove(temp_file)
        if log_handler:
            if count_try < try_times:
                log_handler.error("Downloading %s failed. Try again .." % os.path.basename(output_file))
            else:
                log_handler.error("Downloading %s failed, cost %.2f s" %
                                  (os.path.basename(output_file), time.time() - time_0))
        else:
            if count_try < try_times:
                sys.stdout.write("Downloading %s failed. Try again ..\n" % os.path.basename(output_file))
            else:
                sys.stdout.write("Downloading %s failed, cost %.2f s\n" %
                                 (os.path.basename(output_file), time.time() - time_0))
        return {"status": False, "info": ",".join(info_list), "content": ""}


def cal_f_sha256(file_name):
    hash_class = hashlib.sha256()
    hash_class.update(open(file_name, "rb").read())
    return hash_class.hexdigest()


def get_current_db_versions(db_type, seq_db_path, lbl_db_path, clean_mode=True,
                            check_hash=False, silent=False, log_handler=None):
    """
    :param db_type: seed, label, both
    :param seq_db_path:
    :param lbl_db_path:
    :param check_hash:
    :param silent:
    :param log_handler:
    :return:
    """
    existing_seed_db = {}
    if not silent:
        if check_hash:
            if log_handler:
                log_handler.info("Checking existing database(s) ..")
            else:
                sys.stdout.write("Checking existing database(s) .. \n")
        else:
            if not log_handler:
                sys.stdout.write("Existing databases(s): \n")
    if db_type in ("seed", "both"):
        seed_version_file = os.path.join(seq_db_path, "VERSION")
        if os.path.isfile(seed_version_file):
            with open(seed_version_file) as open_version:
                for line in open_version:
                    og_type, db_version, db_hash = line.strip().split("\t")
                    seed_f = os.path.join(seq_db_path, og_type + ".fasta")
                    if not os.path.exists(seed_f):
                        if not clean_mode:
                            if log_handler:
                                log_handler.error(seed_f + " not available! "
                                                           "Please run `get_organelle_config.py --clean` to reset!")
                            else:
                                sys.stdout.write(seed_f + " not available! "
                                                          "Please run `get_organelle_config.py --clean` to reset!\n")
                            sys.exit()
                        else:
                            existing_seed_db[og_type] = {"version": db_version, "sha256": db_hash}
                    else:
                        existing_seed_db[og_type] = {"version": db_version, "sha256": db_hash}
                        if check_hash:
                            seed_f_hash = cal_f_sha256(seed_f)
                            if db_version not in SEED_DB_HASH or \
                                seed_f_hash != SEED_DB_HASH[db_version][og_type]["sha256"] or \
                                seed_f_hash != db_hash:
                                if log_handler:
                                    log_handler.error(seed_f + " not damaged! "
                                                      "Please run `get_organelle_config.py --clean` to reset!")
                                else:
                                    sys.stdout.write(seed_f + " not damaged! "
                                                     "Please run `get_organelle_config.py --clean` to reset!\n")
                                sys.exit()
        if not clean_mode:
            for og_type in sorted(existing_seed_db):
                db_version = existing_seed_db[og_type]["version"]
                db_hash = existing_seed_db[og_type]["sha256"]
                if not silent:
                    if log_handler:
                        log_handler.info(og_type + " Seed Database:\t" + db_version +
                                         str("\t" + db_hash) * int(bool(db_version == "customized")))
                    else:
                        sys.stdout.write(og_type + " Seed Database:\t" + db_version +
                                         str("\t" + db_hash) * int(bool(db_version == "customized")) + "\n")

    existing_label_db = {}
    if db_type in ("label", "both"):
        label_version_file = os.path.join(lbl_db_path, "VERSION")
        if os.path.isfile(label_version_file):
            with open(label_version_file) as open_version:
                for line in open_version:
                    og_type, db_version, db_hash = line.strip().split("\t")
                    label_f = os.path.join(lbl_db_path, og_type + ".fasta")
                    if not os.path.isfile(label_f):
                        if not clean_mode:
                            if log_handler:
                                log_handler.error(label_f + " not available! "
                                                  "Please run `get_organelle_config.py --clean` to reset!")
                            else:
                                sys.stdout.write(label_f + " not available! "
                                                 "Please run `get_organelle_config.py --clean` to reset!\n")
                            sys.exit()
                        else:
                            existing_label_db[og_type] = {"version": db_version, "sha256": db_hash}
                    else:
                        existing_label_db[og_type] = {"version": db_version, "sha256": db_hash}
                        if check_hash:
                            label_f_hash = cal_f_sha256(label_f)
                            if db_version not in LABEL_DB_HASH or \
                                label_f_hash != LABEL_DB_HASH[db_version][og_type]["sha256"] or \
                                label_f_hash != db_hash:
                                if log_handler:
                                    log_handler.error(label_f + " not damaged! "
                                                      "Please run `get_organelle_config.py --clean` to reset!")
                                else:
                                    sys.stdout.write(label_f + " not damaged! "
                                                     "Please run `get_organelle_config.py --clean` to reset!\n")
                                sys.exit()
        if not clean_mode:
            for og_type in sorted(existing_label_db):
                db_version = existing_label_db[og_type]["version"]
                db_hash = existing_label_db[og_type]["sha256"]
                if not silent:
                    if log_handler:
                        log_handler.info(og_type + " Label Database:\t" + db_version +
                                         str("\t" + db_hash) * int(bool(db_version == "customized")))
                    else:
                        sys.stdout.write(og_type + " Label Database:\t" + db_version +
                                         str("\t" + db_hash) * int(bool(db_version == "customized")) + "\n")
    if not silent and not log_handler:
        sys.stdout.write("\n")
    return existing_seed_db, existing_label_db


def single_line_db_versions(existing_db, check_types, default_value="customized"):
    default_dict = {"version": default_value}
    return "; ".join([check_type + " " + existing_db.get(check_type, default_dict)["version"]
                      for check_type in check_types])



LABEL_DB_HASH = \
    {
    "0.0.0":
        {
        "embplant_nr":
            {"sha256": "603033541683b7c53fb63970c188eb2891844f6419eca342fa648f6ff5e29d71", "size": 16500},
        "embplant_pt":
            {"sha256": "cca977f68f5764bdd9e4b82e711adcfc12d752e78f088cea05c8283dea74fa7e", "size": 370580},
        "animal_mt":
            {"sha256": "af9541621107623ed13f668b4744ddbdcc1e6a8157aa53132717d90e7ac63c4e", "size": 6674918},
        "fungus_mt":
            {"sha256": "e0c6968372ba3fec892532a9231d1e8de789323156b3d23fb5216915a2811b27", "size": 659823},
        "embplant_mt":
            {"sha256": "f3621444441fbf4fa98835999b20e7ac1673a8866cb119707200412c46263570", "size": 148338},
        "other_pt":
            {"sha256": "84346c42ff7e85e5d0859ecdee1de4f06f45fd072851d57ba31c894bcbee3cd8", "size": 412866}
        },
    "0.0.1":
        {
        "embplant_nr":
            {"sha256": "603033541683b7c53fb63970c188eb2891844f6419eca342fa648f6ff5e29d71", "size": 16500},
        "embplant_pt":
            {"sha256": "cca977f68f5764bdd9e4b82e711adcfc12d752e78f088cea05c8283dea74fa7e", "size": 370580},
        "animal_mt":
            {"sha256": "af9541621107623ed13f668b4744ddbdcc1e6a8157aa53132717d90e7ac63c4e", "size": 6674918},
        "fungus_mt":
            {"sha256": "e0c6968372ba3fec892532a9231d1e8de789323156b3d23fb5216915a2811b27", "size": 659823},
        "fungus_nr":
            {"sha256": "c794d09fae34ebd10d97ee7d9075c3ba35c2560dae52786302e59081e2326aac", "size": 19786},
        "embplant_mt":
            {"sha256": "f3621444441fbf4fa98835999b20e7ac1673a8866cb119707200412c46263570", "size": 148338},
        "other_pt":
            {"sha256": "84346c42ff7e85e5d0859ecdee1de4f06f45fd072851d57ba31c894bcbee3cd8", "size": 412866}
        },
    "0.0.1.minima":
        {
        "embplant_nr":
            {"sha256": "603033541683b7c53fb63970c188eb2891844f6419eca342fa648f6ff5e29d71", "size": 16500},
        "embplant_pt":
            {"sha256": "a38f3d65009c75aa00a13c8a737fb8390843409671d0e061d56d26b6b9c7ed14", "size": 88006},
        "animal_mt":
            {"sha256": "2d4e8f441a531cbee64d1514542edf3d158b2e1fcedccc15d32b7a3e9fa0a5a3", "size": 14483},
        "fungus_mt":
            {"sha256": "903e3a3c82aaedece3033218e0623caccabce3d58866d1577ecbfdf0d7115cca", "size": 20281},
        "embplant_mt":
            {"sha256": "f3621444441fbf4fa98835999b20e7ac1673a8866cb119707200412c46263570", "size": 148338},
        "other_pt":
            {"sha256": "62fb528edcd956ce605c9292063d3d7a0f168d36b6f4f114bb1050236622cc41", "size": 141269},
        "fungus_nr":
            {"sha256": "4532be1e3f1cc2627c876f435dc72e6a9def6060f322e41f84c24912017f879a", "size": 2984},
        }
     }

SEED_DB_HASH = \
    {
    "0.0.0":
        {
        "embplant_nr":
            {"sha256": "e19365f85b3bda29aabb5cf1ceb5c814e667ba251b08388d805b52b1f1fe1445", "size": 18309},
        "embplant_pt":
            {"sha256": "e8e5f461d6e5fe67c5314e46c61265e6ff1079886af3b959609dde2be97d870d", "size": 15342405},
        "animal_mt":
            {"sha256": "a293c02e0c0496beb29383927cc7c643b13e09f8cfa03f7688b352315a43f898", "size": 30285897},
        "fungus_mt":
            {"sha256": "abbc7658c9431d11454f2fec75b7b5f4deeb12bc3590351dc93883655c7c194e", "size": 7977301},
        "embplant_mt":
            {"sha256": "2f28612e7c2280a7273738eded0dd2fcfb59c6153c1cf3bac15e4e7ed1bf4e89", "size": 407052},
        "other_pt":
            {"sha256": "a548538ef6560ededefb7e0d9c41f1dcb8585dccb1d06124ffb95e6770df2c6b", "size": 14667508}
        },
    "0.0.1":
        {
        "embplant_nr":
            {"sha256": "e19365f85b3bda29aabb5cf1ceb5c814e667ba251b08388d805b52b1f1fe1445", "size": 18309},
        "embplant_pt":
            {"sha256": "e8e5f461d6e5fe67c5314e46c61265e6ff1079886af3b959609dde2be97d870d", "size": 15342405},
        "animal_mt":
            {"sha256": "a293c02e0c0496beb29383927cc7c643b13e09f8cfa03f7688b352315a43f898", "size": 30285897},
        "fungus_mt":
            {"sha256": "abbc7658c9431d11454f2fec75b7b5f4deeb12bc3590351dc93883655c7c194e", "size": 7977301},
        "fungus_nr":
            {"sha256": "73bec014ebbbd88f7e51aa5575ea0cd05b8336a2b41808ff6c7cfb825271e5a9", "size": 375273},
        "embplant_mt":
            {"sha256": "2f28612e7c2280a7273738eded0dd2fcfb59c6153c1cf3bac15e4e7ed1bf4e89", "size": 407052},
        "other_pt":
            {"sha256": "a548538ef6560ededefb7e0d9c41f1dcb8585dccb1d06124ffb95e6770df2c6b", "size": 14667508}
        },
    "0.0.1.minima": {
        "embplant_nr":
            {"sha256": "e19365f85b3bda29aabb5cf1ceb5c814e667ba251b08388d805b52b1f1fe1445", "size": 18309},
        "embplant_pt":
            {"sha256": "9e988da116df42107cd760eaac247e76ab2b155d1072112790d16c5a25e5188d", "size": 156097},
        "animal_mt":
            {"sha256": "be089dbad2bfeb75ecfde1a48e01101c156f197d7874bf41d1a20c2c99edabb3", "size": 17135},
        "fungus_mt":
            {"sha256": "6e5434bbd880c5063e6932913fc76ccfca4c7038d12d7aea36eb06207915ed75", "size": 64609},
        "embplant_mt":
            {"sha256": "2f28612e7c2280a7273738eded0dd2fcfb59c6153c1cf3bac15e4e7ed1bf4e89", "size": 407052},
        "other_pt":
            {"sha256": "50311a0c96798b2745efbbe2d908f2d26da021fa2e8469fae82b5f5e752bb2db", "size": 192036},
        "fungus_nr":
            {"sha256": "1258941c417c3544e79b08d71be3f7d0590563b24c7bc7403358304357b20cc8", "size": 1582},
        }
    }


