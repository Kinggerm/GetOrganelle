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
DEAD_CODE = {"2.7+": 32512, "3.5+": 127}[python_version]

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(PATH_OF_THIS_SCRIPT, ".."))
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.assembly_parser import Assembly
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]


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


def set_time_limit(num, flag_str="'--time-limit'"):
    def wrap(func):
        def handle(sig_num, interrupted_stack_frame):
            raise RuntimeError("\n\nIncrease " + flag_str + " to meet the specific need of data!")

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
    return True if os.access(test_this, os.X_OK) or getstatusoutput(test_this)[0] != DEAD_CODE else False


def draw_assembly_graph_using_bandage(input_graph_file, output_image_file, assembly_graph_ob=None,
                                      resume=False, log_handler=None, verbose_log=False, which_bandage=""):
    if resume and os.path.exists(output_image_file):
        return True
    elif executable(os.path.join(which_bandage, "Bandage -v")):
        assert output_image_file.split(".")[-1] in {"png", "jpg", "svg"}, "Output image format must be png/jpg/svg!"
        temp_file = os.path.join(os.path.split(output_image_file)[0], os.path.split(output_image_file)[1])
        # preparing for color
        if assembly_graph_ob is None:
            assembly_graph_ob = Assembly(input_graph_file)
            assembly_graph_ob.estimate_copy_and_depth_by_cov(mode="all", verbose=verbose_log, log_handler=log_handler)
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
        if "error" in output.decode("utf8") or "Abort trap" in output.decode("utf8") or "Error:" in output.decode("utf8"):
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
    if "error" in output.decode("utf8") or "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
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
    if "(ERR)" in str(output) or "Error:" in str(output) or "error" in str(output):
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
        if "unrecognized option" in output.decode("utf8"):
            this_command = os.path.join(which_bowtie2, "bowtie2-build") + " --seed " + str(random_seed) + \
                           " " + seed_file + " " + seed_index_base
            building = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            if not silent and verbose_log:
                if log_handler:
                    log_handler.info(this_command)
                else:
                    sys.stdout.write(this_command + "\n")
            output, err = building.communicate()
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
        make_seed_bowtie2 = subprocess.Popen(this_command,
                                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            if not silent:
                if log_handler:
                    log_handler.error("\n" + output.decode("utf8"))
                else:
                    sys.stdout.write("\nError: " + output.decode("utf8") + "\n")
            raise Exception("")
        elif "No such file" in output.decode("utf8") or "not found" in output.decode("utf8"):
            if not silent:
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
        if os.path.exists(os.path.join(sample_out_dir, prefix + "filtered_spades", "spades.log")):
            with open(os.path.join(sample_out_dir, prefix + "filtered_spades", "spades.log"), "r") as spades_log:
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
            "Error" not in output.decode("utf8") and \
            "error" not in output.decode("utf8"):
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


def unzip(source, target, line_limit=1000000000000000, verbose_log=False, log_handler=None):
    target_temp = target + ".Temp"
    try_commands = [
        "tar -x -f " + source + " -O | head -n " + str(line_limit) + " > " + target_temp,
        "gunzip -c " + source + " | head -n " + str(line_limit) + " > " + target_temp]
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
                "Error" not in output.decode("utf8") and \
                "error" not in output.decode("utf8"):
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
