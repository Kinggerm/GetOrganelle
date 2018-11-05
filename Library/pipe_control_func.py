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

