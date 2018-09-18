import signal
import time
import logging
import sys
import os
from multiprocessing import Pool


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
            raise RuntimeError

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
