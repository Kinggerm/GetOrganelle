#!/usr/bin/env python

import sympy
import scipy
try:
    from math import inf
except ImportError:
    inf = float("inf")
from optparse import OptionParser
import GetOrganelleLib
from GetOrganelleLib.versions import get_versions
from GetOrganelleLib.pipe_control_func import *
from GetOrganelleLib.assembly_parser import *
import time
import random
import subprocess
import sys
import os
from shutil import copyfile
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
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
NOT_DB_PATH = os.path.realpath(os.path.join(GO_LIB_PATH, "LabelDatabase"))
SEQ_DB_PATH = os.path.realpath(os.path.join(GO_LIB_PATH, "SeedDatabase"))
UTILITY_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "Utilities")


MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    PYTHON_VERSION = "2.7+"
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    PYTHON_VERSION = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)


def get_options(description, version):
    usage = str(os.path.basename(__file__)+" -g assembly_graph_file -F embplant_pt -o output --min-depth 10")
    parser = OptionParser(usage=usage, version=version, description=description, add_help_option=False)
    parser.add_option("-F", dest="organelle_type",
                      help="This flag should be followed with embplant_pt (embryophyta plant plastome), "
                           "other_pt (non-embryophyta plant plastome), embplant_mt "
                           "(plant mitochondrion), embplant_nr (plant nuclear ribosomal RNA), animal_mt "
                           "(animal mitochondrion), fungus_mt (fungus mitochondrion), "
                           "or embplant_mt,other_pt,fungus_mt "
                           "(the combination of any of above organelle genomes split by comma(s), "
                           "which might be computationally more intensive than separate runs), "
                           "or anonym (uncertain organelle genome type, with customized gene database "
                           "('--genes'), which is suggested only when the above database is genetically distant "
                           "from your sample) or raw (disentangle the raw graph directly without tagging).")
    parser.add_option("-g", dest="input_graph", help="Input assembly graph (fastg/gfa) file.")
    parser.add_option("-o", dest="output_base", help="Output directory. Overwriting files if directory exists.")
    parser.add_option('--min-depth', dest='min_depth', default=0., type=float,
                      help='Input a float or integer number. Filter graph file by a minimum depth. Default: %default.')
    parser.add_option('--max-depth', dest='max_depth', default=inf, type=float,
                      help='Input a float or integer number. filter graph file by a maximum depth. Default: %default.')
    parser.add_option("--genes", dest="genes_fasta",
                      help="Followed with a customized database (a fasta file or the base name of a "
                           "blast database) containing or made of ONE set of protein coding genes "
                           "and ribosomal RNAs extracted from ONE reference genome that you want to assemble. "
                           "Should be a list of databases split by comma(s) on a multi-organelle mode, "
                           "with the same list length to organelle_type (followed by '-F'). "
                           "This is optional for any organelle mentioned in '-F' but required for 'anonym'. "
                           "By default, certain database(s) in " + str(NOT_DB_PATH) + " would be used "
                           "contingent on the organelle types chosen (-F). "
                           "The default value no longer holds when '--genes' or '--ex-genes' is used.")
    parser.add_option("--ex-genes", dest="exclude_genes",
                      help="This is optional and Not suggested, since non-target contigs could contribute "
                           "information for better downstream coverage-based clustering. "
                           "Followed with a customized database (a fasta file or the base name of a "
                           "blast database) containing or made of protein coding genes "
                           "and ribosomal RNAs extracted from reference genome(s) that you want to exclude. "
                           "Could be a list of databases split by comma(s) but "
                           "NOT required to have the same list length to organelle_type (followed by '-F'). "
                           "The default value no longer holds when '--genes' or '--ex-genes' is used.")
    parser.add_option("--no-slim", dest="no_slim", default=False, action="store_true",
                      help="Disable slimming process and directly disentangle the original assembly graph. "
                           "Default: %default")
    parser.add_option("--slim-options", dest="slim_options", default="",
                      help="Other options for calling slim_fastg.py")
    parser.add_option("--depth-factor", dest="depth_factor", default=10.0, type=float,
                      help="Depth factor for differentiate genome type of contigs. "
                           "The genome type of contigs are determined by blast. "
                           "Default: %default")
    parser.add_option("--type-f", dest="type_factor", type=float, default=3.,
                      help="Type factor for identifying contig type tag when multiple tags exist in one contig. "
                           "Default:%default")
    parser.add_option("--contamination-depth", dest="contamination_depth", default=3., type=float,
                      help="Depth factor for confirming contamination in parallel contigs. Default: %default")
    parser.add_option("--contamination-similarity", dest="contamination_similarity", default=0.9, type=float,
                      help="Similarity threshold for confirming contaminating contigs. Default: %default")
    parser.add_option("--no-degenerate", dest="degenerate", default=True, action="store_false",
                      help="Disable making consensus from parallel contig based on nucleotide degenerate table.")
    parser.add_option("--degenerate-depth", dest="degenerate_depth", default=1.5, type=float,
                      help="Depth factor for confirming parallel contigs. Default: %default")
    parser.add_option("--degenerate-similarity", dest="degenerate_similarity", default=0.98, type=float,
                      help="Similarity threshold for confirming parallel contigs. Default: %default")
    parser.add_option("--disentangle-time-limit", dest="disentangle_time_limit", default=3600, type=int,
                      help="Time limit (second) for each try of disentangling a graph file as a circular "
                           "genome. Disentangling a graph as contigs is not limited. Default: %default")
    parser.add_option("--expected-max-size", dest="expected_max_size", default='200000', type=str,
                      help="Expected maximum target genome size(s) for disentangling. "
                           "Should be a list of INTEGER numbers split by comma(s) on a multi-organelle mode, "
                           "with the same list length to organelle_type (followed by '-F'). "
                           "Default: 250000 (-F embplant_pt/fungus_mt), 25000 (-F embplant_nr/animal_mt), "
                           "1000000 (-F embplant_mt/other_pt), "
                           "1000000,1000000,250000 (-F other_pt,embplant_mt,fungus_mt)")
    parser.add_option("--expected-min-size", dest="expected_min_size", default=10000, type=int,
                      help="Expected minimum target genome size(s) for disentangling. "
                           "Should be a list of INTEGER numbers split by comma(s) on a multi-organelle mode, "
                           "with the same list length to organelle_type (followed by '-F'). "
                           "Default: %default for all.")
    parser.add_option("--keep-all-polymorphic", dest="only_keep_max_cov", default=True, action="store_false",
                      help="By default, this script would pick the contig with highest coverage among all parallel "
                           "(polymorphic) contigs when degenerating was not applicable. "
                           "Choose this flag to export all combinations.")
    parser.add_option("--min-sigma", dest="min_sigma_factor", type=float, default=0.1,
                      help="Minimum deviation factor for excluding non-target contigs. Default:%default")
    parser.add_option("-t", dest="threads", type=int, default=1,
                      help="Maximum threads to use.")
    parser.add_option("--prefix", dest="prefix", default="",
                      help="Add extra prefix to resulting files under the output directory.")
    parser.add_option("--which-blast", dest="which_blast", default="",
                      help="Assign the path to BLAST binary files if not added to the path. "
                           "Default: try GetOrganelleDep/" + SYSTEM_NAME + "/ncbi-blast first, then $PATH")
    parser.add_option("--which-bandage", dest="which_bandage", default="",
                      help="Assign the path to bandage binary file if not added to the path. Default: try $PATH")
    parser.add_option("--keep-temp", dest="keep_temp_files", action="store_true", default=False,
                      help="Choose to keep the running temp/index files.")
    parser.add_option("--continue", dest="script_resume", default=False, action="store_true",
                      help="Several check points based on files produced, rather than on the log file, "
                           "so keep in mind that this script will not detect the difference "
                           "between this input parameters and the previous ones.")
    parser.add_option("--random-seed", dest="random_seed", default=12345, type=int,
                      help="Default: %default")
    parser.add_option("--verbose", dest="verbose_log", action="store_true", default=False,
                      help="Verbose output. Choose to enable verbose running log_handler.")
    parser.add_option("-h", dest="simple_help", default=False, action="store_true",
                      help="print brief introduction for frequently-used options.")
    parser.add_option("--help", dest="verbose_help", default=False, action="store_true",
                      help="print verbose introduction for all options.")
    
    if "--help" in sys.argv:
        parser.print_help()
        exit()
    elif "-h" in sys.argv:
        for not_often_used in ("--genes", "--ex-genes", "--slim-options", "--depth-factor", "--type-f",
                               "--contamination-depth", "--contamination-similarity", "--no-degenerate",
                               "--degenerate-depth", "--degenerate-similarity", "--disentangle-time-limit",
                               "--expected-max-size", "--expected-min-size", "--keep-all-polymorphic", "--min-sigma",
                               "--prefix", "--which-blast", "--which-bandage",
                               "--keep-temp", "--random-seed", "--verbose"):
            parser.remove_option(not_often_used)
        parser.remove_option("-F")
        parser.add_option("-F", dest="organelle_type",
                          help="Target organelle genome type: "
                               "embplant_pt/other_pt/embplant_mt/embplant_nr/animal_mt/fungus_mt/anonym/"
                               "embplant_pt,embplant_mt/other_pt,embplant_mt,fungus_mt")
        parser.remove_option("-g")
        parser.add_option("-g", dest="input_graph", help="Input assembly graph (fastg/gfa) file.")
        parser.remove_option("-o")
        parser.add_option("-o", dest="output_base", help="Output directory.")
        parser.remove_option("--min-depth")
        parser.add_option("--min-depth", dest="min_depth", default=0., type=float,
                          help="Minimum depth threshold of contigs. Default: %default.")
        parser.remove_option("--max-depth")
        parser.add_option("--max-depth", dest="max_depth", default=inf, type=float,
                          help="Maximum depth threshold of contigs. Default: %default.")
        parser.remove_option("--no-slim")
        parser.add_option("--no-slim", dest="no_slim",
                          help="Disable the slimming process and directly disentangle the assembly graph.")
        parser.remove_option("-t")
        parser.add_option("-t", dest="threads", type=int, default=1, help="Maximum threads to use. Default: %default.")
        parser.remove_option("--continue")
        parser.add_option("--continue", dest="script_resume", default=False, action="store_true",
                          help="Resume a previous run. Default: %default.")
        parser.print_help()
        sys.stdout.write("\n")
        exit()
    else:
        pass
    
    # redirect organelle types before parsing arguments
    redirect_organelle_types = {"plant_cp": "embplant_pt",
                                "plant_pt": "embplant_pt",
                                "plant_mt": "embplant_mt",
                                "plant_nr": "embplant_nr"}
    for go_arg, candidate_arg in enumerate(sys.argv):
        if go_arg > 1 and sys.argv[go_arg - 1] in {"-F", "-E"}:
            if candidate_arg in redirect_organelle_types:
                sys.argv[go_arg] = redirect_organelle_types[candidate_arg]
            elif "," in candidate_arg:
                new_arg = []
                for sub_arg in candidate_arg.split(","):
                    if sub_arg in redirect_organelle_types:
                        new_arg.append(redirect_organelle_types[sub_arg])
                    else:
                        new_arg.append(sub_arg)
                sys.argv[go_arg] = ",".join(new_arg)
    #
    try:
        (options, args) = parser.parse_args()
    except Exception as e:
        sys.stdout.write("\n############################################################################" + str(e))
        sys.stdout.write("\n\"-h\" for more usage")
        exit()
    else:
        if not options.output_base:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\n")
            sys.stdout.write("Missing option: output directory (followed after '-o')!\n")
            exit()
        if not options.input_graph:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\n")
            sys.stdout.write("Missing option: input assembly graph (followed after '-g')!\n")
            exit()
        if not options.organelle_type:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\n")
            sys.stdout.write("Missing option: organelle type (followed after '-F')!\n")
            exit()
        else:
            options.organelle_type = options.organelle_type.split(",")
        for sub_organelle_t in options.organelle_type:
            if sub_organelle_t not in {"embplant_pt", "other_pt", "embplant_mt", "embplant_nr", "animal_mt",
                                       "fungus_mt", "anonym"}:
                sys.stdout.write("\n############################################################################"
                                 "\nERROR: \"-F\" MUST be one of 'embplant_pt', 'other_pt', 'embplant_mt', "
                                 "'embplant_nr', 'animal_mt', 'fungus_mt', 'anonym', or a "
                                 "combination of above split by comma(s)!\n\n")
                exit()
            elif sub_organelle_t == "anonym":
                if not options.genes_fasta:
                    sys.stdout.write("\n############################################################################"
                                     "\nERROR: \"--genes\" must be specified when \"-F anonym\"!\n\n")
                    exit()
        assert options.threads > 0
        organelle_type_len = len(options.organelle_type)
        options.prefix = os.path.basename(options.prefix)
        if not os.path.isdir(options.output_base):
            os.mkdir(options.output_base)
            options.script_resume = False
        log_handler = simple_log(logging.getLogger(), options.output_base, options.prefix + "get_org.")
        log_handler.info("")
        log_handler.info(description)
        log_handler.info("Python " + str(sys.version).replace("\n", " "))
        # log versions of dependencies
        lib_versions_info = list()
        lib_versions_info.append("numpy " + np.__version__)
        lib_versions_info.append("sympy " + sympy.__version__)
        lib_versions_info.append("scipy " + scipy.__version__)
        log_handler.info("Python libs: " + "; ".join(lib_versions_info))
        dep_versions_info = []

        if not options.no_slim:
            if not options.which_blast:
                try_this_bin = os.path.join(GO_DEP_PATH, "ncbi-blast", "blastn")
                if os.path.isfile(try_this_bin) and executable(try_this_bin):
                    output, err = subprocess.Popen(
                        try_this_bin + " -version", stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, shell=True).communicate()
                    if "not found" in output.decode("utf8"):
                        log_handler.warning(output.decode("utf8"))
                    else:
                        options.which_blast = os.path.split(try_this_bin)[0]
            if executable(os.path.join(options.which_blast, "blastn")):
                output, err = subprocess.Popen(
                    os.path.join(options.which_blast, "blastn") + " -version", stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, shell=True).communicate()
                this_lines = output.decode("utf8").split("\n")[:2]
                dep_versions_info.append("Blast " + this_lines[1].strip().split()[2].replace(",", "").strip())
            else:
                log_handler.error(os.path.join(options.which_blast, "blastn") + " not accessible!")
                exit()
        if options.genes_fasta and not executable(os.path.join(options.which_blast, "makeblastdb")):
            log_handler.error(os.path.join(options.which_blast, "makeblastdb") + " not accessible!")
            exit()
        if executable(os.path.join(options.which_bandage, "Bandage -v")):
            output, err = subprocess.Popen(
                os.path.join(options.which_bandage, "Bandage") + " -v", stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, shell=True).communicate()
            dep_versions_info.append("Bandage " + output.decode("utf8").strip().split()[-1])
        if dep_versions_info:
            log_handler.info("Dependencies: " + "; ".join(dep_versions_info))
        log_handler.info(" ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n")

        # for sub_organelle_t in options.organelle_type:
        #     if sub_organelle_t in ("animal_mt", "fungus_mt"):
        #         log_handler.info("Animal & fungus samples have much higher substitution rates, therefore no guarantee "
        #                          "for success rate in animal & fungus samples with default label database!")
        #         log_handler.info("If you failed after this run, a customized fasta-format label database made "
        #                          "of protein coding & ribosomal genes extracted from a close-related species "
        #                          "(--genes) is suggested for another run!\n")
        #         break

        log_handler = timed_log(log_handler, options.output_base, options.prefix + "get_org.")

        # using the default
        if "--genes" not in sys.argv:
            options.genes_fasta = []  # None] * organelle_type_len
        else:
            temp_val_len = len(str(options.genes_fasta).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --genes is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            temp_vals = []
            for sub_genes in str(options.genes_fasta).split(","):
                # if sub_genes == "":
                #     temp_vals.append(sub_genes)
                if not (os.path.exists(sub_genes) or os.path.exists(remove_db_postfix(sub_genes) + ".nhr")):
                    log_handler.error(sub_genes + " not found!")
                    exit()
                else:
                    temp_vals.append(sub_genes)
            options.genes_fasta = temp_vals
        if "--ex-genes" not in sys.argv:
            options.exclude_genes = []
        else:
            temp_vals = []
            for sub_genes in str(options.exclude_genes).split(","):
                if not (os.path.exists(sub_genes) or os.path.exists(remove_db_postfix(sub_genes) + ".nhr")):
                    log_handler.error(sub_genes + " not found!")
                    exit()
                else:
                    temp_vals.append(sub_genes)
            options.exclude_genes = temp_vals

        if "--expected-max-size" not in sys.argv:
            raw_default_value = int(str(options.expected_max_size))
            options.expected_max_size = []
            for got_t, sub_organelle_t in enumerate(options.organelle_type):
                if sub_organelle_t == "embplant_pt":
                    options.expected_max_size.append(raw_default_value)
                elif sub_organelle_t == "embplant_mt":
                    options.expected_max_size.append(int(raw_default_value * 4))
                elif sub_organelle_t == "fungus_mt":
                    options.expected_max_size.append(raw_default_value)
                elif sub_organelle_t in ("embplant_nr", "animal_mt"):
                    options.expected_max_size.append(int(raw_default_value / 10))
                elif sub_organelle_t == "anonym":
                    ref_seqs = read_fasta(options.genes_fasta[got_t])[1]
                    options.expected_max_size.append(10 * sum([len(this_seq) for this_seq in ref_seqs]))
                    log_handler.info("Setting '--expected-max-size " + str(options.expected_max_size) +
                                     "' for estimating the word size value for anonym type.")
        else:
            temp_val_len = len(str(options.expected_max_size).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --expected-max-size is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            try:
                options.expected_max_size = [int(sub_size) for sub_size in str(options.expected_max_size).split(",")]
            except ValueError:
                log_handler.error("Invalid --expected-max-size value(s): " + str(options.expected_max_size))
                exit()

        if "--expected-min-size" not in sys.argv:
            raw_default_value = int(str(options.expected_min_size))
            options.expected_min_size = []
            for sub_organelle_t in options.organelle_type:
                options.expected_min_size.append(raw_default_value)
        else:
            temp_val_len = len(str(options.expected_min_size).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --expected-min-size is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            try:
                options.expected_min_size = [int(sub_size) for sub_size in str(options.expected_min_size).split(",")]
            except ValueError:
                log_handler.error("Invalid --expected-min-size value(s): " + str(options.expected_min_size))
                exit()
        # check slim options
        if options.slim_options.strip():
            slim_op_parts = options.slim_options.split()
            ban_argv = {"-F": 2, "-E": 2, "--min-depth": 2, "--max-depth": 2, "--merge": 1,
                        "--include": 2, "--include-priority": 2, "--exclude": 2, "--exclude-priority": 2,
                        "--no-hits-labeled-tab": 1, "-o": 2, "--prefix": 2, "--out-base": 2, "--log": 1,
                        "--verbose": 1, "--continue": 1, "--no-overwrite": 1, "--which-blast": 2, "-t": 2,
                        "--threads": 2}
            remove_ops = []
            go_op = 0
            while go_op < len(slim_op_parts):
                if slim_op_parts[go_op] in ban_argv:
                    for repeat_time in range(ban_argv[slim_op_parts[go_op]]):
                        remove_ops.append(slim_op_parts.pop(go_op))
                else:
                    go_op += 1
            if remove_ops:
                log_handler.info("Options \"" + " ".join(remove_ops) +
                                 "\" taken/invalid for wrapped slim_fastg.py, removed.")
                options.slim_options = " ".join(slim_op_parts)
        random.seed(options.random_seed)
        np.random.seed(options.random_seed)
        return options, log_handler


def slim_spades_result(organelle_types, in_custom, ex_custom, graph_in, graph_out_base,
                       verbose_log, log_handler, threads, which_slim, which_blast="", other_options="",
                       resume=False, keep_temp=False):
    include_priority_db = []
    exclude_db = []
    if in_custom or ex_custom:
        include_priority_db = in_custom
        exclude_db = ex_custom
    else:
        if organelle_types == ["embplant_pt"]:
            include_priority_db = [os.path.join(NOT_DB_PATH, "embplant_pt.fasta"),
                                   os.path.join(NOT_DB_PATH, "embplant_mt.fasta")]
        elif organelle_types == ["embplant_mt"]:
            include_priority_db = [os.path.join(NOT_DB_PATH, "embplant_mt.fasta"),
                                   os.path.join(NOT_DB_PATH, "embplant_pt.fasta")]
        # elif organelle_types == ["embplant_mt"]:
        #     include_priority_db = [os.path.join(NOT_DB_PATH, "embplant_mt.fasta")]
        #     exclude_db = [os.path.join(NOT_DB_PATH, "embplant_pt.fasta")]
        else:
            include_priority_db = [os.path.join(NOT_DB_PATH, sub_organelle_t + ".fasta")
                                   for sub_organelle_t in organelle_types]
    if resume:
        if os.path.exists(graph_out_base + ".fastg") and os.path.exists(graph_out_base + ".csv"):
            if log_handler:
                log_handler.info("Slimming " + graph_in + " ... skipped.")
            return 0
    run_command = ""
    if include_priority_db:
        run_command += " --include-priority " + ",".join(include_priority_db)
    if exclude_db:
        run_command += " --exclude " + ",".join(exclude_db)
    which_bl_str = " --which-blast " + which_blast if which_blast else ""
    run_command = os.path.join(which_slim, "slim_fastg.py") + " -t " + str(threads) + which_bl_str + \
                  " " + graph_in + " --out-base " + graph_out_base + " " + run_command + " --log " + \
                  "--verbose " * int(bool(verbose_log)) + "--continue " * int(bool(resume)) + other_options + \
                  " --keep-temp" * int(bool(keep_temp))
    # + ' -o ' + out_base + (' --prefix ' + prefix if prefix else "")
    slim_spades = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log:
        log_handler.info(run_command)
    output, err = slim_spades.communicate()
    # if "not recognized" in output.decode("utf8") or "command not found" in output.decode("utf8") \
    #         or "no such file" in output.decode("utf8"):
    #     if log_handler:
    #         log_handler.error("Slimming " + graph_in + " failed.")
    #     return os.path.join(which_slim, "slim_fastg.py") + " not accessible!\n" + output.decode("utf8").strip()
    if "failed" in output.decode("utf8") or "error" in output.decode("utf8") or "Error" in output.decode("utf8"):
        if log_handler:
            log_handler.error("Slimming " + graph_in + " failed.")
        return "\n" + output.decode("utf8").strip()
    else:
        if log_handler:
            if verbose_log:
                log_handler.info(output.decode("utf8"))
            log_handler.info("Slimming " + graph_in + " finished!")
        return 0


def extract_organelle_genome(out_base, slim_out_fg, slim_out_csv, organelle_prefix, organelle_type, blast_db,
                             verbose, log_handler, expected_maximum_size, expected_minimum_size, no_slim, options):
    def disentangle_assembly(fastg_file, tab_file, output, weight_factor, log_dis, time_limit, type_factor=3.,
                             mode="embplant_pt", blast_db_base="embplant_pt", contamination_depth=3.,
                             contamination_similarity=0.95, degenerate=True,
                             degenerate_depth=1.5, degenerate_similarity=0.98,
                             expected_max_size=inf, expected_min_size=0, hard_cov_threshold=10.,
                             min_sigma_factor=0.1, here_only_max_c=True, here_acyclic_allowed=False,
                             here_verbose=False, timeout_flag_str="'--disentangle-time-limit'", temp_graph=None):
        @set_time_limit(time_limit, flag_str=timeout_flag_str)
        def disentangle_inside(fastg_f, tab_f, o_p, w_f, log_in, type_f=3., mode_in="embplant_pt",
                               in_db_n="embplant_pt", c_d=3., c_s=0.95,
                               deg=True, deg_dep=1.5, deg_sim=0.98, hard_c_t=10., min_s_f=0.1, max_c_in=True,
                               max_s=inf, min_s=0, acyclic_allowed_in=False, verbose_in=False, in_temp_graph=None):
            if acyclic_allowed_in:
                log_in.info("Disentangling " + fastg_f + " as contig(s) ... ")
            else:
                log_in.info("Disentangling " + fastg_f + " as a circular genome ... ")
            image_produced = False
            input_graph = Assembly(fastg_f)
            if no_slim:
                input_graph.estimate_copy_and_depth_by_cov(mode=mode_in, log_handler=log_in, verbose=verbose_in)
                target_results = input_graph.estimate_copy_and_depth_precisely(
                    broken_graph_allowed=acyclic_allowed_in, return_new_graphs=True, verbose=verbose_in,
                    log_handler=log_in)
            else:
                target_results = input_graph.find_target_graph(tab_f,
                                                               mode=mode_in, database_name=in_db_n, type_factor=type_f,
                                                               log_hard_cov_threshold=hard_c_t,
                                                               contamination_depth=c_d,
                                                               contamination_similarity=c_s,
                                                               degenerate=deg, degenerate_depth=deg_dep,
                                                               degenerate_similarity=deg_sim,
                                                               expected_max_size=max_s, expected_min_size=min_s,
                                                               only_keep_max_cov=max_c_in,
                                                               min_sigma_factor=min_s_f,
                                                               weight_factor=w_f,
                                                               broken_graph_allowed=acyclic_allowed_in,
                                                               log_handler=log_in, verbose=verbose_in,
                                                               temp_graph=in_temp_graph)
            if len(target_results) > 1:
                log_in.warning(str(len(target_results)) + " sets of graph detected!")
            # log_in.info("Slimming and disentangling graph finished!")

            log_in.info("Writing output ...")
            degenerate_base_used = False
            if acyclic_allowed_in:
                contig_num = set()
                still_complete = []
                for go_res, res in enumerate(target_results):
                    broken_graph = res["graph"]
                    count_path = 0
                    these_paths = broken_graph.get_all_paths(mode=mode_in, log_handler=log_in)
                    for this_paths, other_tag in these_paths:
                        count_path += 1
                        all_contig_str = []
                        contig_num.add(len(this_paths))
                        contigs_are_circular = []
                        for go_contig, this_p_part in enumerate(this_paths):
                            this_contig = broken_graph.export_path(this_p_part)
                            if DEGENERATE_BASES & set(this_contig.seq):
                                degenerate_base_used = True
                            if this_contig.label.endswith("(circular)"):
                                contigs_are_circular.append(True)
                            else:
                                contigs_are_circular.append(False)
                            if len(this_paths) == 1 and contigs_are_circular[-1]:
                                all_contig_str.append(this_contig.fasta_str())
                            else:
                                all_contig_str.append(">contig_" + str(go_contig + 1) + "--" + this_contig.label +
                                                      "\n" + this_contig.seq + "\n")
                        if len(all_contig_str) == 1 and set(contigs_are_circular) == {True}:
                            still_complete.append(True)
                        else:
                            still_complete.append(False)
                        if still_complete[-1]:
                            out_n = o_p + ".complete.graph" + str(go_res) + "." + \
                                    str(count_path) + other_tag + ".path_sequence.fasta"
                            log_in.info("Writing PATH" + str(count_path) + " of complete " + mode_in + " to " + out_n)
                        else:
                            out_n = o_p + ".contigs.graph" + str(go_res + 1) + other_tag + "." + \
                                    str(count_path) + ".path_sequence.fasta"
                            log_in.info("Writing PATH" + str(count_path) + " of " + mode_in + " contig(s) to " + out_n)
                        open(out_n, "w").write("\n".join(all_contig_str))
                    if set(still_complete[-len(these_paths):]) == {True}:
                        log_in.info(
                            "Writing GRAPH to " + o_p + ".complete.graph" + str(go_res + 1) + ".selected_graph.gfa")
                        broken_graph.write_to_gfa(o_p + ".complete.graph" + str(go_res + 1) + ".selected_graph.gfa")
                        image_produced = draw_assembly_graph_using_bandage(
                            input_graph_file=o_p + ".complete.graph" + str(go_res + 1) + ".selected_graph.gfa",
                            output_image_file=o_p + ".complete.graph" + str(go_res + 1) + ".selected_graph.png",
                            assembly_graph_ob=broken_graph,
                            log_handler=log_handler, verbose_log=verbose_in, which_bandage=options.which_bandage)
                    else:
                        log_in.info(
                            "Writing GRAPH to " + o_p + ".contigs.graph" + str(go_res + 1) + ".selected_graph.gfa")
                        broken_graph.write_to_gfa(o_p + ".contigs.graph" + str(go_res + 1) + ".selected_graph.gfa")
                        image_produced = draw_assembly_graph_using_bandage(
                            input_graph_file=o_p + ".contigs.graph" + str(go_res + 1) + ".selected_graph.gfa",
                            output_image_file=o_p + ".contigs.graph" + str(go_res + 1) + ".selected_graph.png",
                            assembly_graph_ob=broken_graph,
                            log_handler=log_handler, verbose_log=verbose_in, which_bandage=options.which_bandage)
                if set(still_complete) == {True}:
                    log_in.info("Result status of " + mode_in + ": circular genome")
                else:
                    log_in.info("Result status of " + mode_in + ": " +
                                ",".join(sorted([str(c_n) for c_n in contig_num])) + " contig(s)")
            else:
                for go_res, res in enumerate(target_results):
                    go_res += 1
                    idealized_graph = res["graph"]
                    count_path = 0
                    for this_path, other_tag in idealized_graph.get_all_circular_paths(mode=mode_in,
                                                                                       log_handler=log_in):
                        count_path += 1
                        out_n = o_p + ".complete.graph" + str(go_res) + "." + str(
                            count_path) + other_tag + ".path_sequence.fasta"
                        this_seq_obj = idealized_graph.export_path(this_path)
                        if DEGENERATE_BASES & set(this_seq_obj.seq):
                            degenerate_base_used = True
                        open(out_n, "w").write(this_seq_obj.fasta_str())
                        log_in.info("Writing PATH" + str(count_path) + " of complete " + mode_in + " to " + out_n)
                    log_in.info("Writing GRAPH to " + o_p + ".complete.graph" + str(go_res) + ".selected_graph.gfa")
                    idealized_graph.write_to_gfa(o_p + ".complete.graph" + str(go_res) + ".selected_graph.gfa")
                    image_produced = draw_assembly_graph_using_bandage(
                        input_graph_file=o_p + ".complete.graph" + str(go_res) + ".selected_graph.gfa",
                        output_image_file=o_p + ".complete.graph" + str(go_res) + ".selected_graph.png",
                        assembly_graph_ob=idealized_graph,
                        log_handler=log_handler, verbose_log=verbose_in, which_bandage=options.which_bandage)
                log_in.info("Result status of " + mode_in + ": circular genome")
            if degenerate_base_used:
                log_in.warning("Degenerate base(s) used!")
            if not acyclic_allowed_in:
                if slim_out_csv:
                    if image_produced:
                        log_in.info("Please check the produced assembly image"
                                    " or manually visualize " + slim_out_fg + " and load " + slim_out_csv +
                                    " using Bandage to confirm the final result.")
                    else:
                        log_in.info("Please visualize " + slim_out_fg + " and load " + slim_out_csv +
                                    " using Bandage to confirm the final result.")
            log_in.info("Writing output finished.")

        disentangle_inside(fastg_f=fastg_file, tab_f=tab_file, o_p=output, w_f=weight_factor, log_in=log_dis,
                           type_f=type_factor, mode_in=mode, in_db_n=blast_db_base,
                           c_d=contamination_depth, c_s=contamination_similarity,
                           deg=degenerate, deg_dep=degenerate_depth, deg_sim=degenerate_similarity,
                           hard_c_t=hard_cov_threshold, min_s_f=min_sigma_factor, max_c_in=here_only_max_c,
                           max_s=expected_max_size, min_s=expected_min_size,
                           acyclic_allowed_in=here_acyclic_allowed, verbose_in=here_verbose, in_temp_graph=temp_graph)

    # start
    timeout_flag = "'--disentangle-time-limit'"
    export_succeeded = False
    path_prefix = os.path.join(out_base, organelle_prefix)
    graph_temp_file = path_prefix + ".temp.gfa" if options.keep_temp_files else None
    try:
        """disentangle"""
        # if it is the first round (the largest kmer), copy the slimmed result to the main spades output
        # if go_k == 0:
        #     main_spades_folder = os.path.split(kmer_dir)[0]
        #     os.system("cp " + out_fastg + " " + main_spades_folder)
        #     os.system("cp " + out_csv + " " + main_spades_folder)
        disentangle_assembly(fastg_file=slim_out_fg, blast_db_base=blast_db, mode=organelle_type,
                             tab_file=slim_out_csv, output=path_prefix,
                             weight_factor=100, type_factor=options.type_factor, 
                             hard_cov_threshold=options.depth_factor,
                             contamination_depth=options.contamination_depth,
                             contamination_similarity=options.contamination_similarity,
                             degenerate=options.degenerate, degenerate_depth=options.degenerate_depth,
                             degenerate_similarity=options.degenerate_similarity,
                             expected_max_size=expected_maximum_size,
                             expected_min_size=expected_minimum_size,
                             here_only_max_c=options.only_keep_max_cov,
                             min_sigma_factor=options.min_sigma_factor,
                             here_acyclic_allowed=False, here_verbose=verbose, log_dis=log_handler,
                             time_limit=options.disentangle_time_limit, timeout_flag_str=timeout_flag,
                             temp_graph=graph_temp_file)
        # currently time is not limited for exporting contigs
    except ImportError as e:
        log_handler.warning("Disentangling failed: numpy/scipy/sympy not installed!")
        if verbose:
            log_handler.error(str(e))
        return False
    except AttributeError as e:
        if verbose:
            raise e
    except RuntimeError:
        log_handler.info("Disentangling timeout. (see " + timeout_flag + " for more)")
    except ProcessingGraphFailed as e:
        log_handler.info("Disentangling failed: " + str(e).strip())
    except Exception as e:
        raise e
    else:
        export_succeeded = True

    if not export_succeeded:
        try:
            """disentangle the graph as contig(s)"""
            disentangle_assembly(fastg_file=slim_out_fg, blast_db_base=blast_db, mode=organelle_type,
                                 tab_file=slim_out_csv, output=path_prefix,
                                 weight_factor=100, type_factor=options.type_factor,
                                 here_verbose=verbose, log_dis=log_handler,
                                 hard_cov_threshold=options.depth_factor * 0.8,
                                 contamination_depth=options.contamination_depth,
                                 contamination_similarity=options.contamination_similarity,
                                 degenerate=options.degenerate,
                                 degenerate_depth=options.degenerate_depth,
                                 degenerate_similarity=options.degenerate_similarity,
                                 expected_max_size=expected_maximum_size,
                                 expected_min_size=expected_minimum_size,
                                 min_sigma_factor=options.min_sigma_factor,
                                 here_only_max_c=options.only_keep_max_cov, here_acyclic_allowed=True,
                                 time_limit=3600, timeout_flag_str=timeout_flag,
                                 temp_graph=graph_temp_file)
        except (ImportError, AttributeError) as e:
            if verbose:
                log_handler.error(str(e))
        except RuntimeError:
            log_handler.info("Disentangling timeout. (see " + timeout_flag + " for more)")
        except ProcessingGraphFailed as e:
            log_handler.info("Disentangling failed: " + str(e).strip())
        except Exception as e:
            raise e
        else:
            export_succeeded = True
            if slim_out_csv:
                log_handler.info(
                    "Please visualize " + slim_out_fg + " and load " +
                    slim_out_csv + " to confirm the incomplete result.")
            # log.info("-------------------------------------------------------")
            log_handler.info("If the result is nearly complete, you can try join_spades_fastg_by_blast.py to fill "
                             "N-gaps in-between contigs with a closely-related reference.")
            log_handler.info("If you have questions for us, "
                             "please provide us with the get_org.log.txt file "
                             "and the graph in the format you like!")
            # log.info("-------------------------------------------------------")
        if not export_succeeded:
            if slim_out_csv:
                log_handler.info(
                    "Please visualize " + slim_out_fg + " and load " +
                    slim_out_csv + " to confirm the incomplete result.")
            # log.info("-------------------------------------------------------")
            log_handler.info("If you have questions for us, "
                             "please provide us with the get_org.log.txt file "
                             "and the graph in the format you like!")
    return export_succeeded


def main():
    time0 = time.time()
    title = "GetOrganelle v" + str(get_versions()) + \
            "\n" \
            "\nget_organelle_from_assembly.py isolates organelle genomes from assembly graph." \
            "\nFind updates in https://github.com/Kinggerm/GetOrganelle and see README.md for more information." \
            "\n"
    options, log_handler = get_options(description=title, version=get_versions())
    try:
        if executable(os.path.join(UTILITY_PATH, "slim_fastg.py")):
            which_slim = UTILITY_PATH
        elif executable(os.path.join(PATH_OF_THIS_SCRIPT, "slim_fastg.py")):
            which_slim = PATH_OF_THIS_SCRIPT
        elif executable("slim_fastg.py"):
            which_slim = ""
        else:
            raise Exception("slim_fastg.py not found!")
        log_handler.info("Processing assembly graph ...")
        if os.path.getsize(options.input_graph) == 0:
            raise Exception("No vertices found in " + options.input_graph + "!")
        processed_graph_file = os.path.join(options.output_base, options.prefix + "initial_assembly_graph.fastg")
        if options.max_depth != inf or options.min_depth != 0. or options.input_graph.endswith(".gfa"):
            this_graph = Assembly(options.input_graph, max_cov=options.max_depth, min_cov=options.min_depth)
            this_graph.write_to_fastg(
                out_file=processed_graph_file, check_postfix=False, rename_if_needed=True,
                out_renaming_table=os.path.join(
                    options.output_base, options.prefix + "initial_assembly_graph.vertex_trans.tab"),
                echo_rename_warning=True, log_handler=log_handler)
        else:
            copyfile(options.input_graph, processed_graph_file)
        log_handler.info("Processing assembly graph finished.\n")
        if os.path.getsize(processed_graph_file) == 0:
            raise Exception("No vertices left in " + processed_graph_file + "! Please adjust the depth range!")
        if options.no_slim:
            slimmed_fastg_file = processed_graph_file
            slimmed_csv_file = None
        else:
            log_handler.info("Slimming assembly graph ...")
            slimmed_graph_file_base = os.path.join(options.output_base, options.prefix + "slimmed_assembly_graph")
            slimmed_fastg_file = os.path.join(options.output_base, options.prefix + "slimmed_assembly_graph.fastg")
            slimmed_csv_file = os.path.join(options.output_base, options.prefix + "slimmed_assembly_graph.csv")
            run_stat = slim_spades_result(organelle_types=options.organelle_type, in_custom=options.genes_fasta,
                                          ex_custom=options.exclude_genes, graph_in=processed_graph_file,
                                          graph_out_base=slimmed_graph_file_base, verbose_log=options.verbose_log,
                                          log_handler=log_handler, threads=options.threads,
                                          which_blast=options.which_blast,
                                          which_slim=which_slim, other_options=options.slim_options,
                                          resume=options.script_resume, keep_temp=options.keep_temp_files)
            if run_stat:
                log_handler.error(run_stat + "\n")
                exit()
            else:
                if os.path.getsize(slimmed_fastg_file) == 0:
                    return "Slimming " + processed_graph_file + " finished with no target organelle contigs found!"
                log_handler.info("Slimming assembly graph finished.\n")

        organelle_type_prefix = []
        duplicated_o_types = {o_type: 1
                              for o_type in options.organelle_type if options.organelle_type.count(o_type) > 1}
        for here_type in options.organelle_type:
            if here_type in duplicated_o_types:
                organelle_type_prefix.append(here_type + "-" + str(duplicated_o_types[here_type]))
                duplicated_o_types[here_type] += 1
            else:
                organelle_type_prefix.append(here_type)
        for go_t, sub_organelle_type in enumerate(options.organelle_type):
            og_prefix = options.prefix + organelle_type_prefix[go_t]
            graph_existed = bool([gfa_f for gfa_f in os.listdir(options.output_base)
                                  if gfa_f.startswith(og_prefix) and gfa_f.endswith(".selected_graph.gfa")])
            fasta_existed = bool([fas_f for fas_f in os.listdir(options.output_base)
                                  if fas_f.startswith(og_prefix) and fas_f.endswith(".path_sequence.fasta")])
            if options.script_resume and graph_existed and fasta_existed:
                log_handler.info("Extracting " + sub_organelle_type + " from the assemblies ... skipped.\n")
            else:
                # log_handler.info("Parsing assembly graph and outputting ...")
                log_handler.info("Extracting " + sub_organelle_type + " from the assemblies ...")
                if options.genes_fasta:
                    db_base_name = remove_db_postfix(os.path.basename(options.genes_fasta[go_t]))
                else:
                    db_base_name = sub_organelle_type
                ext_res = extract_organelle_genome(out_base=options.output_base,
                                                   slim_out_fg=slimmed_fastg_file, slim_out_csv=slimmed_csv_file,
                                                   organelle_prefix=og_prefix,
                                                   organelle_type=sub_organelle_type,
                                                   blast_db=db_base_name,
                                                   verbose=options.verbose_log, log_handler=log_handler,
                                                   expected_minimum_size=options.expected_min_size[go_t],
                                                   expected_maximum_size=options.expected_max_size[go_t],
                                                   no_slim=options.no_slim, options=options)
                if ext_res:
                    log_handler.info("Extracting " + sub_organelle_type + " from the assemblies finished.\n")
                else:
                    log_handler.info("Extracting " + sub_organelle_type + " from the assemblies failed.\n")
        log_handler = simple_log(log_handler, options.output_base, prefix=options.prefix + "get_org.")
        log_handler.info("\nTotal cost " + "%.2f" % (time.time() - time0) + " s")
        log_handler.info("Thank you!")
    except:
        log_handler.exception("")
        log_handler = simple_log(log_handler, options.output_base, prefix=options.prefix + "get_org.")
        log_handler.info("\nTotal cost " + "%.2f" % (time.time() - time0) + " s")
        log_handler.info("Please email jinjianjun@mail.kib.ac.cn or phylojin@163.com if you find bugs!")
        log_handler.info("Please provide me with the get_org.log.txt file!")
    logging.shutdown()


if __name__ == '__main__':
    main()


"""Copyright 2016 Jianjun Jin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License."""