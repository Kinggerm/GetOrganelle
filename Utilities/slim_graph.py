#!/usr/bin/env python
# coding: utf8
import time
import os
import sys
import subprocess
try:
    # python2
    import commands
except:
    pass
inf = float("inf")
from optparse import OptionParser, OptionConflictError
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
import GetOrganelleLib
from GetOrganelleLib.versions import get_versions
from GetOrganelleLib.pipe_control_func import *
import math
import copy
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
VALID_PATH_CHARS = set("1234567890qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM,./ _-+")

_GO_PATH = GO_PATH
_LBL_DB_PATH = LBL_DB_PATH
_SEQ_DB_PATH = SEQ_DB_PATH


def check_path_chars(path_or_file_name):
    for this_char in path_or_file_name:
        if this_char not in VALID_PATH_CHARS:
            raise Exception("Invalid char \"" + this_char + "\" from the input file/path " + path_or_file_name + "!")


def get_options(print_title):
    usage = "python "+str(os.path.basename(__file__)+" your_fastg_files -F embplant_pt -E embplant_mt")
    parser = OptionParser(usage=usage)
    # parser.add_option("-o", dest="out_fastg_file", help="Output file")
    # filters
    parser.add_option("-F", dest="organelle_types",
                      help="followed with mode embplant_pt, other_pt, embplant_mt, embplant_nr, animal_mt, fungus_mt, "
                           "fungus_nr "
                           "(which means embryophyta plastome, non-embryophyta plastome, "
                           "plant mitochondrion, plant nuclear ribosomal RNA, animal mitochondrion, "
                           "fungus mitochondrion, fungus nuclear ribosomal RNA separately), "
                           "or a combination of above split by comma(s) "
                           "(corresponds to certain arguments as following listed). "
                           "\t"
                           " ------------------------------------------------------ "
                           "\nembplant_pt \t \" --include-priority " + os.path.join(LBL_DB_PATH, "embplant_pt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nother_pt \t \" --include-priority " + os.path.join(LBL_DB_PATH, "other_pt.fasta") + "\""
                           " ------------------------------------------------------ "                                                                                                                                        
                           "\nembplant_mt \t \" --include-priority " + os.path.join(LBL_DB_PATH, "embplant_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nembplant_nr \t \" --include-priority " + os.path.join(LBL_DB_PATH, "embplant_nr.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nanimal_mt \t \" --include-priority " + os.path.join(LBL_DB_PATH, "animal_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nfungus_mt \t \" --include-priority " + os.path.join(LBL_DB_PATH, "fungus_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nfungus_nr \t \" --include-priority " + os.path.join(LBL_DB_PATH, "fungus_nr.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nother_pt,embplant_mt,fungus_mt \t \" --include-priority " + os.path.join(LBL_DB_PATH, "other_pt.fasta") + "," + os.path.join(LBL_DB_PATH, "embplant_mt.fasta") + "," + os.path.join(LBL_DB_PATH, "fungus_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "For easy usage and compatibility of old versions, following redirection "
                           "would be automatically fulfilled without warning:\t"
                           "\nplant_cp->embplant_pt; plant_pt->embplant_pt; "
                           "\nplant_mt->embplant_mt; plant_nr->embplant_nr"
                      )
    parser.add_option("-E", dest="exclude_organelle_types",
                      help="followed with mode embplant_pt, other_pt, embplant_mt, embplant_nr, animal_mt, fungus_mt,"
                           "fungus_nr "
                           "(which means embryophyta plastome, non-embryophyta plastome, "
                           "plant mitochondrion, plant nuclear ribosomal RNA, animal mitochondrion, "
                           "fungus mitochondrion, fungus nuclear ribosomal RNA separately), "
                           "or a combination of above split by comma(s) "
                           "(be similar to -F and corresponds to certain arguments as following listed). "
                           "\t"
                           " ------------------------------------------------------ "
                           "\nembplant_pt \t \" --exclude " + os.path.join(LBL_DB_PATH, "embplant_pt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nembplant_mt \t \" --exclude " + os.path.join(LBL_DB_PATH, "embplant_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nembplant_nr \t \" --exclude " + os.path.join(LBL_DB_PATH, "embplant_nr.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nanimal_mt \t \" --exclude " + os.path.join(LBL_DB_PATH, "animal_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nfungus_mt \t \" --exclude " + os.path.join(LBL_DB_PATH, "fungus_mt.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nfungus_nr \t \" --exclude " + os.path.join(LBL_DB_PATH, "fungus_nr.fasta") + "\""
                           " ------------------------------------------------------ "
                           "\nembplant_mt,embplant_nr \t \" --exclude " + os.path.join(LBL_DB_PATH, "embplant_mt.fasta") + "," + os.path.join(LBL_DB_PATH, "embplant_nr.fasta") + "\""
                           " ------------------------------------------------------ "
                           "For easy usage and compatibility of old versions, following redirection "
                           "would be automatically fulfilled without warning:\t"
                           "\nplant_cp->embplant_pt; plant_pt->embplant_pt; "
                           "\nplant_mt->embplant_mt; plant_nr->embplant_nr")
    parser.add_option("--no-hits", dest="treat_no_hits", default="ex_no_con",
                      help="Provide treatment for non-hitting contigs.\t"
                           " ------------------------------------------------------ "
                           "\nex_no_con \t keep those connect with hitting-include contigs. (Default)"
                           " ------------------------------------------------------ "
                           "\nex_no_hit \t exclude all."
                           " ------------------------------------------------------ "
                           "\nkeep_all \t keep all"
                           " ------------------------------------------------------ ")
    parser.add_option("--max-slim-extending-len", dest="max_slim_extending_len",
                      default=MAX_SLIM_EXTENDING_LENS["anonym"],
                      type=float,
                      help="This is used to limit the extending length, below which a \"non-hit contig\" is allowed "
                           "to be distant from a \"hit contig\" to be kept. This distance is measured by the shortest "
                           "distance connecting those two contigs, weighted by the depth of the \"hit contig\". "
                           "This is used only when \"--no-hits ex_no_con\" was chosen. "
                           "Should be a single INTEGER number or inf (meaning infinite). "
                           "It is supposed to be half of the maximum expected "
                           "genome size to be safe, but could be much smaller if the LabelDatabse is closely related. "
                           "Default: " +
                           str(MAX_SLIM_EXTENDING_LENS["embplant_pt"]) + " (-F embplant_pt), " +
                           str(MAX_SLIM_EXTENDING_LENS["embplant_mt"]) + " (-F embplant_mt/fungus_mt/other_pt), " +
                           str(MAX_SLIM_EXTENDING_LENS["embplant_nr"]) + " (-F embplant_nr/fungus_nr/animal_mt), "
                           "maximum_of_type1_type2 (-F type1,type2), %default (cases without using -F)")
    parser.add_option("--significant", dest="significant", default=3.0, type=float,
                      help="Within a contig, if the query-score of hitting A is more than given times (Default: 3.0) "
                           "of the query-score of hitting B, this contig would be treated as only A related, "
                           "rather than both.")
    parser.add_option("--depth-cutoff", dest="depth_cutoff", default=10000.0, type=float,
                      help="After detection for target coverage, those beyond certain times (depth cutoff) of the"
                           " detected coverage would be excluded. Default: %default")
    parser.add_option("--min-depth", dest="min_depth", default=0., type=float,
                      help="Input a float or integer number. Filter fastg file by a minimum depth. Default: %default.")
    parser.add_option("--max-depth", dest="max_depth", default=inf, type=float,
                      help="Input a float or integer number. filter fastg file by a maximum depth. Default: %default.")
    parser.add_option("--merge", dest="merge_contigs", default=False, action="store_true",
                      help="Merge all possible contigs. ")
    parser.add_option("--include", dest="include",
                      help="followed by Blastn database(s)")
    parser.add_option("--include-priority", dest="include_priority",
                      help="followed by Blastn database(s).")
    parser.add_option("--exclude", dest="exclude",
                      help="followed by Blastn database(s).")
    parser.add_option("--exclude-priority", dest="exclude_priority",
                      help="followed by Blastn database(s)")
    parser.add_option("--no-hits-labeled-tab", dest="no_tab", default=False, action="store_true",
                      help="Choose to disable producing tab file")
    parser.add_option("--keep-temp", dest="keep_temp", default=False, action="store_true",
                      help="Choose to disable deleting temp files produced by blast and this script")
    parser.add_option("-o", "--out-dir", dest="out_dir",
                      help="By default the output would be along with the input fastg file. "
                           "But you could assign a new directory with this option.")
    parser.add_option("-e", "--evalue", dest="evalue", default=1e-25, type=float,
                      help="blastn evalue threshold. Default: %default")
    parser.add_option("--prefix", dest="prefix", default="",
                      help="Add prefix to the output basename. Conflict with \"--out-base\".")
    parser.add_option("--out-base", dest="out_base", default="",
                      help="By default the output basename would be modified based on the input fastg file. "
                           "But you could assign a new basename with this option. Conflict with \"--prefix\". "
                           "Conflict with multiple input files!")
    parser.add_option("--log", dest="generate_log", default=False, action="store_true",
                      help="Generate log file.")
    parser.add_option("--wrapper", dest="wrapper_mode", default=False, action="store_true",
                      help="Wrapper mode logging when called by get_organelle*.py. Default: %default")
    parser.add_option("--verbose", dest="verbose_log", default=False, action="store_true",
                      help="For debug usage.")
    parser.add_option("--continue", dest="resume", default=False, action="store_true",
                      help="Specified for calling from get_organelle_from_reads.py")
    parser.add_option("--no-overwrite", dest="overwrite", default=True, action="store_false",
                      help="Overwrite existing output result.")
    parser.add_option("--which-blast", dest="which_blast", default="",
                      help="Assign the path to BLAST binary files if not added to the path. "
                           "Default: try \"" + os.path.realpath("GetOrganelleDep") + "/" + SYSTEM_NAME +
                           "/ncbi-blast\" first, then $PATH")
    parser.add_option("--config-dir", dest="get_organelle_path", default=None,
                      help="The directory where the default databases were placed. "
                           "The default value also can be changed by adding 'export GETORG_PATH=your_favor' "
                           "to the shell script (e.g. ~/.bash_profile or ~/.bashrc) "
                           "Default: " + GO_PATH)
    parser.add_option("-t", "--threads", dest="threads", default=4, type=int,
                      help="Threads for blastn.")

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
        options, args = parser.parse_args()
    except OptionConflictError as e:
        sys.stdout.write('\n\n######################################'+str(e))
        sys.stdout.write('\n\n"-h" for more usage\n')
        exit()
    else:
        if not options.which_blast:
            try_this_bin = os.path.join(GO_DEP_PATH, "ncbi-blast", "blastn")
            if os.path.isfile(try_this_bin) and executable(try_this_bin):
                output, err = subprocess.Popen(
                    try_this_bin + " -version", stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, shell=True).communicate()
                if "not found" in output.decode("utf8"):
                    sys.stdout.write(output.decode("utf8") + "\n")
                else:
                    options.which_blast = os.path.split(try_this_bin)[0]
        if not executable(os.path.join(options.which_blast, "blastn")):
            sys.stdout.write(os.path.join(options.which_blast, "blastn") + " not accessible!")
            exit()
        if not executable(os.path.join(options.which_blast, "makeblastdb")):
            sys.stdout.write(os.path.join(options.which_blast, "makeblastdb") + " not accessible!")
            exit()
        if options.treat_no_hits not in ["ex_no_con", "ex_no_hit", "keep_all"]:
            sys.stdout.write('\n\nOption Error: you should choose assign one of "ex_no_con", "ex_no_hit"'
                             ' and "keep_all" to variable treat_no_hits\n')
            exit()
        priority_chosen = int(bool(options.include_priority)) + int(bool(options.exclude_priority))
        secondary_chosen = int(bool(options.include)) + int(bool(options.exclude))
        if priority_chosen + secondary_chosen > 0:
            # sys.stdout.write("\norganelle_types is disabled since you assign the customized index/indices.\n")
            if priority_chosen > 1:
                sys.stdout.write('\n\nLogical Error: only one option with "-priority" allowed!\n')
                exit()
            in_chosen = int(bool(options.include_priority)) + int(bool(options.include))
            if in_chosen > 1:
                sys.stdout.write('\n\nOption Error: you can not simultaneously choose two "--include-*" options!\n')
                exit()
            ex_chosen = int(bool(options.exclude_priority)) + int(bool(options.exclude))
            if ex_chosen > 1:
                sys.stdout.write('\n\nOption Error: you can not simultaneously choose two "--exclude-*" options!\n')
                exit()
            if in_chosen == 1 and ex_chosen == 1 and priority_chosen == 0:
                sys.stdout.write('\n\nOption Error: since you have include and exclude chosen, '
                                 'one of them should be assigned priority!\n')
                exit()
            if ex_chosen == 1 and in_chosen == 0 and (options.treat_no_hits in ["ex_no_con", "ex_no_hit"]):
                sys.stdout.write('\n\nOption Error: no contigs survive according to you choice!\n')
                exit()
            if options.include_priority:
                include_priority_str = str(options.include_priority)
                options.include_priority = []
                for sub_i_p in include_priority_str.split(","):
                    if not os.path.exists(sub_i_p):
                        sys.stdout.write("Error: " + sub_i_p + " not found!\n")
                        exit()
                    # elif not os.path.exists(remove_db_postfix(sub_i_p) + ".nhr"):
                    #     sys.stdout.write("Error: " + remove_db_postfix(sub_i_p) + ".nhr not found!\n")
                    #     exit()
                    else:
                        options.include_priority.append(sub_i_p)
            else:
                options.include_priority = []
            if options.exclude_priority:
                exclude_priority_str = str(options.exclude_priority)
                options.exclude_priority = []
                for sub_e_p in exclude_priority_str.split(","):
                    if not os.path.exists(sub_e_p):
                        sys.stdout.write("Error: " + sub_e_p + " not found!\n")
                        exit()
                    # elif not os.path.exists(remove_db_postfix(sub_e_p) + ".nhr"):
                    #     sys.stdout.write("Error: " + remove_db_postfix(sub_e_p) + ".nhr not found!\n")
                    #     exit()
                    else:
                        options.exclude_priority.append(sub_e_p)
            else:
                options.exclude_priority = []
            if options.include:
                include_str = str(options.include)
                options.include = []
                for sub_i in include_str.split(","):
                    if not os.path.exists(sub_i):
                        sys.stdout.write("Error: " + sub_i + " not found!\n")
                        exit()
                    # elif not os.path.exists(remove_db_postfix(sub_i) + ".nhr"):
                    #     sys.stdout.write("Error: " + remove_db_postfix(sub_i) + ".nhr not found!\n")
                    #     exit()
                    else:
                        options.include.append(sub_i)
            else:
                options.include = []
            if options.exclude:
                exclude_str = str(options.exclude)
                options.exclude = []
                for sub_e in exclude_str.split(","):
                    if not os.path.exists(sub_e):
                        sys.stdout.write("Error: " + sub_e + " not found!\n")
                        exit()
                    # elif not os.path.exists(remove_db_postfix(sub_e) + ".nhr"):
                    #     sys.stdout.write("Error: " + remove_db_postfix(sub_e) + ".nhr not found!\n")
                    #     exit()
                    else:
                        options.exclude.append(sub_e)
            else:
                options.exclude = []
        # elif options.organelle_types == 'embplant_pt':
        #     options.include_priority = [os.path.join(LBL_DB_PATH, 'embplant_pt.fasta')]
        #     options.exclude = [os.path.join(LBL_DB_PATH, 'embplant_mt.fasta')]
        # elif options.organelle_types == 'embplant_mt':
        #     options.include_priority = [os.path.join(LBL_DB_PATH, 'embplant_mt.fasta')]
        #     options.exclude = [os.path.join(LBL_DB_PATH, 'embplant_pt.fasta')]
        else:

            global _GO_PATH, _LBL_DB_PATH, _SEQ_DB_PATH
            if options.get_organelle_path:
                _GO_PATH = os.path.expanduser(options.get_organelle_path)
                if os.path.isdir(_GO_PATH):
                    _LBL_DB_PATH = os.path.join(_GO_PATH, LBL_NAME)
                    _SEQ_DB_PATH = os.path.join(_GO_PATH, SEQ_NAME)
                else:
                    sys.stdout.write(
                        "\n############################################################################"
                        "\nERROR: path " + _GO_PATH + " invalid!\n")

            if options.organelle_types:

                def _check_default_db(this_sub_organelle, extra_type=""):
                    if not ((os.path.isfile(os.path.join(_LBL_DB_PATH, this_sub_organelle + ".fasta")) and
                             os.path.isfile(os.path.join(_SEQ_DB_PATH, this_sub_organelle + ".fasta")))):
                        sys.stdout.write(
                            "\n############################################################################"
                            "\nERROR: default " + this_sub_organelle + "," * int(bool(extra_type)) + extra_type +
                            " database has not been added yet!\n"
                            "\nInstall it by: get_organelle_config.py -a " + this_sub_organelle +
                            "," * int(bool(extra_type)) + extra_type +
                            "\nor\nInstall all types by: get_organelle_config.py -a all\n")
                        exit()
                options.include_priority = []
                options.organelle_types = options.organelle_types.split(",")
                for sub_organelle_t in options.organelle_types:
                    if sub_organelle_t not in \
                            ("embplant_pt", "other_pt", "embplant_mt", "embplant_nr",
                             "animal_mt", "fungus_mt", "fungus_nr"):
                        sys.stdout.write(
                            "\n############################################################################"
                            "\nERROR: \"-F\" MUST be one of 'embplant_pt', 'other_pt', 'embplant_mt', 'embplant_nr', "
                            "'animal_mt', 'fungus_mt', 'fungus_nr', or a combination of above split by comma(s)!\n\n")
                        exit()
                    else:
                        if not (os.path.isfile(os.path.join(_LBL_DB_PATH, sub_organelle_t + ".fasta")) and
                                os.path.isfile(os.path.join(_SEQ_DB_PATH, sub_organelle_t + ".fasta"))):
                            if sub_organelle_t in ("embplant_pt", "embplant_mt"):
                                for go_t, check_sub in enumerate(["embplant_pt", "embplant_mt"]):
                                    _check_default_db(check_sub, ["embplant_pt", "embplant_mt"][not go_t])
                            else:
                                _check_default_db(sub_organelle_t)
                        else:
                            options.include_priority.append(os.path.join(_LBL_DB_PATH, sub_organelle_t + ".fasta"))
            else:
                options.organelle_types = []
            if options.exclude_organelle_types:
                options.exclude = []
                options.exclude_organelle_types = options.exclude_organelle_types.split(",")
                for sub_organelle_t in options.exclude_organelle_types:
                    if sub_organelle_t not in \
                            ("embplant_pt", "other_pt", "embplant_mt", "embplant_nr",
                             "animal_mt", "fungus_mt", "fungus_nr"):
                        sys.stdout.write(
                            "\n############################################################################"
                            "\nERROR: \"-E\" MUST be one of 'embplant_pt', 'other_pt', 'embplant_mt', 'embplant_nr', "
                            "'animal_mt', 'fungus_mt', 'fungus_nr', or a combination of above split by comma(s)!\n\n")
                        exit()
                    else:
                        if not (os.path.isfile(os.path.join(_LBL_DB_PATH, sub_organelle_t + ".fasta")) and
                                 os.path.isfile(os.path.join(_SEQ_DB_PATH, sub_organelle_t + ".fasta"))):
                            sys.stdout.write(
                                "\n############################################################################"
                                "\nERROR: default " + sub_organelle_t + " database has not been added yet!"
                                "\nInstall it by: get_organelle_config.py -a " + sub_organelle_t +
                                "\nor install all types by: get_organelle_config.py -a all\n")
                            exit()
                        else:
                            options.exclude.append(os.path.join(_LBL_DB_PATH, sub_organelle_t + ".fasta"))
        if "--max-slim-extending-len" not in sys.argv:
            default = options.max_slim_extending_len
            if options.organelle_types:
                options.max_slim_extending_len = 0
                for sub_organelle_t in options.organelle_types:
                    if sub_organelle_t in ("embplant_mt", "other_pt", "fungus_mt"):
                        options.max_slim_extending_len = max(MAX_SLIM_EXTENDING_LENS["embplant_mt"], options.max_slim_extending_len)
                    elif sub_organelle_t == "embplant_pt":
                        options.max_slim_extending_len = max(MAX_SLIM_EXTENDING_LENS["embplant_pt"], options.max_slim_extending_len)
                    elif sub_organelle_t in ("embplant_nr", "fungus_nr", "animal_mt"):
                        options.max_slim_extending_len = max(MAX_SLIM_EXTENDING_LENS["embplant_nr"], options.max_slim_extending_len)
                    else:
                        options.max_slim_extending_len = default
        else:
            if options.max_slim_extending_len != inf:
                options.max_slim_extending_len = int(options.max_slim_extending_len)
        if not len(args):
            sys.stdout.write('\n\nInput Error: you must choose one input fasta or fastg file!\n')
            exit()
        if options.out_dir:
            check_path_chars(options.out_dir)
        log_output_dir = options.out_dir if options.out_dir else os.path.realpath(os.path.split(args[0])[0])
        if not os.path.isdir(log_output_dir):
            os.mkdir(log_output_dir)
        assert not (options.out_base and options.prefix), "\"--out-base\" conflicts with \"--prefix\"!"
        assert not (options.out_base and len(args) > 1), "\"--out-base\" conflicts with multiple input files!"
        if options.out_base:
            # Replace illegal characters for blastn
            options.out_base = os.path.basename(options.out_base.replace("'", "_"))
        if options.prefix:
            # Replace illegal characters for blastn
            options.prefix = os.path.basename(options.prefix.replace("'", "_"))
        if options.generate_log:
            log_handler = simple_log(logging.getLogger(), log_output_dir,
                                     str((options.prefix + options.out_base + ".") * int(
                                         bool(options.out_base))) + "slim.")
            if not options.wrapper_mode:
                log_handler.info(print_title)
                log_handler.info(" ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n")
            log_handler = timed_log(log_handler, log_output_dir,
                                    str((options.prefix + options.out_base + ".") * int(
                                         bool(options.out_base))) + "slim.")
        else:
            log_handler = None
            sys.stdout.write(print_title + "\n")
            sys.stdout.write("\n" + " ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n\n")
        return options, args, log_handler


def check_db(include_priority_f, include_f, exclude_priority_f, exclude_f, new_db_dir,
             which_blast="", log_handler=None, verbose_log=False):
    databases_made_for_this_run = []
    in_index = []
    ex_index = []
    if include_priority_f:
        in_index = include_priority_f
    elif include_f:
        in_index = include_f
    if exclude_priority_f:
        ex_index = exclude_priority_f
    elif exclude_f:
        ex_index = exclude_f
    for check_blast_databases in [in_index, ex_index]:
        for go_db, check_blast_db in enumerate(check_blast_databases):
            if sum([os.path.exists(remove_db_postfix(check_blast_db) + postfix)
                    for postfix in (".nsq", ".nin", ".nhr")]) != 3:
                if log_handler:
                    log_handler.info("Making BLAST database ... ")
                else:
                    sys.stdout.write("\nMaking BLAST database ... \n")
                # make output name
                output_file_base = os.path.join(new_db_dir, os.path.basename(check_blast_db))
                output_file_base = remove_db_postfix(output_file_base)
                # make blast
                make_blast_db(input_file=check_blast_db, output_base=output_file_base, which_blast=which_blast,
                              log_handler=log_handler, verbose_log=verbose_log)
                databases_made_for_this_run.extend([output_file_base + postfix for postfix in (".nsq", ".nin", ".nhr")])
                # modifiy the in_index and ex_index
                check_blast_databases[go_db] = output_file_base
                if log_handler:
                    log_handler.info("Making BLAST database finished.")
                else:
                    sys.stdout.write("Making BLAST database finished.")
            else:
                check_blast_databases[go_db] = remove_db_postfix(check_blast_db)
    return in_index, ex_index, databases_made_for_this_run


def make_new_matrix_with_names(names, old_matrix):
    i = 0
    while i < len(old_matrix[0]):
        if old_matrix[0][i] in names:
            i += 1
        else:
            del old_matrix[0][i]
            del old_matrix[1][i]
    return old_matrix


def blast_and_call_names(fasta_file, index_files, out_file, threads, e_value=1e-25, which_blast="", log_handler=None):
    names = {}
    if index_files:
        time0 = time.time()
        for index_f in index_files:
            if log_handler:
                log_handler.info('Executing BLAST to ' + index_f + ' ...')
            else:
                sys.stdout.write('\nExecuting BLAST to ' + index_f + '...')
            execute_blast(query=fasta_file, blast_db=index_f, output=out_file, threads=threads,
                          outfmt=6, e_value=e_value, which_blast=which_blast, log_handler=log_handler)
            time1 = time.time()
            if log_handler:
                log_handler.info('Executing BLAST to ' + index_f + ' finished.')
            else:
                sys.stdout.write('\nExecuting BLAST to '+os.path.split(index_f)[-1]+
                                 ' cost: '+str(round(time1-time0, 2)) + "\n")
            try:
                blast_out_lines = open(out_file, 'r')
            except IOError:
                if log_handler:
                    log_handler.error(os.path.join(which_blast, "blastn") + " was not properly configured.")
                else:
                    sys.stdout.write("\n" + os.path.join(which_blast, "blastn") + " was not properly configured.")
                raise EnvironmentError
            this_database = os.path.split(index_f)[-1].split(" ")[0]
            for line in blast_out_lines:
                line_split = line.strip().split('\t')
                query, hit = line_split[0], line_split[1]
                q_start, q_end, q_score = int(line_split[6]), int(line_split[7]), float(line_split[2])
                q_min, q_max = min(q_start, q_end), max(q_start, q_end)
                # q_score = abs(q_max - q_min + 1)*q_score
                if query in names:
                    if this_database not in names[query]:
                        names[query][this_database] = {}
                    if hit not in names[query][this_database]:
                        names[query][this_database][hit] = [(q_min, q_max, q_score)]
                    else:
                        # if overlap, then merge
                        i = 0
                        while i < len(names[query][this_database][hit]):
                            this_min, this_max, previous_q_score = names[query][this_database][hit][i]
                            if q_max < this_min:
                                break
                            elif q_min > this_max:
                                i += 1
                                continue
                            else:
                                q_min = min(q_min, this_min)
                                q_max = max(q_max, this_max)
                                q_score = max(previous_q_score, q_score)
                                del names[query][this_database][hit][i]
                        names[query][this_database][hit].insert(i, (q_min, q_max, q_score))
                else:
                    names[query] = {this_database: {}}
                    names[query][this_database][hit] = [(q_min, q_max, q_score)]
            if log_handler:
                log_handler.info("Parsing blast result finished.")
            else:
                sys.stdout.write('Parsing blast result cost: '+str(round(time.time()-time1, 2)) + "\n")
    return names


def summarize_q_score(info_of_a_query):
    all_loci = []
    for this_db in info_of_a_query:
        for hits in info_of_a_query[this_db].values():
            for locus in hits:
                all_loci.append(list(locus))
    all_loci.sort(key=lambda x: (x[0], x[1]))
    i = 0
    while i < len(all_loci) - 1:
        if all_loci[i][1] >= all_loci[i+1][0]:
            if all_loci[i][2] > all_loci[i+1][2]:
                all_loci[i+1][0] = all_loci[i][1] + 1
                if all_loci[i+1][0] > all_loci[i+1][1]:
                    del all_loci[i+1]
                else:
                    i += 1
            else:
                all_loci[i][1] = all_loci[i+1][0] - 1
                if all_loci[i][0] > all_loci[i][1]:
                    del all_loci[i]
                else:
                    i += 1
        else:
            i += 1
    return sum([this_locus[2]*(this_locus[1]-this_locus[0]+1) for this_locus in all_loci])


def modify_in_ex_according_to_depth(in_names, ex_names, significant, assembly_graph, depth_cutoff, log_handler=None):
    here_in_names, here_ex_names = copy.deepcopy(in_names), copy.deepcopy(ex_names)
    if assembly_graph:
        training_in = set()
        for query_name in list(in_names):
            if query_name not in ex_names:
                training_in.add(query_name)
        training_ex = set()
        for query_name in list(ex_names):
            if query_name not in in_names:
                training_ex.add(query_name)

        for query_name in list(in_names):
            if query_name in ex_names:
                in_name_score = summarize_q_score(in_names[query_name])
                ex_name_score = summarize_q_score(ex_names[query_name])
                if in_name_score / float(ex_name_score) > significant:
                    training_in.add(query_name)
                elif ex_name_score / float(in_name_score) > significant:
                    training_ex.add(query_name)

        def get_average_coverage(here_q_list, here_info):
            total_in_base = 0
            total_in_len = 0
            for here_q_name in here_q_list:
                this_coverage = assembly_graph.vertex_info[here_q_name].cov
                for here_db in here_info[here_q_name]:
                    for hits in here_info[here_q_name][here_db].values():
                        for q_min, q_max, q_score in hits:
                            this_len = abs(q_max - q_min) + 1
                            total_in_len += this_len
                            total_in_base += this_len * this_coverage
            return total_in_base / float(total_in_len)

        def combine_coverage_to_check():
            aver_in_coverage = get_average_coverage(training_in, in_names)
            aver_ex_coverage = get_average_coverage(training_ex, ex_names)
            if log_handler:
                log_handler.info("average in coverage: " + str(aver_in_coverage))
                log_handler.info("average ex coverage: " + str(aver_ex_coverage))
            else:
                sys.stdout.write("\naverage in coverage: " + str(aver_in_coverage))
                sys.stdout.write("\naverage ex coverage: " + str(aver_ex_coverage))
            for q_name in list(here_in_names):
                if q_name in here_ex_names:
                    in_score = summarize_q_score(here_in_names[q_name]) / \
                               (1 + abs(math.log(assembly_graph.vertex_info[q_name].cov / aver_in_coverage)))
                    ex_score = summarize_q_score(here_ex_names[q_name]) / \
                               (1 + abs(math.log(assembly_graph.vertex_info[q_name].cov / aver_ex_coverage)))
                    if in_score / float(ex_score) > significant:
                        del here_ex_names[q_name]
                        training_in.add(q_name)
                        training_ex.discard(q_name)
                        # sys.stdout.write('\n' + ' ' * 4 + q_name + ' included: '
                        #                  + str(round(in_name_score, 2)) + '~' + str(round(ex_name_score, 2)))
                    elif ex_score / float(in_score) > significant:
                        del here_in_names[q_name]
                        training_in.discard(q_name)
                        training_ex.add(q_name)
                        if 2 ** abs(math.log(assembly_graph.vertex_info[q_name].cov / aver_in_coverage, 2)) <= 5:
                            if log_handler:
                                log_handler.info(" " * 4 + q_name + " excluded: "
                                                 + str(round(in_score, 2)) + '~' + str(round(ex_score, 2)))
                            else:
                                sys.stdout.write('\n' + ' ' * 4 + q_name + ' excluded: '
                                                 + str(round(in_score, 2)) + '~' + str(round(ex_score, 2)) + "\n")
                    else:
                        pass
                        # sys.stdout.write('\n' + ' ' * 4 + query_name + ' with both hits: '
                        #                  + str(round(in_name_score, 2)) + '~' + str(round(ex_name_score, 2)))
            aver_in_coverage = get_average_coverage(training_in, in_names)
            aver_ex_coverage = get_average_coverage(training_ex, ex_names)
            if log_handler:
                log_handler.info("average in coverage: " + str(aver_in_coverage))
                log_handler.info("average ex coverage: " + str(aver_ex_coverage))
            else:
                sys.stdout.write("\naverage in coverage: " + str(aver_in_coverage))
                sys.stdout.write("\naverage ex coverage: " + str(aver_ex_coverage))
            for q_name in list(here_in_names):
                if 2 ** abs(math.log(assembly_graph.vertex_info[q_name].cov / aver_in_coverage, 2)) >= depth_cutoff:
                    del here_in_names[q_name]
            return here_in_names, here_ex_names, aver_in_coverage

        if training_in and training_ex:
            combine_coverage_to_check()
        else:
            for query_name in list(here_in_names):
                if query_name in here_ex_names:
                    in_name_score = summarize_q_score(here_in_names[query_name])
                    ex_name_score = summarize_q_score(here_ex_names[query_name])
                    if in_name_score / float(ex_name_score) > significant:
                        training_in.add(query_name)
                        training_ex.discard(query_name)
                    elif ex_name_score / float(in_name_score) > significant:
                        training_in.discard(query_name)
                        training_ex.add(query_name)
            if training_in and training_ex:
                combine_coverage_to_check()
            else:
                if log_handler:
                    log_handler.info("No enough coverage information found.")
                else:
                    sys.stdout.write("No enough coverage information found.\n")
    else:
        if log_handler:
            log_handler.info("No coverage information found.")
        else:
            sys.stdout.write("No coverage information found.\n")
    for query_name in list(here_in_names):
        if query_name in here_ex_names:
            in_name_score = summarize_q_score(here_in_names[query_name])
            ex_name_score = summarize_q_score(here_ex_names[query_name])
            if in_name_score/float(ex_name_score) > significant:
                del here_ex_names[query_name]
                # sys.stdout.write('\n'+' ' * 4 + query_name + ' included: '
                #                  + str(round(in_name_score, 2))+'~'+str(round(ex_name_score, 2)))
            elif ex_name_score/float(in_name_score) > significant:
                del here_in_names[query_name]
                # sys.stdout.write('\n'+' '*4+query_name+' excluded: '
                #                  + str(round(in_name_score, 2))+'~'+str(round(ex_name_score, 2)))
            else:
                pass
                # sys.stdout.write('\n'+' '*4+query_name+' with both hits: '
                #                  + str(round(in_name_score, 2))+'~'+str(round(ex_name_score, 2)))
    return here_in_names, here_ex_names, None


def generate_baits_offsets(in_names, databases, assembly_graph):
    vertex_trimming = {}
    # structure: names[query][this_database][label] = [(q_min, q_max, q_score)]
    for vertex_name in in_names:
        these_hit_list = []
        for include_db_file in databases:
            include_db_file = os.path.split(include_db_file)[-1].split(" ")[0]
            if include_db_file in in_names[vertex_name]:
                for hit_info_list in in_names[vertex_name][include_db_file].values():
                    these_hit_list.extend(hit_info_list)
        all_loc_values = [x[0] for x in these_hit_list] + [x[1] for x in these_hit_list]
        min_loc = min(all_loc_values)
        max_loc = min(all_loc_values)
        vertex_trimming[(vertex_name, False)] = min_loc - 1
        vertex_trimming[(vertex_name, True)] = assembly_graph.vertex_info[vertex_name].len - max_loc
    return vertex_trimming


def reduce_matrix(in_names, ex_names, seq_matrix, assembly_graph, max_slim_extending_len, bait_offsets,
                  aver_in_dep, depth_cutoff, treat_no_hits,
                  include_priority_assigned, include_assigned, exclude_priority_assigned, exclude_assigned,
                  log_handler=None):
    if log_handler:
        log_handler.info("Mapping names ...")
    time0 = time.time()
    # candidate_short_to_2fulls = {}
    if assembly_graph:
        if aver_in_dep:
            rm_contigs = [this_v.v_name
                          for this_v in assembly_graph
                          if 2 ** abs(math.log(this_v.cov / aver_in_dep, 2)) < depth_cutoff]
            assembly_graph.remove_vertex(rm_contigs)
        if exclude_priority_assigned:
            assembly_graph.remove_vertex(ex_names)
            if treat_no_hits == "ex_no_con":
                assembly_graph.reduce_to_subgraph(bait_vertices=in_names,
                                                  bait_offsets=bait_offsets,
                                                  limit_extending_len=max_slim_extending_len,
                                                  extending_len_weighted_by_depth=True)
            elif treat_no_hits == "ex_no_hit":
                assembly_graph.remove_vertex([rm_c for rm_c in assembly_graph.vertex_info if rm_c not in in_names])
            else:
                pass
        elif include_priority_assigned or include_assigned:
            if treat_no_hits == "ex_no_con":
                assembly_graph.remove_vertex([rm_c for rm_c in ex_names if rm_c not in in_names])
                assembly_graph.reduce_to_subgraph(bait_vertices=in_names,
                                                  bait_offsets=bait_offsets,
                                                  limit_extending_len=max_slim_extending_len,
                                                  extending_len_weighted_by_depth=True)
            elif treat_no_hits == "ex_no_hit":
                assembly_graph.remove_vertex([rm_c for rm_c in assembly_graph.vertex_info if rm_c not in in_names])
            else:
                pass
        elif exclude_assigned:
            assembly_graph.remove_vertex(ex_names)
        else:
            pass
    else:
        # accepted = set()
        if exclude_priority_assigned:
            seq_matrix.remove(ex_names)
            if treat_no_hits in ("ex_no_con", "ex_no_hit"):
                seq_matrix.remove([rm_c.label for rm_c in seq_matrix if rm_c.label not in in_names])
            else:
                pass
        elif include_priority_assigned or include_assigned:
            if treat_no_hits in ["ex_no_con", "ex_no_hit"]:
                seq_matrix.remove([rm_c.label for rm_c in seq_matrix if rm_c.label not in in_names])
            else:
                seq_matrix.remove(ex_names)
        elif exclude_assigned:
            seq_matrix.remove(ex_names)
        else:
            pass
    if log_handler:
        log_handler.info("Mapping names finished.")
    else:
        sys.stdout.write('\nMapping names cost: '+str(round(time.time()-time0, 2)) + "\n")
    return assembly_graph, seq_matrix


def write_hits_tab_for_bandage(in_names, ex_names, out_file, overwrite, assembly_graph, log_handler=None):
    if not overwrite:
        while os.path.exists(out_file+'.csv'):
            out_file += '_'
    out_file += '.csv'
    out_lines = []
    # if include_file:
    #     in_database = os.path.split(include_file)[-1]
    # if exclude_file:
    #     ex_database = os.path.split(exclude_file)[-1]
    edges_set = set(in_names)
    for edge in ex_names:
        edges_set.add(edge)
    edges = list(edges_set)
    for edge in edges:
        this_string = edge
        if edge in in_names and edge in ex_names:
            databases = sorted(set(in_names[edge]) | set(ex_names[edge]))
            this_string += '\t' + ";".join(databases)
            this_db_hits = {db_name: set() for db_name in databases}
            for this_db in in_names[edge]:
                for hit_n in in_names[edge][this_db]:
                    this_db_hits[this_db].add(hit_n)
            for this_db in ex_names[edge]:
                for hit_n in ex_names[edge][this_db]:
                    this_db_hits[this_db].add(hit_n)
            this_string += "\t" + ";".join([",".join(sorted(this_db_hits[in_db])) for in_db in databases])
            # this_string += '\t' + ','.join(
            #     [hit_n for this_db in in_names[edge] for hit_n in in_names[edge][this_db]]) + ';' + ','.join(
            #     [hit_n for this_db in ex_names[edge] for hit_n in ex_names[edge][this_db]])
            loci = []
            for in_database in in_names[edge]:
                for locus in in_names[edge][in_database]:
                    for region in in_names[edge][in_database][locus]:
                        loci.append([region[0], region[1], locus, in_database])
            for ex_database in ex_names[edge]:
                for locus in ex_names[edge][ex_database]:
                    for region in ex_names[edge][ex_database][locus]:
                        loci.append([region[0], region[1], locus, ex_database])
        elif edge in in_names:
            databases = sorted(set(in_names[edge]))
            this_string += '\t' + ";".join(databases)
            this_db_hits = {db_name: set() for db_name in databases}
            for this_db in in_names[edge]:
                for hit_n in in_names[edge][this_db]:
                    this_db_hits[this_db].add(hit_n)
            this_string += "\t" + ";".join([",".join(sorted(this_db_hits[in_db])) for in_db in databases])
            # this_string += '\t' + ','.join(
            #     [hit_n for this_db in in_names[edge] for hit_n in in_names[edge][this_db]])
            loci = []
            for in_database in in_names[edge]:
                for locus in in_names[edge][in_database]:
                    for region in in_names[edge][in_database][locus]:
                        loci.append([region[0], region[1], locus, in_database])
        else:
            databases = sorted(set(ex_names[edge]))
            this_string += '\t' + ";".join(databases)
            this_db_hits = {db_name: set() for db_name in databases}
            for this_db in ex_names[edge]:
                for hit_n in ex_names[edge][this_db]:
                    this_db_hits[this_db].add(hit_n)
            this_string += "\t" + ";".join([",".join(sorted(this_db_hits[in_db])) for in_db in databases])
            # this_string += '\t' + ','.join(
            #     [hit_n for this_db in ex_names[edge] for hit_n in ex_names[edge][this_db]])
            loci = []
            for ex_database in ex_names[edge]:
                for locus in ex_names[edge][ex_database]:
                    for region in ex_names[edge][ex_database][locus]:
                        loci.append([region[0], region[1], locus, ex_database])
        loci.sort(key=lambda x: x[0])
        i = 1
        while i < len(loci):
            if loci[i-1][2] == loci[i][2]:
                loci[i-1][1] = loci[i][1]
                del loci[i]
            else:
                i += 1
        postfix = ''
        if assembly_graph:
            if edge in assembly_graph.vertex_info:
                next_edges = [next_v
                              for next_v, next_e in assembly_graph.vertex_info[edge].connections[True]
                              if next_v in assembly_graph.vertex_info]
                if next_edges:
                    next_edges.sort(
                        key=lambda x: (x not in in_names, x not in ex_names, -assembly_graph.vertex_info[x].cov))
                    postfix = ">>" + next_edges[0]
        this_string += '\t'+'>>'.join([x[2] for x in loci if x[2] != 'noncoding'])+postfix
        this_string += '\t'+'>>'.join([x[2] for x in loci])+postfix
        this_string += '\t'+'>>'.join([x[2]+'('+str(x[0])+'-'+str(x[1])+','+x[3]+')' for x in loci])+postfix
        this_string += '\n'
        out_lines.append(this_string)
    out_lines.sort()
    out_lines = ['EDGE\tdatabase\tloci\tloci_gene_sequential\tloci_sequential\tdetails\n'] + out_lines
    open(out_file, 'w').writelines(out_lines)
    if log_handler:
        log_handler.info("Generating notation table to " + out_file)
    else:
        sys.stdout.write("Generating notation table to " + out_file + "\n")


def remove_temp_files(fastg_base, keep_temp, other_files_to_remove):
    if not keep_temp:
        try:
            os.remove(fastg_base + '.blast_in')
        except OSError:
            pass
        try:
            os.remove(fastg_base + '.blast_ex')
        except OSError:
            pass
        for rm_f in other_files_to_remove:
            try:
                os.remove(rm_f)
            except OSError:
                pass


def main():
    time0 = time.time()
    print_title = "GetOrganelle v" + str(get_versions()) + \
                  "\n" \
                  "\nThis is a script for excluding certain contigs " \
                  "from assembly graph file (*.fastg/*.fasta) by blast\n" \
                  "By jinjianjun@mail.kib.ac.cn\n"
    options, args, log_handler = get_options(print_title)
    from GetOrganelleLib.assembly_parser import Assembly
    from GetOrganelleLib.seq_parser import SequenceList
    log_output_dir = options.out_dir if options.out_dir else os.path.split(args[0])[0]
    log_output_name = str((options.prefix + options.out_base + ".") * int(bool(options.out_base))) + "slim."
    try:
        if not log_handler:
            sys.stdout.write('\n' + '=' * 100)
        if options.out_dir and not os.path.exists(options.out_dir):
            os.mkdir(options.out_dir)
        # make blast database if not made
        other_files_to_remove = []
        include_indices, exclude_indices, database_made = \
            check_db(options.include_priority, options.include,
                     options.exclude_priority, options.exclude,
                     new_db_dir=options.out_dir if options.out_dir else os.path.split(args[0])[0],
                     which_blast=options.which_blast, log_handler=log_handler, verbose_log=options.verbose_log)
        other_files_to_remove.extend(database_made)

        for i in range(len(args)):
            fas_file = args[i]
            if log_handler:
                log_handler.info('Slimming file '+str(i+1)+'/'+str(len(args))+': '+fas_file)
            else:
                sys.stdout.write('\nSlimming file '+str(i+1)+'/'+str(len(args))+': '+fas_file+'\n')
            #######################################
            # 0. initialization: define variables #
            #######################################
            is_fastg = fas_file.endswith('.fastg')
            is_gfa = fas_file.endswith('.gfa')
            is_graph = is_fastg or is_gfa
            output_dir = options.out_dir if options.out_dir else os.path.split(os.path.realpath(fas_file))[0]
            check_path_chars(output_dir)
            if options.out_base:
                # ignore other modifications
                constrained_full_out_base = os.path.join(output_dir, options.out_base)
                full_out_fasta = constrained_full_out_base
            else:
                constrained_full_out_base = ""
                # Replace illegal characters for blastn
                full_out_fasta = os.path.join(output_dir,
                                              options.prefix + os.path.basename(fas_file.replace("'", "_")))
            check_path_chars(full_out_fasta)
            blast_fas = full_out_fasta + ".Temp"
            if constrained_full_out_base:
                out_fas = constrained_full_out_base + '.' + str(fas_file).split('.')[-1]
                out_csv = constrained_full_out_base
            else:  # if there's no out_base requirements, add modifications to the name
                in_ex_info = generate_in_ex_info_name(include_indices=include_indices,
                                                      exclude_indices=exclude_indices,
                                                      exclude_no_con=options.treat_no_hits == 'ex_no_con',
                                                      exclude_no_hit=options.treat_no_hits == 'ex_no_hit')
                out_fas = full_out_fasta + in_ex_info + '.' + str(fas_file).split('.')[-1]
                out_csv = full_out_fasta + in_ex_info
            #######################
            # 1. parsing sequence #
            #######################
            time_r0 = time.time()
            this_assembly = Assembly()
            this_matrix = None
            if is_graph:
                this_assembly = Assembly(graph_file=fas_file, min_cov=options.min_depth, max_cov=options.max_depth)
                # merge contigs
                if options.merge_contigs:
                    this_assembly.merge_all_possible_vertices()
            else:
                this_matrix = SequenceList(fas_file)
            if log_handler:
                log_handler.info("Parsing input finished.")
            else:
                sys.stdout.write("\nParsing input cost: " + str(round(time.time() - time_r0, 2)) + "\n")
            ###################################
            # 2. prepare fasta file for blast #
            ###################################
            if bool(include_indices) or bool(exclude_indices):
                t_dc = time.time()
                if is_graph:
                    # del complementary to avoid using both the forward and reverse sequence for blast
                    this_assembly.write_to_fasta(blast_fas, interleaved=70, check_postfix=False)
                else:
                    # remove space or other characters from the seq names
                    changed_name_dict = {}
                    for seq_record in this_matrix:
                        new_name = seq_record.label.replace(" ", "_").replace("&", "_").replace("#", "_")
                        if new_name in changed_name_dict:
                            raise Exception("Conflict names:\n" + changed_name_dict[new_name] + " -> " + new_name +
                                            "\n" + seq_record.label + " -> " + new_name)
                        else:
                            changed_name_dict[new_name] = seq_record.label
                        seq_record.label = new_name
                    this_matrix.write_fasta(blast_fas, interleaved=70)
                other_files_to_remove.append(blast_fas)
                if log_handler:
                    log_handler.info("Preparing fasta file finished.")
                else:
                    sys.stdout.write('\nPreparing fasta file cost: %.2f\n' % (time.time() - t_dc))
            ############
            # 3. BLAST #
            ############
            try:
                # structure: names[query][this_database][label] = [(q_min, q_max, q_score)]
                in_names = blast_and_call_names(fasta_file=blast_fas, index_files=include_indices,
                                                out_file=blast_fas+'.blast_in', threads=options.threads,
                                                e_value=options.evalue,
                                                which_blast=options.which_blast, log_handler=log_handler)
                ex_names = blast_and_call_names(fasta_file=blast_fas, index_files=exclude_indices,
                                                out_file=blast_fas+'.blast_ex', threads=options.threads,
                                                e_value=options.evalue,
                                                which_blast=options.which_blast, log_handler=log_handler)
                if bool(include_indices) or bool(exclude_indices):
                    in_names_r, ex_names_r, aver_dep = modify_in_ex_according_to_depth(
                        in_names=in_names, ex_names=ex_names, significant=options.significant,
                        assembly_graph=this_assembly, depth_cutoff=options.depth_cutoff, log_handler=log_handler)
                else:
                    in_names_r, ex_names_r, aver_dep = in_names, ex_names, None
                # prepare bait_offsets: trim unlabeled terminal regions from bait vertices, for more accurate
                # control of "maximum slimming extending length"
                if this_assembly and options.treat_no_hits == "ex_no_con" and \
                        options.max_slim_extending_len not in (None, inf):
                    bait_offsets = generate_baits_offsets(
                        in_names=in_names, databases=include_indices, assembly_graph=this_assembly)
                else:
                    bait_offsets = {}
                new_assembly, new_matrix = \
                    reduce_matrix(in_names=in_names_r, ex_names=ex_names_r, seq_matrix=this_matrix,
                                  assembly_graph=this_assembly, max_slim_extending_len=options.max_slim_extending_len,
                                  bait_offsets=bait_offsets, aver_in_dep=aver_dep,
                                  depth_cutoff=options.depth_cutoff, treat_no_hits=options.treat_no_hits,
                                  include_priority_assigned=options.include_priority,
                                  include_assigned=options.include,
                                  exclude_priority_assigned=options.exclude_priority,
                                  exclude_assigned=options.exclude, log_handler=log_handler)
                if log_handler:
                    log_handler.info("Generating slimmed file to " + out_fas)
                else:
                    sys.stdout.write("Generating slimmed file to " + out_fas + "\n")
                if new_assembly:
                    if is_fastg:
                        new_assembly.write_to_fastg(out_fas, check_postfix=False)
                    else:
                        new_assembly.write_to_gfa(out_fas, check_postfix=False)
                else:
                    new_matrix.write_fasta(out_fas, overwrite=options.overwrite)
                # write out hits tab according to blast
                if not options.no_tab:
                    write_hits_tab_for_bandage(
                        in_names=in_names, ex_names=ex_names, out_file=out_csv, overwrite=options.overwrite,
                        assembly_graph=this_assembly, log_handler=log_handler)
            except Exception as e:
                if log_handler:
                    log_handler.error("\n" + str(e).strip())
                    log_handler.error(
                        'Slimming file ' + str(i + 1) + '/' + str(len(args)) + ': ' + args[i] + ' failed!\n')
                    if options.verbose_log:
                        raise e

                else:
                    sys.stdout.write("\n" + str(e).strip() + "\n")
                    sys.stdout.write(
                        '\nSlimming file ' + str(i + 1) + '/' + str(len(args)) + ': ' + args[i] + ' failed!'
                        '\n' + '=' * 100 + "\n")
                    if options.verbose_log:
                        raise e
            else:
                if log_handler:
                    log_handler.info('Slimming file ' + str(i+1) + '/' + str(len(args)) + ': ' + args[i] + ' finished!\n')
                else:
                    sys.stdout.write('\nSlimming file ' + str(i+1) + '/' + str(len(args)) + ': ' + args[i] + ' finished!'
                                     '\n' + '=' * 100 + "\n")
            #######################
            # 4. clean temp files #
            #######################
            if i == len(args) - 1:
                remove_temp_files(
                    fastg_base=blast_fas, keep_temp=options.keep_temp, other_files_to_remove=other_files_to_remove)
            else:
                remove_temp_files(
                    fastg_base=blast_fas, keep_temp=options.keep_temp, other_files_to_remove=[])
        if not options.wrapper_mode:
            if log_handler:
                log_handler = simple_log(log_handler, log_output_dir, log_output_name)
                log_handler.info("\nTotal cost " + "%.2f" % (time.time() - time0) + " s")
                log_handler.info("Thank you!")
            else:
                sys.stdout.write('\nTotal cost: ' + str(round(time.time() - time0, 2)) + '\nThank you!\n')
    except Exception as e:
        if log_handler:
            log_handler.exception("")
            if not options.wrapper_mode:
                log_handler = simple_log(log_handler, log_output_dir, log_output_name)
                log_handler.info("\nTotal cost " + "%.2f" % (time.time() - time0) + " s")
                log_handler.info("Thank you!")
        else:
            raise e
    if log_handler:
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