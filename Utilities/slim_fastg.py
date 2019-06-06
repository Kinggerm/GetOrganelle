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
from optparse import OptionParser
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(PATH_OF_THIS_SCRIPT, ".."))
import GetOrganelleLib
from GetOrganelleLib.versions import get_versions
from GetOrganelleLib.assembly_parser import filter_fastg_by_depth_simple, get_graph_coverage_dict_simple
from GetOrganelleLib.pipe_control_func import *
from GetOrganelleLib.seq_parser import *
import optparse
import copy
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


def get_options(print_title):
    usage = 'python '+str(os.path.basename(__file__)+' your_fastg_files -F embplant_pt')
    parser = OptionParser(usage=usage)
    # parser.add_option('-o', dest='out_fastg_file', help='Output file')
    # filters
    parser.add_option('-F', dest='organelle_types',
                      help='followed with mode embplant_pt, other_pt, embplant_mt, embplant_nr, animal_mt, fungus_mt '
                           '(which means embryophyta plastome, non-embryophyta plastome, '
                           'plant mitochondrion, plant nrDNA, animal mitochondrion, fungus mitochondrion separately), '
                           'or a combination of above split by comma(s) '
                           '(corresponding to certain arguments as following listed). '
                           '\t'
                           ' ------------------------------------------------------ '
                           '\nembplant_pt \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'embplant_pt.fasta') + ' --exclude ' + os.path.join(NOT_DB_PATH, 'embplant_mt.fasta') + '"'
                           ' ------------------------------------------------------ '
                           '\nother_pt \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'other_pt.fasta') + '"'
                           ' ------------------------------------------------------ '                                                                                                                                        
                           '\nembplant_mt \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'embplant_mt.fasta') + ' --exclude ' + os.path.join(NOT_DB_PATH, 'embplant_pt.fasta') + '"'
                           ' ------------------------------------------------------ '
                           '\nembplant_nr \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'embplant_nr.fasta') + '"'
                           ' ------------------------------------------------------ '
                           '\nanimal_mt \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'animal_mt.fasta') + '"'
                           ' ------------------------------------------------------ '
                           '\nfungus_mt \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'fungus_mt.fasta') + '"'
                           ' ------------------------------------------------------ '
                           '\nother_pt,embplant_mt,fungus_mt \t " --include-priority ' + os.path.join(NOT_DB_PATH, 'other_pt.fasta') + ',' + os.path.join(NOT_DB_PATH, 'embplant_mt.fasta') + ',' + os.path.join(NOT_DB_PATH, 'fungus_mt.fasta') +
                           ' ------------------------------------------------------ '
                           "For easy usage and compatibility of old versions, following redirection "
                           "would be automatically fulfilled without warning:\t"
                           "\nplant_cp->embplant_pt; plant_pt->embplant_pt; "
                           "\nplant_mt->embplant_mt; plant_nr->embplant_nr"
                      )
    parser.add_option('--no-hits', dest='treat_no_hits', default='ex_no_con',
                      help='Provide treatment for non-hitting contigs.\t'
                           ' ------------------------------------------------------ '
                           '\nex_no_con \t keep those connect with hitting-include contigs. (Default)'
                           ' ------------------------------------------------------ '
                           '\nex_no_hit \t exclude all.'
                           ' ------------------------------------------------------ '
                           '\nkeep_all \t keep all'
                           ' ------------------------------------------------------ ')
    parser.add_option('--significant', dest='significant', default=3.0, type=float,
                      help='Within a contig, if the query-score of hitting A is more than given times (Default: 3.0) '
                           'of the query-score of hitting B, this contig would be treated as only A related, '
                           'rather than both.')
    parser.add_option('--remove-low-dup', dest='remove_low_duplicated', default=False, action='store_true',
                      help='Remove redundant low-coverage contigs that largely overlap some high-coverage contigs.')
    parser.add_option('--depth-cutoff', dest='depth_cutoff', default=10000.0, type=float,
                      help='After detection for target coverage, those beyond certain times (depth cutoff) of the'
                           ' detected coverage would be excluded. Default: 10000.0')
    parser.add_option('--min-depth', dest='min_depth', default=0., type=float,
                      help='Input a float or integer number. Filter fastg file by a minimum depth. Default: %default.')
    parser.add_option('--max-depth', dest='max_depth', default=inf, type=float,
                      help='Input a float or integer number. filter fastg file by a maximum depth. Default: %default.')
    parser.add_option('--no-merge', dest='merge_contigs', default=True, action="store_false",
                      help="Merge all possible contigs. ")
    parser.add_option('--include', dest='include',
                      help='followed by Blastn database(s)')
    parser.add_option('--include-priority', dest='include_priority',
                      help='followed by Blastn database(s).')
    parser.add_option('--exclude', dest='exclude',
                      help='followed by Blastn database(s).')
    parser.add_option('--exclude-priority', dest='exclude_priority',
                      help='followed by Blastn database(s)')
    parser.add_option('--no-hits-labeled-tab', dest='no_tab', default=False, action='store_true',
                      help='Choose to disable producing tab file')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true',
                      help='Choose to disable deleting temp files produced by blast and this script')
    parser.add_option('-o', '--out-dir', dest="out_dir",
                      help="By default the output would be along with the input fastg file. "
                           "But you could assign a new directory with this option.")
    parser.add_option("--prefix", dest="prefix", default="",
                      help="Add prefix to the output basename. Conflict with '--out-base'.")
    parser.add_option("--out-base", dest="out_base", default="",
                      help="By default the output basename would be modified based on the input fastg file. "
                           "But you could assign a new basename with this option. Conflict with '--prefix'.")
    parser.add_option("--log", dest="generate_log", default=False, action="store_true",
                      help="Generate log file.")
    parser.add_option("--verbose", dest="verbose_log", default=False, action="store_true",
                      help="For debug usage.")
    parser.add_option('--continue', dest='resume', default=False, action='store_true',
                      help='Specified for calling from get_organelle_from_reads.py')
    parser.add_option("--no-overwrite", dest="overwrite", default=True, action="store_false",
                      help="Overwrite existing output result.")
    parser.add_option("--which-blast", dest="which_blast", default="",
                      help="Assign the path to BLAST binary files if not added to the path. "
                           "Default: try GetOrganelleDep/" + SYSTEM_NAME + "/ncbi-blast first, then $PATH")
    parser.add_option('-t', '--threads', dest="threads", default=4, type=int,
                      help="Threads for blastn.")

    # redirect organelle types before parsing arguments
    redirect_organelle_types = {"plant_cp": "embplant_pt",
                                "plant_pt": "embplant_pt",
                                "plant_mt": "embplant_mt",
                                "plant_nr": "embplant_nr"}
    for go_arg, candidate_arg in enumerate(sys.argv):
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
    except optparse.OptionConflictError as e:
        sys.stdout.write('\n\n######################################'+str(e))
        sys.stdout.write('\n\n"-h" for more usage\n')
        exit()
    else:
        if not options.which_blast:
            try_this_bin = os.path.join(GO_DEP_PATH, "ncbi-blast", "blastn")
            if os.path.isfile(try_this_bin) and executable(try_this_bin):
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
            sys.stdout.write("\norganelle_types is disabled since you assign the customized index/indices.\n")
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
                    if not (os.path.exists(sub_i_p) or os.path.exists(remove_db_postfix(sub_i_p) + ".nhr")):
                        sys.stdout.write("Error: " + sub_i_p + " not found!\n")
                        exit()
                    else:
                        options.include_priority.append(sub_i_p)
            else:
                options.include_priority = []
            if options.exclude_priority:
                exclude_priority_str = str(options.exclude_priority)
                options.exclude_priority = []
                for sub_e_p in exclude_priority_str.split(","):
                    if not (os.path.exists(sub_e_p) or os.path.exists(remove_db_postfix(sub_e_p) + ".nhr")):
                        sys.stdout.write("Error: " + sub_e_p + " not found!\n")
                        exit()
                    else:
                        options.exclude_priority.append(sub_e_p)
            else:
                options.exclude_priority = []
            if options.include:
                include_str = str(options.include)
                options.include = []
                for sub_i in include_str.split(","):
                    if not (os.path.exists(sub_i) or os.path.exists(remove_db_postfix(sub_i) + ".nhr")):
                        sys.stdout.write("Error: " + sub_i + " not found!\n")
                        exit()
                    else:
                        options.include.append(sub_i)
            else:
                options.include = []
            if options.exclude:
                exclude_str = str(options.exclude)
                options.exclude = []
                for sub_e in exclude_str.split(","):
                    if not (os.path.exists(sub_e) or os.path.exists(remove_db_postfix(sub_e) + ".nhr")):
                        sys.stdout.write("Error: " + sub_e + " not found!\n")
                        exit()
                    else:
                        options.exclude.append(sub_e)
            else:
                options.exclude = []
        elif options.organelle_types == 'embplant_pt':
            options.include_priority = [os.path.join(NOT_DB_PATH, 'embplant_pt.fasta')]
            options.exclude = [os.path.join(NOT_DB_PATH, 'embplant_mt.fasta')]
        elif options.organelle_types == 'embplant_mt':
            options.include_priority = [os.path.join(NOT_DB_PATH, 'embplant_mt.fasta')]
            options.exclude = [os.path.join(NOT_DB_PATH, 'embplant_pt.fasta')]
        # elif options.organelle_types == 'embplant_nr':
        #     options.include_priority = [os.path.join(NOT_DB_PATH, 'embplant_nr')]
        # elif options.organelle_types == 'animal_mt':
        #     options.include_priority = [os.path.join(NOT_DB_PATH, 'animal_mt')]
        # elif options.organelle_types == 'fungus_mt':
        #     options.include_priority = [os.path.join(NOT_DB_PATH, 'fungus_mt')]
        else:
            if options.organelle_types:
                options.include_priority = []
                options.organelle_types = options.organelle_types.split(",")
                for sub_organelle_t in options.organelle_types:
                    if sub_organelle_t not in {"embplant_pt", "other_pt", "embplant_mt", "embplant_nr", "animal_mt", "fungus_mt"}:
                        sys.stdout.write(
                            "\n############################################################################"
                            "\nERROR: \"-F\" MUST be one of 'embplant_pt', 'other_pt', 'embplant_mt', 'embplant_nr', "
                            "'animal_mt', 'fungus_mt', or a combination of above split by comma(s)!\n\n")
                        exit()
                    else:
                        options.include_priority.append(os.path.join(NOT_DB_PATH, sub_organelle_t + ".fasta"))
        if not len(args):
            sys.stdout.write('\n\nInput Error: you must choose one input fasta or fastg file!\n')
            exit()
        output_dir = options.out_dir if options.out_dir else os.path.split(args[0])[0]
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        if options.out_base:
            options.out_base = os.path.basename(options.out_base)
        if options.generate_log:
            log_handler = simple_log(logging.getLogger(), output_dir,
                                     str((options.prefix + options.out_base + ".") * int(
                                          bool(options.out_base))) + "slim.")
            log_handler.info(print_title)
            log_handler.info(" ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n")
            log_handler = timed_log(log_handler, output_dir,
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


def blast_and_call_names(fasta_file, index_files, out_file, is_fastg, threads, which_blast="", log_handler=None):
    names = {}
    if index_files:
        time0 = time.time()
        if is_fastg:
            fasta_file += '.Temp'
        for index_f in index_files:
            if log_handler:
                log_handler.info('Executing BLAST to ' + index_f + ' ...')
            else:
                sys.stdout.write('\nExecuting BLAST to ' + index_f + '...')
            execute_blast(query=fasta_file, blast_db=index_f, output=out_file, threads=threads,
                          outfmt=6, e_value="1e-25", which_blast=which_blast, log_handler=log_handler)
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
                if is_fastg:
                    query, hit = line_split[0].split('_')[1], line_split[1]
                else:
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


def modify_in_ex(in_names, ex_names, significant, coverages, depth_cutoff, log_handler=None):
    here_in_names, here_ex_names = copy.deepcopy(in_names), copy.deepcopy(ex_names)
    if coverages:
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
                this_coverage = coverages[here_q_name]
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
                               (1+abs(math.log(coverages[q_name]/aver_in_coverage)))
                    ex_score = summarize_q_score(here_ex_names[q_name]) / \
                               (1+abs(math.log(coverages[q_name]/aver_ex_coverage)))
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
                        if 2 ** abs(math.log(coverages[q_name]/aver_in_coverage, 2)) <= 5:
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
                if 2 ** abs(math.log(coverages[q_name] / aver_in_coverage, 2)) >= depth_cutoff:
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


def map_names(in_names, ex_names, candidates, is_fastg, aver_in_dep, coverages, depth_cutoff, treat_no_hits,
              include_priority_assigned, include_assigned, exclude_priority_assigned, exclude_assigned,
              log_handler=None):
    if log_handler:
        log_handler.info("Mapping names ...")
    time0 = time.time()
    short_candidates = {}
    if is_fastg:
        if treat_no_hits == 'ex_no_con':
            short_connections = {}
            for candidate in candidates:
                this_short = candidate.split('_')[1]
                if aver_in_dep and 2 ** abs(math.log(coverages[this_short] / aver_in_dep, 2)) >= depth_cutoff:
                    continue
                #
                if this_short in short_candidates:
                    short_candidates[this_short].append(candidate)
                else:
                    short_candidates[this_short] = [candidate]
                #
            for candidate in candidates:
                this_short = candidate.split('_')[1]
                if this_short not in short_connections:
                    short_connections[this_short] = set()
                if ':' in candidate:
                    for node in candidate.split(':')[1].strip().split(','):
                        if node.split('_')[1] in short_candidates:
                            short_connections[this_short].add(node.split('_')[1])
            # expand include
            new_nodes = 1
            old_nodes = copy.deepcopy(in_names)
            while new_nodes:
                new_nodes = set()
                for in_node in list(old_nodes):
                    if in_node in short_connections:
                        for new_node in short_connections[in_node]:
                            if new_node not in in_names and new_node not in ex_names:
                                in_names.add(new_node)
                                new_nodes.add(new_node)
                old_nodes = copy.deepcopy(new_nodes)
        else:
            for candidate in candidates:
                this_short = candidate.split('_')[1]
                if this_short in short_candidates:
                    short_candidates[this_short].append(candidate)
                else:
                    short_candidates[this_short] = [candidate]
    else:
        for candidate in candidates:
            this_short = candidate.split(' ')[0]
            if this_short in short_candidates:
                short_candidates[this_short].append(candidate)
            else:
                short_candidates[this_short] = [candidate]
    accepted = set()
    if exclude_priority_assigned:
        if treat_no_hits in ["ex_no_con", "ex_no_hit"]:
            for this_short, full_name in short_candidates.items():
                if this_short in ex_names:
                    pass
                elif this_short in in_names:
                    for name in full_name:
                        accepted.add(name)
                else:
                    pass
        else:
            for this_short, full_name in short_candidates.items():
                if this_short in ex_names:
                    pass
                else:
                    for name in full_name:
                        accepted.add(name)
    elif include_priority_assigned or include_assigned:
        if treat_no_hits in ["ex_no_con", "ex_no_hit"]:
            for this_short, full_name in short_candidates.items():
                if this_short in in_names:
                    for name in full_name:
                        accepted.add(name)
        else:
            for this_short, full_name in short_candidates.items():
                if this_short in in_names:
                    for name in full_name:
                        accepted.add(name)
                elif this_short in ex_names:
                    pass
                else:
                    for name in full_name:
                        accepted.add(name)
    elif exclude_assigned:
        for this_short, full_name in short_candidates.items():
            if this_short in ex_names:
                pass
            else:
                for name in full_name:
                    accepted.add(name)
    else:
        for this_short, full_name in short_candidates.items():
            for name in full_name:
                accepted.add(name)
    if log_handler:
        log_handler.info("Mapping names finished.")
    else:
        sys.stdout.write('\nMapping names cost: '+str(round(time.time()-time0, 2)) + "\n")
    return accepted, short_candidates


def del_complementary(fastg_file, is_fastg, log_handler=None):
    if is_fastg:
        if log_handler:
            log_handler.info("Deleting complementary ...")
        time0 = time.time()
        temp_matrix = read_fasta(fasta_dir=fastg_file)
        i = 0
        while i < len(temp_matrix[0]):
            if temp_matrix[0][i].rstrip(';').split(':')[0].endswith('\''):
                del temp_matrix[0][i]
                del temp_matrix[1][i]
            else:
                i += 1
        write_fasta(out_file=fastg_file + '.Temp', matrix=temp_matrix, overwrite=True)
        if log_handler:
            log_handler.info("Deleting complementary finished.")
        else:
            sys.stdout.write('\nDeleting complementary cost: '+str(round(time.time()-time0, 2)) + "\n")


def write_hits_tab_for_bandage(in_names, ex_names, out_file, overwrite, is_fastg,
                               short_candidates, log_handler=None):
    time0 = time.time()
    if not overwrite:
        while os.path.exists(out_file+'.csv'):
            out_file += '_'
    out_file += '.csv'
    out_lines = []
    # if include_file:
    #     in_database = os.path.split(include_file)[-1]
    # if exclude_file:
    #     ex_database = os.path.split(exclude_file)[-1]
    edges = set(in_names)
    for edge in ex_names:
        edges.add(edge)
    edges = list(edges)
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
        if is_fastg:
            try:
                full_name = short_candidates[edge]
                full_name.sort()
                if ':' in full_name[1]:
                    next_edges = [x.rstrip('\'').split('_') for x in full_name[1].split(':')[1].split(',')]
                    next_edges.sort(key=lambda x: -float(x[5].rstrip(';').rstrip('\'')))
                    postfix = '>>' + next_edges[0][1]
            except (IndexError, KeyError):
                pass
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
    return short_candidates


def remove_temp_files(fastg_file, keep_temp, is_fastg, other_files_to_remove):
    if not keep_temp:
        if is_fastg and os.path.exists(fastg_file+'.Temp'):
            os.remove(fastg_file+'.Temp')
        try:
            os.remove(fastg_file+'.blast_in')
        except OSError:
            pass
        try:
            os.remove(fastg_file+'.blast_ex')
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
    output_dir = options.out_dir if options.out_dir else os.path.split(args[0])[0]
    log_output_name = str((options.prefix + options.out_base + ".") * int(bool(options.out_base))) + "slim."
    out_base = os.path.join(output_dir, options.prefix + options.out_base) * int(bool( options.out_base))
    try:
        if not log_handler:
            sys.stdout.write('\n' + '=' * 100)
        if options.out_dir and not os.path.exists(options.out_dir):
            os.mkdir(options.out_dir)
        for i in range(len(args)):
            if log_handler:
                log_handler.info('Slimming file '+str(i+1)+'/'+str(len(args))+': '+args[i])
            else:
                sys.stdout.write('\nSlimming file '+str(i+1)+'/'+str(len(args))+': '+args[i]+'\n')
            fas_file = args[i]
            is_fastg = fas_file.endswith('.fastg')
            if not options.out_dir:
                options.out_dir = os.path.split(os.path.realpath(fas_file))[0]
            # 0. check file name. Replace illegal characters for blastn
            if "'" in fas_file:
                copyfile(fas_file, fas_file.replace("'", "_"))
                fas_file = fas_file.replace("'", "_")
            # 1. filter by depth
            fas_file = filter_fastg_by_depth_simple(fas_file, options.max_depth, options.min_depth, out_dir=options.out_dir)
            # 2. rm low coverage duplicated contigs
            if options.remove_low_duplicated and is_fastg:
                blastn_str = " --which-blast " + options.which_blast if options.which_blast else ""
                os.system(os.path.join(PATH_OF_THIS_SCRIPT, 'rm_low_coverage_duplicated_contigs.py') +
                          ' -t ' + str(options.threads) + blastn_str + " " + fas_file + " -o " + options.out_dir)
                fas_file = os.path.join(options.out_dir, os.path.basename(fas_file) + '.purified.fastg')
            # 3. merge contigs
            if options.merge_contigs and is_fastg:
                this_assembly = Assembly(graph_file=fas_file)
                time_m0 = time.time()
                merged = this_assembly.merge_all_possible_vertices()
                if merged:
                    fas_file = os.path.join(options.out_dir, os.path.basename(fas_file)[:-6] + ".merged.fastg")
                    this_assembly.write_to_fastg(fas_file, rename_if_needed=True)
                    if log_handler:
                        log_handler.info("Merging contigs finished.")
                    else:
                        sys.stdout.write("\nMerging contigs cost: " + str(round(time.time() - time_m0, 2)) + "\n")
            # if 1,2&3 neither work, and the output directory is different from the original file,
            # copy the original file to the destination for downstream 4, 5
            other_files_to_remove = []
            if os.path.realpath(fas_file) != os.path.realpath(os.path.join(options.out_dir, os.path.basename(fas_file))):
                os.system("cp " + fas_file + " " + os.path.join(options.out_dir, os.path.basename(fas_file)))
                fas_file = os.path.join(options.out_dir, os.path.basename(fas_file))
                other_files_to_remove.append(fas_file)
            # make blast database if not made
            include_indices, exclude_indices, database_made = \
                check_db(options.include_priority, options.include,
                         options.exclude_priority, options.exclude, options.out_dir,
                         which_blast=options.which_blast, log_handler=log_handler, verbose_log=options.verbose_log)
            other_files_to_remove.extend(database_made)
            # 4. del complementary to avoid using both the forward and reverse sequence for blast
            if bool(include_indices) or bool(exclude_indices):
                del_complementary(fas_file, is_fastg, log_handler=log_handler)
            in_ex_info = generate_in_ex_info_name(include_indices=include_indices, exclude_indices=exclude_indices,
                                                  exclude_no_con=options.treat_no_hits == 'ex_no_con',
                                                  exclude_no_hit=options.treat_no_hits == 'ex_no_hit')
            # 5.
            try:
                in_names = blast_and_call_names(fasta_file=fas_file, index_files=include_indices,
                                                out_file=fas_file+'.blast_in', is_fastg=is_fastg,
                                                threads=options.threads, which_blast=options.which_blast,
                                                log_handler=log_handler)
                ex_names = blast_and_call_names(fasta_file=fas_file, index_files=exclude_indices,
                                                out_file=fas_file+'.blast_ex', is_fastg=is_fastg,
                                                threads=options.threads, which_blast=options.which_blast,
                                                log_handler=log_handler)
                # write out fasta according to blast
                fasta_matrix = read_fasta(fasta_dir=fas_file)
                coverages = get_graph_coverage_dict_simple(fasta_matrix=fasta_matrix, is_fastg=is_fastg)
                if bool(include_indices) or bool(exclude_indices):
                    in_names_r, ex_names_r, aver_dep = modify_in_ex(in_names=in_names, ex_names=ex_names,
                                                                    significant=options.significant,
                                                                    coverages=coverages,
                                                                    depth_cutoff=options.depth_cutoff,
                                                                    log_handler=log_handler)
                else:
                    in_names_r, ex_names_r, aver_dep = in_names, ex_names, None
                accept_names, short_candidates = \
                    map_names(in_names=set(in_names_r), ex_names=set(ex_names_r), candidates=fasta_matrix[0],
                              is_fastg=is_fastg, aver_in_dep=aver_dep, coverages=coverages,
                              depth_cutoff=options.depth_cutoff, treat_no_hits=options.treat_no_hits,
                              include_priority_assigned=options.include_priority,
                              include_assigned=options.include,
                              exclude_priority_assigned=options.exclude_priority,
                              exclude_assigned=options.exclude, log_handler=log_handler)
                fasta_matrix = make_new_matrix_with_names(names=accept_names, old_matrix=fasta_matrix)
                if out_base:
                    out_fas = out_base + '.' + str(fas_file).split('.')[-1]
                    out_csv = out_base
                else:
                    out_fas = os.path.join(options.out_dir, options.prefix + os.path.basename(fas_file)) + in_ex_info\
                              + '.' + str(fas_file).split('.')[-1]
                    out_csv = os.path.join(options.out_dir, options.prefix + os.path.basename(fas_file)) + in_ex_info
                if log_handler:
                    log_handler.info("Generating slimmed file to " + out_fas)
                else:
                    sys.stdout.write("Generating slimmed file to " + out_fas + "\n")
                write_fasta(out_file=out_fas,
                            matrix=fasta_matrix, overwrite=options.overwrite)
                # write out hits tab according to blast
                if not options.no_tab:
                    write_hits_tab_for_bandage(
                        in_names=in_names, ex_names=ex_names, out_file=out_csv, overwrite=options.overwrite,
                        is_fastg=is_fastg, short_candidates=short_candidates, log_handler=log_handler)
            except Exception as e:
                if log_handler:
                    log_handler.error("\n" + str(e).strip())
                    log_handler.error('Slimming file ' + str(i + 1) + '/' + str(len(args)) + ': ' + args[i] + ' failed!\n')
                else:
                    sys.stdout.write("\n" + str(e).strip() + "\n")
                    sys.stdout.write('\nSlimming file ' + str(i + 1) + '/' + str(len(args)) + ': ' + args[i] + ' failed!'
                                     '\n' + '=' * 100 + "\n")
            else:
                if log_handler:
                    log_handler.info('Slimming file ' + str(i+1) + '/' + str(len(args)) + ': ' + args[i] + ' finished!\n')
                else:
                    sys.stdout.write('\nSlimming file ' + str(i+1) + '/' + str(len(args)) + ': ' + args[i] + ' finished!'
                                     '\n' + '=' * 100 + "\n")
            remove_temp_files(fas_file, options.keep_temp, is_fastg, other_files_to_remove)

        if log_handler:
            log_handler = simple_log(log_handler, output_dir, log_output_name)
            log_handler.info("\nTotal cost " + "%.2f" % (time.time() - time0) + " s")
            log_handler.info("Thank you!")
        else:
            sys.stdout.write('\nTotal cost: ' + str(round(time.time() - time0, 2)) + '\nThank you!\n')
    except Exception as e:
        if log_handler:
            log_handler.exception("")
            log_handler = simple_log(log_handler, output_dir, log_output_name)
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