#!/usr/bin/env python
# coding: utf8
import time
import os
import sys
import subprocess
import math
try:
    # python2
    import commands
except:
    pass
from optparse import OptionParser
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
import optparse
import copy



options = ''
args = []
short_candidates = {}


def require_commands():
    global options, args
    try:
        # python3
        blast_in_path = subprocess.getstatusoutput('blastn')
    except AttributeError:
        # python2
        blast_in_path = commands.getstatusoutput('blastn')
    if blast_in_path[0] == 32512:
        sys.stdout.write('\nError: blastn not in the path!')
        exit()
    try:
        # python3
        makeblastdb_in_path = subprocess.getstatusoutput('makeblastdb')
    except AttributeError:
        # python2
        makeblastdb_in_path = commands.getstatusoutput('makeblastdb')
    if makeblastdb_in_path[0] == 32512:
        sys.stdout.write('\nError: makeblastdb not in the path!')
        exit()
    usage = 'python '+str(os.path.basename(__file__)+' your_fastg_files -F cp')
    parser = OptionParser(usage=usage)
    # parser.add_option('-o', dest='out_fastg_file', help='Output file')
    # filters
    parser.add_option('-F', dest='builtin_mode', default='cp',
                      help='followed with mode cp, mt, nr (which means chloroplast, mitochondria, nrDNA'
                           'separately; corresponding to certain arguments as following listed). '
                           'Modify the arguments activated by this flag with your more custom options.'
                           '\t'
                           ' ------------------------------------------------------ '
                           '\ncp \t " --include-priority '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'cp')+' --exclude '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'mt')+'"'
                           ' ------------------------------------------------------ '
                           '\nmt \t " --include-priority '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'mt')+' --exclude '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'cp')+'"'
                           ' ------------------------------------------------------ '
                           '\nnr \t " --include-priority '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'nr')+'"'
                           ' ------------------------------------------------------ ')
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
                           ' detected coverage would be excluded. Default: 100.0')
    parser.add_option('--depth-threshold', dest='depth_threshold',
                      help='Input a float or integer number. filter fastg file by depth. Default: no threshold.')
    parser.add_option('--include', dest='include',
                      help='followed by Blast index format')
    parser.add_option('--include-priority', dest='include_priority',
                      help='followed by Blast index format.')
    parser.add_option('--exclude', dest='exclude',
                      help='followed by Blast index format.')
    parser.add_option('--exclude-priority', dest='exclude_priority',
                      help='followed by Blast index format')
    parser.add_option('--no-hits-labeled-tab', dest='no_tab', default=False, action='store_true',
                      help='Choose to disable producing tab file')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true',
                      help='Choose to disable deleting temp files produced by blast and this script')
    parser.add_option('-o', '--out-dir', dest="out_dir",
                      help="By default the output would be along with the input fastg file. "
                           "But you could assign a new directory with this option.")
    parser.add_option("--prefix", dest="prefix", default="",
                      help="Add prefix to the output file.")
    parser.add_option('--continue', dest='resume', default=False, action='store_true',
                      help='Specified for calling from get_organelle_reads.py')
    parser.add_option('-t', '--threads', dest="threads", default=4, type=int,
                      help="Threads for blastn.")
    try:
        options, args = parser.parse_args()
    except optparse.OptionConflictError as e:
        sys.stdout.write('\n\n######################################'+str(e))
        sys.stdout.write('\n\n"-h" for more usage')
        exit()
    else:
        if options.treat_no_hits not in {"ex_no_con", "ex_no_hit", "keep_all"}:
            sys.stdout.write('\n\nOption Error: you should choose assign one of "ex_no_con", "ex_no_hit"'
                             ' and "keep_all" to variable treat_no_hits')
            exit()
        priority_chosen = int(bool(options.include_priority)) + int(bool(options.exclude_priority))
        secondary_chosen = int(bool(options.include)) + int(bool(options.exclude))
        if priority_chosen + secondary_chosen > 0:
            sys.stdout.write("\nbuiltin_mode is disabled since you assign the custom index/indices.")
            if priority_chosen > 1:
                sys.stdout.write('\n\nLogical Error: only one option with "-priority" allowed!')
                exit()
            in_chosen = int(bool(options.include_priority)) + int(bool(options.include))
            if in_chosen > 1:
                sys.stdout.write('\n\nOption Error: you can not simultaneously choose two "--include-*" options!')
                exit()
            ex_chosen = int(bool(options.exclude_priority)) + int(bool(options.exclude))
            if ex_chosen > 1:
                sys.stdout.write('\n\nOption Error: you can not simultaneously choose two "--exclude-*" options!')
                exit()
            if in_chosen == 1 and ex_chosen == 1 and priority_chosen == 0:
                sys.stdout.write('\n\nOption Error: since you have include and exclude chosen, one of them should be assigned priority!')
                exit()
            if ex_chosen == 1 and in_chosen == 0 and (options.treat_no_hits in {"ex_no_con", "ex_no_hit"}):
                sys.stdout.write('\n\nOption Error: no contigs survive according to you choice!')
                exit()
        elif options.builtin_mode == 'cp':
            options.include_priority = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'cp')
            options.exclude = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'mt')
        elif options.builtin_mode == 'mt':
            options.include_priority = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'mt')
            options.exclude = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'cp')
        elif options.builtin_mode == 'nr':
            options.include_priority = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'NotationReference', 'nr')
        else:
            sys.stdout.write('\n\nOption Error: illegal value for builtin mode!\n')
            exit()
        if not len(args):
            sys.stdout.write('\n\nInput Error: you must choose one input fasta or fastg file!\n')
            exit()
        sys.stdout.write('\n'+' '.join(sys.argv)+'\n')


def check_db(options):
    in_index = ""
    ex_index = ""
    if options.include_priority:
        in_index = options.include_priority
    elif options.include:
        in_index = options.include
    if options.exclude_priority:
        ex_index = options.exclude_priority
    elif options.exclude:
        ex_index = options.exclude
    return in_index, ex_index


def make_new_matrix_with_names(names, old_matrix):
    i = 0
    while i < len(old_matrix[0]):
        if old_matrix[0][i] in names:
            i += 1
        else:
            del old_matrix[0][i]
            del old_matrix[1][i]
    return old_matrix


def blast_and_call_names(fasta_file, index_files, out_file, is_fastg, threads):
    if index_files:
        time0 = time.time()
        sys.stdout.write('\nblast ...')
        if is_fastg:
            fasta_file += '.Temp'
        try:
            blast_result = subprocess.getstatusoutput('blastn -num_threads '+ str(threads) +' -query ' + fasta_file + ' -db ' + index_files + ' -out ' + out_file + ' -outfmt 6 -evalue 1e-15')
        except AttributeError:
            blast_result = commands.getstatusoutput('blastn -num_threads '+ str(threads) +' -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
        # blastn -num_threads 4 -query assembly_graph.fastg -db db_f -out out_f -outfmt 6
        if 'Error' in str(blast_result[1]) or 'error' in str(blast_result[1]) or '不是内部或外部命令' in str(blast_result[1]):
            # os.system('blastn -num_threads '+ str(threads) +' -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
            # if not os.path.exists(out_file):
            sys.stdout.write('Blast command: blastn -num_threads ' + str(
                threads) + ' -query ' + fasta_file + ' -db ' + index_files + ' -out ' + out_file +
                ' -outfmt 6 -evalue 1e-15\n')
            sys.stdout.write('\nBlast terminated with following info:\n'+str(blast_result[1]) + "\n")
            raise EnvironmentError
        # windows
        # if not os.path.exists(out_file):
        #     os.system('blastn -num_threads '+ str(threads) +' -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
        time1 = time.time()
        sys.stdout.write('\nblast to '+os.path.split(index_files)[-1]+' cost: '+str(time1-time0))
        names = {}
        try:
            blast_out_lines = open(out_file, 'rU')
        except IOError:
            sys.stdout.write('\nBlast was not properly installed or configurated.')
            raise EnvironmentError
        for line in blast_out_lines:
            line_split = line.split('\t')
            if is_fastg:
                query, hit = line_split[0].split('_')[1], line_split[1]
            else:
                query, hit = line_split[0], line_split[1]
            q_start, q_end, q_score = int(line_split[6]), int(line_split[7]), float(line_split[2])
            q_min, q_max = min(q_start, q_end), max(q_start, q_end)
            # q_score = abs(q_max - q_min + 1)*q_score
            if query in names:
                if hit not in names[query]:
                    names[query][hit] = [(q_min, q_max, q_score)]
                else:
                    # if overlap, then merge
                    i = 0
                    while i < len(names[query][hit]):
                        this_min, this_max, previous_q_score = names[query][hit][i]
                        if q_max < this_min:
                            break
                        elif q_min > this_max:
                            i += 1
                            continue
                        else:
                            q_min = min(q_min, this_min)
                            q_max = max(q_max, this_max)
                            q_score = max(previous_q_score, q_score)
                            del names[query][hit][i]
                    names[query][hit].insert(i, (q_min, q_max, q_score))
            else:
                names[query] = {}
                names[query][hit] = [(q_min, q_max, q_score)]
        sys.stdout.write('\nparse blast result cost'+str(time.time()-time1))
        return names
    else:
        return {}


def get_coverages(matrix, is_fastg):
    if is_fastg:
        coverages = {}
        for fastg_name in matrix[0]:
            this_coverage = float(fastg_name.split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0])
            coverages[fastg_name.split('_')[1]] = this_coverage
        return coverages
    else:
        return []


def summarize_q_score(hits):
    all_loci = []
    for hit in hits:
        for locus in hit:
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


def modify_in_ex(in_names, ex_names, significant, coverages, depth_cutoff):
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
                in_name_score = summarize_q_score(list(in_names[query_name].values()))
                ex_name_score = summarize_q_score(list(ex_names[query_name].values()))
                if in_name_score / float(ex_name_score) > significant:
                    training_in.add(query_name)
                elif ex_name_score / float(in_name_score) > significant:
                    training_ex.add(query_name)

        def get_average_coverage(here_q_list, here_info):
            total_in_base = 0
            total_in_len = 0
            for here_q_name in here_q_list:
                this_coverage = coverages[here_q_name]
                for hit in here_info[here_q_name].values():
                    for q_min, q_max, q_score in hit:
                        this_len = abs(q_max - q_min) + 1
                        total_in_len += this_len
                        total_in_base += this_len * this_coverage
            return total_in_base / float(total_in_len)

        def combine_coverage_to_check():
            aver_in_coverage = get_average_coverage(training_in, in_names)
            aver_ex_coverage = get_average_coverage(training_ex, ex_names)
            sys.stdout.write('\naverage in coverage: '+str(aver_in_coverage))
            sys.stdout.write('\naverage ex coverage: '+str(aver_ex_coverage))
            for q_name in list(here_in_names):
                if q_name in here_ex_names:
                    in_score = summarize_q_score(list(here_in_names[q_name].values())) / \
                               (1+abs(math.log(coverages[q_name]/aver_in_coverage)))
                    ex_score = summarize_q_score(list(here_ex_names[q_name].values())) / \
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
                            sys.stdout.write('\n' + ' ' * 4 + q_name + ' excluded: '
                                             + str(round(in_score, 2)) + '~' + str(round(ex_score, 2)))
                    else:
                        pass
                        # sys.stdout.write('\n' + ' ' * 4 + query_name + ' with both hits: '
                        #                  + str(round(in_name_score, 2)) + '~' + str(round(ex_name_score, 2)))
            aver_in_coverage = get_average_coverage(training_in, in_names)
            aver_ex_coverage = get_average_coverage(training_ex, ex_names)
            sys.stdout.write('\naverage in coverage: ' + str(aver_in_coverage))
            sys.stdout.write('\naverage ex coverage: ' + str(aver_ex_coverage))
            for q_name in list(here_in_names):
                if 2 ** abs(math.log(coverages[q_name] / aver_in_coverage, 2)) >= depth_cutoff:
                    del here_in_names[q_name]
            return here_in_names, here_ex_names, aver_in_coverage

        if training_in and training_ex:
            combine_coverage_to_check()
        else:
            for query_name in list(here_in_names):
                if query_name in here_ex_names:
                    in_name_score = summarize_q_score(list(here_in_names[query_name].values()))
                    ex_name_score = summarize_q_score(list(here_ex_names[query_name].values()))
                    if in_name_score / float(ex_name_score) > significant:
                        training_in.add(query_name)
                        training_ex.discard(query_name)
                    elif ex_name_score / float(in_name_score) > significant:
                        training_in.discard(query_name)
                        training_ex.add(query_name)
            if training_in and training_ex:
                combine_coverage_to_check()
            else:
                sys.stdout.write('\nNot enough coverage information found.')
    else:
        sys.stdout.write('\nNo coverage information found.')
    for query_name in list(here_in_names):
        if query_name in here_ex_names:
            in_name_score = summarize_q_score(list(here_in_names[query_name].values()))
            ex_name_score = summarize_q_score(list(here_ex_names[query_name].values()))
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


def map_names(in_names, ex_names, candidates, is_fastg, aver_in_dep, coverages, depth_cutoff):
    global options, short_candidates
    time0 = time.time()
    if is_fastg:
        if options.treat_no_hits == 'ex_no_con':
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
    if options.exclude_priority:
        if options.treat_no_hits in {"ex_no_con", "ex_no_hit"}:
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
    elif options.include_priority or options.include:
        if options.treat_no_hits in {"ex_no_con", "ex_no_hit"}:
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
    elif options.exclude:
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
    sys.stdout.write('\nmap names cost: '+str(time.time()-time0))
    return accepted


def filter_fastg_by_depth(fas_file, depth):
    if fas_file.endswith('.fastg') and depth and float(depth):
        depth = float(depth)
        time0 = time.time()
        fastg_matrix = read_fasta(fas_file)
        new_fastg_matrix = [[], [], fastg_matrix[2]]
        for i in range(len(fastg_matrix[0])):
            if float(fastg_matrix[0][i].split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0]) >= depth:
                new_fastg_matrix[0].append(fastg_matrix[0][i])
                new_fastg_matrix[1].append(fastg_matrix[1][i])
        out_fasta = '.'.join(fas_file.split('.')[:-1]) + '.depth' + str(depth) + '.' + fas_file.split('.')[-1]
        write_fasta(out_dir=out_fasta, matrix=new_fastg_matrix, overwrite=True)
        sys.stdout.write('\nfilter by depth cost: '+str(time.time()-time0))
        return out_fasta
    else:
        return fas_file


def del_complementary(fastg_file, is_fastg):
    if is_fastg:
        time0 = time.time()
        temp_matrix = read_fasta(fasta_dir=fastg_file)
        i = 0
        while i < len(temp_matrix[0]):
            if temp_matrix[0][i].rstrip(';').split(':')[0].endswith('\''):
                del temp_matrix[0][i]
                del temp_matrix[1][i]
            else:
                i += 1
        write_fasta(out_dir=fastg_file + '.Temp', matrix=temp_matrix, overwrite=True)
        sys.stdout.write('\ndel complementary cost: '+str(time.time()-time0))


def write_hits_tab_for_bandage(in_names, include_file, ex_names, exclude_file, out_file, overwrite, is_fastg):
    global options, short_candidates
    if options.no_tab:
        return ''
    else:
        time0 = time.time()
        if not overwrite:
            while os.path.exists(out_file+'.csv'):
                out_file += '_'
        out_file += '.csv'
        out_lines = []
        if include_file:
            in_database = os.path.split(include_file)[-1]
        if exclude_file:
            ex_database = os.path.split(exclude_file)[-1]
        edges = set(in_names)
        for edge in ex_names:
            edges.add(edge)
        edges = list(edges)
        for edge in edges:
            this_string = edge
            if edge in in_names and edge in ex_names:
                this_string += '\t'+in_database+';'+ex_database
                this_string += '\t'+','.join(list(in_names[edge]))+';'+','.join(list(ex_names[edge]))
                loci = []
                for locus in in_names[edge]:
                    for region in in_names[edge][locus]:
                        loci.append([region[0], region[1], locus, in_database])
                for locus in ex_names[edge]:
                    for region in ex_names[edge][locus]:
                        loci.append([region[0], region[1], locus, ex_database])
            elif edge in in_names:
                this_string += '\t'+in_database
                this_string += '\t'+','.join(list(in_names[edge]))
                loci = []
                for locus in in_names[edge]:
                    for region in in_names[edge][locus]:
                        loci.append([region[0], region[1], locus, in_database])
            else:
                this_string += '\t'+ex_database
                this_string += '\t'+','.join(list(ex_names[edge]))
                loci = []
                for locus in ex_names[edge]:
                    for region in ex_names[edge][locus]:
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
        sys.stdout.write('\ncreate tab cost: '+str(time.time()-time0))


def remove_temp_files(fastg_file, keep_temp, is_fastg):
    if not keep_temp:
        if is_fastg:
            os.remove(fastg_file+'.Temp')
        try:
            os.remove(fastg_file+'.blast_in')
        except OSError:
            pass
        try:
            os.remove(fastg_file+'.blast_ex')
        except OSError:
            pass


__version__ = "1.8.0"


def main():
    time0 = time.time()
    sys.stdout.write("\nThis is a script for excluding certain contigs from "
                     "assembly graph file (*.fastg/*.fasta) by blast\n"
                     "By jinjianjun@mail.kib.ac.cn\n")
    require_commands()
    global options, args
    for i in range(len(args)):
        sys.stdout.write('\nRound '+str(i+1)+'/'+str(len(args))+': '+args[i]+'\n')
        fas_file = args[i]
        is_fastg = fas_file.endswith('.fastg')
        # prepare fasta file
        fas_file = filter_fastg_by_depth(fas_file, options.depth_threshold)
        # rm low coverage duplicated contigs
        if options.remove_low_duplicated and is_fastg:
            os.system('rm_low_coverage_duplicated_contigs.py -t '+str(options.threads)+" "+fas_file)
            fas_file += '.purified.fastg'
        del_complementary(fas_file, is_fastg)
        # make blast database if not made
        include_index, exclude_index = check_db(options)
        in_ex_info = 'only' * int(options.treat_no_hits == 'ex_no_hit') + 'extend' * int(
            options.treat_no_hits == 'ex_no_con') + '+' + os.path.split(include_index)[-1] + '-' + \
            os.path.split(exclude_index)[-1]
        # make blast
        try:
            in_names = blast_and_call_names(fasta_file=fas_file, index_files=include_index,
                                            out_file=fas_file+'.blast_in', is_fastg=is_fastg, threads=options.threads)
            ex_names = blast_and_call_names(fasta_file=fas_file, index_files=exclude_index,
                                            out_file=fas_file+'.blast_ex', is_fastg=is_fastg, threads=options.threads)
            # write out fasta according to blast
            fasta_matrix = read_fasta(fasta_dir=fas_file)
            coverages = get_coverages(matrix=fasta_matrix, is_fastg=is_fastg)
            in_names_r, ex_names_r, aver_dep = modify_in_ex(in_names=in_names, ex_names=ex_names,
                                                            significant=options.significant, coverages=coverages,
                                                            depth_cutoff=options.depth_cutoff)
            accept_names = map_names(in_names=set(in_names_r), ex_names=set(ex_names_r), candidates=fasta_matrix[0],
                                     is_fastg=is_fastg, aver_in_dep=aver_dep, coverages=coverages,
                                     depth_cutoff=options.depth_cutoff)
            fasta_matrix = make_new_matrix_with_names(names=accept_names, old_matrix=fasta_matrix)
            if options.out_dir:
                if not os.path.exists(options.out_dir):
                    os.mkdir(options.out_dir)
                out_fas = os.path.join(options.out_dir, options.prefix + os.path.basename(fas_file)+'.'+in_ex_info+'.'+fas_file.split('.')[-1])
                out_csv = os.path.join(options.out_dir, options.prefix + os.path.basename(fas_file) + '.' + in_ex_info)
            else:
                out_fas = os.path.join(os.path.split(fas_file)[0], options.prefix + os.path.basename(fas_file)) + '.'+ in_ex_info + '.' + fas_file.split('.')[-1]
                out_csv = os.path.join(os.path.split(fas_file)[0], options.prefix + os.path.basename(fas_file)) + '.' + in_ex_info
            write_fasta(out_dir=out_fas,
                        matrix=fasta_matrix, overwrite=False)
            # write out hits tab according to blast
            write_hits_tab_for_bandage(in_names=in_names, include_file=include_index, ex_names=ex_names,
                                       exclude_file=exclude_index, out_file=out_csv, overwrite=False,
                                       is_fastg=is_fastg)
        except EnvironmentError:
            sys.stdout.write('\nRound ' + str(i + 1) + '/' + str(len(args)) + ': ' + args[i] + ' failed!\n')
        else:
            sys.stdout.write('\nRound ' + str(i+1) + '/' + str(len(args)) + ': ' + args[i] + ' finished!\n')
        remove_temp_files(fas_file, options.keep_temp, is_fastg)
    sys.stdout.write('\nTotal cost: '+str(time.time()-time0)+'\n\n')


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