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
from optparse import OptionParser
import optparse
import copy


path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
options = ''
short_candidates = {}


def require_commands():
    global options
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
    usage = 'python '+str(os.path.basename(__file__)+' -g your_fastg_file -m cp')
    parser = OptionParser(usage=usage)
    parser.add_option('-g', dest='in_fastg_file', help='followed by your input fastg file')
    parser.add_option('-f', dest='in_fasta_file', help='followed by your input fasta file')
    # parser.add_option('-o', dest='out_fastg_file', help='Output file')
    # filters
    parser.add_option('-m', dest='builtin_mode', default='cp',
                      help='followed with mode cp, mt, nr (which means chloroplast, mitochondria, nrDNA'
                           'separately; corresponding to certain arguments as following listed). '
                           'Modify the arguments activated by this flag with your more custom options.'
                           '\t'
                           ' ------------------------------------------------------ '
                           '\ncp \t " --include-index-priority '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'cp')+' --exclude-index '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'mt')+'"'
                           ' ------------------------------------------------------ '
                           '\nmt \t " --include-index-priority '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'mt')+' --exclude-index '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'cp')+'"'
                           ' ------------------------------------------------------ '
                           '\nnr \t " --include-index-priority '+os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'nr')+'"'
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
    parser.add_option('--depth-threshold', dest='depth_threshold',
                      help='Input a float or integer number. filter fastg file by depth. Default: no threshold.')
    parser.add_option('--include-index', dest='in_nc_base',
                      help='followed by Blast index format')
    parser.add_option('--include-index-priority', dest='in_nc_base_priority',
                      help='followed by Blast index format.')
    parser.add_option('--exclude-index', dest='ex_nc_base',
                      help='followed by Blast index format.')
    parser.add_option('--exclude-index-priority', dest='ex_nc_base_priority',
                      help='followed by Blast index format')
    parser.add_option('--include-fasta', dest='in_fa_base',
                      help='followed by Fasta index format')
    parser.add_option('--include-fasta-priority', dest='in_fa_base_priority',
                      help='followed by Fasta index format')
    parser.add_option('--exclude-fasta', dest='ex_fa_base',
                      help='followed by Fasta index format')
    parser.add_option('--exclude-fasta-priority', dest='ex_fa_base_priority',
                      help='followed by Fasta index format')
    parser.add_option('--no-hits-labeled-csv', dest='no_csv', default=False, action='store_true',
                      help='Choose to disable producing csv file')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true',
                      help='Choose to disable deleting temp files produced by blast and this script')
    try:
        (options, args) = parser.parse_args()
    except optparse.OptionConflictError as e:
        sys.stdout.write('\n\n######################################'+str(e))
        sys.stdout.write('\n\n"-h" for more usage')
        exit()
    else:
        if options.treat_no_hits not in {"ex_no_con", "ex_no_hit", "keep_all"}:
            sys.stdout.write('\n\nOption Error: you should choose assign one of "ex_no_con", "ex_no_hit"'
                             ' and "keep_all" to variable treat_no_hits')
            exit()
        priority_chosen = int(bool(options.in_nc_base_priority)) + int(bool(options.ex_nc_base_priority)) + int(
            bool(options.in_fa_base_priority)) + int(bool(options.ex_fa_base_priority))
        secondary_chosen = int(bool(options.in_nc_base)) + int(bool(options.ex_nc_base)) + int(
            bool(options.in_fa_base)) + int(bool(options.ex_fa_base))
        if priority_chosen + secondary_chosen > 0:
            sys.stdout.write("\nbuiltin_mode is disabled since you assign the custom index/indices.")
            if priority_chosen > 1:
                sys.stdout.write('\n\nLogical Error: only one option with "-priority" allowed!')
                exit()
            in_chosen = int(bool(options.in_nc_base_priority)) + int(bool(options.in_nc_base)) + int(
                bool(options.in_fa_base_priority)) + int(bool(options.in_fa_base))
            if in_chosen > 1:
                sys.stdout.write('\n\nOption Error: you can not simultaneously choose two "--include-*" options!')
                exit()
            ex_chosen = int(bool(options.ex_nc_base_priority)) + int(bool(options.ex_nc_base)) + int(
                bool(options.ex_fa_base_priority)) + int(bool(options.ex_fa_base))
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
            options.in_nc_base_priority = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'cp')
            options.ex_nc_base = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'mt')
        elif options.builtin_mode == 'mt':
            options.in_nc_base_priority = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'mt')
            options.ex_nc_base = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'cp')
        elif options.builtin_mode == 'nr':
            options.in_nc_base_priority = os.path.join(os.path.split(path_of_this_script)[0], 'Library', 'Reference', 'nr')
        else:
            sys.stdout.write('\n\nOption Error: illegal value for builtin mode!')
            exit()
        if not options.in_fastg_file and (options.treat_no_hits == "ex_no_con"):
            sys.stdout.write('\n\nOption Error: ex_no_connect is available only when you input a fastg file!')
            exit()
        if int(bool(options.in_fastg_file))+int(bool(options.in_fasta_file)) != 1:
            sys.stdout.write('\n\nInput Error: you must choose one input fasta or fastg file!')
            exit()
        if int(bool(options.depth_threshold))+int(bool(options.in_fasta_file)) > 1:
            sys.stdout.write('\n\nOption Warning: you can use depth threshold only when you input a fastg file!'
                             '\nDepth threshold disabled.')
            options.depth_threshold = None
        sys.stdout.write('\n'+' '.join(sys.argv)+'\n')


def check_db():
    global options
    in_index = ''
    if options.in_fa_base_priority or options.in_fa_base:
        time0 = time.time()
        if options.in_fa_base_priority:
            fasta_file = options.in_fa_base_priority
        else:
            fasta_file = options.in_fa_base
        try:
            makedb_result = subprocess.getstatusoutput('makeblastdb -dbtype nucl -in '+fasta_file+' -out '+fasta_file+'.index')
        except AttributeError:
            makedb_result = commands.getstatusoutput('makeblastdb -dbtype nucl -in ' + fasta_file + ' -out ' + fasta_file + '.index')
        if 'Error' in str(makedb_result[1]) or 'error' in str(makedb_result[1]) or '不是内部或外部命令' in str(makedb_result[1]):
            os.system('makeblastdb -dbtype nucl -in '+fasta_file+' -out '+fasta_file+'.index')
            if not os.path.exists(fasta_file+'.index.nhr'):
                sys.stdout.write('\nBlast terminated with following info:\n'+str(makedb_result[1]))
                exit()
        in_index = fasta_file+'.index'
        sys.stdout.write('\nmake blast databse cost '+str(time.time()-time0))
    elif options.in_nc_base_priority:
        in_index = options.in_nc_base_priority
    elif options.in_nc_base:
        in_index = options.in_nc_base
    ex_index = ''
    if options.ex_fa_base_priority or options.ex_fa_base:
        time0 = time.time()
        if options.ex_fa_base_priority:
            fasta_file = options.ex_fa_base_priority
        else:
            fasta_file = options.ex_fa_base
        try:
            makedb_result = subprocess.getstatusoutput('makeblastdb -dbtype nucl -in '+fasta_file+' -out '+fasta_file+'.index')
        except AttributeError:
            makedb_result = commands.getstatusoutput('makeblastdb -dbtype nucl -in ' + fasta_file + ' -out ' + fasta_file + '.index')
        if 'Error' in str(makedb_result[1]) or 'error' in str(makedb_result[1]):
            sys.stdout.write('\nBlast terminated with following info:\n'+str(makedb_result[1]))
            exit()
        ex_index = fasta_file+'.index'
        sys.stdout.write('\nmake blast databse cost '+str(time.time()-time0))
    elif options.ex_nc_base_priority:
        ex_index = options.ex_nc_base_priority
    elif options.ex_nc_base:
        ex_index = options.ex_nc_base
    return in_index, ex_index


def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    interleaved = 0
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = fasta_file.readline()
            seq_line_count = 0
            while this_line and not this_line.startswith('>'):
                if seq_line_count == 1:
                    interleaved = len(this_seq)
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
                seq_line_count += 1
            seqs.append(this_seq)
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs, interleaved]


def write_fasta(out_dir, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_dir):
            out_dir = '.'.join(out_dir.split('.')[:-1])+'_.'+out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
                j += matrix[2]
            fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
    else:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            fasta_file.write(matrix[1][i]+'\n')
    fasta_file.close()


def make_new_matrix_with_names(names, old_matrix):
    i = 0
    while i < len(old_matrix[0]):
        if old_matrix[0][i] in names:
            i += 1
        else:
            del old_matrix[0][i]
            del old_matrix[1][i]
    return old_matrix


def blast_and_call_names(fasta_file, index_files, out_file):
    global options
    if index_files:
        time0 = time.time()
        sys.stdout.write('\nblast ...')
        if options.in_fastg_file:
            fasta_file += '.Temp'
        try:
            blast_result = subprocess.getstatusoutput('blastn -num_threads 4 -query ' + fasta_file + ' -db ' + index_files + ' -out ' + out_file + ' -outfmt 6 -evalue 1e-15')
        except AttributeError:
            blast_result = commands.getstatusoutput('blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
        # blastn -num_threads 4 -query assembly_graph.fastg -db db_f -out out_f -outfmt 6
        if 'Error' in str(blast_result[1]) or 'error' in str(blast_result[1]) or '不是内部或外部命令' in str(blast_result[1]):
            os.system('blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
            if not os.path.exists(out_file):
                sys.stdout.write('\nBlast terminated with following info:\n'+str(blast_result[1]))
                exit()
        # windows
        if not os.path.exists(out_file):
            os.system('blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
        time1 = time.time()
        sys.stdout.write('\nblast to '+os.path.split(index_files)[-1]+' cost '+str(time1-time0))
        names = {}
        try:
            blast_out_lines = open(out_file, 'rU')
        except IOError:
            sys.stdout.write('\nBlast was not properly installed or configurated.')
            exit()
        for line in blast_out_lines:
            line_split = line.split('\t')
            if options.in_fastg_file:
                query, template = line_split[0].split('_')[1], line_split[1]
            else:
                query, template = line_split[0], line_split[1]
            q_start, q_end = int(line_split[6]), int(line_split[7])
            q_min, q_max = min(q_start, q_end), max(q_start, q_end)
            if query in names:
                if template not in names[query]:
                    names[query][template] = [(q_min, q_max)]
                else:
                    i = 0
                    while i < len(names[query][template]):
                        this_min, this_max = names[query][template][i]
                        if q_max < this_min:
                            break
                        elif q_min > this_max:
                            i += 1
                            continue
                        else:
                            q_min = min(q_min, this_min)
                            q_max = min(q_max, this_max)
                            del names[query][template][i]
                    names[query][template].insert(i, (q_min, q_max))
            else:
                names[query] = {}
                names[query][template] = [(q_min, q_max)]
        sys.stdout.write('\nparse blast result cost'+str(time.time()-time1))
        return names
    else:
        return {}


def map_names(come_include, come_exclude, candidates):
    global options, short_candidates
    time0 = time.time()
    here_include = copy.deepcopy(come_include)
    here_exclude = copy.deepcopy(come_exclude)
    if options.in_fastg_file:
        if options.treat_no_hits == 'ex_no_con':
            short_connections = {}
            for candidate in candidates:
                this_short = candidate.split('_')[1]
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
            # expand here_include
            new_nodes = 1
            old_nodes = copy.deepcopy(here_include)
            while new_nodes:
                new_nodes = set()
                for in_node in list(old_nodes):
                    for new_node in short_connections[in_node]:
                        if new_node not in here_include and new_node not in here_exclude:
                            here_include.add(new_node)
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
    if options.ex_nc_base_priority or options.ex_fa_base_priority:
        if options.treat_no_hits in {"ex_no_con", "ex_no_hit"}:
            for this_short, full_name in short_candidates.items():
                if this_short in here_exclude:
                    pass
                elif this_short in here_include:
                    for name in full_name:
                        accepted.add(name)
                else:
                    pass
        else:
            for this_short, full_name in short_candidates.items():
                if this_short in here_exclude:
                    pass
                else:
                    for name in full_name:
                        accepted.add(name)
    elif options.in_nc_base_priority or options.in_fa_base_priority or options.in_nc_base or options.in_fa_base:
        if options.treat_no_hits in {"ex_no_con", "ex_no_hit"}:
            for this_short, full_name in short_candidates.items():
                if this_short in here_include:
                    for name in full_name:
                        accepted.add(name)
        else:
            for this_short, full_name in short_candidates.items():
                if this_short in here_include:
                    for name in full_name:
                        accepted.add(name)
                elif this_short in here_exclude:
                    pass
                else:
                    for name in full_name:
                        accepted.add(name)
    elif options.ex_nc_base or options.ex_fa_base:
        for this_short, full_name in short_candidates.items():
            if this_short in here_exclude:
                pass
            else:
                for name in full_name:
                    accepted.add(name)
    else:
        accepted = set(candidates)
    sys.stdout.write('\nmap names cost '+str(time.time()-time0))
    return accepted


def filter_fastg_by_depth():
    global options
    if options.in_fastg_file:
        in_fasta = options.in_fastg_file
    else:
        in_fasta = options.in_fasta_file
    depth = options.depth_threshold
    if depth and float(depth):
        depth = float(depth)
        time0 = time.time()
        fastg_matrix = read_fasta(in_fasta)
        new_fastg_matrix = [[], [], fastg_matrix[2]]
        for i in range(len(fastg_matrix[0])):
            if float(fastg_matrix[0][i].split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0]) >= depth:
                new_fastg_matrix[0].append(fastg_matrix[0][i])
                new_fastg_matrix[1].append(fastg_matrix[1][i])
        out_fasta = ''.join(in_fasta.split('.')[:-1]) + '.depth' + str(depth) + '.' + in_fasta.split('.')[-1]
        write_fasta(out_dir=out_fasta, matrix=new_fastg_matrix, overwrite=True)
        sys.stdout.write('\nfilter by depth cost '+str(time.time()-time0))
        return out_fasta
    else:
        return in_fasta


def del_complementary(fastg_file):
    global options
    if options.in_fastg_file:
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
        sys.stdout.write('\ndel complementary cost '+str(time.time()-time0))


def write_hits_csv_for_bandage(in_names, include_file, ex_names, exclude_file, out_file, overwrite):
    global options, short_candidates
    if options.no_csv:
        return ''
    else:
        time0 = time.time()
        if not overwrite:
            while os.path.exists(out_file+'.csv'):
                out_file += '_'
        out_file += '.csv'
        out_lines = ['EDGE\tdatabase\tloci\tloci_gene_sequential\tloci_sequential\tdetails\n']
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
            if options.in_fastg_file:
                full_name = short_candidates[edge]
                full_name.sort()
                if ':' in full_name[1]:
                    next_edges = [x.rstrip('\'').split('_') for x in full_name[1].split(':')[1].split(',')]
                    next_edges.sort(key=lambda x: -float(x[5].rstrip(';').rstrip('\'')))
                    postfix = '>>' + next_edges[0][1]
            this_string += '\t'+'>>'.join([x[2] for x in loci if x[2] != 'noncoding'])+postfix
            this_string += '\t'+'>>'.join([x[2] for x in loci])+postfix
            this_string += '\t'+'>>'.join([x[2]+'('+str(x[0])+'-'+str(x[1])+','+x[3]+')' for x in loci])+postfix
            this_string += '\n'
            out_lines.append(this_string)
        open(out_file, 'w').writelines(out_lines)
        sys.stdout.write('\ncreate csv cost '+str(time.time()-time0))


def remove_temp_files(fastg_file):
    global options
    if not options.keep_temp:
        if options.in_fastg_file:
            os.remove(fastg_file+'.Temp')
        try:
            os.remove(fastg_file+'.blast_in')
        except OSError:
            pass
        try:
            os.remove(fastg_file+'.blast_ex')
        except OSError:
            pass


__version__ = "1.7.01"


def main():
    time0 = time.time()
    sys.stdout.write("\nThis is a script for filtering spades assembly_graph.fastg contigs by blast\n")
    require_commands()
    global options
    # prepare fasta file
    fasta_file = filter_fastg_by_depth()
    del_complementary(fasta_file)
    # make blast database if not made
    include_index, exclude_index = check_db()
    # make blast
    in_names = blast_and_call_names(fasta_file=fasta_file, index_files=include_index, out_file=fasta_file+'.blast_in')
    ex_names = blast_and_call_names(fasta_file=fasta_file, index_files=exclude_index, out_file=fasta_file+'.blast_ex')
    # write out fasta according to blast
    fasta_matrix = read_fasta(fasta_dir=fasta_file)
    accept_names = map_names(come_include=set(in_names), come_exclude=set(ex_names), candidates=fasta_matrix[0])
    in_ex_info = 'only'*int(options.treat_no_hits=='ex_no_hit')+'extend'*int(options.treat_no_hits=='ex_no_con')+'+'+os.path.split(include_index)[-1]+'-'+os.path.split(exclude_index)[-1]
    fasta_matrix = make_new_matrix_with_names(names=accept_names, old_matrix=fasta_matrix)
    write_fasta(out_dir=fasta_file+'.'+in_ex_info+'.'+fasta_file.split('.')[-1], matrix=fasta_matrix, overwrite=False)
    # write out hits csv according to blast
    write_hits_csv_for_bandage(in_names=in_names, include_file=include_index, ex_names=ex_names, exclude_file=exclude_index, out_file=fasta_file+'.'+in_ex_info+'.', overwrite=False)
    remove_temp_files(fasta_file)
    sys.stdout.write('\nTotal cost: '+str(time.time()-time0)+'\n\n')

if __name__ == '__main__':
    main()
