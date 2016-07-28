#!/usr/bin/env python
# coding: utf8
import time
import os
import sys
import commands
from optparse import OptionParser, OptionGroup
import copy


options = ''
short_candidates = {}


def require_commands():
    global options
    blast_in_path = commands.getstatusoutput('blastn')
    if 'command not found' in blast_in_path[1]:
        print 'Error: blast not in the path!'
        os._exit(0)
    makedb_in_path = commands.getstatusoutput('makeblastdb')
    if 'command not found' in makedb_in_path[1]:
        print 'Error: makeblastdb not in the path!'
        os._exit(0)
    usage = 'python '+str(os.path.basename(__file__))
    parser = OptionParser(usage=usage)
    parser.add_option('-g', dest='in_fastg_file', help='followed by your input fastg file')
    parser.add_option('-f', dest='in_fasta_file', help='followed by your input fasta file')
    # parser.add_option('-o', dest='out_fastg_file', help='Output file')
    # filters
    parser.add_option('--depth-threshold', dest='depth_threshold',
                      help='Input a float or integer number. filter fastg file by depth. Default: no threshold.')
    parser.add_option('--include-index', dest='in_nc_base',
                      help='followed by Blast index format')
    parser.add_option('--include-index-priority', dest='in_nc_base_priority',
                      help='followed by Blast index format')
    parser.add_option('--exclude-index', dest='ex_nc_base',
                      help='followed by Blast index format')
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
    parser.add_option('--exclude-no-hits', dest='ex_no_hits', default=False, action='store_true', help='Choose to exclude all contigs that mismatch the "include index"')
    parser.add_option('--exclude-no-con', dest='ex_no_connect', default=False, action='store_true', help='Choose to exclude all contigs that disconnect with "include index" contigs')
    parser.add_option('--no-hits-labeled-csv', dest='no_csv', default=False, action='store_true', help='Choose to disable producing csv file')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true', help='Choose to disable deleting temp files produced by blast and this script')
    try:
        (options, args) = parser.parse_args()
    except Exception as e:
        print '\n######################################', e
        print '\n"-h" for more usage'
        os._exit(0)
    else:
        priority_chosen = int(bool(options.in_nc_base_priority))+int(bool(options.ex_nc_base_priority))+int(bool(options.in_fa_base_priority))+int(bool(options.ex_fa_base_priority))
        if priority_chosen > 1:
            print '\nLogical Error: only one option with "-priority" allowed!'
            os._exit(0)
        in_chosen = int(bool(options.in_nc_base_priority))+int(bool(options.in_nc_base))+int(bool(options.in_fa_base_priority))+int(bool(options.in_fa_base))
        if in_chosen > 1:
            print '\nOption Error: you can not simultaneously choose two "--include-*" options!'
            os._exit(0)
        ex_chosen = int(bool(options.ex_nc_base_priority))+int(bool(options.ex_nc_base))+int(bool(options.ex_fa_base_priority))+int(bool(options.ex_fa_base))
        if ex_chosen > 1:
            print '\nOption Error: you can not simultaneously choose two "--exclude-*" options!'
            os._exit(0)
        if in_chosen == 1 and ex_chosen == 1 and priority_chosen == 0:
            print '\nOption Error: since you have include and exclude chosen, one of them should be assigned priority!'
            os._exit(0)
        if ex_chosen == 1 and in_chosen == 0 and (options.ex_no_hits or options.ex_no_connect):
            print '\nOption Error: no contigs survive according to you choice!'
            os._exit(0)
        if options.ex_no_hits and options.ex_no_connect:
            print '\nOption Error: you can not simultaneously choose --exclude-no-hits and --exclude-no-con!'
            os._exit(0)
        if not options.in_fastg_file and options.ex_no_connect:
            print '\nOption Error: ex_no_connect is available only when you input a fastg file!'
            os._exit(0)
        if int(bool(options.in_fastg_file))+int(bool(options.in_fasta_file)) != 1:
            print '\nInput Error: you must choose one input fasta or fastg file!'
            os._exit(0)
        if int(bool(options.depth_threshold))+int(bool(options.in_fasta_file)) > 1:
            print '\nOption Warning: you can use depth threshold only when you input a fastg file!' \
                  '\nDepth threshold disabled.'
            options.depth_threshold = None
        print ' '.join(sys.argv)+'\n'


def check_db():
    global options
    in_index = ''
    if options.in_fa_base_priority or options.in_fa_base:
        time0 = time.time()
        if options.in_fa_base_priority:
            fasta_file = options.in_fa_base_priority
        else:
            fasta_file = options.in_fa_base
        makedb_result = commands.getstatusoutput('makeblastdb -dbtype nucl -in '+fasta_file+' -out '+fasta_file+'.index')
        if 'Error' in makedb_result[1] or 'error' in makedb_result[1] or '不是内部或外部命令' in makedb_result[1]:
            os.system('makeblastdb -dbtype nucl -in '+fasta_file+' -out '+fasta_file+'.index')
            if not os.path.exists(fasta_file+'.index.nhr'):
                print 'Blast terminated with following info:\n'+makedb_result[1]
                os._exit(0)
        in_index = fasta_file+'.index'
        print 'make blast databse cost', time.time()-time0
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
        makedb_result = commands.getstatusoutput('makeblastdb -dbtype nucl -in '+fasta_file+' -out '+fasta_file+'.index')
        if 'Error' in makedb_result[1] or 'error' in makedb_result[1]:
            print 'Blast terminated with following info:\n'+makedb_result[1]
            os._exit(0)
        ex_index = fasta_file+'.index'
        print 'make blast databse cost', time.time()-time0
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
    fasta_file = open(out_dir, 'wb')
    if matrix[2]:
        for i in xrange(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
                j += matrix[2]
            fasta_file.write(matrix[1][i][(j-matrix[2]):j]+'\n')
    else:
        for i in xrange(len(matrix[0])):
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
        print 'blast ...'
        if options.in_fastg_file:
            fasta_file += '.Temp'
        blast_result = commands.getstatusoutput(
            'blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
        # blastn -num_threads 4 -query assembly_graph.fastg -db db_f -out out_f -outfmt 6
        if 'Error' in blast_result[1] or 'error' in blast_result[1] or '不是内部或外部命令' in blast_result[1]:
            os.system('blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
            if not os.path.exists(out_file):
                print 'Blast terminated with following info:\n'+blast_result[1]
                os._exit(0)
        # windows
        if not os.path.exists(out_file):
            os.system('blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6 -evalue 1e-15')
        time1 = time.time()
        print 'blast to '+os.path.split(index_files)[-1]+' cost', time1-time0
        names = {}
        try:
            blast_out_lines = open(out_file, 'rU')
        except IOError:
            print 'Blast was not properly installed or configurated.'
            os._exit(0)
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
        print 'parse blast result cost', time.time()-time1
        return names
    else:
        return {}


def map_names(come_include, come_exclude, candidates):
    global options, short_candidates
    time0 = time.time()
    here_include = copy.deepcopy(come_include)
    here_exclude = copy.deepcopy(come_exclude)
    if options.in_fastg_file:
        if options.ex_no_connect:
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
        if options.ex_no_hits or options.ex_no_connect:
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
        if options.ex_no_hits or options.ex_no_connect:
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
    print 'map names cost', time.time()-time0
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
        for i in xrange(len(fastg_matrix[0])):
            if float(fastg_matrix[0][i].split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0]) >= depth:
                new_fastg_matrix[0].append(fastg_matrix[0][i])
                new_fastg_matrix[1].append(fastg_matrix[1][i])
        out_fasta = ''.join(in_fasta.split('.')[:-1]) + '.depth' + str(depth) + '.' + in_fasta.split('.')[-1]
        write_fasta(out_dir=out_fasta, matrix=new_fastg_matrix, overwrite=True)
        print 'filter by depth cost', time.time()-time0
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
        print 'del complementary cost', time.time()-time0


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
        open(out_file, 'wb').writelines(out_lines)
        print 'create csv cost', time.time()-time0


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


def main():
    time0 = time.time()
    print "\nThis is a script for filtering spades assembly_graph.fastg contigs by blast\n"
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
    in_ex_info = 'only'*int(options.ex_no_hits)+'extend'*int(options.ex_no_connect)+'+'+os.path.split(include_index)[-1]+'-'+os.path.split(exclude_index)[-1]
    fasta_matrix = make_new_matrix_with_names(names=accept_names, old_matrix=fasta_matrix)
    write_fasta(out_dir=fasta_file+'.'+in_ex_info+'.'+fasta_file.split('.')[-1], matrix=fasta_matrix, overwrite=False)
    # write out hits csv according to blast
    write_hits_csv_for_bandage(in_names=in_names, include_file=include_index, ex_names=ex_names, exclude_file=exclude_index, out_file=fasta_file+'.'+in_ex_info+'.', overwrite=False)
    remove_temp_files(fasta_file)
    print 'Total cost:', time.time()-time0

if __name__ == '__main__':
    main()
