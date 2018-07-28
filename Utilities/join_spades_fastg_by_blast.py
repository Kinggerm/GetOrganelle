#!/usr/bin/env python
# coding:utf8
import time
import os
import sys
import platform
import subprocess
try:
    # python2
    import commands
except:
    pass
from optparse import OptionParser, OptionGroup
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *


# V1_4

this_dir_split = '/'
if 'Win' in platform.architecture()[1]:
    this_dir_split = '\\'
options = ''
short_candidates = {}
try:
    # python2
    import string
    translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")

    def complementary_seq(input_seq):
        return string.translate(input_seq, translator)[::-1]
except AttributeError:
    # python3
    translator = str.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")

    def complementary_seq(input_seq):
        return str.translate(input_seq, translator)[::-1]


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
    usage = 'python '+str(os.path.basename(__file__))+' -g input.fastg -f refernce.fasta'
    parser = OptionParser(usage=usage)
    parser.add_option('-g', dest='in_fastg_file', help='followed by your input fastg file')
    parser.add_option('-f', dest='reference_fa_base', help='followed by Fasta index format')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true', help='Choose to disable deleting temp files produced by blast and this script')
    parser.add_option('--bt', dest='blast_hits_threshold', default=0.60, help='Default: 0.60', type=float)
    parser.add_option('--max-gap', dest='max_gap_to_add', default=1500, help='Default: 1500', type=int)
    parser.add_option('--con-all', dest='connect_inner_contig', default=False, action='store_true', help='Choose to activate connecting all possible contigs. Default: False')
    parser.add_option('--depth', dest='depth_to_connect', default=1.0, help='Default: 1.0', type=float)
    # parser.add_option('--merge-overlaps', default=False, action='store_true', help='Choose to activate automatically merging overlapping contigs')
    # parser.add_option('--min-os', dest='min_overlap_similarity', default=0.9, help='The similarity threshold to merge overlapping contigs. Default: 0.9', type=float)
    # parser.add_option('--min-ol', dest='min_overlap_length', default=15, help='The length threshold to merge overlapping contigs. Default: 15', type=int)
    try:
        (options, args) = parser.parse_args()
    except Exception as e:
        sys.stdout.write('\n######################################'+str(e))
        sys.stdout.write('\n"-h" for more usage')
        exit()


def check_db():
    global options
    if options.reference_fa_base:
        time0 = time.time()
        ref_fasta = read_fasta(options.reference_fa_base)
        if len(ref_fasta[0]) > 1:
            options.reference_fa_base += '.1st.fasta'
            write_fasta(out_dir=options.reference_fa_base, matrix=[[ref_fasta[0][0]], [ref_fasta[1][0]], ref_fasta[2]], overwrite=True)
            sys.stdout.write('\nWarning: multi-seqs in reference file, only use the 1st sequence.')
        elif len(ref_fasta[0]) == 0:
            sys.stdout.write('\nError: illegal reference file!')
            exit()
        try:
            # python2
            makedb_result = subprocess.getstatusoutput('makeblastdb -dbtype nucl -in '+options.reference_fa_base+' -out '+options.reference_fa_base+'.index')
        except AttributeError:
            # python3
            makedb_result = commands.getstatusoutput('makeblastdb -dbtype nucl -in ' + options.reference_fa_base + ' -out ' + options.reference_fa_base + '.index')
        if 'Error' in str(makedb_result[1]) or 'error' in str(makedb_result[1]) or '不是内部或外部命令' in str(makedb_result[1]):
            os.system('makeblastdb -dbtype nucl -in '+options.reference_fa_base+' -out '+options.reference_fa_base+'.index')
            if not os.path.exists(options.reference_fa_base+'.index.nhr'):
                sys.stdout.write('Blast terminated with following info:\n'+str(makedb_result[1]))
                exit()
        in_index = options.reference_fa_base+'.index'
        sys.stdout.write('\nMaking BLAST db cost '+str(time.time()-time0))
    else:
        sys.stdout.write('\nError: No reference input!')
        exit()
    return in_index


def blast_and_call_new_matrix(fasta_file, index_files, out_file, len_db):
    global options
    time0 = time.time()
    sys.stdout.write('\nMaking BLAST ...')
    fasta_file += '.Temp'
    try:
        blast_result = subprocess.getstatusoutput(
        'blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6')
    except AttributeError:
        blast_result = commands.getstatusoutput(
            'blastn -num_threads 4 -query ' + fasta_file + ' -db ' + index_files + ' -out ' + out_file + ' -outfmt 6')
    if 'Error' in str(blast_result[1]) or 'error' in str(blast_result[1]) or '不是内部或外部命令' in str(blast_result[1]):
        sys.stdout.write('\nBlast terminated with following info:\n'+str(blast_result[1]))
        exit()
    # windows
    if not os.path.exists(out_file):
        os.system('blastn -num_threads 4 -query '+fasta_file+' -db '+index_files+' -out '+out_file+' -outfmt 6')
    time1 = time.time()
    sys.stdout.write('\nBLAST to '+index_files.split(this_dir_split)[-1]+' cost '+str(time1-time0))
    # ----------------------------------------
    # find start and end points of query
    # initialize candidates: fastq topologies and sequences
    query_matrix = read_fasta(options.in_fastg_file)
    len_fastg = len(query_matrix[0])
    hits_candidates = {}
    short_names = []
    for i in range(len_fastg):
        full_name = query_matrix[0][i]
        short_name = '_'.join(full_name.split()[0].split('_')[1:]).split('_length')[0]
        coverage = float(full_name.split('cov_')[1].split(';')[0].split('\'')[0].split(':')[0])
        hits_candidates[short_name] = {False: set(), True: set(), 'coverage': coverage}
        short_names.append(short_name)
    for i in range(len_fastg):
        full_name = query_matrix[0][i]
        short_name = short_names[i]
        connected_edges = set()
        if ':' in full_name:
            for edge in full_name.rstrip(';').split(':')[1].split(','):
                edge_short_name = '_'.join(edge.split('_')[1:]).split('_length')[0]
                if edge_short_name in hits_candidates:
                    if edge.endswith('\''):
                        connected_edges.add((edge_short_name, False))
                    else:
                        connected_edges.add((edge_short_name, True))
        if full_name.split(';')[0].split(':')[0].endswith('\''):
            sequence = query_matrix[1][i]
            len_seq = len(sequence)
            new_items = {'identity': [0 for j in range(len_seq)],
                         ('index', False): i,
                         ('seq', False): sequence,
                         ('seq', True): complementary_seq(sequence),
                         'len_seq': len_seq,
                         False: connected_edges}
            hits_candidates[short_name].update(new_items)
        else:
            sequence = query_matrix[1][i]
            len_seq = len(sequence)
            new_items = {'identity': [0 for j in range(len_seq)],
                         'start_block': {'q': (len_seq, len_seq), 'r':[]},
                         'end_block': {'q': (0, 0), 'r': []},
                         ('index', True): i,
                         ('seq', True): sequence,
                         ('seq', False): complementary_seq(sequence),
                         'len_seq': len_seq,
                         True: connected_edges}
            hits_candidates[short_name].update(new_items)

    # -----------------------------------
    # detect k-mer
    k_mer = 0
    try:
        for short_name in hits_candidates:
            for direction in [True, False]:
                for next_edge_info in hits_candidates[short_name][direction]:
                    if k_mer:
                        if hits_candidates[short_name][('seq', direction)][-k_mer:] != hits_candidates[next_edge_info[0]][('seq', next_edge_info[1])][:k_mer]:
                            raise ValueError
                    else:
                        for k_mer in range(127, 19, -2):
                            if hits_candidates[short_name][('seq', direction)][-k_mer:] == hits_candidates[next_edge_info[0]][('seq', next_edge_info[1])][:k_mer]:
                                break
                        else:
                            raise ValueError
    except ValueError:
        k_mer = 0
        pass
    sys.stdout.write('\nDetected k-mer:'+str(k_mer))

    # calculate edge connections according to hits_candidates and max_gap
    #
    # wait to improve:
    # miss the directions for jointed edges!
    def get_jointed_edges_within_distance(all_infos, this_edge, this_direction, length_left, jointed_edges, k_mer, recurse_depth=0):
        for this_next_edge in all_infos[this_edge][this_direction]:
            this_length_left = length_left - all_infos[this_next_edge[0]]['len_seq'] + k_mer
            if this_length_left >= 0 and this_next_edge not in jointed_edges:
                # try:
                # arbitrarily set recurse_depth to 20
                if recurse_depth < 20:
                    jointed_edges = get_jointed_edges_within_distance(all_infos, this_next_edge[0], this_direction==this_next_edge[1], this_length_left, jointed_edges, k_mer, recurse_depth+1)
                # except RuntimeError:
                #     sys.stdout.write('\nWarning: RuntimeError!')
                #     pass
            jointed_edges.add(this_next_edge)
        return jointed_edges
    edge_connections = {}
    for edge in hits_candidates:
        for direction in [False, True]:
            edge_connections[(edge, direction)] = get_jointed_edges_within_distance(hits_candidates, edge, direction, options.max_gap_to_add+k_mer, set(), k_mer)

    # compare candidates with blast results
    blast_out_lines = open(out_file, 'rU')
    for line in blast_out_lines:
        line_split = line.split('\t')
        query = '_'.join(line_split[0].split('_')[1:]).split('_length')[0]
        q_start, q_end = int(line_split[6]), int(line_split[7])
        r_start, r_end = int(line_split[8]), int(line_split[9])
        identity = float(line_split[2])
        for i in range(q_start-1, q_end):
            hits_candidates[query]['identity'][i] = max(identity, hits_candidates[query]['identity'][i])
        if q_start < hits_candidates[query]['start_block']['q'][0]:
            hits_candidates[query]['start_block']['q'] = (q_start, q_end)
            hits_candidates[query]['start_block']['r'] = [(r_start, r_end)]
        elif q_start == hits_candidates[query]['start_block']['q'][0]:
            if q_end > hits_candidates[query]['start_block']['q'][1]:
                hits_candidates[query]['start_block']['q'] = (q_start, q_end)
                hits_candidates[query]['start_block']['r'] = [(r_start, r_end)]
            elif q_end == hits_candidates[query]['start_block']['q'][1]:
                hits_candidates[query]['start_block']['r'].append((r_start, r_end))
        if q_end > hits_candidates[query]['end_block']['q'][1]:
            hits_candidates[query]['end_block']['q'] = (q_start, q_end)
            hits_candidates[query]['end_block']['r'] = [(r_start, r_end)]
        elif q_end == hits_candidates[query]['end_block']['q'][1]:
            if q_start < hits_candidates[query]['end_block']['q'][0]:
                hits_candidates[query]['end_block']['q'] = (q_start, q_end)
                hits_candidates[query]['end_block']['r'] = [(r_start, r_end)]
            elif q_start == hits_candidates[query]['end_block']['q'][0]:
                hits_candidates[query]['end_block']['r'].append((r_start, r_end))
    blast_out_lines.close()
    time2 = time.time()
    sys.stdout.write('\nParsing BLAST result cost '+str(time2-time1))
    # ------------------------------------
    # map terminal blocks of candidates to reference bases
    # workout points to connect
    # {base: [(query name, query identity, is_start_of_query, direction_in_reference)]}
    ref_bases_dict = {}
    for hit in hits_candidates.keys():
        average_identity = sum(hits_candidates[hit]['identity'])/float(len(hits_candidates[hit]['identity']))
        hits_candidates[hit]['identity'] = average_identity
        if average_identity >= options.blast_hits_threshold:
            for block in ['start_block', 'end_block']:
                is_start_of_query = bool(block == 'start_block')
                if options.connect_inner_contig or not bool(hits_candidates[hit][not is_start_of_query]):
                    if hits_candidates[hit]['coverage'] >= options.depth_to_connect:
                        query_loci = hits_candidates[hit][block]['q']
                        if is_start_of_query:
                            length_to_terminal = query_loci[0] - 1
                        else:
                            length_to_terminal = hits_candidates[hit]['len_seq'] - query_loci[1]
                        for reference_block in hits_candidates[hit][block]['r']:
                            direction_in_ref = bool(bool(reference_block[0] <= reference_block[1]) == is_start_of_query)
                            ref_block_to_mark = int(not is_start_of_query)
                            if reference_block[ref_block_to_mark] in ref_bases_dict:
                                ref_bases_dict[reference_block[ref_block_to_mark]].append((hit, length_to_terminal, is_start_of_query, direction_in_ref))
                            else:
                                ref_bases_dict[reference_block[ref_block_to_mark]] = [(hit, length_to_terminal, is_start_of_query, direction_in_ref)]

    # ------------------------------------
    # search for new connections
    used_edge_numbers = []
    for crazy_string in list(hits_candidates):
        for numbers in ''.join(filter(lambda ch: ch in '0123456789-_', crazy_string)).split('_'):
            for num in numbers.split('-'):
                used_edge_numbers.append(int(num))
    used_edge_numbers.sort()
    variances_to_pass = {'edge': used_edge_numbers[-1]+1, 'index': len_fastg}

    def make_connections(edge1, base1, edge2, base2, k_mer):
        # if end to end and disable self-connection
        if edge1[3] != edge2[3] and edge1[0] != edge2[0]:
            # if not connected
            if (edge2[0], edge2[2]) not in edge_connections[(edge1[0], not edge1[2])]:
                # if Overlaps
                if edge1[3] or base1 == base2:
                    overlap_or_gap_length = (base2-base1)%len_db+1 + edge1[1] + edge2[1]
                    edge_name = str(variances_to_pass['edge'])+'overlap'+str(overlap_or_gap_length)
                    new_full_name = 'EDGE_'+edge_name+'_length_'+str(overlap_or_gap_length+2*k_mer)+'_cov_80'
                    forward_edge_sequence = hits_candidates[edge1[0]][('seq', not edge1[2])][-k_mer:] + '?'*overlap_or_gap_length + hits_candidates[edge2[0]][('seq', edge2[2])][:k_mer]
                    reverse_edge_sequence = hits_candidates[edge2[0]][('seq', not edge2[2])][-k_mer:] + '?'*overlap_or_gap_length + hits_candidates[edge1[0]][('seq', edge1[2])][:k_mer]
                else:
                    overlap_or_gap_length = (base2-base1)%len_db-1 - edge1[1] - edge2[1]
                    # if still overlaps
                    if overlap_or_gap_length < 0:
                        overlap_or_gap_length = -overlap_or_gap_length
                        edge_name = str(variances_to_pass['edge'])+'overlap'+str(overlap_or_gap_length)
                        new_full_name = 'EDGE_'+edge_name+'_length_'+str(overlap_or_gap_length+2*k_mer)+'_cov_20'
                        forward_edge_sequence = hits_candidates[edge1[0]][('seq', not edge1[2])][-k_mer:] + '?'*overlap_or_gap_length + hits_candidates[edge2[0]][('seq', edge2[2])][:k_mer]
                        reverse_edge_sequence = hits_candidates[edge2[0]][('seq', not edge2[2])][-k_mer:] + '?'*overlap_or_gap_length + hits_candidates[edge1[0]][('seq', edge1[2])][:k_mer]
                    # if Gaps
                    else:
                        edge_name = str(variances_to_pass['edge'])+'gap'+str(overlap_or_gap_length)
                        new_full_name = 'EDGE_'+edge_name+'_length_'+str(overlap_or_gap_length+2*k_mer)+'_cov_5'
                        forward_edge_sequence = hits_candidates[edge1[0]][('seq', not edge1[2])][-k_mer:] + 'N'*overlap_or_gap_length + hits_candidates[edge2[0]][('seq', edge2[2])][:k_mer]
                        reverse_edge_sequence = hits_candidates[edge2[0]][('seq', not edge2[2])][-k_mer:] + 'N'*overlap_or_gap_length + hits_candidates[edge1[0]][('seq', edge1[2])][:k_mer]
                variances_to_pass['edge'] += 1
                # these_directions = {'to_edge1':False, 'edge1':not edge1[2],'to_edge2':True, 'edge2':edge2[2]}
                # add new edge to matrix
                query_matrix[0].append(new_full_name+':'+query_matrix[0][hits_candidates[edge2[0]][('index', edge2[2])]].split(';')[0].split(':')[0]+';')
                edge2_full_name = query_matrix[0][hits_candidates[edge2[0]][('index', not edge2[2])]]
                if ':' in edge2_full_name:
                    query_matrix[0][hits_candidates[edge2[0]][('index', not edge2[2])]] = edge2_full_name.rstrip(';')+','+new_full_name+'\';'
                else:
                    query_matrix[0][hits_candidates[edge2[0]][('index', not edge2[2])]] = edge2_full_name.rstrip(';')+':'+new_full_name+'\';'
                query_matrix[0].append(new_full_name+'\':'+query_matrix[0][hits_candidates[edge1[0]][('index', edge1[2])]].split(';')[0].split(':')[0]+';')
                edge1_full_name = query_matrix[0][hits_candidates[edge1[0]][('index', not edge1[2])]]
                if ':' in edge1_full_name:
                    query_matrix[0][hits_candidates[edge1[0]][('index', not edge1[2])]] = edge1_full_name.rstrip(';')+','+new_full_name+';'
                else:
                    query_matrix[0][hits_candidates[edge1[0]][('index', not edge1[2])]] = edge1_full_name.rstrip(';')+':'+new_full_name+';'
                query_matrix[1].append(forward_edge_sequence)
                query_matrix[1].append(reverse_edge_sequence)
                # add new edge to hits_candidates
                hits_candidates[edge_name] = {('index', True): variances_to_pass['index'],
                                              ('index', False): variances_to_pass['index']+1,
                                              ('seq', True): forward_edge_sequence,
                                              ('seq', False): forward_edge_sequence,
                                              'len_seq': overlap_or_gap_length+2*k_mer,
                                              True: [(edge2[0], edge2[2])],
                                              False: [(edge1[0], edge1[2])]}
                variances_to_pass['index'] += 2
                hits_candidates[edge1[0]][not edge1[2]].add((edge_name, True))
                hits_candidates[edge2[0]][not edge2[2]].add((edge_name, False))
                # add new edge to edge_connections (update)
                edge_connections[(edge1[0], not edge1[2])] = get_jointed_edges_within_distance(hits_candidates, edge1[0], not edge1[2], options.max_gap_to_add+k_mer, set(), k_mer)
                edge_connections[(edge2[0], not edge2[2])] = get_jointed_edges_within_distance(hits_candidates, edge2[0], not edge2[2], options.max_gap_to_add+k_mer, set(), k_mer)
                edge_connections[(edge_name, True)] = get_jointed_edges_within_distance(hits_candidates, edge_name, True, options.max_gap_to_add+k_mer, set(), k_mer)
                edge_connections[(edge_name, False)] = get_jointed_edges_within_distance(hits_candidates, edge_name, False, options.max_gap_to_add+k_mer, set(), k_mer)

    ref_bases_list = sorted(list(ref_bases_dict))
    len_ref_base = len(ref_bases_list)
    for i in range(len_ref_base):
        candidates = ref_bases_dict[ref_bases_list[i]]
        # the same base
        len_candidates = len(candidates)
        if len_candidates >= 2:
            for k in range(len_candidates):
                for l in range(1, len_candidates):
                    make_connections(candidates[k], ref_bases_list[i], candidates[l], ref_bases_list[i], k_mer)
        # next bases
        for candidate_infos in candidates:
            i_plus = i + 1
            base = ref_bases_list[i_plus % len_ref_base]
            while i_plus-i < len_ref_base and (base - ref_bases_list[i]) % len_db <= options.max_gap_to_add:
                for hit_infos in ref_bases_dict[base]:
                    make_connections(candidate_infos, ref_bases_list[i], hit_infos, base, k_mer)
                i_plus += 1
                base = ref_bases_list[i_plus%len_ref_base]
    sys.stdout.write('\nRedirecting contig path cost '+str(time.time()-time2))
    return query_matrix


def del_complementary(fastg_file):
    global options

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
    sys.stdout.write('\nDel complementary cost'+str(time.time()-time0))


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
            os.remove(options.reference_fa_base+'.index.nhr')
            os.remove(options.reference_fa_base+'.index.nin')
            os.remove(options.reference_fa_base+'.index.nsq')
        except OSError:
            pass


def main():
    time0 = time.time()
    sys.stdout.write(
        "\nThis script would join the spades fastg contigs according to the reference."
        "\nIt would add extra gap nodes and/or overlap nodes in between the worth connecting nodes in a fastg file."
        "\n"
        "\nHowever, this is a BETA version:"
        "\nAlthough it will not produce error connections, it usually replicates the same right connection."
        "\nDon't be surprised if you find any other bugs.\n")
    require_commands()
    global options
    # fastg to fasta
    fasta_file = options.in_fastg_file
    del_complementary(fasta_file)
    # make blast database if not made
    include_index = check_db()
    len_db = len(read_fasta(options.reference_fa_base)[1][0])
    # make blast
    new_fasta_matrix = blast_and_call_new_matrix(fasta_file=fasta_file, index_files=include_index, out_file=fasta_file + '.blast_in', len_db=len_db)
    # write out fastg
    write_fasta(out_dir=fasta_file+'.Ncontigs_added.'+fasta_file.split('.')[-1], matrix=new_fasta_matrix, overwrite=False)
    remove_temp_files(fasta_file)
    sys.stdout.write('\n\nTotal cost: '+str(time.time()-time0)+'\n\n')


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