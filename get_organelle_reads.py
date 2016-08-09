#!/usr/bin/env python
import time, datetime
import sys
import os
import subprocess
try:
    import commands
except:
    pass
import logging
from optparse import OptionParser, OptionGroup

path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
word_size = int
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

try:
    import psutil
except ImportError:
    this_process = None
else:
    this_process = psutil.Process(os.getpid())


def chop_seqs(seq_generator_or_list):
    return_words = set()
    for seed in seq_generator_or_list:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            for i in range(0, this_seq_len-word_size+1):
                forward = seed[i:i+word_size]
                return_words.add(forward)
                reverse = cpt_seed[i:i+word_size]
                return_words.add(reverse)
    return return_words


def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = fasta_file.readline()
            while this_line and not this_line.startswith('>'):
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
            seqs.append(this_seq)
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs]


def read_self_fq_seq_generator(fq_dir_list):
    if not ((type(fq_dir_list) is list) or (type(fq_dir_list) is tuple)):
        fq_dir_list = [fq_dir_list]
    for fq_dir in fq_dir_list:
        count = 0
        for fq_line in open(fq_dir, 'rU'):
            if count % 4 == 1:
                yield fq_line
            count += 1


def write_fq_results(original_fq_dir, accepted_contig_id, pair_end_out, out_file_name, temp2_clusters_dir, temp3_pair_dir, fastq_indices_in_memory, verbose, options, log):

    if verbose:
        sys.stdout.write(' ' * 100 + '\b' * 100)
        sys.stdout.flush()
        log.info("Producing output ...")
        log.info("reading indices ...")
    accepted_lines = []
    if options.index_in_memory:
        # read cluster indices
        for this_index in accepted_contig_id:
            accepted_lines += fastq_indices_in_memory[1][this_index]
        # produce the pair-end output
        accepted_lines = set(accepted_lines)
        if pair_end_out:
            dump_candidates = set()
            for candidate in accepted_lines:
                if fastq_indices_in_memory[2][candidate] not in accepted_lines:
                    dump_candidates.add(candidate)
            accepted_lines -= dump_candidates

    else:
        # read cluster indices
        temp2_indices_file_in = open(temp2_clusters_dir, 'rU')
        this_index = 0
        for line in temp2_indices_file_in:
            if this_index in accepted_contig_id:
                accepted_lines += [int(x) for x in line.strip().split('\t')]
            this_index += 1
        accepted_lines = set(accepted_lines)

        # produce the pair-end output
        if pair_end_out:
            pair_to_each_other = dict()
            temp3_indices_file_in = open(temp3_pair_dir, 'rU')
            line_count = 0
            for line in temp3_indices_file_in:
                if line_count in accepted_lines:
                    pair_to_each_other[line_count] = int(line.strip())
                line_count += 4
            temp3_indices_file_in.close()
            dump_candidates = set()
            for candidate in accepted_lines:
                if pair_to_each_other[candidate] not in accepted_lines:
                    dump_candidates.add(candidate)
            accepted_lines -= dump_candidates
            del pair_to_each_other

    # write by line
    if verbose:
        log.info("writing fastq lines ...")
    post_reading = [open(original_fq_dir[0], 'rU'), open(original_fq_dir[1], 'rU')]
    file_out = [open(out_file_name+'_1.fq', 'w'), open(out_file_name+'_2.fq', 'w')]
    line_count = 0
    for i in range(2):
        line = post_reading[i].readline()
        while line:
            if line_count in accepted_lines:
                file_out[i].write(line)
                for j in range(3):
                    file_out[i].write(post_reading[i].readline())
                    line_count += 1
                line = post_reading[i].readline()
                line_count += 1
            else:
                for j in range(4):
                    line = post_reading[i].readline()
                    line_count += 1
        file_out[i].close()
        post_reading[i].close()
    del accepted_lines
    if verbose:
        log.info("writing fastq lines finished.")


def read_fq_and_pair_infos(original_fq_dir, pair_end_out, rm_duplicates, output_base, anti_lines, options, log):
    # read original reads
    # pair_to_each_other
    # line_cluster (list) ~ forward_reverse_reads
    line_clusters = []
    seq_duplicates = {}
    forward_reverse_reads = []
    line_count = 0
    this_index = 0
    #
    pair_to_each_other = {}
    name_to_line = {}
    #
    temp1_contig_dir = [os.path.join(output_base, k+'temp.indices.1') for k in ("_", "")]
    temp2_clusters_dir = [os.path.join(output_base, k+'temp.indices.2') for k in ("_", "")]
    temp3_pair_dir = [os.path.join(output_base, k+'temp.indices.3') for k in ("_", "")]
    if options.script_resume and os.path.exists(temp1_contig_dir[1]) and os.path.exists(temp2_clusters_dir[1]) and (os.path.exists(temp3_pair_dir[1]) or not options.pair_end_out):
        if options.pseudo_assembled or options.index_in_memory:
            log.info("Reading existed indices for fastq ...")
            #
            forward_reverse_reads = [x.strip() for x in open(temp1_contig_dir[1], 'rU')]
            #
            line_clusters = [[int(x) for x in y.split('\t')] for y in open(temp2_clusters_dir[1], 'rU')]
            # read pair infos
            if options.pair_end_out:
                line_count = 0
                for line in open(temp3_pair_dir[1], 'rU'):
                    pair_to_each_other[line_count] = int(line)
                    line_count += 4
            else:
                line_count = sum([len(x) for x in line_clusters])*4
            # log
            len_indices = len(line_clusters)
            if this_process:
                memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + ", "
            else:
                memory_usage = ''
            log.info(memory_usage + str(len_indices) + " unique in all " + str(line_count // 4) + " reads")
        else:
            log.info("indices for fastq existed!")
    else:
        if not options.index_in_memory:
            temp1_contig_out = open(temp1_contig_dir[0], 'w')
        pre_reading = [open(original_fq_dir[0], 'rU'), open(original_fq_dir[1], 'rU')]
        for file_in in pre_reading:
            line = file_in.readline()
            if options.bowtie2_anti_seed or options.anti_seed:
                while line:
                    if line.startswith("@"):
                        try:
                            if ' ' in line:
                                this_head = line[1:].split(' ')
                                this_name, direction = this_head[0], int(this_head[1][0])
                            elif '#' in line:
                                this_head = line[1:].split('#')
                                this_name, direction = this_head[0], int(this_head[1][0])
                            else:
                                this_name, direction = line[1:].strip(), 1
                        except (ValueError, IndexError):
                            log.error('Unrecognized fq format in '+str(line_count)+' '+str(line))
                            exit()
                        else:
                            if pair_end_out:
                                try:
                                    if this_name in name_to_line:
                                        pair_to_each_other[line_count] = name_to_line[this_name]
                                        pair_to_each_other[name_to_line[this_name]] = line_count
                                    else:
                                        name_to_line[this_name] = line_count
                                    if this_name in anti_lines:
                                        line_count += 4
                                        for i in range(4):
                                            line = file_in.readline()
                                        continue
                                except KeyError:
                                    log.error('Unknown error in format. Please try to disable "--pair_end_out".')
                                    exit()
                            else:
                                if (this_name, direction) in anti_lines:
                                    line_count += 4
                                    for i in range(4):
                                        line = file_in.readline()
                                    continue
                        this_seq = file_in.readline().strip()
                        # drop illegal reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue
                        this_c_seq = complementary_seq(this_seq)
                        if rm_duplicates:
                            if this_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_seq]].append(line_count)
                            elif this_c_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_c_seq]].append(line_count)
                            else:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                                if not options.index_in_memory:
                                    temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                                seq_duplicates[this_seq] = this_index
                                line_clusters.append([line_count])
                                this_index += 1
                        else:
                            line_clusters.append([line_count])
                            if options.index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            else:
                                temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                    else:
                        log.error("Illegal fq format in line "+str(line_count)+' '+str(line))
                        exit()
                    if line_count % 54321 == 0:
                        to_print = str("%s"%datetime.datetime.now())[:23].replace('.', ',')+" - INFO: "+str((line_count+4)//4)+" reads"
                        sys.stdout.write(to_print+'\b'*len(to_print))
                        sys.stdout.flush()
                    line_count += 1
                    for i in range(3):
                        line = file_in.readline()
                        line_count += 1
            else:
                while line:
                    if line.startswith("@"):
                        if pair_end_out:
                            try:
                                if ' ' in line:
                                    this_head = line[1:].split(' ')
                                    this_name, direction = this_head[0], int(this_head[1][0])
                                elif '#' in line:
                                    this_head = line[1:].split('#')
                                    this_name, direction = this_head[0], int(this_head[1][0])
                                else:
                                    this_name, direction = line[1:].strip(), 1
                            except (ValueError, IndexError):
                                log.error('Unrecognized fq format.')
                                exit()
                            else:
                                try:
                                    if this_name in name_to_line:
                                        pair_to_each_other[line_count] = name_to_line[this_name]
                                        pair_to_each_other[name_to_line[this_name]] = line_count
                                    else:
                                        name_to_line[this_name] = line_count
                                except KeyError:
                                    log.error('Unknown error in format. Please try to disable "--pair_end_out".')
                                    exit()
                        this_seq = file_in.readline().strip()
                        # drop illegal reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue
                        this_c_seq = complementary_seq(this_seq)
                        if rm_duplicates:
                            if this_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_seq]].append(line_count)
                            elif this_c_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_c_seq]].append(line_count)
                            else:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                                if not options.index_in_memory:
                                    temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                                seq_duplicates[this_seq] = this_index
                                line_clusters.append([line_count])
                                this_index += 1
                        else:
                            line_clusters.append([line_count])
                            if options.index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            else:
                                temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                    else:
                        log.error("Illegal fq format in line "+str(line_count)+' '+str(line))
                        exit()
                    if line_count % 54321 == 0:
                        to_print = str("%s"%datetime.datetime.now())[:23].replace('.', ',')+" - INFO: "+str((line_count+4)//4)+" reads"
                        sys.stdout.write(to_print+'\b'*len(to_print))
                        sys.stdout.flush()
                    line_count += 1
                    for i in range(3):
                        line = file_in.readline()
                        line_count += 1
            file_in.close()
        if not options.index_in_memory:
            temp1_contig_out.close()
            os.rename(temp1_contig_dir[0], temp1_contig_dir[1])

        if this_process:
            memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + ", "
        else:
            memory_usage = ''

        del name_to_line

        if not options.index_in_memory:
            # dump line clusters
            len_indices = len(line_clusters)
            temp2_indices_file_out = open(temp2_clusters_dir[0], 'w')
            for this_index in range(len_indices):
                temp2_indices_file_out.write('\t'.join([str(x) for x in line_clusters[this_index]]))
                temp2_indices_file_out.write('\n')
            temp2_indices_file_out.close()
            os.rename(temp2_clusters_dir[0], temp2_clusters_dir[1])

            if pair_end_out:
                # dump pair_to_each_other to file to save memory
                # read in again in the last
                temp3_pair_file_out = open(temp3_pair_dir[0], 'w')
                for line_count_2 in range(0, line_count, 4):
                    temp3_pair_file_out.write(str(pair_to_each_other[line_count_2])+'\n')
                temp3_pair_file_out.close()
                os.rename(temp3_pair_dir[0], temp3_pair_dir[1])
                pair_to_each_other = None

        del seq_duplicates
        del pre_reading
        len_indices = len(line_clusters)
        log.info(memory_usage+str(len_indices)+" unique in all "+str(line_count//4)+" reads")

    return forward_reverse_reads, line_clusters, pair_to_each_other


def pseudo_assembly(fastq_indices_in_memory, dupli_threshold, log):
    global word_size
    forward_and_reverse_reads, line_clusters, pairs = fastq_indices_in_memory
    log.info("Pseudo-assembly reads...")
    lines_with_duplicates = {}
    count_dupli = 0
    for j in range(len(line_clusters)):
        if len(line_clusters[j]) >= 2:
            if count_dupli < dupli_threshold:
                lines_with_duplicates[j] = int
            count_dupli += 1
    if this_process:
        memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + ", "
    else:
        memory_usage = ''
    log.info(memory_usage+str(len(lines_with_duplicates))+"/"+str(count_dupli)+" used/duplicated")

    groups_of_duplicate_lines = {}
    count_groups = 0
    these_words = {}
    for this_unique_read_id in list(lines_with_duplicates):
        this_seq = forward_and_reverse_reads[2*this_unique_read_id]
        this_c_seq = forward_and_reverse_reads[2*this_unique_read_id+1]
        seq_len = len(this_seq)
        temp_length = seq_len - word_size
        these_group_id = set()
        this_words = []
        for i in range(0, temp_length+1):
            forward = this_seq[i:i+word_size]
            reverse = this_c_seq[temp_length-i:seq_len-i]
            if forward in these_words:
                these_group_id.add(these_words[forward])
            else:
                this_words.append(forward)
                this_words.append(reverse)
        len_groups = len(these_group_id)
        # create a new group
        if len_groups == 0:
            new_group_id = count_groups
            groups_of_duplicate_lines[new_group_id] = [{this_unique_read_id}, set(this_words)]
            for this_word in this_words:
                these_words[this_word] = new_group_id
            lines_with_duplicates[this_unique_read_id] = new_group_id
            count_groups += 1
        # belongs to one group
        elif len_groups == 1:
            this_group_id = these_group_id.pop()
            groups_of_duplicate_lines[this_group_id][0].add(this_unique_read_id)
            for this_word in this_words:
                groups_of_duplicate_lines[this_group_id][1].add(this_word)
                these_words[this_word] = this_group_id
            lines_with_duplicates[this_unique_read_id] = this_group_id
        # connect different groups
        else:
            these_group_id = list(these_group_id)
            these_group_id.sort()
            this_group_to_keep = these_group_id[0]
            # for related group to merge
            for to_merge in range(1, len_groups):
                this_group_to_merge = these_group_id[to_merge]
                lines_to_merge, words_to_merge = groups_of_duplicate_lines[this_group_to_merge]
                for line_to_merge in lines_to_merge:
                    groups_of_duplicate_lines[this_group_to_keep][0].add(line_to_merge)
                    lines_with_duplicates[line_to_merge] = this_group_to_keep
                for word_to_merge in words_to_merge:
                    groups_of_duplicate_lines[this_group_to_keep][1].add(word_to_merge)
                    these_words[word_to_merge] = this_group_to_keep
                    # print 'words to merge', words_to_merge
                del groups_of_duplicate_lines[this_group_to_merge]
            # for the remain group to grow
            for this_word in this_words:
                groups_of_duplicate_lines[this_group_to_keep][1].add(this_word)
                these_words[this_word] = this_group_to_keep
            groups_of_duplicate_lines[this_group_to_keep][0].add(this_unique_read_id)
            lines_with_duplicates[this_unique_read_id] = this_group_to_keep
    for del_words in groups_of_duplicate_lines:
        groups_of_duplicate_lines[del_words] = groups_of_duplicate_lines[del_words][0]
    count_del_single = 0
    for del_words in list(groups_of_duplicate_lines):
        if len(groups_of_duplicate_lines[del_words]) == 1:
            del_line = groups_of_duplicate_lines[del_words].pop()
            del lines_with_duplicates[del_line]
            del groups_of_duplicate_lines[del_words]
            count_del_single += 1
    # print 'deleting', count_del_single, 'single line groups'
    # print lines_with_duplicates
    if this_process:
        memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + ", "
    else:
        memory_usage = ''
    del these_words
    log.info(memory_usage+str(len(groups_of_duplicate_lines))+" pseudo-contigs assembled.\n")
    return groups_of_duplicate_lines, lines_with_duplicates


class RoundLimitException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NoMoreReads(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def extending_reads(accepted_words, original_fq_dir, len_indices, fg_out_per_round, pseudo_assembled, groups_of_duplicate_lines, lines_with_duplicates, fastq_indices_in_memory, round_limit, output_base, jump_step, mesh_size, verbose, resume, options, log):
    global word_size
    accepted_contig_id = set()
    accepted_contig_id_this_round = set()
    line_to_accept = set()
    round_count = 1
    previous_aw_count = 0
    if fg_out_per_round:
        round_dir = os.path.join(output_base, "Reads_per_round")
        if not os.path.exists(round_dir):
            os.mkdir(round_dir)
    if this_process and verbose:
        log.warning("Package psutil is not installed, so that memory usage will not be logged\n"
                    "Don't worry. This will not affect the result.")
    try:
        def summarise_round(acc_words, acc_contig_id_this_round, pre_aw, r_count, unique_id):
            len_aw = len(acc_words)
            len_al = len(acc_contig_id_this_round)
            if this_process:
                memory_usage = " Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
            else:
                memory_usage = ''
            if fg_out_per_round:
                write_fq_results(original_fq_dir, acc_contig_id_this_round, False, os.path.join(round_dir, "Round."+str(r_count)), os.path.join(output_base, 'temp.indices.2'), os.path.join(output_base, 'temp.indices.2'), fastq_indices_in_memory, verbose, options, log)
                # clear former accepted words from memory
                del acc_words
                # then add new accepted words into memory
                acc_words = chop_seqs(read_self_fq_seq_generator([os.path.join(round_dir, "Round."+str(r_count)+'_'+str(x+1)+'.fq') for x in range(2)]))
                acc_contig_id_this_round = set()
            log.info("Round " + str(r_count) + ': ' + str(unique_id + 1) + '/' + str(len_indices) + " AI " + str(len_al) + " AW " + str(len_aw) + memory_usage)
            #
            if len_aw == pre_aw:
                raise NoMoreReads('')
            pre_aw = len(acc_words)
            #
            if r_count == round_limit:
                raise RoundLimitException(r_count)
            r_count += 1
            return acc_words, acc_contig_id_this_round, pre_aw, r_count

        while True:
            if verbose:
                log.info("Round "+str(round_count)+": Start ...")

            if options.index_in_memory:
                reads_generator = (this_read for this_read in fastq_indices_in_memory[0])
            else:
                reads_generator = (this_read.strip() for this_read in open(os.path.join(output_base, 'temp.indices.1'), 'rU'))
            if pseudo_assembled and groups_of_duplicate_lines:
                for unique_read_id in range(len_indices):
                    this_seq = next(reads_generator)
                    this_c_seq = next(reads_generator)
                    if unique_read_id not in accepted_contig_id:
                        # if not universal_len:
                        seq_len = len(this_seq)
                        temp_length = seq_len - word_size
                        if unique_read_id in line_to_accept:
                            accepted_contig_id.add(unique_read_id)
                            accepted_contig_id_this_round.add(unique_read_id)
                            line_to_accept.remove(unique_read_id)
                            for i in range(0, temp_length+1, mesh_size):
                                # add forward
                                accepted_words.add(this_seq[i:i + word_size])
                                # add reverse
                                accepted_words.add(this_c_seq[temp_length - i:seq_len - i])
                        else:
                            accepted = False
                            for i in range(0, (temp_length+1)//2, jump_step):
                                # from first kmer
                                if this_seq[i:i + word_size] in accepted_words:
                                    accepted = True
                                    break
                                # from last kmer
                                if this_seq[temp_length-i:seq_len-i] in accepted_words:
                                    accepted = True
                                    break
                            if accepted:
                                for i in range(0, temp_length + 1, mesh_size):
                                    # add forward
                                    accepted_words.add(this_seq[i:i + word_size])
                                    # add reverse
                                    accepted_words.add(this_c_seq[temp_length - i:seq_len - i])
                                accepted_contig_id.add(unique_read_id)
                                accepted_contig_id_this_round.add(unique_read_id)
                                if unique_read_id in lines_with_duplicates:
                                    which_group = lines_with_duplicates[unique_read_id]
                                    for id_to_accept in groups_of_duplicate_lines[which_group]:
                                        line_to_accept.add(id_to_accept)
                                        del lines_with_duplicates[id_to_accept]
                                    line_to_accept.remove(unique_read_id)
                                    del groups_of_duplicate_lines[which_group]
                    if unique_read_id % 14321 == 0:
                        this_print = str("%s"%datetime.datetime.now())[:23].replace('.', ',')+" - INFO: Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + " AI " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words))
                        sys.stdout.write(this_print + '\b' * len(this_print))
                        sys.stdout.flush()
                accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count = summarise_round(accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count, unique_read_id)
            else:
                for unique_read_id in range(len_indices):
                    this_seq = next(reads_generator)
                    this_c_seq = next(reads_generator)
                    if unique_read_id not in accepted_contig_id:
                        accepted = False
                        seq_len = len(this_seq)
                        temp_length = seq_len - word_size
                        for i in range(0, (temp_length + 1) // 2, jump_step):
                            # from first kmer
                            if this_seq[i:i + word_size] in accepted_words:
                                accepted = True
                                break
                            # from last kmer
                            if this_seq[temp_length - i:seq_len - i] in accepted_words:
                                accepted = True
                                break
                        if accepted:
                            for i in range(0, temp_length + 1, mesh_size):
                                accepted_words.add(this_seq[i:i + word_size])
                                accepted_words.add(this_c_seq[temp_length - i:seq_len - i])
                            accepted_contig_id.add(unique_read_id)
                            accepted_contig_id_this_round.add(unique_read_id)
                    if unique_read_id % 14321 == 0:
                        this_print = str("%s"%datetime.datetime.now())[:23].replace('.', ',')+" - INFO: Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + " AI " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words))
                        sys.stdout.write(this_print + '\b' * len(this_print))
                        sys.stdout.flush()
                accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count = summarise_round(accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count, unique_read_id)
            reads_generator.close()
    except KeyboardInterrupt:
        reads_generator.close()
        sys.stdout.write(' ' * 100 + '\b' * 100)
        sys.stdout.flush()
        log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + " AI " + str(
                len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words)))
        log.info("KeyboardInterrupt")
    except NoMoreReads:
        reads_generator.close()
        log.info("No more reads found and terminated ...")
    except RoundLimitException as r_lim:
        reads_generator.close()
        log.info("Hit the round limit "+str(r_lim)+" and terminated ...")
    del reads_generator
    del accepted_words
    del accepted_contig_id_this_round
    del lines_with_duplicates
    return accepted_contig_id


def get_anti_lines_via_built_in_mapping(anti_words, options, log):
    anti_lines = set()
    pre_reading = [open(options.fastq_file_1, 'rU'), open(options.fastq_file_2, 'rU')]
    line_count = 0

    def add_to_anti_lines(here_head):
        try:
            if ' ' in here_head:
                here_head_split = here_head.split(' ')
                this_name, direction = here_head_split[0], int(here_head_split[1][0])
            elif '#' in here_head:
                here_head_split = here_head.split('#')
                this_name, direction = here_head_split[0], int(here_head_split[1][0])
            else:
                this_name, direction = here_head, 1
        except (ValueError, IndexError):
            log.error('Unrecognized fq format in '+str(line_count))
            exit()
        if options.pair_end_out:
            anti_lines.add(this_name)
        else:
            anti_lines.add((this_name, direction))

    for file_in in pre_reading:
        line = file_in.readline()
        if options.bowtie2_anti_seed or options.anti_seed:
            while line:
                if line.startswith("@"):
                    this_head = line[1:].strip()
                    this_seq = file_in.readline().strip()
                    # drop illegal reads
                    seq_len = len(this_seq)
                    if seq_len < word_size:
                        line_count += 4
                        for i in range(3):
                            line = file_in.readline()
                        add_to_anti_lines(this_head)
                        continue
                    this_c_seq = complementary_seq(this_seq)
                    temp_length = seq_len - word_size
                    for i in range(0, temp_length+1):
                        if this_seq[i:i + word_size] in anti_words:
                            add_to_anti_lines(this_head)
                            break
                        if this_c_seq[i:i + word_size] in anti_words:
                            add_to_anti_lines(this_head)
                            break
                else:
                    log.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                    exit()
                line_count += 1
                for i in range(3):
                    line = file_in.readline()
                    line_count += 1
    return anti_lines


def get_heads_from_sam(bowtie_sam_file, pair_end_out):
    hit_heads = set()
    if pair_end_out:
        for line in open(bowtie_sam_file, 'rU'):
            if line.strip() and not line.startswith('@'):
                hit_heads.add(line.strip().split('\t')[0])
    else:
        for line in open(bowtie_sam_file, 'rU'):
            if line.strip() and not line.startswith('@'):
                line_split = line.strip().split('\t')
                this_name = line_split[0]
                flag = int(line_split[1])
                direction = 1 if flag % 32 < 16 else 2
                hit_heads.add((this_name, direction))
    return hit_heads


def mapping_with_bowtie2(options, log):
    if options.seed_file:
        if os.path.exists(options.seed_file+'.index.1.bt2l'):
            log.info("seed bowtie2 index existed!")
        else:
            log.info("Making seed bowtie2 index ...")
            build_seed_index = subprocess.Popen("bowtie2-build --large-index "+options.seed_file+" " + options.seed_file+'.index',
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = build_seed_index.communicate()
            if "unrecognized option" in str(output):
                build_seed_index = subprocess.Popen("bowtie2-build " + options.seed_file + " " + options.seed_file + '.index',
                                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                output, err = build_seed_index.communicate()
            if options.verbose_log:
                log.info("\n"+str(output).strip())
            log.info("Making seed bowtie2 index finished.")
        seed_index_base = options.seed_file + '.index'
    else:
        seed_index_base = options.bowtie2_seed
    total_seed_file = os.path.join(options.output_base, "Initial.mapped.fq")
    total_seed_sam = [os.path.join(options.output_base, x+"seed_bowtie.sam") for x in ("temp.", "")]
    if options.script_resume and os.path.exists(total_seed_sam[1]):
        log.info("Initial seeds existed!")
    else:
        log.info("Mapping reads to seed bowtie2 index ...")
        make_seed_bowtie2 = subprocess.Popen(
            "bowtie2 -p 4 --very-fast-local --al " + total_seed_file + " -x " + seed_index_base + " -U " +
            options.fastq_file_1 + "," + options.fastq_file_2 + " -S " + total_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if options.verbose_log:
            log.info("\n"+str(output).strip())
        if os.path.exists(total_seed_sam[0]):
            os.rename(total_seed_sam[0], total_seed_sam[1])
            log.info("Mapping finished.")
        else:
            log.error("Cannot find bowtie2 result!")
            exit()
    # coverage_statistics = get_coverage(os.path.join(options.output_base, "bowtie.sam"))
    if options.anti_seed or options.bowtie2_anti_seed:
        if options.anti_seed:
            if os.path.exists(options.anti_seed + '.index.1.bt2l'):
                log.info("anti-seed bowtie2 index existed!")
            else:
                log.info("Making anti-seed bowtie2 index ...")
                build_anti_index = subprocess.Popen("bowtie2-build --large-index " + options.anti_seed + " " + options.anti_seed + '.index',
                                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                output, err = build_anti_index.communicate()
                if "unrecognized option" in str(output):
                    build_anti_index = subprocess.Popen("bowtie2-build " + options.anti_seed + " " + options.anti_seed + '.index',
                                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                    output, err = build_anti_index.communicate()
                if options.verbose_log:
                    log.info("\n" + str(output).strip())
                log.info("Making anti-seed bowtie2 index finished.")
            anti_index_base = options.anti_seed + '.index'
        else:
            anti_index_base = options.bowtie2_anti_seed

        anti_seed_sam = [os.path.join(options.output_base, x+"anti_seed_bowtie.sam") for x in ("temp.", "")]
        if options.script_resume and os.path.exists(anti_seed_sam[1]):
            log.info("Anti-seed mapping information existed!")
        else:
            log.info("Mapping reads to anti-seed bowtie2 index ...")
            make_anti_seed_bowtie2 = subprocess.Popen("bowtie2 -p 4 --very-fast-local -x " + anti_index_base + " -1 " +
                                                      options.fastq_file_1 + " -2 " + options.fastq_file_2 + " -S " +
                                                      anti_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
                                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = make_anti_seed_bowtie2.communicate()
            if options.verbose_log:
                log.info("\n" + str(output).strip())
            if os.path.exists(anti_seed_sam[0]):
                os.rename(anti_seed_sam[0], anti_seed_sam[1])
                log.info("Mapping finished.")
            else:
                log.error("Cannot find bowtie2 result!")
                exit()
        log.info("Parsing bowtie2 result ...")
        anti_lines = get_heads_from_sam(anti_seed_sam[1], options.pair_end_out) - \
                     get_heads_from_sam(total_seed_sam[1], options.pair_end_out)
        log.info("Parsing bowtie2 result finished ...")
    else:
        anti_lines = set()
    if False and not options.keep_temp_files:
        os.remove(total_seed_sam[1])
        if options.anti_seed or options.bowtie2_anti_seed:
            os.remove(anti_seed_sam[1])
    return total_seed_file, anti_lines


def assembly_with_spades(options, log):
    parameters = options.other_spades_options.strip('"')
    if '-k' in parameters:
        kmer = ''
    else:
        kmer = '-k '+options.spades_kmer
    if '-o' in parameters:
        spades_out_put = ''
    else:
        spades_out_put = '-o '+os.path.join(options.output_base, "filtered_spades")
    spades_command = ' '.join(['spades.py', parameters, '-1', os.path.join(options.output_base, "filtered_1.fq"), '-2', os.path.join(options.output_base, "filtered_2.fq"), kmer, spades_out_put]).strip()
    spades_running = subprocess.Popen(spades_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = spades_running.communicate()
    if "not recognized" in str(output):
        if options.verbose_log:
            log.warning('Problem with running SPAdes:')
            log.warning(str(output))
        log.warning("WAINING in SPAdes: unrecognized option in "+options.other_spades_options)
        log.warning('Assembling failed.')
    else:
        if options.verbose_log:
            log.info(str(output))
        log.info('Assembling finished.\n')


def slim_spades_result(options, log, depth_threshold=0):
    try:
        # python3
        blast_in_path = subprocess.getstatusoutput('blastn')
    except AttributeError:
        # python2
        blast_in_path = commands.getstatusoutput('blastn')
    if blast_in_path[0] == 32512:
        log.warning('blastn not in the path!\nSkip slimming assembly result ...')
        return
    try:
        # python3
        makeblastdb_in_path = subprocess.getstatusoutput('makeblastdb')
    except AttributeError:
        # python2
        makeblastdb_in_path = commands.getstatusoutput('makeblastdb')
    if makeblastdb_in_path[0] == 32512:
        log.warning('makeblastdb not in the path!\nSkip slimming assembly result ...')
        return
    scheme_tranlation = {'cp': ' --include-index-priority ' + os.path.join(path_of_this_script, 'Library', 'Reference', 'cp') + ' --exclude-index ' + os.path.join(path_of_this_script, 'Library', 'Reference', 'mt'),
                         'mt': ' --include-index-priority ' + os.path.join(path_of_this_script, 'Library', 'Reference', 'mt') + ' --exclude-index ' + os.path.join(path_of_this_script, 'Library', 'Reference', 'cp'),
                         'nr': ' --include-index-priority ' + os.path.join(path_of_this_script, 'Library', 'Reference', 'nr')}
    if options.scheme_for_slimming_spades_result in scheme_tranlation:
        run_command = scheme_tranlation[options.scheme_for_slimming_spades_result]
    else:
        run_command = options.scheme_for_slimming_spades_result
    if '-o' in options.other_spades_options:
        graph_file = os.path.join(options.other_spades_options.split('-o')[1].strip().split(' ')[0], "assembly_graph.fastg")
    else:
        graph_file = os.path.join(options.output_base, "filtered_spades", "assembly_graph.fastg")
    run_command = os.path.join(path_of_this_script, 'Utilities', 'slim_spades_fastg_by_blast.py')+' -g '+ graph_file + run_command+' --depth-threshold '+str(depth_threshold)
    slim_spades = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = slim_spades.communicate()
    if "not recognized" in str(output) or "command not found" in str(output):
        if options.verbose_log:
            log.warning(os.path.join(path_of_this_script, "Utilities", "slim_spades_fastg_by_blast.py")+' not found!')
            log.warning(str(output))
        log.warning('Processing assembly result failed.')
    else:
        if options.verbose_log:
            log.info(str(output))
        log.info('Processing assembly result finished!')


def require_commands(print_title, version):
    version = version
    usage = "\n"+str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq {-s reference.fasta|--bs reference.fasta.index} -w 103 -o chloroplast " \
            "[-a mitochondria.fasta -P 500000 -R 20 -k 75,85,95,105 -J 2]"
    description = print_title
    parser = OptionParser(usage=usage, version=version, description=description)
    # group1
    group_need = OptionGroup(parser, "COMMON OPTIONS", "All these arguments are required unless alternations provided")
    group_need.add_option('-1', dest='fastq_file_1', help='Input 1st fastq format file as pool.')
    group_need.add_option('-2', dest='fastq_file_2', help='Input 2nd fastq format file as pool.')
    group_need.add_option('-s', dest='seed_file', help='Reference. Input fasta format file as initial seed '
                                                       'or input bowtie index base name as pre-seed (see flag "--bs")')
    group_need.add_option('-w', dest='word_size', type=int,
                          help='Word size for extension, which have to be smaller or equal to (read length - 1)')
    group_need.add_option('-o', dest='output_base', help='Out put directory. Overwriting files if directory exists.')
    # group2
    group_result = OptionGroup(parser, "INFLUENTIAI OPTIONS", "These option will affect the final results"
                                                              " or serve as alternations of the required options")
    group_result.add_option('-R', dest='round_limit', type=int,
                            help='Limit running rounds (>=2).')
    group_result.add_option('--bs', dest='bowtie2_seed',
                            help='Input bowtie2 index base name as pre-seed. '
                                 'This flag serves as an alternation of flag "-s".')
    group_result.add_option('-a', dest='anti_seed', help='Anti-reference. Input fasta format file as anti-seed, '
                                                         'where the extension process stop. Typically serves as '
                                                         'excluding chloroplast reads when extending mitochondrial '
                                                         'reads, or the other way around.')
    group_result.add_option('--ba', dest='bowtie2_anti_seed',
                            help='Input bowtie2 index base name as pre-anti-seed. '
                                 'This flag serves as an alternation of flag "-a".')
    group_result.add_option('-J', dest='jump_step', type=int, default=1,
                            help='The wide of steps of checking words in reads during extension (integer >= 1). '
                                 'When you have reads of high quality, the larger the number is, '
                                 'the faster the extension will be, '
                                 'the more risk of missing reads in low coverage area. '
                                 'Choose 1 to choose the slowest but safest extension strategy. Default: 1')
    group_result.add_option('-M', dest='mesh_size', type=int, default=1,
                            help='(Beta parameter) '
                                 'The wide of steps of building words from seeds during extension (integer >= 1). '
                                 'When you have reads of high quality, the larger the number is, '
                                 'the faster the extension will be, '
                                 'the more risk of missing reads in low coverage area. '
                                 'Another usage of this mesh size is to choose a larger mesh size coupled with '
                                 'a smaller word size, which makes smaller word size feasible when memory is limited.'
                                 'Choose 1 to choose the slowest but safest extension strategy. Default: 1')
    group_result.add_option('-F', dest='scheme_for_slimming_spades_result', default='cp',
                            help='followed with cp, mt, nr, 0 (corresponding to certain arguments as following listed) '
                                 'or custom arguments with double quotation marks. By default, this script calls '
                                 'slim_spades_fastg_by_blast.py (should be under the same directory) '
                                 'with "cp". You can also make the index by your self.\t'
                                 ' ------------------------------------------------------ '
                                 '\ncp \t " --include-index-priority '+os.path.join(path_of_this_script, 'Library', 'Reference', 'cp')+' --exclude-index '+os.path.join(path_of_this_script, 'Library', 'Reference', 'mt')+'"'
                                 ' ------------------------------------------------------ '
                                 '\nmt \t " --include-index-priority '+os.path.join(path_of_this_script, 'Library', 'Reference', 'mt')+' --exclude-index '+os.path.join(path_of_this_script, 'Library', 'Reference', 'cp')+'"'
                                 ' ------------------------------------------------------ '
                                 '\nnr \t " --include-index-priority '+os.path.join(path_of_this_script, 'Library', 'Reference', 'nr')+'"'
                                 ' ------------------------------------------------------ '
                                 '\n0 \t disable this slimming function'
                                 ' ------------------------------------------------------ ')
    group_result.add_option('-k', dest='spades_kmer', default='65,75,85',
                            help='SPAdes kmer settings. Use the same format as in SPAdes. Default=65,75,85')
    group_result.add_option('--spades-options', dest='other_spades_options', default='--careful',
                            help='Other SPAdes options. Use double quotation marks to include all the arguments '
                                 'and parameters, such as "--careful -o test"')
    group_result.add_option('--pair-end-out', dest='pair_end_out', action="store_true", default=False,
                            help='Only accept reads with both sides of a pair been extended.')
    group_result.add_option('--no-spades', dest='run_spades', action="store_false", default=True,
                            help='Disable SPAdes.')
    group_result.add_option('--no-mapping', dest='utilize_mapping', action="store_false", default=True,
                            help='Choose to disable mapping process.'
                                 'By default, this script would map original reads to pre-seed (or bowtie2 index) '
                                 'to acquire the initial seed. This requires bowtie2 to be installed '
                                 '(and thus Linux/Unix only). It is suggested to keep mapping as default '
                                 'when the seed is too diversed from target.')
    # group3
    group_computational = OptionGroup(parser, "ADDITIONAL OPTIONS", "These options will NOT affect the final results. "
                                                                    "Take easy to pick some according your computer's flavour")
    group_computational.add_option('-P', dest='pseudo_assembled', type=int, default=200000,
                                   help='The maximum potential-organ reads to be pseudo-assembled before extension.'
                                        ' The default value is 200000.'
                                        ' For personal computer with 8G memory, we suggest no more than 300000.'
                                        ' A larger number (ex. 600000) would run faster but exhaust memory'
                                        ' in the first few minutes. Choose 0 to disable this process.')
    group_computational.add_option('--continue', dest='script_resume', default=False, action="store_true",
                                   help='limited function in this version. Only help you skip mapping process.')
    group_computational.add_option('--index-in-memory', dest='index_in_memory', action="store_true", default=False,
                                   help="Keep index in memory. Choose save index in memory than disk.")
    group_computational.add_option('--out-per-round', dest='fg_out_per_round', action="store_true", default=False,
                                   help='Enable output per round. Choose to save memory but cost more time per round.')
    group_computational.add_option('--no-remove-duplicates', dest='rm_duplicates', action="store_false", default=True,
                                   help='By default this script use unique reads to extend. Choose to disable'
                                        ' removing duplicates, which save memory but run more slow.'
                                        ' Note that whether choose or not will not disable'
                                        ' the calling of replicate reads.')
    group_computational.add_option('--keep-temp', dest='keep_temp_files', action="store_true", default=False,
                                   help="Choose to keep the running temp/index files.")
    group_computational.add_option('--verbose', dest='verbose_log', action="store_true", default=False,
                                   help='Verbose output. Choose to enable verbose running log.')
    parser.add_option_group(group_need)
    parser.add_option_group(group_result)
    parser.add_option_group(group_computational)
    try:
        (options, args) = parser.parse_args()
    except Exception as e:
        sys.stdout.write('\n############################################################################'+str(e))
        sys.stdout.write('\n"-h" for more usage')
        exit()
    else:
        if not ((options.seed_file or options.bowtie2_seed) and options.fastq_file_1 and options.fastq_file_2 and options.word_size and options.output_base):
            parser.print_help()
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Insufficient arguments!\nUsage:')
            sys.stdout.write(usage+'\n\n')
            exit()
        if options.jump_step < 1:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Jump step MUST be an integer that >= 1')
            exit()
        if options.mesh_size < 1:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Mesh size MUST be an integer that >= 1')
            exit()
        if options.fastq_file_1 == options.fastq_file_2:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: 1st fastq file is the same with 2nd fastq file!')
            exit()
        if not os.path.isdir(options.output_base):
            os.mkdir(options.output_base)
        log = simple_log(logging.getLogger(), options.output_base)
        log.info(print_title)
        log.info(' '.join(sys.argv)+'\n')
        log = complicated_log(log, options.output_base)
        if options.seed_file and options.bowtie2_seed:
            log.error('Simultaneously "-s" and "--bs" is not allowed!')
            exit()
        if options.anti_seed and options.bowtie2_anti_seed:
            log.error('Simultaneously "-a" and "--as" is not allowed!')
            exit()
        if options.bowtie2_seed or options.bowtie2_anti_seed:
            options.utilize_mapping = True
        if options.utilize_mapping:
            try:
                # python3
                bowtie2_in_path = subprocess.getstatusoutput('bowtie2')
            except:
                # python2
                bowtie2_in_path = commands.getstatusoutput('bowtie2')
            if bowtie2_in_path[0] == 32512:
                options.utilize_mapping = False
                if options.seed_file:
                    log.warning('bowtie2 not in the path! Take the seed file as initial seed.')
                else:
                    log.error('bowtie2 not in the path!')
                    exit()
                if options.anti_seed:
                    log.warning('bowtie2 not in the path! Take the anti-seed file as initial anti-seed.')
                else:
                    log.warning('bowtie2 not in the path! Anti-seed disabled!')
        if options.run_spades:
            try:
                # python3
                spades_in_path = subprocess.getstatusoutput('spades.py -h')
            except:
                # python2
                spades_in_path = commands.getstatusoutput('spades.py -h')
            if spades_in_path[0] == 32512:
                options.run_spades = False
        if not options.rm_duplicates and options.pseudo_assembled:
            log.warning("remove duplicates was inactive, so that the pseudo-assembly was disabled.")
            options.pseudo_assembled = False
        if options.round_limit and options.round_limit < 2:
            log.warning("illegal limit for rounds! Been set to default: unlimited.")
            options.round_limit = None
        return options, log


def simple_log(log, output_base):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'get_org.log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_simple.addHandler(console)
    log_simple.addHandler(logfile)
    return log_simple


def complicated_log(log, output_base):
    log_timed = log
    for handler in list(log_timed.handlers):
        log_timed.removeHandler(handler)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'get_org.log.txt'), mode='a')
    logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_timed.addHandler(console)
    log_timed.addHandler(logfile)
    return log_timed


__version__ = "1.9.79"


def main():
    time0 = time.time()
    title = "\nGetOrganelle released by Jianjun Jin on Aug 8 2016." \
            "\n" \
            "\nThis pipeline get organelle reads and genomes from genome skimming data by extending." \
            "\nFind updates in https://github.com/Kinggerm/GetOrganelle and see README.md for more information." \
            "\n"
    options, log = require_commands(print_title=title, version=__version__)
    try:
        global word_size
        word_size = options.word_size
        original_fq_dir = (options.fastq_file_1, options.fastq_file_2)

        """reading seeds"""
        log.info("Reading seeds...")
        if not options.utilize_mapping:
            initial_accepted_words = chop_seqs(read_fasta(options.seed_file)[1])
            anti_lines = get_anti_lines_via_built_in_mapping(chop_seqs(read_fasta(options.anti_seed)[1]), options, log) if options.anti_seed else set()
        else:
            seed_fastq, anti_lines = mapping_with_bowtie2(options, log)
            initial_accepted_words = chop_seqs(read_self_fq_seq_generator(seed_fastq))
        log.info("Reading seeds finished.\n")

        """reading fastq files"""
        log.info("Pre-reading fastq ...")
        fastq_indices_in_memory = read_fq_and_pair_infos(original_fq_dir, options.pair_end_out,
                                                         options.rm_duplicates, options.output_base,
                                                         anti_lines, options, log)
        len_indices = len(fastq_indices_in_memory[1])
        log.info("Pre-reading fastq finished.\n")

        """pseudo-assembly if asked"""
        # pseudo_assembly is suggested when the whole genome coverage is shallow but the organ genome coverage is deep
        if options.pseudo_assembled:
            groups_of_duplicate_lines, lines_with_duplicates = pseudo_assembly(fastq_indices_in_memory,
                                                                               options.pseudo_assembled, log)
        else:
            groups_of_duplicate_lines = None
            lines_with_duplicates = None
        if not options.index_in_memory:
            fastq_indices_in_memory = None

        """extending process"""
        log.info("Extending ...")
        accepted_contig_id = extending_reads(initial_accepted_words, original_fq_dir, len_indices, options.fg_out_per_round,
                                             options.pseudo_assembled, groups_of_duplicate_lines, lines_with_duplicates,
                                             fastq_indices_in_memory, options.round_limit, options.output_base, options.jump_step,
                                             options.mesh_size, options.verbose_log, options.script_resume, options, log)
        write_fq_results(original_fq_dir, accepted_contig_id, options.pair_end_out,
                         os.path.join(options.output_base, "filtered"), os.path.join(options.output_base, 'temp.indices.2'),
                         os.path.join(options.output_base, 'temp.indices.3'), fastq_indices_in_memory, options.verbose_log, options, log)
        del accepted_contig_id, fastq_indices_in_memory, groups_of_duplicate_lines, \
            anti_lines, initial_accepted_words, lines_with_duplicates
        if not options.keep_temp_files:
            try:
                os.remove(os.path.join(options.output_base, 'temp.indices.1'))
                os.remove(os.path.join(options.output_base, 'temp.indices.2'))
                os.remove(os.path.join(options.output_base, 'temp.indices.3'))
            except OSError:
                pass
        log.info("Extending finished.\n")

        """assembly"""
        if options.run_spades:
            log.info('Assembling with SPAdes ...')
            assembly_with_spades(options, log)
            if options.scheme_for_slimming_spades_result != '0':
                slim_spades_result(options, log)

        log = simple_log(log, options.output_base)
        log.info("\nTotal Calc-cost "+str(time.time() - time0))
        log.info("Thanks you!")
    except:
        log.exception("")
        log = simple_log(log, options.output_base)
        log.info("\nTotal Calc-cost " + str(time.time() - time0))
        log.info("Please email me if you find bugs!")
    logging.shutdown()


if __name__ == '__main__':
    main()

