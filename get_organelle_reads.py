#!/usr/bin/env python
import time
import datetime
import sys
import os
import logging
from optparse import OptionParser, OptionGroup


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
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
word_size = int


# test whether an external binary is executable
def executable(test_this):
    return True if os.access(test_this, os.X_OK) or getstatusoutput(test_this)[0] != dead_code else False


if python_version == "2.7+":
    # python2
    import string
    translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")

    def complementary_seq(input_seq):
        return string.translate(input_seq, translator)[::-1]
else:
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


def read_self_fq_seq_generator(fq_dir_list, this_trim_values):
    if not ((type(fq_dir_list) is list) or (type(fq_dir_list) is tuple)):
        fq_dir_list = [fq_dir_list]
    for fq_dir in fq_dir_list:
        count = 0
        if this_trim_values:
            trim1, trim2 = [int(trim_value) for trim_value in this_trim_values.split(',')]
            if trim2:
                for fq_line in open(fq_dir, 'rU'):
                    if count % 4 == 1:
                        yield fq_line[trim1:(len(fq_line)-trim2)]
                    count += 1
            elif trim1:
                for fq_line in open(fq_dir, 'rU'):
                    if count % 4 == 1:
                        yield fq_line[trim1:]
                    count += 1
            else:
                for fq_line in open(fq_dir, 'rU'):
                    if count % 4 == 1:
                        yield fq_line
                    count += 1
        else:
            for fq_line in open(fq_dir, 'rU'):
                if count % 4 == 1:
                    yield fq_line
                count += 1


def write_fq_results(original_fq_files, accepted_contig_id, out_file_name, temp2_clusters_dir, fastq_indices_in_memory, verbose, index_in_memory, log):

    if verbose:
        sys.stdout.write(' ' * 100 + '\b' * 100)
        sys.stdout.flush()
        log.info("Producing output ...")
        log.info("reading indices ...")
    accepted_lines = []
    if index_in_memory:
        # read cluster indices
        for this_index in accepted_contig_id:
            accepted_lines += fastq_indices_in_memory[1][this_index]
        # produce the pair-end output
        accepted_lines = set(accepted_lines)
    else:
        # read cluster indices
        temp2_indices_file_in = open(temp2_clusters_dir, 'rU')
        this_index = 0
        for line in temp2_indices_file_in:
            if this_index in accepted_contig_id:
                accepted_lines += [int(x) for x in line.strip().split('\t')]
            this_index += 1
        accepted_lines = set(accepted_lines)

    # write by line
    if verbose:
        log.info("writing fastq lines ...")
    post_reading = [open(fq_file, 'rU') for fq_file in original_fq_files]
    files_out = [open(out_file_name+'_'+str(i+1)+'.temp', 'w') for i in range(len(original_fq_files))]
    line_count = 0
    for i in range(len(original_fq_files)):
        line = post_reading[i].readline()
        while line:
            if line_count in accepted_lines:
                files_out[i].write(line)
                for j in range(3):
                    files_out[i].write(post_reading[i].readline())
                    line_count += 1
                line = post_reading[i].readline()
                line_count += 1
            else:
                for j in range(4):
                    line = post_reading[i].readline()
                    line_count += 1
        files_out[i].close()
        post_reading[i].close()
    del accepted_lines
    for i in range(len(original_fq_files)):
        os.rename(out_file_name + '_' + str(i + 1) + '.temp', out_file_name + '_' + str(i + 1) + '.fq')
    if verbose:
        log.info("writing fastq lines finished.")


def read_fq_infos(original_fq_files, rm_duplicates, output_base, anti_lines, pseudo_assembled, index_in_memory, bowtie2_anti_seed, anti_seed, trim_values, resume, log):
    if trim_values:
        trim1, trim2 = [int(trim_value) for trim_value in trim_values.split(',')]
    # read original reads
    # line_cluster (list) ~ forward_reverse_reads
    line_clusters = []
    seq_duplicates = {}
    forward_reverse_reads = []
    line_count = 0
    this_index = 0
    #
    name_to_line = {}
    #
    temp1_contig_dir = [os.path.join(output_base, k+'temp.indices.1') for k in ("_", "")]
    temp2_clusters_dir = [os.path.join(output_base, k+'temp.indices.2') for k in ("_", "")]
    if resume and os.path.exists(temp1_contig_dir[1]) and os.path.exists(temp2_clusters_dir[1]):
        if pseudo_assembled or index_in_memory:
            log.info("Reading existed indices for fastq ...")
            #
            forward_reverse_reads = [x.strip() for x in open(temp1_contig_dir[1], 'rU')]
            #
            line_clusters = [[int(x) for x in y.split('\t')] for y in open(temp2_clusters_dir[1], 'rU')]
            if rm_duplicates:
                line_count = sum([len(x) for x in line_clusters])*4
            # log
            len_indices = len(line_clusters)
            if this_process:
                memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + ", "
            else:
                memory_usage = ''
            if rm_duplicates:
                log.info(memory_usage + str(len_indices) + " candidates in all " + str(line_count // 4) + " reads")
            else:
                log.info(memory_usage + str(len_indices) + " reads")
        else:
            log.info("indices for fastq existed!")
            len_indices = len([x for x in open(temp2_clusters_dir[1], 'rU')])
    else:
        if not index_in_memory:
            temp1_contig_out = open(temp1_contig_dir[0], 'w')
        pre_reading = [open(fq_file, 'rU') for fq_file in original_fq_files]
        for file_in in pre_reading:
            line = file_in.readline()
            if bowtie2_anti_seed or anti_seed:
                while line:
                    if line.startswith("@"):
                        try:
                            if ' ' in line:
                                this_head = line[1:].split(' ')
                                this_name, direction = this_head[0], int(this_head[1][0])
                            elif '#' in line:
                                this_head = line[1:].split('#')
                                this_name, direction = this_head[0], int(this_head[1].strip("/")[0])
                            else:
                                this_name, direction = line[1:].strip(), 1
                        except (ValueError, IndexError):
                            log.error('Unrecognized fq format in '+str(line_count)+' '+str(line))
                            exit()
                        else:
                            if (this_name, direction) in anti_lines:
                                line_count += 4
                                for i in range(4):
                                    line = file_in.readline()
                                continue
                        this_seq = file_in.readline().strip()
                        if trim_values:
                            this_seq = this_seq[trim1:(len(this_seq)-trim2)].strip("N")
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
                                if not index_in_memory:
                                    temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                                seq_duplicates[this_seq] = this_index
                                line_clusters.append([line_count])
                                this_index += 1
                        else:
                            line_clusters.append([line_count])
                            if index_in_memory:
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
                                if not index_in_memory:
                                    temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                                seq_duplicates[this_seq] = this_index
                                line_clusters.append([line_count])
                                this_index += 1
                        else:
                            line_clusters.append([line_count])
                            if index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            else:
                                temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                    else:
                        log.error("Illegal fq format in line "+str(line_count)+' '+str(line))
                        exit()
                    if line_count % 54321 == 0:
                        to_print = str("%s"%datetime.datetime.now())[:23].replace('.', ',')+" - INFO: "+str((line_count + 4) // 4)+" reads"
                        sys.stdout.write(to_print+'\b'*len(to_print))
                        sys.stdout.flush()
                    line_count += 1
                    for i in range(3):
                        line = file_in.readline()
                        line_count += 1
            file_in.close()
        if not index_in_memory:
            temp1_contig_out.close()
            os.rename(temp1_contig_dir[0], temp1_contig_dir[1])

        if this_process:
            memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + ", "
        else:
            memory_usage = ''

        del name_to_line

        if not index_in_memory:
            # dump line clusters
            len_indices = len(line_clusters)
            temp2_indices_file_out = open(temp2_clusters_dir[0], 'w')
            for this_index in range(len_indices):
                temp2_indices_file_out.write('\t'.join([str(x) for x in line_clusters[this_index]]))
                temp2_indices_file_out.write('\n')
            temp2_indices_file_out.close()
            os.rename(temp2_clusters_dir[0], temp2_clusters_dir[1])

        del seq_duplicates
        del pre_reading
        len_indices = len(line_clusters)
        if rm_duplicates:
            log.info(memory_usage + str(len_indices) + " candidates in all " + str(line_count // 4) + " reads")
        else:
            log.info(memory_usage + str(len_indices) + " reads")
    return forward_reverse_reads, line_clusters, len_indices


def pseudo_assembly(fastq_indices_in_memory, dupli_threshold, log):
    global word_size
    forward_and_reverse_reads, line_clusters, len_indices = fastq_indices_in_memory
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


def extending_reads(accepted_words, original_fq_dir, len_indices, pseudo_assembled, groups_of_duplicate_lines, lines_with_duplicates, fastq_indices_in_memory, output_base, round_limit, fg_out_per_round, jump_step, mesh_size, verbose, resume, trim_values, log):
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
                write_fq_results(original_fq_dir, acc_contig_id_this_round, os.path.join(round_dir, "Round."+str(r_count)), os.path.join(output_base, 'temp.indices.2'), fastq_indices_in_memory, verbose, bool(fastq_indices_in_memory), log)
                # clear former accepted words from memory
                del acc_words
                # then add new accepted words into memory
                acc_words = chop_seqs(read_self_fq_seq_generator([os.path.join(round_dir, "Round."+str(r_count)+'_'+str(x+1)+'.fq') for x in range(len(original_fq_dir))], trim_values))
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

            if fastq_indices_in_memory:
                reads_generator = (this_read for this_read in fastq_indices_in_memory[0])
            else:
                reads_generator = (this_read.strip() for this_read in open(os.path.join(output_base, 'temp.indices.1'), 'rU'))
            unique_read_id = 0
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


def get_anti_built_in_mapping(anti_words, anti_input, original_fq_files, log):
    anti_lines = set()
    pre_reading = [open(fq_file, 'rU') for fq_file in original_fq_files]
    line_count = 0

    def add_to_anti_lines(here_head):
        try:
            if ' ' in here_head:
                here_head_split = here_head.split(' ')
                this_name, direction = here_head_split[0], int(here_head_split[1][0])
            elif '#' in here_head:
                here_head_split = here_head.split('#')
                this_name, direction = here_head_split[0], int(here_head_split[1].strip("/")[0])
            else:
                this_name, direction = here_head, 1
        except (ValueError, IndexError):
            log.error('Unrecognized fq format in '+str(line_count))
            exit()
        else:
            anti_lines.add((this_name, direction))

    for file_in in pre_reading:
        line = file_in.readline()
        if anti_input:
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
        file_in.close()
    return anti_lines


def get_heads_from_sam(bowtie_sam_file):
    hit_heads = set()
    for line in open(bowtie_sam_file, 'rU'):
        if line.strip() and not line.startswith('@'):
            line_split = line.strip().split('\t')
            this_name = line_split[0]
            flag = int(line_split[1])
            direction = 1 if flag % 32 < 16 else 2
            hit_heads.add((this_name, direction))
    return hit_heads


def mapping_with_bowtie2(seed_file, bowtie2_seed, anti_seed, bowtie2_anti_seed, original_fq_files, out_base, resume, verbose_log, threads, log):
    if seed_file:
        if os.path.exists(seed_file+'.index.1.bt2l'):
            log.info("seed bowtie2 index existed!")
        else:
            log.info("Making seed bowtie2 index ...")
            build_seed_index = subprocess.Popen("bowtie2-build --large-index "+seed_file+" " + seed_file+'.index',
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = build_seed_index.communicate()
            if "unrecognized option" in str(output):
                build_seed_index = subprocess.Popen("bowtie2-build " + seed_file + " " + seed_file + '.index',
                                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                output, err = build_seed_index.communicate()
            if "(ERR)" in str(output) or "Error:" in str(output):
                log.error('\n'+str(output))
                exit()
            if verbose_log:
                log.info("\n"+str(output).strip())
            log.info("Making seed bowtie2 index finished.")
        seed_index_base = seed_file + '.index'
    else:
        seed_index_base = bowtie2_seed
    total_seed_file = [os.path.join(out_base, x + "Initial.mapped.fq") for x in ("temp.", "")]
    total_seed_sam = [os.path.join(out_base, x + "seed_bowtie.sam") for x in ("temp.", "")]
    if resume and os.path.exists(total_seed_file[1]):
        log.info("Initial seeds existed!")
    else:
        log.info("Mapping reads to seed bowtie2 index ...")
        make_seed_bowtie2 = subprocess.Popen(
            "bowtie2 -p "+str(threads)+" --very-fast-local --al " + total_seed_file[0] + " -x " + seed_index_base + " -U " +
            ",".join(original_fq_files) + " -S " + total_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in str(output) or "Error:" in str(output):
            log.error('\n' + str(output))
            exit()
        if verbose_log:
            log.info("\n"+str(output).strip())
        if os.path.exists(total_seed_sam[0]):
            os.rename(total_seed_sam[0], total_seed_sam[1])
            os.rename(total_seed_file[0], total_seed_file[1])
            log.info("Mapping finished.")
        else:
            log.error("Cannot find bowtie2 result!")
            exit()
    if anti_seed or bowtie2_anti_seed:
        if anti_seed:
            if os.path.exists(anti_seed + '.index.1.bt2l'):
                log.info("anti-seed bowtie2 index existed!")
            else:
                log.info("Making anti-seed bowtie2 index ...")
                build_anti_index = subprocess.Popen("bowtie2-build --large-index " + anti_seed + " " + anti_seed + '.index',
                                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                output, err = build_anti_index.communicate()
                if "unrecognized option" in str(output):
                    build_anti_index = subprocess.Popen("bowtie2-build " + anti_seed + " " + anti_seed + '.index',
                                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                    output, err = build_anti_index.communicate()
                if "(ERR)" in str(output) or "Error:" in str(output):
                    log.error('\n' + str(output))
                    exit()
                if verbose_log:
                    log.info("\n" + str(output).strip())
                log.info("Making anti-seed bowtie2 index finished.")
            anti_index_base = anti_seed + '.index'
        else:
            anti_index_base = bowtie2_anti_seed

        anti_seed_sam = [os.path.join(out_base, x+"anti_seed_bowtie.sam") for x in ("temp.", "")]
        if resume and os.path.exists(anti_seed_sam[1]):
            log.info("Anti-seed mapping information existed!")
        else:
            log.info("Mapping reads to anti-seed bowtie2 index ...")
            make_anti_seed_bowtie2 = subprocess.Popen("bowtie2 -p "+str(threads)+" --very-fast-local -x " + anti_index_base + " -U " +
                                                      ",".join(original_fq_files) + " -S " +
                                                      anti_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
                                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = make_anti_seed_bowtie2.communicate()
            if "(ERR)" in str(output) or "Error:" in str(output):
                log.error('\n'+str(output))
                exit()
            if verbose_log:
                log.info("\n" + str(output).strip())
            if os.path.exists(anti_seed_sam[0]):
                os.rename(anti_seed_sam[0], anti_seed_sam[1])
                log.info("Mapping finished.")
            else:
                log.error("Cannot find bowtie2 result!")
                exit()
        log.info("Parsing bowtie2 result ...")
        anti_lines = get_heads_from_sam(anti_seed_sam[1]) - get_heads_from_sam(total_seed_sam[1])
        log.info("Parsing bowtie2 result finished ...")
    else:
        anti_lines = set()
    return total_seed_file[1], anti_lines


def assembly_with_spades(spades_kmer, spades_out_put, parameters, out_base, original_fq_files, reads_paired, verbose_log, resume, threads, log):
    if '-k' in parameters:
        kmer = ''
    else:
        kmer = '-k '+spades_kmer
    if resume:
        continue_command = '--continue'
    else:
        continue_command = ''
    spades_out_put = '-o '+spades_out_put
    if reads_paired['input'] and reads_paired['pair_out']:
        spades_command = ' '.join(
            ['spades.py', '-t', str(threads), continue_command, parameters, '-1', os.path.join(out_base, "filtered_1_paired.fq"), '-2',
             os.path.join(out_base, "filtered_2_paired.fq"), '--s1', os.path.join(out_base, "filtered_1_unpaired.fq")] +
            ['--s' + str(i + 3) + ' ' + str(os.path.join(out_base, "filtered_" + str(i + 3) + ".fq")) for i in
             range(len(original_fq_files)-2)] +
            ['--s2', os.path.join(out_base, "filtered_2_unpaired.fq"), kmer, spades_out_put]).strip()
    else:
        spades_command = ' '.join(
            ['spades.py', '-t', str(threads), continue_command, parameters] +
            ['--s' + str(i + 1) + ' ' + str(os.path.join(out_base, "filtered_" + str(i + 1) + ".fq")) for i in
             range(len(original_fq_files))] +
            [kmer, spades_out_put]).strip()
    spades_running = subprocess.Popen(spades_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = spades_running.communicate()
    if "not recognized" in str(output):
        if verbose_log:
            log.error('Problem with running SPAdes:')
            log.error(str(output))
        log.error("WAINING in SPAdes: unrecognized option in "+parameters)
        log.error('Assembling failed.')
        return False
    elif "== Error ==" in str(output):
        log.error("Error in SPAdes: \n== Error =="+str(output).split("== Error ==")[-1].split("In case you")[0])
        log.error('Assembling failed.')
        return False
    else:
        if verbose_log:
            log.info(str(output))
        log.info('Assembling finished.\n')
        return True


def slim_spades_result(scheme, spades_output, verbose_log, log, threads, depth_threshold=0):
    if not executable("blastn"):
        log.warning('blastn not in the path!\nSkip slimming assembly result ...')
        return
    if not executable("makeblastdb"):
        log.warning('makeblastdb not in the path!\nSkip slimming assembly result ...')
        return
    scheme_tranlation = {'cp': ' --include-priority ' + os.path.join(path_of_this_script, 'Library', 'NotationReference', 'cp') + ' --exclude ' + os.path.join(path_of_this_script, 'Library', 'NotationReference', 'mt'),
                         'mt': ' --include-priority ' + os.path.join(path_of_this_script, 'Library', 'NotationReference', 'mt') + ' --exclude ' + os.path.join(path_of_this_script, 'Library', 'NotationReference', 'cp'),
                         'nr': ' --include-priority ' + os.path.join(path_of_this_script, 'Library', 'NotationReference', 'nr')}
    if scheme in scheme_tranlation:
        run_command = scheme_tranlation[scheme]
    else:
        run_command = scheme
    graph_file = os.path.join(spades_output, "assembly_graph.fastg")
    run_command = os.path.join(path_of_this_script, 'Utilities', 'slim_fastg_by_blast.py') + ' -t ' + str(threads) +' ' + graph_file + run_command+' --depth-threshold '+str(depth_threshold)
    slim_spades = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = slim_spades.communicate()
    if "not recognized" in str(output) or "command not found" in str(output):
        if verbose_log:
            log.warning(os.path.join(path_of_this_script, "Utilities", "slim_spades_fastg_by_blast.py")+' not found!')
            log.warning(str(output))
        log.warning("Processing assembly result failed.")
    else:
        if verbose_log:
            log.info(str(output))
        log.info("Processing assembly result finished!")


def separate_fq_by_pair(out_base, verbose_log, log):
    log.info("Separating filtered fastq file ... ")
    run_command = os.path.join(path_of_this_script, "Utilities", "get_pair_reads.py")\
                  + ' ' + os.path.join(out_base, "filtered_1.fq")\
                  + ' ' + os.path.join(out_base, "filtered_2.fq")
    separating = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = separating.communicate()
    if "not recognized" in str(output) or "command not found" in str(output):
        if verbose_log:
            log.warning(os.path.join(path_of_this_script, "Utilities", "get_pair_reads.py") + "not found!")
            log.warning(str(output))
        log.warning("Separating filtered fastq file failed.")
        return False
    elif not os.path.exists(os.path.join(out_base, "filtered_2_paired.fq")):
        log.warning("Separating filtered fastq file failed with unknown error: " + run_command)
        return False
    else:
        if verbose_log:
            log.info(str(output))
        log.info("Separating filtered fastq file finished!")
        return True


def require_commands(print_title, version):
    version = version
    usage = "\n###  Chloroplast, Normal, 2G raw data, 150 bp reads\n"+str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output " \
            " -R 10 -k 75,85,95,105 -P 300000\n" \
            "###  Chloroplast, Fast, Memory-consuming\n"+str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output " \
            " -R 5 -k 75,85,95,105 -P 1000000 -a mitochondria.fasta -J 3 -M 5\n" \
            "###  Chloroplast, Slow, Memory-economic\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output " \
            " -R 10 -k 75,85,95,105 -P 0 --out-per-round --no-remove-duplicates\n" \
            "###  Mitochondria\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s mt_reference.fasta -w 93 -o mitochondria_output " \
            " -R 30 -k 65,75,85,95 -P 1000000 -F mt\n" \
            "###  Nuclear Ribosomal RNA (18S-ITS1-5.8S-ITS2-26S)\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s nr_reference.fasta -w 115 -o nr_output " \
            " -R 7 -k 95,105,115 -P 0 -F nr\n"
    description = print_title
    parser = OptionParser(usage=usage, version=version, description=description)
    # group1
    group_need = OptionGroup(parser, "COMMON OPTIONS", "All these arguments are required unless alternations provided")
    group_need.add_option('-1', dest='fastq_file_1', help='Input file with forward paired-end reads as pool.')
    group_need.add_option('-2', dest='fastq_file_2', help='Input file with reverse paired-end reads as pool.')
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
    group_result.add_option('-u', dest='unpaired_fastq_files',
                            help='Input file(s) with unpaired reads as pool. '
                                 'files could be comma-separated lists such as "seq1,seq2".')
    group_result.add_option('--bs', dest='bowtie2_seed',
                            help='Input bowtie2 index base name as pre-seed. '
                                 'This flag serves as an alternation of flag "-s".')
    group_result.add_option('-a', dest='anti_seed', help='Anti-reference. Input fasta format file as anti-seed, '
                                                         'where the extension process stop. Typically serves as '
                                                         'excluding chloroplast reads when extending mitochondrial '
                                                         'reads, or the other way around. You should be cautious about '
                                                         'using this option, because if the anti-seed includes '
                                                         'some word in the target but not in the seed, the result '
                                                         'would have gaps. Typically, use the mt and cp from the same '
                                                         'species as seed and anti-seed.')
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
                            help='This is designed to capture target associated contigs from total graph, by calling '
                                 '"Utilities/slim_spades_fastg_by_blast.py". This flag should be followed with '
                                 'cp (if you want get chloroplast), mt (mitochondria), nr (nuclear ribosomal RNA), '
                                 '0 (disable this) or custom arguments with double quotation marks. Default: cp. '
                                 'You can also make the index by your self and add those index to ' +
                                 os.path.join(path_of_this_script, 'Library', '/NotationReference') + '')
    group_result.add_option('--trim', dest='trim_values',
                            help='Assign the number of bases in the ends to trim in extending process. '
                                 'This function will not change the length of the out put reads. '
                                 'Input format: int,int (Example: 4,4). Default: 0,0')
    group_result.add_option('-k', dest='spades_kmer', default='65,75,85',
                            help='SPAdes kmer settings. Use the same format as in SPAdes. Default=65,75,85')
    group_result.add_option('--spades-options', dest='other_spades_options', default='--careful',
                            help='Other SPAdes options. Use double quotation marks to include all the arguments '
                                 'and parameters, such as "--careful -o test"')
    group_result.add_option('--no-spades', dest='run_spades', action="store_false", default=True,
                            help='Disable SPAdes.')
    group_result.add_option('--no-bowtie2', dest='utilize_mapping', action="store_false", default=True,
                            help='Choose to disable mapping process.'
                                 'By default, this script would map original reads to pre-seed (or bowtie2 index) '
                                 'to acquire the initial seed. This requires bowtie2 to be installed '
                                 '(and thus Linux/Unix only). It is suggested to keep mapping as default '
                                 'when the seed is too diversed from target.')
    # group3
    group_computational = OptionGroup(parser, "ADDITIONAL OPTIONS", "These options will NOT affect the final results. "
                                                                    "Take easy to pick some according your computer's flavour")
    group_computational.add_option('-t', dest='threads', type=int, default=4,
                                   help="Threads for third-party tools (bowtie2, blastn, SPAdes).")
    group_computational.add_option('-P', dest='pseudo_assembled', type=int, default=200000,
                                   help='The maximum potential-organ reads to be pseudo-assembled before extension. '
                                        'pseudo_assembly is suggested when the whole genome coverage is shallow but '
                                        'the organ genome coverage is deep.'
                                        'The default value is 200000. '
                                        'For personal computer with 8G memory, we suggest no more than 300000. '
                                        'A larger number (ex. 600000) would run faster but exhaust memory '
                                        'in the first few minutes. Choose 0 to disable this process.')
    group_computational.add_option('--continue', dest='script_resume', default=False, action="store_true",
                                   help='Several check point based on files produced, rather than log, '
                                        'so keep in mind that this script will not detect the difference '
                                        'between this input parameters and the previous ones.')
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
        if not ((options.seed_file or options.bowtie2_seed) and ((options.fastq_file_1 and options.fastq_file_2) or options.unpaired_fastq_files) and options.word_size and options.output_base):
            parser.print_help()
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Insufficient arguments!\nUsage:')
            sys.stdout.write(usage+'\n\n')
            exit()
        if int(bool(options.fastq_file_1)) + int(bool(options.fastq_file_2)) == 1:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: unbalanced paired reads!\n\n')
            exit()
        for check_file in (options.fastq_file_1, options.fastq_file_2, options.seed_file, options.anti_seed):
            if check_file:
                if not os.path.exists(check_file):
                    sys.stdout.write('\n############################################################################'
                                     '\nERROR: ' + check_file + ' not found!\n\n')
                    exit()
        if options.unpaired_fastq_files:
            options.unpaired_fastq_files = options.unpaired_fastq_files.split(',')
            for fastq_file in options.unpaired_fastq_files:
                if not os.path.exists(fastq_file):
                    sys.stdout.write('\n############################################################################'
                                     '\nERROR: ' + fastq_file + ' not found!\n\n')
                    exit()
        else:
            options.unpaired_fastq_files = []
        if options.jump_step < 1:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Jump step MUST be an integer that >= 1')
            exit()
        if options.mesh_size < 1:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Mesh size MUST be an integer that >= 1')
            exit()
        if options.fastq_file_1 == options.fastq_file_2 and options.fastq_file_1:
            sys.stdout.write('\n############################################################################'
                             '\nERROR: 1st fastq file is the same with 2nd fastq file!')
            exit()
        if not os.path.isdir(options.output_base):
            os.mkdir(options.output_base)
        log = simple_log(logging.getLogger(), options.output_base)
        log.info(print_title)
        log.info(' '.join(sys.argv)+'\n')
        log = timed_log(log, options.output_base)
        if options.seed_file and options.bowtie2_seed:
            log.error('Simultaneously "-s" and "--bs" is not allowed!')
            exit()
        if options.anti_seed and options.bowtie2_anti_seed:
            log.error('Simultaneously "-a" and "--as" is not allowed!')
            exit()
        if options.bowtie2_seed or options.bowtie2_anti_seed:
            options.utilize_mapping = True
        if options.utilize_mapping:
            if not executable("bowtie2"):
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
            if not executable("spades.py -h"):
                log.warning("spades.py not found in the path. Only get the reads and skip assembly.")
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


def timed_log(log, output_base):
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


__version__ = "1.9.81"


def main():
    time0 = time.time()
    title = "\nGetOrganelle released by Jianjun Jin on Oct 19 2016." \
            "\n" \
            "\nThis pipeline get organelle reads and genomes from genome skimming data by extending." \
            "\nFind updates in https://github.com/Kinggerm/GetOrganelle and see README.md for more information." \
            "\n"
    options, log = require_commands(print_title=title, version=__version__)
    try:
        """ initialization """
        global word_size
        word_size = options.word_size
        if options.fastq_file_1 and options.fastq_file_2:
            reads_paired = {'input': True, 'pair_out': bool}
            original_fq_files = [options.fastq_file_1, options.fastq_file_2] + \
                                [fastq_file for fastq_file in options.unpaired_fastq_files]
        else:
            reads_paired = {'input': False, 'pair_out': False}
            original_fq_files = [fastq_file for fastq_file in options.unpaired_fastq_files]

        resume = options.script_resume
        verb_out = options.verbose_log
        out_base = options.output_base
        other_options = options.other_spades_options.split(' ')
        if '-o' in other_options:
            which_out = other_options.index('-o')
            spades_output = other_options[which_out+1]
            del other_options[which_out: which_out+2]
        else:
            spades_output = os.path.join(out_base, "filtered_spades")
        other_options = ' '.join(other_options)

        """get reads"""
        filtered_files_exist = max(
            min([os.path.exists(str(os.path.join(out_base, "filtered")) + '_' + str(i + 1) + '_unpaired.fq') for i in
                 range(2)] +
                [os.path.exists(str(os.path.join(out_base, "filtered")) + '_' + str(i + 1) + '.fq') for i in
                 range(2, len(original_fq_files))]),
            min([os.path.exists(str(os.path.join(out_base, "filtered")) + '_' + str(i + 1) + '.fq') for i in
                range(len(original_fq_files))]))
        if not (resume and filtered_files_exist):

            seed_file = options.seed_file
            bowt_seed = options.bowtie2_seed
            anti_seed = options.anti_seed
            b_at_seed = options.bowtie2_anti_seed
            pseudoasm = options.pseudo_assembled
            trim_ends = options.trim_values
            in_memory = options.index_in_memory

            """reading seeds"""
            log.info("Reading seeds...")
            if not options.utilize_mapping:
                initial_accepted_words = chop_seqs(read_fasta(seed_file)[1])
                anti_lines = get_anti_built_in_mapping(chop_seqs(read_fasta(anti_seed)[1]),
                                                       (anti_seed or b_at_seed),
                                                       original_fq_files, log) if anti_seed else set()
            else:
                seed_fastq, anti_lines = mapping_with_bowtie2(seed_file, bowt_seed, anti_seed,
                                                              b_at_seed, original_fq_files,
                                                              out_base, resume,
                                                              verb_out, options.threads, log)
                initial_accepted_words = chop_seqs(read_self_fq_seq_generator(seed_fastq, trim_ends))
            log.info("Reading seeds finished.\n")

            """reading fastq files"""
            log.info("Pre-reading fastq ...")
            fastq_indices_in_memory = read_fq_infos(original_fq_files, options.rm_duplicates, out_base,
                                                    anti_lines, pseudoasm, in_memory,
                                                    b_at_seed, anti_seed, trim_ends,
                                                    resume, log)
            len_indices = fastq_indices_in_memory[2]
            log.info("Pre-reading fastq finished.\n")

            """pseudo-assembly if asked"""
            if pseudoasm:
                groups_of_duplicate_lines, lines_with_duplicates = pseudo_assembly(fastq_indices_in_memory,
                                                                                   pseudoasm, log)
            else:
                groups_of_duplicate_lines = None
                lines_with_duplicates = None
            if not in_memory:
                fastq_indices_in_memory = None

            """extending process"""
            log.info("Extending ...")
            accepted_contig_id = extending_reads(initial_accepted_words, original_fq_files, len_indices,
                                                 pseudoasm, groups_of_duplicate_lines, lines_with_duplicates,
                                                 fastq_indices_in_memory, out_base, options.round_limit,
                                                 options.fg_out_per_round,
                                                 options.jump_step,
                                                 options.mesh_size, verb_out, resume,
                                                 trim_ends, log)
            write_fq_results(original_fq_files, accepted_contig_id,
                             os.path.join(out_base, "filtered"), os.path.join(out_base, 'temp.indices.2'),
                             fastq_indices_in_memory, verb_out, in_memory, log)
            del accepted_contig_id, fastq_indices_in_memory, groups_of_duplicate_lines, \
                anti_lines, initial_accepted_words, lines_with_duplicates

            if not options.keep_temp_files:
                try:
                    os.remove(os.path.join(out_base, 'temp.indices.1'))
                    os.remove(os.path.join(out_base, 'temp.indices.2'))
                except OSError:
                    pass

            log.info("Extending finished.\n")
        else:
            log.info("Extending ... skipped.")

        if reads_paired['input']:
            if not (resume and min([os.path.exists(x) for x in (os.path.join(out_base, "filtered_"+y+"_"+z+"paired.fq") for y in ('1', '2') for z in ('', 'un'))])):
                resume = False
                reads_paired['pair_out'] = separate_fq_by_pair(out_base, verb_out, log)
                if reads_paired['pair_out'] and not options.keep_temp_files:
                    os.remove(os.path.join(out_base, "filtered_1.fq"))
                    os.remove(os.path.join(out_base, "filtered_2.fq"))
            else:
                log.info("Separating filtered fastq file ... skipped.")

        """assembly"""
        if options.run_spades:
            if not (resume and os.path.exists(os.path.join(spades_output, 'assembly_graph.fastg'))):
                # resume = False
                log.info('Assembling with SPAdes ...')
                assembly_with_spades(options.spades_kmer, spades_output, other_options, out_base,
                                     original_fq_files, reads_paired, options.verbose_log, resume, options.threads, log)
            else:
                log.info('Assembling with SPAdes ... skipped.')

        """slim"""
        if os.path.exists(os.path.join(spades_output, 'assembly_graph.fastg'))\
                and options.scheme_for_slimming_spades_result != '0':
            slim_spades_result(options.scheme_for_slimming_spades_result, spades_output,
                               options.verbose_log, log, options.threads)

        log = simple_log(log, out_base)
        log.info("\nTotal Calc-cost "+str(time.time() - time0))
        log.info("Thanks you!")
    except:
        log.exception("")
        log = simple_log(log, out_base)
        log.info("\nTotal Calc-cost " + str(time.time() - time0))
        log.info("Please email me if you find bugs!")
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