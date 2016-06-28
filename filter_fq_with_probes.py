#!/usr/bin/env python
import time
import sys
import os
import string
import commands, subprocess
from optparse import OptionParser, OptionGroup


word_size = int
translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


def complementary_seq(input_seq):
    return string.translate(input_seq, translator)[::-1]


def chop_seqs(seq_generator_or_list, word_size):
    return_words = set()
    for seed in seq_generator_or_list:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            for i in xrange(0, this_seq_len-word_size+1):
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
    for fq_dir in fq_dir_list:
        count = 0
        for fq_line in open(fq_dir, 'rU'):
            if count % 4 == 1:
                yield fq_line
            count += 1


def write_fq_results(original_fq_dir, accepted_contig_id, pair_end_out, print_out, out_file_name):
    global temp3_pair_dir, temp2_clusters_dir
    if print_out:
        print "\nWriting ..."
        time_write = time.time()
        print "reading indices ..."
        time_indices = time.time()
    temp2_indices_file_in = open(temp2_clusters_dir, 'rU')
    accepted_lines = []
    this_index = 0
    for line in temp2_indices_file_in:
        if this_index in accepted_contig_id:
            accepted_lines += [int(x) for x in line.strip().split('\t')]
        this_index += 1
    accepted_lines = set(accepted_lines)
    if print_out:
        print time.time()-time_indices
    # produce the pair-end output
    if pair_end_out:
        if print_out:
            print "reading pair infos ..."
            time_pair = time.time()
        pair_to_each_other = dict()
        temp1_indices_file_in = open(temp3_pair_dir, 'rU')
        line_count = 0
        for line in temp1_indices_file_in:
            if line_count in accepted_lines:
                pair_to_each_other[line_count] = int(line.strip())
            line_count += 4
        temp1_indices_file_in.close()
        # os.remove(temp1_indices_file_dir)
        dump_candidates = set()
        for candidate in accepted_lines:
            if pair_to_each_other[candidate] not in accepted_lines:
                dump_candidates.add(candidate)
        accepted_lines -= dump_candidates
        del pair_to_each_other
        if print_out:
            print time.time() - time_pair

    # write by line
    if print_out:
        print "writing lines ..."
        time_line = time.time()
    post_reading = [open(original_fq_dir[0], 'rU'), open(original_fq_dir[1], 'rU')]
    file_out = [open(out_file_name+'_1.fq', 'wb'), open(out_file_name+'_2.fq', 'wb')]
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
    if print_out:
        print time.time() - time_line
        print "Writing cost", time.time()-time_write


def read_fq_and_pair_infos(original_fq_dir, pair_end_out, rm_duplicates, output_base):

    global temp3_pair_dir, temp2_clusters_dir
    ##
    # read original reads
    # create pair list
    # pair_to_each_other
    # line_cluster (list) ~ contig_seqs ~ contig_c_seqs
    pre_reading = [open(original_fq_dir[0], 'rU'), open(original_fq_dir[1], 'rU')]
    line_clusters = []
    seq_duplicates = {}
    contig_seqs = []
    contig_c_seqs = []
    line_count = 0
    this_index = 0
    #
    pair_to_each_other = {}
    name_to_line = {}
    #
    temp1_contig_dir = os.path.join(output_base, 'temp.indices.1')
    temp1_contig_out = open(temp1_contig_dir, 'wb')
    for file_in in pre_reading:
        line = file_in.readline()
        while line:
            if line.startswith("@"):
                if pair_end_out:
                    try:
                        # take the line startswith @ as the seq name, skip the last char
                        this_name = line[1:].split(' ')[0].split('#')[0]
                        if this_name in name_to_line:
                            pair_to_each_other[line_count] = name_to_line[this_name]
                            pair_to_each_other[name_to_line[this_name]] = line_count
                        else:
                            name_to_line[this_name] = line_count
                    except KeyError:
                        print 'Unknown error in format. Please try again with "--no_pair_end_out".'
                        exit()
                #
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
                        contig_seqs.append(this_seq)
                        contig_c_seqs.append(this_c_seq)
                        temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
                        seq_duplicates[this_seq] = this_index
                        line_clusters.append([line_count])
                        this_index += 1
                else:
                    line_clusters.append([line_count])
                    temp1_contig_out.write(this_seq+'\n'+this_c_seq+'\n')
            else:
                print "\nWARNING: illegal fq format in line", line_count, line
                exit()
            if line_count % 3571 == 0:
                back = len(str(line_count+4))
                sys.stdout.write(str(line_count+4)+'\b'*back)
                sys.stdout.flush()
            line_count += 1
            for i in range(3):
                line = file_in.readline()
                line_count += 1
        file_in.close()
    temp1_contig_out.close()
    del name_to_line

    # dump line clusters
    len_indices = len(line_clusters)
    temp2_clusters_dir = os.path.join(output_base, 'temp.indices.2')
    temp2_indices_file_out = open(temp2_clusters_dir, 'wb')
    for this_index in xrange(len_indices):
        temp2_indices_file_out.write('\t'.join([str(x) for x in line_clusters[this_index]]))
        temp2_indices_file_out.write('\n')
    temp2_indices_file_out.close()

    if pair_end_out:
        # dump pair_to_each_other to file to save memory
        # read in again in the last
        temp3_pair_dir = os.path.join(output_base, 'temp.indices.3')
        temp3_pair_file_out = open(temp3_pair_dir, 'wb')
        for line_count_2 in xrange(0, line_count, 4):
            temp3_pair_file_out.write(str(pair_to_each_other[line_count_2])+'\n')
        temp3_pair_file_out.close()
        del pair_to_each_other

    del seq_duplicates
    del pre_reading
    # if universal_len:
    #     seq_len = len(contig_seqs[0])
    len_before_assembly = len(line_clusters)
    print len_before_assembly, "identical of", line_count/4

    return line_clusters, contig_seqs, contig_c_seqs


def pre_assembly(line_clusters, dupli_threshold, contig_seqs, contig_c_seqs):
    global word_size
    print "\nPre-assembly reads..."
    time_assembly = time.time()
    check_time = time.time()
    lines_with_duplicates = {}
    count_dupli = 0
    for j in xrange(len(line_clusters)):
        if len(line_clusters[j]) >= 2:
            lines_with_duplicates[j] = int
            count_dupli += 1
            if count_dupli == dupli_threshold:
                break
    groups_of_duplicate_lines = {}
    count_groups = 0
    print len(lines_with_duplicates), "duplicated"
    these_words = {}
    list_to_pre_assemble = list(lines_with_duplicates)
    len_pre_assemble = len(list_to_pre_assemble)
    for list_id in xrange(len(list_to_pre_assemble)):
        this_identical_read_id = list_to_pre_assemble[list_id]
        this_seq = contig_seqs[this_identical_read_id]
        this_c_seq = contig_c_seqs[this_identical_read_id]
        # if not universal_len:
        seq_len = len(this_seq)
        temp_length = seq_len - word_size
        these_group_id = set()
        this_words = []
        for i in xrange(0, temp_length+1):
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
            groups_of_duplicate_lines[new_group_id] = [{this_identical_read_id}, set(this_words)]
            for this_word in this_words:
                these_words[this_word] = new_group_id
            lines_with_duplicates[this_identical_read_id] = new_group_id
            count_groups += 1
        # belongs to one group
        elif len_groups == 1:
            this_group_id = these_group_id.pop()
            groups_of_duplicate_lines[this_group_id][0].add(this_identical_read_id)
            for this_word in this_words:
                groups_of_duplicate_lines[this_group_id][1].add(this_word)
                these_words[this_word] = this_group_id
            lines_with_duplicates[this_identical_read_id] = this_group_id
        # connect different groups
        else:
            these_group_id = list(these_group_id)
            these_group_id.sort()
            this_group_to_keep = these_group_id[0]
            # for related group to merge
            for to_merge in xrange(1, len_groups):
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
            groups_of_duplicate_lines[this_group_to_keep][0].add(this_identical_read_id)
            lines_with_duplicates[this_identical_read_id] = this_group_to_keep
        if list_id % 31 == 0:
            to_print = str(list_id+1)+'/'+str(len_pre_assemble)
            sys.stdout.write(to_print+'\b'*len(to_print))
            sys.stdout.flush()
    # print str(list_id+1)+'/'+str(len_pre_assemble)
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
    del these_words
    print "Assembly", len(groups_of_duplicate_lines), "pseudo-contigs cost", time.time()-time_assembly
    return groups_of_duplicate_lines, lines_with_duplicates


class RoundLimitException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def extending_reads(accepted_words, original_fq_dir, len_indices, fg_out_per_round, pre_assembled, groups_of_duplicate_lines, lines_with_duplicates, round_limit, output_base):
    global word_size
    accepted_contig_id = set()
    accepted_contig_id_this_round = set()
    line_to_accept = set()
    round_count = 1
    line_to_stop = len_indices - 1
    time_last_round = time.time()
    if fg_out_per_round:
        round_dir = os.path.join(output_base, "Reads_per_round")
        if not os.path.exists(round_dir):
            os.mkdir(round_dir)
    try:
        def summarise_round(acc_words, acc_contig_id_this_round, r_count, time_last_point):
            if fg_out_per_round:
                time_w_out = time.time()
                len_aw = len(acc_words)
                len_al = len(acc_contig_id_this_round)
                write_fq_results(original_fq_dir, acc_contig_id_this_round, False, False, os.path.join(round_dir, "Round."+str(r_count)))
                # clear former accepted words from memory
                del acc_words
                # then add new accepted words into memory
                acc_words = chop_seqs(read_self_fq_seq_generator([os.path.join(round_dir, "Round."+str(r_count)+'_'+str(x+1)+'.fq') for x in range(2)]), word_size)
                acc_contig_id_this_round = set()
                print "Round " + str(r_count) + ': ' + str(identical_read_id + 1) + '/' + str(len_indices) + " AL " + str(len_al) + " AW " + str(len_aw) + ' write_out cost', time.time()-time_w_out, 'round cost', time.time() - time_last_point
            else:
                print "Round " + str(r_count) + ': ' + str(identical_read_id + 1) + '/' + str(len_indices) + " AL " + str(len(acc_contig_id_this_round)) + " AW " + str(len(acc_words)) + ' round cost', time.time() - time_last_point
            if r_count == round_limit:
                raise RoundLimitException(r_count)
            r_count += 1
            return acc_words, acc_contig_id_this_round, r_count, time.time()
        while True:
            identical_reads_file = open(os.path.join(output_base, 'temp.indices.1'), 'rU')
            if pre_assembled and groups_of_duplicate_lines:
                for identical_read_id in xrange(len_indices):
                    if identical_read_id == line_to_stop:
                        raise KeyboardInterrupt
                    this_seq = identical_reads_file.readline().strip()
                    this_c_seq = identical_reads_file.readline().strip()
                    if identical_read_id not in accepted_contig_id:
                        # if not universal_len:
                        seq_len = len(this_seq)
                        temp_length = seq_len - word_size
                        if identical_read_id in line_to_accept:
                            accepted_contig_id.add(identical_read_id)
                            accepted_contig_id_this_round.add(identical_read_id)
                            line_to_accept.remove(identical_read_id)
                            for i in xrange(0, temp_length+1):
                                forward = this_seq[i:i+word_size]
                                reverse = this_c_seq[temp_length-i:seq_len-i]
                                if forward not in accepted_words:
                                    accepted_words.add(forward)
                                    accepted_words.add(reverse)
                                    line_to_stop = identical_read_id
                        else:
                            this_words = []
                            accepted = False
                            for i in xrange(0, temp_length+1):
                                forward = this_seq[i:i+word_size]
                                reverse = this_c_seq[temp_length-i:seq_len-i]
                                if forward in accepted_words:
                                    accepted = True
                                else:
                                    this_words.append(forward)
                                    this_words.append(reverse)
                            if accepted:
                                if this_words:
                                    for this_word in this_words:
                                        accepted_words.add(this_word)
                                    line_to_stop = identical_read_id
                                accepted_contig_id.add(identical_read_id)
                                accepted_contig_id_this_round.add(identical_read_id)
                                if identical_read_id in lines_with_duplicates:
                                    which_group = lines_with_duplicates[identical_read_id]
                                    for id_to_accept in groups_of_duplicate_lines[which_group]:
                                        line_to_accept.add(id_to_accept)
                                        del lines_with_duplicates[id_to_accept]
                                    line_to_accept.remove(identical_read_id)
                                    del groups_of_duplicate_lines[which_group]
                    if identical_read_id % 1213 == 0:
                        this_print = "Round " + str(round_count) + ': ' + str(identical_read_id + 1) + '/' + str(len_indices) + " AL " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words)) + ' ST ' + str(line_to_stop) + ' cost ' + str(time.time() - time_last_round)
                        sys.stdout.write(this_print + '\b' * len(this_print))
                        sys.stdout.flush()
                accepted_words, accepted_contig_id_this_round, round_count, time_last_round = summarise_round(accepted_words, accepted_contig_id_this_round, round_count, time_last_round)
            else:
                for identical_read_id in xrange(len_indices):
                    if identical_read_id == line_to_stop:
                        raise KeyboardInterrupt
                    this_seq = identical_reads_file.readline().strip()
                    this_c_seq = identical_reads_file.readline().strip()
                    if identical_read_id not in accepted_contig_id:
                        accepted = False
                        this_words = []
                        # if not universal_len:
                        seq_len = len(this_seq)
                        temp_length = seq_len - word_size
                        for i in xrange(0, temp_length+1):
                            forward = this_seq[i:i+word_size]
                            reverse = this_c_seq[temp_length-i:seq_len-i]
                            if forward in accepted_words:
                                accepted = True
                            else:
                                this_words.append(forward)
                                this_words.append(reverse)
                        if accepted:
                            if this_words:
                                for this_word in this_words:
                                    accepted_words.add(this_word)
                                line_to_stop = identical_read_id
                            accepted_contig_id.add(identical_read_id)
                            accepted_contig_id_this_round.add(identical_read_id)
                    if identical_read_id % 1213 == 0:
                        this_print = "Round " + str(round_count) + ': ' + str(identical_read_id + 1) + '/' + str(len_indices) + " AL " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words)) + ' ST ' + str(line_to_stop) + ' cost ' + str(time.time() - time_last_round)
                        sys.stdout.write(this_print + '\b' * len(this_print))
                        sys.stdout.flush()
                accepted_words, accepted_contig_id_this_round, round_count, time_last_round = summarise_round(accepted_words, accepted_contig_id_this_round, round_count, time_last_round)
            identical_reads_file.close()
    except KeyboardInterrupt:
        identical_reads_file.close()
        print "Round " + str(round_count) + ': ' + str(identical_read_id + 1) + '/' + str(len_indices) + " AL " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words)) + ' round cost', time.time() - time_last_round
        pass
    except RoundLimitException as r_lim:
        print "\nHit the round limit "+str(r_lim)+" and terminated ..."
    time_del_words = time.time()
    del accepted_words
    print "del words cost", time.time()-time_del_words
    return accepted_contig_id


def mapping_with_bowtie2(options):
    if options.seed_dir:
        build_index = subprocess.Popen("bowtie2-build --large-index "+options.seed_dir+" "+options.seed_dir+'.index', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = build_index.communicate()
        if "unrecognized option" in output:
            build_index = subprocess.Popen("bowtie2-build " + options.seed_dir + " " + options.seed_dir + '.index', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = build_index.communicate()
        index_base = options.seed_dir+'.index'
        if options.verbose_log:
            print output
    else:
        index_base = options.bowtie2_index
    print "Mapping reads to bowtie2 index ..."
    make_bowtie2 = subprocess.Popen("bowtie2 -p 4 --very-fast-local --al-conc "+os.path.join(options.output_base, "%.initial.mapped.fq")+" -q -x "+index_base+" -1 "+options.fastq_file_1+" -2 "+options.fastq_file_2+" -S "+os.path.join(options.output_base, "bowtie.sam")+" --no-unal -t", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = make_bowtie2.communicate()
    if options.verbose_log:
        print output
    print "Mapping finished."
    os.remove(os.path.join(options.output_base,"bowtie.sam"))
    return (os.path.join(options.output_base, str(x+1)+".initial.mapped.fq") for x in range(2))


def require_commands():
    usage = "python this_script.py -s seed.fasta -1 original1.fq -2 original2.fq -w 90"
    parser = OptionParser(usage=usage)
    # group1
    group_need = OptionGroup(parser, "COMMON OPTIONS", "All these arguments are required unless alternations provided")
    group_need.add_option('-1', dest='fastq_file_1', help='Input 1st fastq format file as pool')
    group_need.add_option('-2', dest='fastq_file_2', help='Input 2nd fastq format file as pool')
    group_need.add_option('-s', dest='seed_dir', help='Input fasta format file as initial seed or pre-seed (see "-m")')
    group_need.add_option('-w', dest='word_size', type=int,
                          help='Word size for extension, which have to be smaller or equal to (read length - 1)')
    group_need.add_option('-o', dest='output_base', help='Out put directory. New directory recommended because this script overwrites existing files.')
    # group2
    group_result = OptionGroup(parser, "INFLUENTIAL OPTIONS", "These option will affect the final results or serve as an alternation of the required options")
    group_result.add_option('--no_pair_end_out', dest='pair_end_out', action="store_false", default=True,
                            help='Disable output result by pairs')
    group_result.add_option('-R', dest='round_limit', type=int,
                            help='Limit running rounds (>=2).')
    group_result.add_option('-m', dest='utilize_mapping', action="store_true", default=False,
                            help='Map original reads to pre-seed (or bowtie2 index) to acquire the initial seed. This option requires bowtie2 to be installed (and thus Linux/Unix only) and is suggested when the seed is too diversed from target.')
    group_result.add_option('-b', dest='bowtie2_index', help='Input bowtie2 index base name as pre-seed. This requires "-m" to activate the mapping process.')
    # group3
    group_computational = OptionGroup(parser, "COMPUTATIONAL OPTIONS", "These options will NOT affect the final results. Take easy to pick some according your computer's flavour")
    group_computational.add_option('-A', dest='pre_assembled', type=int, default=200000,
                                   help='The maximum potential-organ reads to be pre assembled before extension. The default is 300,000. Choose a larger number to run faster but exhaust memory in the first few minutes, and choose 0 to disable.')
    group_computational.add_option('--no_out_per_round', dest='fg_out_per_round', action="store_false", default=True,
                                   help='Disable output per round. Choose to save time but exhaust memory all the whole running time')
    group_computational.add_option('--no_remove_dulplicates', dest='rm_duplicates', action="store_false", default=True,
                                   help='By default this script use unique reads to extend. Choose to disable remove duplicates, which save memory but run more slow. Note that whether choose or not will not disable the calling of replicate reads.')
    group_computational.add_option('--verbose', dest='verbose_log', action="store_true", default=False,
                                   help='Verbose output. Choose to enable verbose running log.')
    parser.add_option_group(group_need)
    parser.add_option_group(group_result)
    parser.add_option_group(group_computational)
    try:
        (options, args) = parser.parse_args()
    except Exception as e:
        print '\n######################################', e
        print '\n"-h" for more usage'
        exit()
    else:
        if not ((options.seed_dir or (options.utilize_mapping and options.bowtie2_index)) and options.fastq_file_1 and options.fastq_file_2 and options.word_size and options.output_base):
            parser.print_help()
            print '\n######################################\nERROR: Insufficient REQUIRED arguments!\n'
            exit()
        if options.seed_dir and options.bowtie2_index:
            print 'Error: Simultaneously "-s" and "-m" is not allowed!'
            exit()
        if options.utilize_mapping:
            bowtie2_in_path = commands.getstatusoutput('bowtie2')
            if bowtie2_in_path[0] == 32512:
                options.utilize_mapping = False
                if options.seed_dir:
                    print 'WARNING: bowtie2 not in the path!\nTake the seed file as initial seed.'
                else:
                    print 'Error: bowtie2 not in the path!'
                    exit()
        if not options.rm_duplicates and options.pre_assembled:
            print "WARNING: remove duplicates was inactive, so that the pre-assembly was disabled."
            options.pre_assembled = False
        if options.round_limit and options.round_limit <= 2:
            print "WARNING: illegal limit for rounds! Been set to default: unlimited."
        if not os.path.isdir(options.output_base):
            os.mkdir(options.output_base)
        return options


def main():
    time0 = time.time()
    global word_size
    print "\nGetOrganelle version 1.9.70 released by Jianjun Jin on June 28 2016." \
          "\n" \
          "\nThis script filters out organelle reads from genome skimming data by extending." \
          "\nSee https://github.com/Kinggerm/GetOrganelle for more informations." \
          "\n"
    options = require_commands()
    original_fq_dir = (options.fastq_file_1, options.fastq_file_2)

    """reading seeds"""
    print "\nReading seeds..."
    time_seed = time.time()
    if not options.utilize_mapping:
        initial_accepted_words = chop_seqs(read_fasta(options.seed_dir)[1], word_size)
    else:
        fastq_1, fastq_2 = mapping_with_bowtie2(options)
        initial_accepted_words = chop_seqs(read_self_fq_seq_generator((fastq_1, fastq_2)), word_size)
    print "reading seeds cost", time.time()-time_seed

    """reading fastq files"""
    print "\nPre-reading fastq ..."
    time_read = time.time()
    # universal_len = True
    word_size = options.word_size
    line_clusters, contig_seqs, contig_c_seqs = read_fq_and_pair_infos(original_fq_dir, options.pair_end_out, options.rm_duplicates, options.output_base)
    len_indices = len(line_clusters)
    print "Reading fastq files cost", time.time()-time_read

    """pre-assembly if asked"""
    # pre_assembly is suggested when the whole genome coverage is shallow but the organ genome coverage is deep
    if options.pre_assembled:
        groups_of_duplicate_lines, lines_with_duplicates = pre_assembly(line_clusters, options.pre_assembled, contig_seqs, contig_c_seqs)
    else:
        groups_of_duplicate_lines = None
        lines_with_duplicates = None
    del line_clusters, contig_seqs, contig_c_seqs

    """extending process"""
    print "\nExtending ..."
    accepted_contig_id = extending_reads(initial_accepted_words, original_fq_dir, len_indices, options.fg_out_per_round, options.pre_assembled, groups_of_duplicate_lines, lines_with_duplicates, options.round_limit, options.output_base)
    write_fq_results(original_fq_dir, accepted_contig_id, options.pair_end_out, True, os.path.join(options.output_base, "filtered"))

    """end"""
    os.remove(os.path.join(options.output_base, 'temp.indices.1'))
    os.remove(os.path.join(options.output_base, 'temp.indices.2'))
    if options.pair_end_out:
        os.remove(os.path.join(options.output_base, 'temp.indices.3'))
    print 'Total Calc-cost', time.time() - time0


if __name__ == '__main__':
    main()

