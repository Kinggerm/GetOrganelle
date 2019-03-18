import os
import sys
import math
import re
major_version, minor_version = sys.version_info[:2]
if major_version == 2 and minor_version >= 7:
    python_version = "2.7+"
elif major_version == 3 and minor_version >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)


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


def complementary_seqs(input_seq_iter):
    return tuple([complementary_seq(seq) for seq in input_seq_iter])


class Sequence:
    def __init__(self, label, seq):
        self.label = label
        self.seq = seq
        
    def __len__(self):
        return len(self.seq)

    def fasta_str(self, interleaved=False):
        out_str = []
        if interleaved:
            out_str.extend(['>', self.label, '\n'])
            j = interleaved
            while j < len(self):
                out_str.append(self.seq[(j - interleaved):j])
                out_str.append('\n')
                j += interleaved
            out_str.append(self.seq[(j - interleaved):j])
            out_str.append('\n')
        else:
            out_str = ['>', self.label, '\n', self.seq, "\n"]
        return "".join(out_str)


class SequenceList:
    def __init__(self, input_fasta_file=None):
        self.sequences = []
        self.interleaved = False
        if input_fasta_file:
            self.read_fasta(input_fasta_file)
        
    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        for seq in self.sequences:
            yield seq
    
    def append(self, sequence):
        self.sequences.append(sequence)
        
    def read_fasta(self, fasta_file):
        fasta_file = open(fasta_file, 'rU')
        this_line = fasta_file.readline()
        interleaved = 0
        while this_line:
            if this_line.startswith('>'):
                this_name = this_line[1:].strip()
                this_seq = ''
                this_line = fasta_file.readline()
                seq_line_count = 0
                while this_line and not this_line.startswith('>'):
                    if seq_line_count == 1:
                        interleaved = len(this_seq)
                    this_seq += this_line.strip()
                    this_line = fasta_file.readline()
                    seq_line_count += 1
                self.append(Sequence(this_name, this_seq))
            else:
                this_line = fasta_file.readline()
        fasta_file.close()
        self.interleaved = interleaved

    def write_fasta(self, fasta_file, overwrite=True):
        if not overwrite:
            while os.path.exists(fasta_file):
                fasta_file = '.'.join(fasta_file.split('.')[:-1]) + '_.' + fasta_file.split('.')[-1]
        fasta_file_handler = open(fasta_file, 'w')
        for seq in self:
            fasta_file_handler.write(seq.fasta_str(self.interleaved))
        fasta_file_handler.close()


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


def read_fasta_as_list(fasta_dir):
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
            seqs.append(list(this_seq))
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs, interleaved]


def write_fasta_with_list(out_dir, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_dir):
            out_dir = '.'.join(out_dir.split('.')[:-1])+'_.'+out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(''.join(matrix[1][i][(j-matrix[2]):j])+'\n')
                j += matrix[2]
            fasta_file.write(''.join(matrix[1][i][(j-matrix[2]):j])+'\n')
    else:
        for i in range(len(matrix[0])):
            fasta_file.write('>'+matrix[0][i]+'\n')
            fasta_file.write(''.join(matrix[1][i])+'\n')
    fasta_file.close()


# from https://github.com/Kinggerm/PersonalUtilities
# Hashing methods.
# I naively wrote this function by myself and latter found the name and description of this algorithm
# in Kurtz et al. 2001. REPuter: the manifold applications of repeat analysis on a genomic scale.
def find_exact_repeats(sequence_string, min_repeat_length, circular,
                       accepted_char=set(list("ATGCRMYKHBDVatgcrmykhbdv"))):
    word_size = min(13, min_repeat_length)
    if len(sequence_string) < min_repeat_length:
        return []
    if circular:
        long_sequence = sequence_string + sequence_string[:word_size - 1]
        here_seq = long_sequence
    else:
        here_seq = sequence_string
    here_seq_r = complementary_seq(here_seq)
    here_seq_length = len(here_seq)
    raw_seq_length = len(sequence_string)
    """create accepted id set"""
    if accepted_char:
        accepted_id = set()
        count_cal = 0
        for base_count in range(here_seq_length):
            if here_seq[base_count] not in accepted_char:
                count_cal = 0
            else:
                count_cal += 1
                if count_cal >= word_size:
                    accepted_id.add(base_count - word_size + 1)
    else:
        accepted_id = set(list(range(here_seq_length - word_size + 1)))
    #
    """initialization"""
    words_to_index = {}
    index_to_words = {}

    def add_to_words(add_index, this_forward, this_reverse):
        if this_forward in words_to_index:
            words_to_index[this_forward].add((add_index, 1))
            words_to_index[this_reverse].add((add_index, -1))
        else:
            words_to_index[this_forward] = {(add_index, 1)}
            if this_reverse == this_forward:
                words_to_index[this_reverse].add((add_index, -1))
            else:
                words_to_index[this_reverse] = {(add_index, -1)}

    for i in range(0, here_seq_length):
        if i in accepted_id:
            forward_s = here_seq[i:i + word_size]
            reverse_s = here_seq_r[here_seq_length - i - word_size: here_seq_length - i]
            add_to_words(i, forward_s, reverse_s)
            index_to_words[i] = forward_s

    """find repeats"""
    repeat_indices = set()
    for repeat_tuples in words_to_index.values():
        if len(repeat_tuples) >= 2:
            for repeat_index in repeat_tuples:
                repeat_indices.add(repeat_index[0])
    repeat_indices = sorted(list(repeat_indices))
    repeats = []
    active_connection_to_repeats = {}
    last_connection = set()
    len_indices = len(repeat_indices)

    if circular:
        if len_indices != raw_seq_length and len(repeat_indices):
            while (repeat_indices[0] - repeat_indices[-1]) % raw_seq_length == 1:
                repeat_indices.insert(0, repeat_indices.pop(-1))
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            this_connection = words_to_index[this_word]
            """part 1: dealing with old connection"""
            # Loop 1: find repeats_to_stop
            # Loop 2: delete the pointers pointing to the stopped repeats
            # Loop 3: update the pointers
            # Loop 4: add new pointers if shorter repeats should be continued with less alias
            # Loop 5: update the repeats according to pointers
            repeats_to_stop = {}
            kinds_del_from_active = set()
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_length, direction_trans)
                if candidate_new_connect not in this_connection:
                    # if one_connection in active_connection_to_repeats:
                    for repeat_kind, repeat_num in active_connection_to_repeats[one_connection]:
                        if repeat_kind in repeats_to_stop:
                            repeats_to_stop[repeat_kind][repeat_num] = one_connection
                        else:
                            repeats_to_stop[repeat_kind] = {repeat_num: one_connection}
                        kinds_del_from_active.add(repeat_kind)
            # print(repeats_to_stop)
            # Loop 2
            for repeat_kind in kinds_del_from_active:
                for now_start, now_go_to, n_direction in repeats[repeat_kind]:
                    connection_del_from_points = ((now_go_to - (word_size - 1) * (n_direction == 1)) % raw_seq_length,
                                                  n_direction)
                    if connection_del_from_points in active_connection_to_repeats:
                        count_this_group = 0
                        while count_this_group < len(active_connection_to_repeats[connection_del_from_points]):
                            if active_connection_to_repeats[connection_del_from_points][count_this_group][0]\
                                    == repeat_kind:
                                del active_connection_to_repeats[connection_del_from_points][count_this_group]
                            else:
                                count_this_group += 1
                        if not len(active_connection_to_repeats[connection_del_from_points]):
                            del active_connection_to_repeats[connection_del_from_points]
            # print("Cleared pointer", active_connection_to_repeats)
            # Loop 3
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_length, direction_trans)
                if candidate_new_connect in this_connection:
                    if one_connection in active_connection_to_repeats:
                        active_connection_to_repeats[candidate_new_connect] = []
                        for one_repeat_id in active_connection_to_repeats[one_connection]:
                            active_connection_to_repeats[candidate_new_connect].append(one_repeat_id)
                        del active_connection_to_repeats[one_connection]
            # print("Updated pointer", active_connection_to_repeats)
            # Loop 4
            for repeat_kind in repeats_to_stop:
                if len(repeats[repeat_kind]) - len(repeats_to_stop[repeat_kind]) >= 2:
                    repeat_to_be_continued = False
                    for repeat_num in range(len(repeats[repeat_kind])):
                        if repeat_num not in repeats_to_stop[repeat_kind]:
                            start_id, go_to_id, gt_direction = repeats[repeat_kind][repeat_num]
                            new_connect = (
                                (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction) % raw_seq_length,
                                gt_direction)
                            if new_connect in this_connection and new_connect not in active_connection_to_repeats:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[repeat_kind])):
                            if inside_repeat_num not in repeats_to_stop[repeat_kind]:
                                start_id, go_to_id, gt_direction = repeats[repeat_kind][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, gt_direction])
                                new_connect = (
                                    (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction) % raw_seq_length,
                                    gt_direction)
                                if new_connect in active_connection_to_repeats:
                                    active_connection_to_repeats[new_connect].append(
                                        (len(repeats) - 1, len(repeats[-1]) - 1))
                                else:
                                    active_connection_to_repeats[new_connect] = [
                                        (len(repeats) - 1, len(repeats[-1]) - 1)]
            # print("Post-add pointer", active_connection_to_repeats)
            # Loop 5
            for one_connection in active_connection_to_repeats:
                for repeat_kind, repeat_num in active_connection_to_repeats[one_connection]:
                    start_id, previous_id, this_direction = repeats[repeat_kind][repeat_num]
                    repeats[repeat_kind][repeat_num][1] += this_direction
                    repeats[repeat_kind][repeat_num][1] %= raw_seq_length
            # print("Repeats", repeats)
            """part 2: dealing with new connection"""
            for one_connection in this_connection:
                here_id, direction_trans = one_connection
                candidate_last_connect = ((here_id - direction_trans) % raw_seq_length, direction_trans)
                if candidate_last_connect not in last_connection:
                    repeats.append([])
                    for inside_connection in this_connection:
                        inside_id, inside_direction = inside_connection
                        repeats[-1].append([(inside_id + (word_size - 1) * (inside_direction == -1)) % raw_seq_length,
                                            (inside_id + (word_size - 1) * (inside_direction == 1)) % raw_seq_length,
                                            inside_direction])
                        if (inside_id, inside_direction) in active_connection_to_repeats:
                            active_connection_to_repeats[inside_connection].append(
                                (len(repeats) - 1, len(repeats[-1]) - 1))
                        else:
                            active_connection_to_repeats[inside_connection] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
                    break

            if i + 1 < len_indices:
                next_index = repeat_indices[i + 1]
            else:
                next_index = None
            if next_index != (this_index + 1) % raw_seq_length:
                active_connection_to_repeats = {}
                last_connection = set()
            else:
                last_connection = this_connection
            # the whole seq is a repeat, problematic?
            if repeats and (repeats[0][0][1] - repeats[0][0][0] + repeats[0][0][2]) % raw_seq_length == 0:
                break
    else:
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            this_connection = words_to_index[this_word]
            """part 1: dealing with old connection"""
            # Loop 1: find repeats_to_stop
            # Loop 2: delete the pointers pointing to the stopped repeats
            # Loop 3: update the pointers
            # Loop 4: add new pointers if shorter repeats should be continued with less alias
            # Loop 5: update the repeats according to pointers
            repeats_to_stop = {}
            kinds_del_from_active = set()
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = (here_id + direction_trans, direction_trans)
                if candidate_new_connect not in this_connection:
                    # if one_connection in active_connection_to_repeats:
                    for repeat_kind, repeat_num in active_connection_to_repeats[one_connection]:
                        if repeat_kind in repeats_to_stop:
                            repeats_to_stop[repeat_kind][repeat_num] = one_connection
                        else:
                            repeats_to_stop[repeat_kind] = {repeat_num: one_connection}
                        kinds_del_from_active.add(repeat_kind)
            # print("kinds_del_from_active", kinds_del_from_active)
            for repeat_kind in kinds_del_from_active:
                for now_start, now_go_to, n_direction in repeats[repeat_kind]:
                    connection_del_from_points = (now_go_to - (word_size - 1) * (n_direction == 1),
                                                  n_direction)
                    if connection_del_from_points in active_connection_to_repeats:
                        count_this_group = 0
                        while count_this_group < len(active_connection_to_repeats[connection_del_from_points]):
                            if active_connection_to_repeats[connection_del_from_points][count_this_group][0] == repeat_kind:
                                del active_connection_to_repeats[connection_del_from_points][count_this_group]
                            else:
                                count_this_group += 1
                        if not len(active_connection_to_repeats[connection_del_from_points]):
                            del active_connection_to_repeats[connection_del_from_points]
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = (here_id + direction_trans, direction_trans)
                if candidate_new_connect in this_connection:
                    if one_connection in active_connection_to_repeats:
                        active_connection_to_repeats[candidate_new_connect] = []
                        for one_repeat_id in active_connection_to_repeats[one_connection]:
                            active_connection_to_repeats[candidate_new_connect].append(one_repeat_id)
                        del active_connection_to_repeats[one_connection]
            for repeat_kind in repeats_to_stop:
                if len(repeats[repeat_kind]) - len(repeats_to_stop[repeat_kind]) >= 2:
                    repeat_to_be_continued = False
                    for repeat_num in range(len(repeats[repeat_kind])):
                        if repeat_num not in repeats_to_stop[repeat_kind]:
                            start_id, go_to_id, gt_direction = repeats[repeat_kind][repeat_num]
                            new_connect = (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction,
                                           gt_direction)
                            if new_connect in this_connection and new_connect not in active_connection_to_repeats:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[repeat_kind])):
                            if inside_repeat_num not in repeats_to_stop[repeat_kind]:
                                start_id, go_to_id, gt_direction = repeats[repeat_kind][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, gt_direction])
                                new_connect = (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction,
                                               gt_direction)
                                if new_connect in active_connection_to_repeats:
                                    active_connection_to_repeats[new_connect].append((len(repeats) - 1, len(repeats[-1]) - 1))
                                else:
                                    active_connection_to_repeats[new_connect] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
            for one_connection in active_connection_to_repeats:
                for repeat_kind, repeat_num in active_connection_to_repeats[one_connection]:
                    start_id, previous_id, this_direction = repeats[repeat_kind][repeat_num]
                    repeats[repeat_kind][repeat_num][1] += this_direction
            """part 2: dealing with new connection"""
            for one_connection in this_connection:
                here_id, direction_trans = one_connection
                candidate_last_connect = (here_id - direction_trans, direction_trans)
                if candidate_last_connect not in last_connection:
                        repeats.append([])
                        for inside_connection in this_connection:
                            inside_id, inside_direction = inside_connection
                            repeats[-1].append([inside_id + (word_size - 1) * (inside_direction == -1),
                                                inside_id + (word_size - 1) * (inside_direction == 1),
                                                inside_direction])
                            if (inside_id, inside_direction) in active_connection_to_repeats:
                                active_connection_to_repeats[inside_connection].append((len(repeats) - 1, len(repeats[-1]) - 1))
                            else:
                                active_connection_to_repeats[inside_connection] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
                        break

            if i + 1 < len_indices:
                next_index = repeat_indices[i + 1]
            else:
                next_index = None
            if next_index != this_index + 1:
                active_connection_to_repeats = {}
                last_connection = set()
            else:
                last_connection = this_connection

    """aftertreatment"""
    # 1.delete repeated repeats
    final_repeat = []
    repeat_dicts = set()
    for repeat_group in repeats:
        if tuple(repeat_group[0]) not in repeat_dicts:
            for single_repeat in repeat_group:
                here_start, here_end, here_direction = single_repeat
                repeat_dicts.add((here_start, here_end, here_direction))
                repeat_dicts.add((here_end, here_start, -here_direction))
            final_repeat.append(repeat_group)
        else:
            continue

    # 2.delete small repeats
    count_group__ = 0
    while count_group__ < len(final_repeat):
        here_start, here_end, here_direction = final_repeat[count_group__][0]
        if not (here_start == 0 and here_end == raw_seq_length - 1) and \
                ((here_end - here_start + here_direction)*here_direction) % raw_seq_length < min_repeat_length:
            del final_repeat[count_group__]
        else:
            count_group__ += 1

    # 3.reorder repeats according to occurrence
    for group_to_sort in range(len(final_repeat)):
        start, end, direction = final_repeat[group_to_sort][0]
        if start == end:
            this_len = 1
        else:
            if (end - start) * direction > 0:
                this_len = (end - start) * direction + 1
            else:
                this_len = (start, end)[direction == 1] + raw_seq_length - (start, end)[direction != 1] + 1
        # transform into dict
        final_repeat[group_to_sort] = [{"start": start, "end": end, "direction": direction, "length": this_len}
                                       for start, end, direction in final_repeat[group_to_sort]]
    # 4.sort according to length: from longest to shortest
    final_repeat.sort(key=lambda x: -x[0]["length"])
    return final_repeat


def reverse_repeats_info(repeats):
    new_repeats = []
    for rep in repeats:
        new_repeats.append({"start": rep["end"], "end": rep["start"],
                            "direction": -rep["direction"], "length": rep["length"]})
    return new_repeats


def re_linear_circular_seqs(sequence, minimum_len_for_flip_flop_recombination=2000, log_handler=None):
    raw_rev = complementary_seq(sequence)
    len_seq = len(sequence)
    min_len = minimum_len_for_flip_flop_recombination
    repeats_set = find_exact_repeats(sequence, min_len, True)
    if repeats_set:
        # currently only for plastomes with inverted repeats, us `Arachis` to generalize this function
        if len(repeats_set) > 1 or len(repeats_set[0]) != 2:
            if log_handler:
                log_handler.error("Currently only for plastomes with inverted/directed repeats! "
                                  "Please inform the author if you want to extend the application.\n"
                                  "Re-linearizing disabled!")
            else:
                sys.stdout.write("Error: Currently only for plastomes with inverted/directed repeats! "
                                 "Please inform the author if you want to extend the application.\n"
                                 "Re-linearizing disabled!\n")
            return sequence
            # raise NotImplementedError("Currently only for plastomes with inverted/directed repeats! "
            #                           "Please inform the author if you want to extend the application.")
        else:
            longest_repeats = repeats_set[0]
            # Sorting makes:
            # direct1==1 and direct2==-1
            # start1 be the smallest forward start
            ir_locations_1 = sorted(longest_repeats, key=lambda x: (-x["direction"], x["start"]))
            ir_locations_2 = sorted(reverse_repeats_info(ir_locations_1), key=lambda x: (-x["direction"], x["start"]))
            ir_locations = sorted([ir_locations_1, ir_locations_2],
                                  key=lambda x: (-max([y["direction"] for y in x]), x[0]["start"]))[0]
            start1, end1, direct1 = ir_locations[0]["start"], ir_locations[0]["end"], ir_locations[0]["direction"]
            start2, end2, direct2 = ir_locations[1]["start"], ir_locations[1]["end"], ir_locations[1]["direction"]
            if ir_locations[0]["direction"] == ir_locations[1]["direction"]:
                # cross the end, meaning site:seq_len in (DR1)
                if end1 < start1:
                    if end2 >= start1:
                        if log_handler:
                            log_handler.error("Currently only for plastomes with inverted/directed repeats! "
                                              "Please inform the author if you want to extend the application.\n"
                                              "Re-linearizing disabled!")
                        else:
                            sys.stdout.write("Error: Currently only for plastomes with inverted/directed repeats! "
                                             "Please inform the author if you want to extend the application.\n"
                                             "Re-linearizing disabled!\n")
                        return sequence
                    else:
                        forward_seq = sequence[end2 + 1:] + sequence[:end2 + 1]
                        reverse_seq = raw_rev[len_seq - start1:] + raw_rev[:len_seq - start1]
                        return sorted([forward_seq, reverse_seq])[0]
                elif end2 < start2:
                    if end2 >= start1:
                        if end1 >= start2:
                            if log_handler:
                                log_handler.error("Currently only for plastomes with inverted/directed repeats! "
                                                  "Please inform the author if you want to extend the application.\n"
                                                  "Re-linearizing disabled!")
                            else:
                                sys.stdout.write("Error: Currently only for plastomes with inverted/directed repeats! "
                                                 "Please inform the author if you want to extend the application.\n"
                                                 "Re-linearizing disabled!\n")
                            return sequence
                        else:
                            forward_seq = sequence[end1 + 1:] + sequence[:end1 + 1]
                            reverse_seq = raw_rev[len_seq - start2:] + raw_rev[:len_seq - start2]
                            return sorted([forward_seq, reverse_seq])[0]
                    elif end1 >= start2:
                        forward_seq = sequence[end2 + 1:] + sequence[:end2 + 1]
                        reverse_seq = raw_rev[len_seq - start1:] + raw_rev[:len_seq - start1]
                        return sorted([forward_seq, reverse_seq])[0]
                    else:
                        if start1 - end2 - 1 > start2 - end1 - 1:
                            forward_seq = sequence[end2 + 1:] + sequence[:end2 + 1]
                            reverse_seq = raw_rev[len_seq - start1:] + raw_rev[:len_seq - start1]
                            return sorted([forward_seq, reverse_seq])[0]
                        elif start1 - end2 - 1 < start2 - end1 - 1:
                            forward_seq = sequence[end1 + 1:] + sequence[:end1 + 1]
                            reverse_seq = raw_rev[len_seq - start2:] + raw_rev[:len_seq - start2]
                            return sorted([forward_seq, reverse_seq])[0]
                        else:
                            forward_seq_1 = sequence[end2 + 1:] + sequence[:end2 + 1]
                            reverse_seq_1 = raw_rev[len_seq - start1:] + raw_rev[:len_seq - start1]
                            forward_seq_2 = sequence[end1 + 1:] + sequence[:end1 + 1]
                            reverse_seq_2 = raw_rev[len_seq - start2:] + raw_rev[:len_seq - start2]
                            return sorted([forward_seq_1, reverse_seq_1, forward_seq_2, reverse_seq_2])[0]
                else:
                    # if start2 - end1 - 1 > len(sequence) + start1 - end2 - 1:
                    forward_seq_1 = [sequence[end1 + 1:], sequence[:end1 + 1]]
                    reverse_seq_1 = [raw_rev[len_seq - start2:], raw_rev[:len_seq - start2]]
                    # elif start2 - end1 - 1 < len(sequence) + start1 - end2 - 1:
                    forward_seq_2 = [sequence[end2 + 1:], sequence[:end2 + 1]]
                    reverse_seq_2 = [raw_rev[len_seq - start1:], raw_rev[:len_seq - start1]]
                    best_res = sorted([forward_seq_1, forward_seq_2, reverse_seq_1, reverse_seq_2],
                                      key=lambda x: (-len(x[0]), x))[0]
                    return "".join(best_res)
            else:
                # cross the end, meaning site:seq_len in (IR1)
                if end1 < start1:
                    # seq_len >= start2 >= start1
                    if start2 >= start1:
                        if log_handler:
                            log_handler.error("Currently only for plastomes with inverted/directed repeats! "
                                              "Please inform the author if you want to extend the application.\n"
                                              "Re-linearizing disabled!")
                        else:
                            sys.stdout.write("Error: Currently only for plastomes with inverted/directed repeats! "
                                             "Please inform the author if you want to extend the application.\n"
                                             "Re-linearizing disabled!\n")
                        return sequence
                    else:
                        forward_seq = sequence[start2 + 1:] + sequence[:start2 + 1]
                        reverse_seq = raw_rev[len_seq - start1:] + raw_rev[:len_seq - start1]
                        return sorted([forward_seq, reverse_seq])[0]
                elif start2 < end2:
                    if start2 >= start1:
                        if end1 >= end2:
                            if log_handler:
                                log_handler.error("Currently only for plastomes with inverted/directed repeats! "
                                                  "Please inform the author if you want to extend the application.\n"
                                                  "Re-linearizing disabled!")
                            else:
                                sys.stdout.write("Error: Currently only for plastomes with inverted/directed repeats! "
                                                 "Please inform the author if you want to extend the application.\n"
                                                 "Re-linearizing disabled!\n")
                            return sequence
                        else:
                            forward_seq = sequence[end1 + 1:] + sequence[:end1 + 1]
                            reverse_seq = raw_rev[len_seq - end2:] + raw_rev[:len_seq - end2]
                            return sorted([forward_seq, reverse_seq])[0]
                    else:
                        # if start1 - start2 - 1 > end2 - end1 - 1:
                        forward_seq_1 = [sequence[start2 + 1:start1],
                                         sequence[start1:end1 + 1],
                                         sequence[end1 + 1:end2],
                                         sequence[end2:] + sequence[:start2 + 1]]
                        forward_seq_2 = [sequence[start2 + 1:start1],
                                         sequence[start1:end1 + 1],
                                         complementary_seq(sequence[end1 + 1:end2]),
                                         sequence[end2:] + sequence[:start2 + 1]]
                        reverse_seq_1 = [raw_rev[len_seq - start1: len_seq - start2 - 1],
                                         raw_rev[len_seq - start2 - 1:] + raw_rev[:len_seq - end2],
                                         raw_rev[len_seq - end2: len_seq - end1 - 1],
                                         raw_rev[len_seq - end1 - 1: len_seq - start1]]
                        reverse_seq_2 = [raw_rev[len_seq - start1: len_seq - start2 - 1],
                                         raw_rev[len_seq - start2 - 1:] + raw_rev[:len_seq - end2],
                                         complementary_seq(raw_rev[len_seq - end2: len_seq - end1 - 1]),
                                         raw_rev[len_seq - end1 - 1: len_seq - start1]]
                        # elif start1 - start2 - 1 < end2 - end1 - 1:
                        forward_seq_3 = [sequence[end1 + 1: end2],
                                         sequence[end2:] + sequence[:start2 + 1],
                                         sequence[start2 + 1:start1],
                                         sequence[start1:end1 + 1]]
                        forward_seq_4 = [sequence[end1 + 1: end2],
                                         sequence[end2:] + sequence[:start2 + 1],
                                         complementary_seq(sequence[start2 + 1:start1]),
                                         sequence[start1:end1 + 1]]
                        reverse_seq_3 = [raw_rev[len_seq - end2: len_seq - end1 - 1],
                                         raw_rev[len_seq - end1 - 1: len_seq - start1],
                                         raw_rev[len_seq - start1: len_seq - start2 - 1],
                                         raw_rev[len_seq - start2 - 1:] + raw_rev[:len_seq - end2]]
                        reverse_seq_4 = [raw_rev[len_seq - end2: len_seq - end1 - 1],
                                         raw_rev[len_seq - end1 - 1: len_seq - start1],
                                         complementary_seq(raw_rev[len_seq - start1: len_seq - start2 - 1]),
                                         raw_rev[len_seq - start2 - 1:] + raw_rev[:len_seq - end2]]
                        best_res = sorted([forward_seq_1, forward_seq_2, forward_seq_3, forward_seq_4,
                                           reverse_seq_1, reverse_seq_2, reverse_seq_3, reverse_seq_4],
                                          key=lambda x: (-len(x[0]), x))[0]
                        return "".join(best_res)
                else:
                    # if len(sequence) + start1 - start2 - 1 > end2 - end1 - 1:
                    forward_seq_1 = [sequence[start2 + 1:] + sequence[:start1],
                                     sequence[start1:end1 + 1],
                                     sequence[end1 + 1:end2],
                                     sequence[end2:start2 + 1]]
                    forward_seq_2 = [sequence[start2 + 1:] + sequence[:start1],
                                     sequence[start1:end1 + 1],
                                     complementary_seq(sequence[end1 + 1:end2]),
                                     sequence[end2:start2 + 1]]
                    reverse_seq_1 = [raw_rev[len_seq - start1:] + raw_rev[:len_seq - start2 - 1],
                                     raw_rev[len_seq - start2 - 1:len_seq - end2],
                                     raw_rev[len_seq - end2: len_seq - end1 - 1],
                                     raw_rev[len_seq - end1 - 1: len_seq - start1]]
                    reverse_seq_2 = [raw_rev[len_seq - start1:] + raw_rev[:len_seq - start2 - 1],
                                     raw_rev[len_seq - start2 - 1:len_seq - end2],
                                     complementary_seq(raw_rev[len_seq - end2: len_seq - end1 - 1]),
                                     raw_rev[len_seq - end1 - 1: len_seq - start1]]
                    # elif len(sequence) + start1 - start2 - 1 < end2 - end1 - 1:
                    forward_seq_3 = [sequence[end1 + 1: end2],
                                     sequence[end2:start2 + 1],
                                     sequence[start2 + 1:] + sequence[:start1],
                                     sequence[start1:end1 + 1]]
                    forward_seq_4 = [sequence[end1 + 1: end2],
                                     sequence[end2:start2 + 1],
                                     complementary_seq(sequence[start2 + 1:] + sequence[:start1]),
                                     sequence[start1:end1 + 1]]
                    reverse_seq_3 = [raw_rev[len_seq - end2: len_seq - end1 - 1],
                                     raw_rev[len_seq - end1 - 1: len_seq - start1],
                                     raw_rev[len_seq - start1:] + raw_rev[:len_seq - start2 - 1],
                                     raw_rev[len_seq - start2 - 1:len_seq - end2]]
                    reverse_seq_4 = [raw_rev[len_seq - end2: len_seq - end1 - 1],
                                     raw_rev[len_seq - end1 - 1: len_seq - start1],
                                     complementary_seq(raw_rev[len_seq - start1:] + raw_rev[:len_seq - start2 - 1]),
                                     raw_rev[len_seq - start2 - 1:len_seq - end2]]
                    best_res = sorted([forward_seq_1, forward_seq_2, forward_seq_3, forward_seq_4,
                                       reverse_seq_1, reverse_seq_2, reverse_seq_3, reverse_seq_4],
                                      key=lambda x: (-len(x[0]), x))[0]
                    return "".join(best_res)
    else:
        extended_template = sequence + sequence[:min_len - 1]
        extended_temp_rev = raw_rev + raw_rev[:min_len - 1]
        words = []
        for go_to_base in range(len_seq):
            words.append([extended_template[go_to_base: go_to_base + min_len], go_to_base, True])
            words.append([extended_temp_rev[go_to_base: go_to_base + min_len], go_to_base, False])
        words.sort()
        initial_word, initial_base, initial_direction = words[0]
        if initial_direction:
            return sequence[initial_base:] + sequence[:initial_base]
        else:
            return raw_rev[initial_base:] + raw_rev[:initial_base]


def find_string_difference(this_string, this_reference, dynamic_span=2.0):
    this_string = this_string.lower()
    this_reference = this_reference.lower()
    len_str = len(this_string)
    len_ref = len(this_reference)
    if dynamic_span == 0:
        difference = sum([this_string[i] != this_reference[i]
                          for i in range(min(len_ref, len_str))]) + abs(len_ref - len_str)
        proper_end = this_string[-1] == this_reference[-1]
        return difference, proper_end
    else:
        dynamic_span = max(abs(len(this_string)-len(this_reference))+1, dynamic_span)
        no_match_penal = this_string[0] != this_reference[0]
        this_matrix = {(0, 0): {"state": no_match_penal}}
        # calculate the first column
        for i in range(1, min(int(math.ceil(dynamic_span))+1, len_str)):
            this_matrix[(i, 0)] = {"right_out": no_match_penal+i, "state": no_match_penal+i}
        # calculate the first line
        for j in range(1, min(int(math.ceil(dynamic_span))+1, len_ref)):
            this_matrix[(0, j)] = {"right_out": no_match_penal+j, "state": no_match_penal+j}
        # calculate iteratively
        start = 0
        for i in range(1, len_str):
            start = max(1, int(i-dynamic_span))
            end = min(len_ref, int(math.ceil(i+dynamic_span)))
            # start: no right_in
            no_match_penal = this_string[i] != this_reference[start]
            this_matrix[(i, start)] = {"diagonal_out": this_matrix[(i-1, start-1)]["state"] + no_match_penal,
                                       "down_out": this_matrix[(i-1, start)]["state"] + 1}
            this_matrix[(i, start)]["state"] = min(this_matrix[(i, start)].values())
            # middle
            for j in range(start+1, end-1):
                no_match_penal = this_string[i] != this_reference[j]
                this_matrix[(i, j)] = {"diagonal_out": this_matrix[(i-1, j-1)]["state"] + no_match_penal,
                                       "down_out": this_matrix[(i-1, j)]["state"] + 1,
                                       "right_out": this_matrix[(i, j-1)]["state"] + 1}
                this_matrix[(i, j)]["state"] = min(this_matrix[(i, j)].values())
            # end
            no_match_penal = this_string[i] != this_reference[end - 1]
            this_matrix[(i, end-1)] = {"diagonal_out": this_matrix[(i-1, end-2)]["state"] + no_match_penal}
            if (i, end-2) in this_matrix:
                this_matrix[(i, end-1)]["right_out"] = this_matrix[(i, end-2)]["state"] + 1
            this_matrix[(i, end-1)]["state"] = min(this_matrix[(i, end-1)].values())
        # print time.time()-time0
        difference = this_matrix[(len_str-1, len_ref-1)]["state"]
        proper_end = True
        for j in range(start, len_ref):
            try:
                if this_matrix[(len_str-1, j)]["state"] < difference:
                    proper_end = False
                    break
            except KeyError:
                pass
        for i in range(max(0, len_str-len_ref+start), len_str):
            try:
                if this_matrix[(i, len_ref-1)]["state"] < difference:
                    proper_end = False
                    break
            except KeyError:
                pass
        return difference, proper_end


DEGENERATE_BASES = {"N", "V", "H", "D", "B", "Y", "R", "K", "M", "S", "W"}


DEGENERATE_DICT = {  # degenerate
    "N": ["A", "C", "G", "T"],
    "V": ["A", "C", "G"], "H": ["A", "C", "T"], "D": ["A", "G", "T"], "B": ["C", "G", "T"],
    "Y": ["C", "T"], "R": ["A", "G"], "K": ["G", "T"], "M": ["A", "C"],
    "S": ["C", "G"], "W": ["A", "T"],
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
    #  consensus
    ('A', 'C', 'G', 'T'): 'N',
    ('A', 'C', 'G'): 'V', ('A', 'C', 'T'): 'H', ('A', 'G', 'T'): 'D', ('C', 'G', 'T'): 'B',
    ('C', 'T'): 'Y', ('A', 'G'): 'R', ('G', 'T'): 'K', ('A', 'C'): 'M',
    ('C', 'G'): 'S', ('A', 'T'): 'W',
    ('A',): 'A', ('C',): 'C', ('G',): 'G', ('T',): 'T'}

DEGENERATE_DICT_DIGIT = {  # degenerate
    "N": [1, 2, 4, 8],
    "V": [1, 2, 4], "H": [1, 2, 8], "D": [1, 4, 8], "B": [2, 4, 8],
    "Y": [2, 8], "R": [1, 4], "K": [4, 8], "M": [1, 2],
    "S": [2, 4], "W": [1, 8],
    "A": [1], "C": [2], "G": [4], "T": [8],
    #  consensus
    15: 'N',
    7: 'V', 11: 'H', 13: 'D', 14: 'B',
    10: 'Y', 5: 'R', 12: 'K', 3: 'M',
    6: 'S', 9: 'W',
    1: "A", 2: "C", 4: "G", 8: "T"}


def generate_consensus(*seq_str):
    consensus_res = []
    seq_num = len(seq_str)
    for go_to_base in range(len(seq_str[0])):
        this_base_set = set()
        for go_to_seq in range(seq_num):
            this_base_set.update(DEGENERATE_DICT_DIGIT.get(seq_str[go_to_seq][go_to_base], []))
        consensus_res.append(DEGENERATE_DICT_DIGIT[sum(this_base_set)])
    return "".join(consensus_res)


def split_seq_by_quality_pattern(sequence, quality_str, low_quality_pattern, min_sub_seq=0):
    seq_list = []
    start = 0
    for splicer in re.finditer(low_quality_pattern, quality_str):
        end, new_start = splicer.span()
        seq_part = sequence[start:end]
        if len(seq_part) > min_sub_seq:
            seq_list.append(seq_part)
        start = new_start
    seq_part = sequence[start:]
    if len(seq_part) > min_sub_seq:
        seq_list.append(seq_part)
    return tuple(seq_list)


def fq_seq_simple_generator(fq_dir_list, go_to_line=1, split_pattern=None, min_sub_seq=0, max_n_reads=float("inf")):
    if not ((type(fq_dir_list) is list) or (type(fq_dir_list) is tuple)):
        fq_dir_list = [fq_dir_list]
    max_n_lines = 4 * max_n_reads
    if split_pattern and len(split_pattern) > 2:
        for fq_dir in fq_dir_list:
            count = 0
            with open(fq_dir, 'rU') as fq_handler:
                seq_line = fq_handler.readline()
                while seq_line:
                    if count % 4 == go_to_line:
                        fq_handler.readline()
                        quality_str = fq_handler.readline()[:-1]
                        count += 2
                        yield split_seq_by_quality_pattern(seq_line[:-1], quality_str, split_pattern, min_sub_seq)
                    count += 1
                    if count >= max_n_lines:
                        break
                    seq_line = fq_handler.readline()
    else:
        for fq_dir in fq_dir_list:
            count = 0
            with open(fq_dir, 'rU') as fq_handler:
                for fq_line in fq_handler:
                    if count % 4 == go_to_line:
                        yield fq_line[:-1]
                    if count >= max_n_lines:
                        break
                    count += 1


def chop_seqs(seq_generator_or_list, word_size, mesh_size=1):
    return_words = set()
    for seed in seq_generator_or_list:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            temp_length = this_seq_len - word_size
            for i in range(0, this_seq_len - word_size + 1, mesh_size):
                return_words.add(seed[i:i + word_size])
                return_words.add(cpt_seed[temp_length - i:this_seq_len - i])
    return return_words


def chop_seqs_as_empty_dict(seq_generator_or_list, word_size, mesh_size=1):
    return_words = dict()
    for seed in seq_generator_or_list:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            temp_length = this_seq_len - word_size
            for i in range(0, this_seq_len - word_size + 1, mesh_size):
                return_words[seed[i:i + word_size]] = 0
                return_words[cpt_seed[temp_length - i:this_seq_len - i]] = 0
    return return_words


def chop_seq_list(seq_generator_or_list, word_size, mesh_size=1):
    return_words = set()
    for seed in seq_generator_or_list:
        for seq_part in seed:
            this_seq_len = len(seq_part)
            if this_seq_len >= word_size:
                cpt_seed = complementary_seq(seq_part)
                temp_length = this_seq_len - word_size
                for i in range(0, this_seq_len - word_size + 1, mesh_size):
                    return_words.add(seq_part[i:i + word_size])
                    return_words.add(cpt_seed[temp_length - i:this_seq_len - i])
    return return_words


def counting_words(seq_generator, words_initial_dict, word_size):
    for seq in seq_generator:
        for i in range(0, len(seq) - word_size + 1):
            this_word = seq[i: i + word_size]
            if this_word in words_initial_dict:
                words_initial_dict[this_word] += 1
    return words_initial_dict


def check_fasta_seq_names(original_fas, log=None):
    fas_matrix = read_fasta(original_fas)
    short_names = [s_n.split(" ")[0] for s_n in fas_matrix[0]]
    if len(short_names) == len(set(short_names)):
        return original_fas
    else:
        new_fas = original_fas + ".modified"
        if os.path.exists(new_fas):
            return check_fasta_seq_names(new_fas)
        else:
            if log:
                log.info("Setting '-s " + new_fas + "'")
            for go_to, name in enumerate(fas_matrix[0]):
                fas_matrix[0][go_to] = name.replace(" ", "_")
            if len(fas_matrix[0]) != len(set(fas_matrix[0])):
                existed = {}
                for go_to, name in enumerate(fas_matrix[0]):
                    if name in existed:
                        existed[name] += 1
                    else:
                        existed[name] = 1
                    fas_matrix[0][go_to] += "-" + str(existed[name])
            write_fasta(new_fas, fas_matrix, True)
            return new_fas


def get_read_len_mean_max_count(fq_files, maximum_n_reads, sampling_percent=1.):
    read_lengths = []
    all_count = 0
    if sampling_percent == 1:
        for fq_f in fq_files:
            count_r = 0
            for seq in fq_seq_simple_generator(fq_f):
                count_r += 1
                read_lengths.append(len(seq.strip("N")))
                if count_r >= maximum_n_reads:
                    break
            all_count += count_r
    else:
        sampling_percent = int(1 / sampling_percent)
        for fq_f in fq_files:
            count_r = 0
            for seq in fq_seq_simple_generator(fq_f):
                count_r += 1
                if count_r % sampling_percent == 0:
                    read_lengths.append(len(seq.strip("N")))
                if count_r >= maximum_n_reads:
                    break
            all_count += count_r
    return sum(read_lengths)/len(read_lengths), max(read_lengths), all_count


def get_read_quality_info(fq_files, maximum_n_reads, min_quality_score, log,
                          maximum_ignore_percent=0.05, sampling_percent=0.1):
    sampling_percent = int(1 / sampling_percent)
    all_quality_chars_list = []
    record_fq_beyond_read_num_limit = []
    for fq_f in fq_files:
        count_r = 0
        record_fq_beyond_read_num_limit.append(False)
        this_fq_generator = fq_seq_simple_generator(fq_f, go_to_line=3)
        for quality_str in this_fq_generator:
            if count_r % sampling_percent == 0:
                all_quality_chars_list.append(quality_str)
            count_r += 1
            if count_r >= maximum_n_reads:
                break
        for quality_str in this_fq_generator:
            if quality_str:
                log.info("Number of reads exceeded " + str(int(maximum_n_reads)) + " in " + os.path.basename(fq_f)
                         + ", only top " + str(int(maximum_n_reads))
                         + " reads are used in downstream analysis.")
                record_fq_beyond_read_num_limit[-1] = True
            break
    all_quality_chars = "".join(all_quality_chars_list)
    len_quality_chars_total = float(len(all_quality_chars))
    max_quality = max(all_quality_chars)
    min_quality = min(all_quality_chars)
    max_quality = ord(max_quality)
    min_quality = ord(min_quality)
    decision_making = []
    for type_name, char_min, char_max, score_min, score_max in [("Sanger", 33, 73, 0, 40),
                                                                ("Solexa", 59, 104, -5, 40),
                                                                ("Illumina 1.3+", 64, 104, 0, 40),
                                                                ("Illumina 1.5+", 67, 105, 3, 41),
                                                                ("Illumina 1.8+", 33, 74, 0, 41)]:
        decision_making.append((type_name, char_min, char_max, score_min, score_max,
                                (max_quality - char_max) ** 2 + (min_quality - char_min) ** 2))
    the_form, the_c_min, the_c_max, the_s_min, the_s_max, deviation = sorted(decision_making, key=lambda x: x[-1])[0]
    log.info("Identified quality encoding format = " + the_form)
    if max_quality > the_c_max:
        log.warning("Max quality score " + repr(chr(max_quality)) +
                    "(" + str(max_quality) + ":" + str(max_quality - the_c_min + the_s_min) +
                    ") in your fastq file exceeds the usual boundary " + str((the_c_min, the_c_max)))
    if min_quality < the_c_min:
        log.warning("Min quality score " + repr(chr(min_quality)) +
                    "(" + str(min_quality) + ":" + str(min_quality - the_c_min + the_s_min) +
                    ") in your fastq file is under the usual lower boundary " + str((the_c_min, the_c_max)))

    # increase --min-quality-score if ignoring too much data.
    low_quality_score = min(the_c_min, min_quality)
    ignore_percent = 0
    trimmed_quality_chars = []

    while low_quality_score < the_c_min + min_quality_score - the_s_min:
        ignore_this_time = all_quality_chars.count(chr(low_quality_score)) / len_quality_chars_total
        ignore_percent += ignore_this_time
        if ignore_percent > maximum_ignore_percent:
            ignore_percent -= ignore_this_time
            break
        else:
            trimmed_quality_chars.append(chr(low_quality_score))
        low_quality_score += 1
    if low_quality_score < the_c_min + min_quality_score - the_s_min:
        log.info("Resetting '--min-quality-score " + str(low_quality_score + the_s_min - the_c_min) + "'")
    if trimmed_quality_chars:
        log.info("Trimming bases with qualities (" + "%.2f" % (ignore_percent * 100) + "%): " +
                 str(ord(min("".join(trimmed_quality_chars)))) + ".." + str(ord(max("".join(trimmed_quality_chars)))) +
                 "  " + "".join(trimmed_quality_chars))
    trimmed_quality_chars = "".join(trimmed_quality_chars)

    # calculate mean error rate
    all_quality_char_dict = {in_quality_char: all_quality_chars.count(in_quality_char)
                             for in_quality_char in set(all_quality_chars)}
    error_prob_func = chose_error_prob_func[the_form]
    mean_error_rate = sum([error_prob_func(in_quality_char) * all_quality_char_dict[in_quality_char]
                           for in_quality_char in all_quality_char_dict]) / len_quality_chars_total
    log.info("Mean error rate = " + str(round(mean_error_rate, 4)))

    return "[" + trimmed_quality_chars + "]", mean_error_rate, record_fq_beyond_read_num_limit  # , post_trimming_mean


def get_cover_range(all_coverages, guessing_percent=0.07):
    # empirical value
    center_p = 1 - guessing_percent
    all_coverages.sort()
    total_len = len(all_coverages)
    start_cov_pos = max(int((center_p - guessing_percent / 2) * total_len), 0)
    end_cov_pos = min(int((center_p + guessing_percent / 2) * total_len), total_len)
    used_len = end_cov_pos-start_cov_pos
    used_points = all_coverages[start_cov_pos:end_cov_pos]
    mean = sum(used_points) / used_len
    up_dev = (sum([(single_p-mean)**2 for single_p in used_points if single_p >= mean])/used_len) ** 0.5
    down_dev = (sum([(single_p-mean)**2 for single_p in used_points if single_p <= mean])/used_len) ** 0.5
    return round(mean - down_dev, 2), round(mean, 2), round(mean + up_dev, 2)


# Q = -10*log10(P)
def score_2_error_prob_sanger_quality(quality_score):
    return 10**(-quality_score / 10.)


# Q = -10*log10(P/(1-P))
def score_2_error_prob_solexa_quality(quality_score):
    this_sanger_prob = score_2_error_prob_sanger_quality(quality_score)
    return this_sanger_prob/(1 + this_sanger_prob)


# solexa_quality_str = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
all_quality_str = set("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
phred_33_quality_trans = {in_trans_char: ord(in_trans_char) - 33
                          for in_trans_char in all_quality_str}
phred_64_quality_trans = {in_trans_char: ord(in_trans_char) - 64
                          for in_trans_char in all_quality_str}


def sanger_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_33_quality_trans.get(quality_char, 0))


def solexa_error_prob(quality_char):
    return score_2_error_prob_solexa_quality(phred_64_quality_trans.get(quality_char, -31))


def illumina_1_3_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_64_quality_trans.get(quality_char, -31))


def illumina_1_5_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_64_quality_trans.get(quality_char, -31))


def illumina_1_8_error_prob(quality_char):
    return score_2_error_prob_sanger_quality(phred_33_quality_trans.get(quality_char, 0))


chose_error_prob_func = {"Sanger": sanger_error_prob,
                         "Solexa": solexa_error_prob,
                         "Illumina 1.3+": illumina_1_3_error_prob,
                         "Illumina 1.5+": illumina_1_5_error_prob,
                         "Illumina 1.8+": illumina_1_8_error_prob}