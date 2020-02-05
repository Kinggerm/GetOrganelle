import os
import sys
import math
import re
import random

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

    PAIRING_TRANSLATOR = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


    def complementary_seq(input_seq):
        return string.translate(input_seq, PAIRING_TRANSLATOR)[::-1]

else:
    # python3
    PAIRING_TRANSLATOR = str.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


    def complementary_seq(input_seq):
        return str.translate(input_seq, PAIRING_TRANSLATOR)[::-1]


def complementary_seqs(input_seq_iter):
    return tuple([complementary_seq(seq) for seq in input_seq_iter])


CLASSIC_START_CODONS = {"ATG", "ATC", "ATA", "ATT", "GTG", "TTG"}
CLASSIC_STOP_CODONS = {"TAA", "TAG", "TGA"}


class Sequence(object):
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


class SequenceList(object):
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
        fasta_file = open(fasta_file, 'r')
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

    def write_fasta(self, fasta_file, overwrite=True, interleaved=None):
        if not overwrite:
            while os.path.exists(fasta_file):
                fasta_file = '.'.join(fasta_file.split('.')[:-1]) + '_.' + fasta_file.split('.')[-1]
        this_interleaved = self.interleaved if interleaved is None else interleaved
        fasta_file_handler = open(fasta_file, 'w')
        for seq in self:
            fasta_file_handler.write(seq.fasta_str(this_interleaved))
        fasta_file_handler.close()


class SeqKmerIndexer(object):
    def __init__(self, sequence, kmer, is_circular=False):
        self.__kmer = kmer
        if is_circular:
            self.__sequence = sequence + sequence[:kmer - 1]
            self.__len = len(sequence)
        else:
            self.__sequence = sequence
            self.__len = len(sequence) - kmer + 1

    def __len__(self):
        return self.__len

    def __getitem__(self, item):
        if -self.__len - 1 < item < -1:
            return self.__sequence[1 - self.__kmer + item: 1 + item]
        elif item == -1:
            return self.__sequence[1 - self.__kmer + item:]
        elif -1 < item < self.__len:
            return self.__sequence[item: item + self.__kmer]
        else:
            raise IndexError("index out of range")

    def __enumerate__(self):
        for i in range(self.__len):
            yield i, self.__getitem__(i)


def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'r')
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


def write_fasta(out_file, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_file):
            out_file = '.'.join(out_file.split('.')[:-1]) + '_.' + out_file.split('.')[-1]
    fasta_file = open(out_file, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            fasta_file.write('>' + matrix[0][i] + '\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(matrix[1][i][(j - matrix[2]):j] + '\n')
                j += matrix[2]
            fasta_file.write(matrix[1][i][(j - matrix[2]):j] + '\n')
    else:
        for i in range(len(matrix[0])):
            fasta_file.write('>' + matrix[0][i] + '\n')
            fasta_file.write(matrix[1][i] + '\n')
    fasta_file.close()


def read_fasta_as_list(fasta_dir):
    fasta_file = open(fasta_dir, 'r')
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
            out_dir = '.'.join(out_dir.split('.')[:-1]) + '_.' + out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            fasta_file.write('>' + matrix[0][i] + '\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                fasta_file.write(''.join(matrix[1][i][(j - matrix[2]):j]) + '\n')
                j += matrix[2]
            fasta_file.write(''.join(matrix[1][i][(j - matrix[2]):j]) + '\n')
    else:
        for i in range(len(matrix[0])):
            fasta_file.write('>' + matrix[0][i] + '\n')
            fasta_file.write(''.join(matrix[1][i]) + '\n')
    fasta_file.close()


# from https://github.com/Kinggerm/PersonalUtilities
# Hashing methods.
def find_exact_repeats(sequence_string, min_repeat_length, circular,
                       accepted_char=set(list("ATGCRMYKHBDVatgcrmykhbdv"))):
    word_size = min(20, min_repeat_length)
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
    # pointer: {active_connection: [(repeat_id_1, sub_repeat_id_1), (repeat_id_2, sub_repeat_id_2)]}
    pointer = {}
    last_connections = set()
    len_indices = len(repeat_indices)

    if circular:
        if len_indices != raw_seq_length and len(repeat_indices):
            # shift start id to the break point
            while (repeat_indices[0] - repeat_indices[-1]) % raw_seq_length == 1:
                repeat_indices.insert(0, repeat_indices.pop(-1))
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            coming_connections = words_to_index[this_word]
            """part 1: dealing with old connection"""
            # Loop 1: find repeats_to_stop
            # Loop 2: delete the pointers pointing to the stopped repeats
            # Loop 3: update the pointers
            # Loop 4: add new pointers if shorter repeats should be continued with less alias
            # Loop 5: update the repeats according to pointers
            repeats_to_stop = {}
            # Loop 1
            for last_con in last_connections:
                here_id, direction_trans = last_con
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_length, direction_trans)
                if candidate_new_connect not in coming_connections:
                    for repeat_id, sub_repeat_id in pointer[last_con]:
                        if repeat_id not in repeats_to_stop:
                            repeats_to_stop[repeat_id] = set()
                        repeats_to_stop[repeat_id].add(sub_repeat_id)

            # Loop 2
            for last_con in list(pointer):
                count_ids = 0
                while count_ids < len(pointer[last_con]):
                    if pointer[last_con][count_ids][0] in repeats_to_stop:
                        del pointer[last_con][count_ids]
                    else:
                        count_ids += 1
                if not len(pointer[last_con]):
                    del pointer[last_con]

            # Loop 3
            new_pointer = {}
            for last_con in pointer:
                here_id, direction_trans = last_con
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_length, direction_trans)
                if candidate_new_connect in coming_connections:  # and last_con in pointer:
                    new_pointer[candidate_new_connect] = list(pointer[last_con])
            pointer = new_pointer

            # Loop 4
            for repeat_id in repeats_to_stop:
                if len(repeats[repeat_id]) - len(repeats_to_stop[repeat_id]) >= 2:
                    repeat_to_be_continued = False
                    for sub_repeat_id in range(len(repeats[repeat_id])):
                        if sub_repeat_id not in repeats_to_stop[repeat_id]:
                            start_id, go_to_id, gt_direction = repeats[repeat_id][sub_repeat_id]
                            new_connect = (
                                (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction) % raw_seq_length,
                                gt_direction)
                            if new_connect in coming_connections and new_connect not in pointer:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[repeat_id])):
                            if inside_repeat_num not in repeats_to_stop[repeat_id]:
                                start_id, go_to_id, gt_direction = repeats[repeat_id][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, gt_direction])
                                new_connect = (
                                    (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction) % raw_seq_length,
                                    gt_direction)
                                if new_connect not in pointer:
                                    pointer[new_connect] = []
                                pointer[new_connect].append((len(repeats) - 1, len(repeats[-1]) - 1))

            # Loop 5
            for last_con in pointer:
                for repeat_id, sub_repeat_id in pointer[last_con]:
                    start_id, previous_id, this_direction = repeats[repeat_id][sub_repeat_id]
                    repeats[repeat_id][sub_repeat_id][1] += this_direction
                    repeats[repeat_id][sub_repeat_id][1] %= raw_seq_length

            """part 2: dealing with new connection"""
            for new_con in coming_connections:
                here_id, direction_trans = new_con
                candidate_last_connect = ((here_id - direction_trans) % raw_seq_length, direction_trans)
                if candidate_last_connect not in last_connections:  # new connection sets
                    repeats.append([])
                    for inside_connection in coming_connections:
                        inside_id, inside_direction = inside_connection
                        repeats[-1].append([(inside_id + (word_size - 1) * (inside_direction == -1)) % raw_seq_length,
                                            (inside_id + (word_size - 1) * (inside_direction == 1)) % raw_seq_length,
                                            inside_direction])
                        if inside_connection not in pointer:
                            pointer[inside_connection] = []
                        pointer[inside_connection].append((len(repeats) - 1, len(repeats[-1]) - 1))
                    break

            if i + 1 < len_indices:
                next_index = repeat_indices[i + 1]
            else:
                next_index = None
            if next_index != (this_index + 1) % raw_seq_length:  # if not continuous, no pointer to update repeats
                pointer = {}
                last_connections = set()
            else:
                last_connections = coming_connections
            # the whole seq is a repeat, problematic?
            if repeats and (repeats[0][0][1] - repeats[0][0][0] + repeats[0][0][2]) % raw_seq_length == 0:
                break
    else:
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            coming_connections = words_to_index[this_word]
            """part 1: dealing with old connection"""
            # Loop 1: find repeats_to_stop
            # Loop 2: delete the pointers pointing to the stopped repeats
            # Loop 3: update the pointers
            # Loop 4: add new pointers if shorter repeats should be continued with less alias
            # Loop 5: update the repeats according to pointers
            repeats_to_stop = {}
            # Loop 1
            for last_con in last_connections:
                here_id, direction_trans = last_con
                candidate_new_connect = (here_id + direction_trans, direction_trans)
                if candidate_new_connect not in coming_connections:
                    # if last_con in pointer:
                    for repeat_id, sub_repeat_id in pointer[last_con]:
                        if repeat_id in repeats_to_stop:
                            repeats_to_stop[repeat_id][sub_repeat_id] = last_con
                        else:
                            repeats_to_stop[repeat_id] = {sub_repeat_id: last_con}
            # Loop 2
            for last_con in list(pointer):
                count_ids = 0
                while count_ids < len(pointer[last_con]):
                    if pointer[last_con][count_ids][0] in repeats_to_stop:
                        del pointer[last_con][count_ids]
                    else:
                        count_ids += 1
                if not len(pointer[last_con]):
                    del pointer[last_con]
            # Loop 3
            new_pointer = {}
            for last_con in pointer:
                here_id, direction_trans = last_con
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_length, direction_trans)
                if candidate_new_connect in coming_connections:  # and last_con in pointer:
                    new_pointer[candidate_new_connect] = list(pointer[last_con])
            pointer = new_pointer
            # Loop 4
            for repeat_id in repeats_to_stop:
                if len(repeats[repeat_id]) - len(repeats_to_stop[repeat_id]) >= 2:
                    repeat_to_be_continued = False
                    for sub_repeat_id in range(len(repeats[repeat_id])):
                        if sub_repeat_id not in repeats_to_stop[repeat_id]:
                            start_id, go_to_id, gt_direction = repeats[repeat_id][sub_repeat_id]
                            new_connect = (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction,
                                           gt_direction)
                            if new_connect in coming_connections and new_connect not in pointer:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[repeat_id])):
                            if inside_repeat_num not in repeats_to_stop[repeat_id]:
                                start_id, go_to_id, gt_direction = repeats[repeat_id][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, gt_direction])
                                new_connect = (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction,
                                               gt_direction)
                                if new_connect not in pointer:
                                    pointer[inside_connection] = []
                                pointer[inside_connection].append((len(repeats) - 1, len(repeats[-1]) - 1))
            for last_con in pointer:
                for repeat_id, sub_repeat_id in pointer[last_con]:
                    start_id, previous_id, this_direction = repeats[repeat_id][sub_repeat_id]
                    repeats[repeat_id][sub_repeat_id][1] += this_direction
            """part 2: dealing with new connection"""
            for last_con in coming_connections:
                here_id, direction_trans = last_con
                candidate_last_connect = (here_id - direction_trans, direction_trans)
                if candidate_last_connect not in last_connections:
                    repeats.append([])
                    for inside_connection in coming_connections:
                        inside_id, inside_direction = inside_connection
                        repeats[-1].append([inside_id + (word_size - 1) * (inside_direction == -1),
                                            inside_id + (word_size - 1) * (inside_direction == 1),
                                            inside_direction])
                        if inside_connection not in pointer:
                            pointer[inside_connection] = []
                        pointer[inside_connection].append((len(repeats) - 1, len(repeats[-1]) - 1))
                    break

            if i + 1 < len_indices:
                next_index = repeat_indices[i + 1]
            else:
                next_index = None
            if next_index != this_index + 1:
                pointer = {}
                last_connections = set()
            else:
                last_connections = coming_connections

    """after-treatment"""
    # 1.delete repeated repeats
    final_repeat = []
    repeat_dicts = set()
    debug_dict = {}
    for repeat_group in repeats:
        if tuple(repeat_group[0]) not in repeat_dicts:
            for single_repeat in repeat_group:
                here_start, here_end, here_direction = single_repeat
                repeat_dicts.add((here_start, here_end, here_direction))
                repeat_dicts.add((here_end, here_start, -here_direction))
                if (here_start, here_end, here_direction) in debug_dict:
                    debug_dict[(here_start, here_end, here_direction)].append(repeat_group)
                else:
                    debug_dict[(here_start, here_end, here_direction)] = [repeat_group]
                if (here_end, here_start, -here_direction) in debug_dict:
                    debug_dict[(here_end, here_start, -here_direction)].append(repeat_group)
                else:
                    debug_dict[(here_end, here_start, -here_direction)] = [repeat_group]
            final_repeat.append(repeat_group)

    # 2.delete small repeats
    count_group__ = 0
    while count_group__ < len(final_repeat):
        here_start, here_end, here_direction = final_repeat[count_group__][0]
        if not (here_start == 0 and here_end == raw_seq_length - 1) and \
                ((here_end - here_start + here_direction) * here_direction) % raw_seq_length < min_repeat_length:
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
                    forward_seq_1 = [sequence[end1 + 1: start2],
                                     sequence[start2: end2 + 1],
                                     sequence[end2 + 1:] + sequence[:start1],
                                     sequence[start1: end1 + 1]]
                    reverse_seq_1 = [raw_rev[len_seq - start2: len_seq - end1 - 1],
                                     raw_rev[len_seq - end1 - 1: len_seq - start1],
                                     raw_rev[len_seq - start1:] + raw_rev[:len_seq - end2 - 1],
                                     raw_rev[len_seq - end2 - 1:len_seq - start2]]
                    # elif start2 - end1 - 1 < len(sequence) + start1 - end2 - 1:
                    forward_seq_2 = [sequence[end2 + 1:] + sequence[:start1],
                                     sequence[start1: end1 + 1],
                                     sequence[end1 + 1: start2],
                                     sequence[start2: end2 + 1]]
                    reverse_seq_2 = [raw_rev[len_seq - start1:] + raw_rev[:len_seq - end2 - 1],
                                     raw_rev[len_seq - end2 - 1:len_seq - start2],
                                     raw_rev[len_seq - start2: len_seq - end1 - 1],
                                     raw_rev[len_seq - end1 - 1: len_seq - start1]]
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
        dynamic_span = max(abs(len(this_string) - len(this_reference)) + 1, dynamic_span)
        no_match_penal = this_string[0] != this_reference[0]
        this_matrix = {(0, 0): {"state": no_match_penal}}
        # calculate the first column
        for i in range(1, min(int(math.ceil(dynamic_span)) + 1, len_str)):
            this_matrix[(i, 0)] = {"right_out": no_match_penal + i, "state": no_match_penal + i}
        # calculate the first line
        for j in range(1, min(int(math.ceil(dynamic_span)) + 1, len_ref)):
            this_matrix[(0, j)] = {"right_out": no_match_penal + j, "state": no_match_penal + j}
        # calculate iteratively
        start = 0
        for i in range(1, len_str):
            start = max(1, int(i - dynamic_span))
            end = min(len_ref, int(math.ceil(i + dynamic_span)))
            # start: no right_in
            no_match_penal = this_string[i] != this_reference[start]
            this_matrix[(i, start)] = {"diagonal_out": this_matrix[(i - 1, start - 1)]["state"] + no_match_penal,
                                       "down_out": this_matrix[(i - 1, start)]["state"] + 1}
            this_matrix[(i, start)]["state"] = min(this_matrix[(i, start)].values())
            # middle
            for j in range(start + 1, end - 1):
                no_match_penal = this_string[i] != this_reference[j]
                this_matrix[(i, j)] = {"diagonal_out": this_matrix[(i - 1, j - 1)]["state"] + no_match_penal,
                                       "down_out": this_matrix[(i - 1, j)]["state"] + 1,
                                       "right_out": this_matrix[(i, j - 1)]["state"] + 1}
                this_matrix[(i, j)]["state"] = min(this_matrix[(i, j)].values())
            # end
            no_match_penal = this_string[i] != this_reference[end - 1]
            this_matrix[(i, end - 1)] = {"diagonal_out": this_matrix[(i - 1, end - 2)]["state"] + no_match_penal}
            if (i, end - 2) in this_matrix:
                this_matrix[(i, end - 1)]["right_out"] = this_matrix[(i, end - 2)]["state"] + 1
            this_matrix[(i, end - 1)]["state"] = min(this_matrix[(i, end - 1)].values())
        # print time.time()-time0
        difference = this_matrix[(len_str - 1, len_ref - 1)]["state"]
        proper_end = True
        for j in range(start, len_ref):
            try:
                if this_matrix[(len_str - 1, j)]["state"] < difference:
                    proper_end = False
                    break
            except KeyError:
                pass
        for i in range(max(0, len_str - len_ref + start), len_str):
            try:
                if this_matrix[(i, len_ref - 1)]["state"] < difference:
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


def fq_simple_generator(fq_dir_list, go_to_line=1, split_pattern=None, min_sub_seq=0, max_n_reads=float("inf")):
    if not ((type(fq_dir_list) is list) or (type(fq_dir_list) is tuple)):
        fq_dir_list = [fq_dir_list]
    max_n_lines = 4 * max_n_reads
    if split_pattern and len(split_pattern) > 2:
        for fq_dir in fq_dir_list:
            count = 0
            with open(fq_dir, 'r') as fq_handler:
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
            with open(fq_dir, 'r') as fq_handler:
                for fq_line in fq_handler:
                    if count % 4 == go_to_line:
                        yield fq_line[:-1]
                    if count >= max_n_lines:
                        break
                    count += 1


def chop_seqs(seq_iter, word_size, mesh_size=1, previous_words=None):
    if not previous_words:
        previous_words = set()
    for seed in seq_iter:
        this_seq_len = len(seed)
        if this_seq_len >= word_size:
            cpt_seed = complementary_seq(seed)
            temp_length = this_seq_len - word_size
            for i in range(0, temp_length + 1, mesh_size):
                if seed[i:i + word_size] not in previous_words:
                    previous_words.add(seed[i:i + word_size])
                if cpt_seed[temp_length - i:this_seq_len - i] not in previous_words:
                    previous_words.add(cpt_seed[temp_length - i:this_seq_len - i])
    return previous_words


def chop_seqs_as_empty_dict(seq_iter, word_size, mesh_size=1, previous_words=None, val_len=float("inf")):
    if previous_words:
        for seed in seq_iter:
            this_seq_len = len(seed)
            if this_seq_len >= word_size:
                cpt_seed = complementary_seq(seed)
                temp_length = this_seq_len - word_size
                for i in range(0, temp_length + 1, mesh_size):
                    forward_word = seed[i:i + word_size]
                    if forward_word in previous_words:
                        previous_words[forward_word] = max(val_len, previous_words[forward_word])
                    else:
                        previous_words[forward_word] = val_len
                    reverse_word = cpt_seed[temp_length - i:this_seq_len - i]
                    if reverse_word in previous_words:
                        previous_words[reverse_word] = max(val_len, previous_words[reverse_word])
                    else:
                        previous_words[reverse_word] = val_len
    else:
        previous_words = dict()
        for seed in seq_iter:
            this_seq_len = len(seed)
            if this_seq_len >= word_size:
                cpt_seed = complementary_seq(seed)
                temp_length = this_seq_len - word_size
                for i in range(0, temp_length + 1, mesh_size):
                    if seed[i:i + word_size] not in previous_words:
                        previous_words[seed[i:i + word_size]] = val_len
                    if cpt_seed[temp_length - i:this_seq_len - i] not in previous_words:
                        previous_words[cpt_seed[temp_length - i:this_seq_len - i]] = val_len
    return previous_words


def chop_seq_list(seq_iter, word_size, mesh_size=1, previous_words=None):
    if not previous_words:
        previous_words = set()
    for seed in seq_iter:
        for seq_part in seed:
            this_seq_len = len(seq_part)
            if this_seq_len >= word_size:
                cpt_seed = complementary_seq(seq_part)
                temp_length = this_seq_len - word_size
                for i in range(0, temp_length + 1, mesh_size):
                    if seq_part[i:i + word_size] not in previous_words:
                        previous_words.add(seq_part[i:i + word_size])
                    if cpt_seed[temp_length - i:this_seq_len - i] not in previous_words:
                        previous_words.add(cpt_seed[temp_length - i:this_seq_len - i])
    return previous_words


def chop_seq_list_as_empty_dict(seq_iter, word_size, mesh_size=1, previous_words=None, val_len=float("inf")):
    if previous_words:
        for seed in seq_iter:
            for seq_part in seed:
                this_seq_len = len(seq_part)
                if this_seq_len >= word_size:
                    cpt_seed = complementary_seq(seq_part)
                    temp_length = this_seq_len - word_size
                    for i in range(0, temp_length + 1, mesh_size):
                        forward_word = seq_part[i:i + word_size]
                        if forward_word in previous_words:
                            previous_words[forward_word] = max(val_len, previous_words[forward_word])
                        else:
                            previous_words[forward_word] = val_len
                        reverse_word = cpt_seed[temp_length - i:this_seq_len - i]
                        if reverse_word in previous_words:
                            previous_words[reverse_word] = max(val_len, previous_words[reverse_word])
                        else:
                            previous_words[reverse_word] = val_len
    else:
        previous_words = dict()
        for seed in seq_iter:
            for seq_part in seed:
                this_seq_len = len(seq_part)
                if this_seq_len >= word_size:
                    cpt_seed = complementary_seq(seq_part)
                    temp_length = this_seq_len - word_size
                    for i in range(0, temp_length + 1, mesh_size):
                        if seq_part[i:i + word_size] not in previous_words:
                            previous_words[seq_part[i:i + word_size]] = val_len
                        if cpt_seed[temp_length - i:this_seq_len - i] not in previous_words:
                            previous_words[cpt_seed[temp_length - i:this_seq_len - i]] = val_len
    return previous_words


def counting_words(seq_generator, words_initial_dict, word_size):
    for seq in seq_generator:
        for i in range(0, len(seq) - word_size + 1):
            this_word = seq[i: i + word_size]
            if this_word in words_initial_dict:
                words_initial_dict[this_word] += 1
    return words_initial_dict


def check_fasta_seq_names(original_fas, new_seed_file, log_handler=None):
    fas_matrix = read_fasta(original_fas)
    short_names = [s_n.split(" ")[0] for s_n in fas_matrix[0]]
    if len(short_names) == len(set(short_names)):
        if os.path.realpath(original_fas) != os.path.realpath(new_seed_file):
            os.system("cp " + original_fas + " " + new_seed_file)
            return True
        else:
            return False
    else:
        if os.path.exists(new_seed_file):
            check_fasta_seq_names(new_seed_file, new_seed_file, log_handler)
            return True
        else:
            if log_handler:
                log_handler.info("Setting '-s " + new_seed_file + "'")
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
            write_fasta(new_seed_file, fas_matrix, True)
            return True


def get_orf_lengths(sequence_string, threshold=200, which_frame=None,
                    here_stop_codons=None, here_start_codons=None):
    """
    :param sequence_string:
    :param threshold: default: 200
    :param which_frame: 1, 2, 3, or None
    :param here_stop_codons: default: CLASSIC_STOP_CODONS
    :param here_start_codons: default: CLASSIC_START_CODONS
    :return: [len_orf1, len_orf2, len_orf3 ...] # longest accumulated orfs among all frame choices
    """
    assert which_frame in {0, 1, 2, None}
    if which_frame is None:
        test_frames = [0, 1, 2]
    else:
        test_frames = [which_frame]
    if here_start_codons is None:
        here_start_codons = CLASSIC_START_CODONS
    if here_stop_codons is None:
        here_stop_codons = CLASSIC_STOP_CODONS
    orf_lengths = {}
    for try_frame in test_frames:
        orf_lengths[try_frame] = []
        this_start = False
        for go in range(try_frame, len(sequence_string), 3):
            if this_start:
                if sequence_string[go:go + 3] not in here_stop_codons:
                    orf_lengths[try_frame][-1] += 3
                else:
                    if orf_lengths[try_frame][-1] < threshold:
                        del orf_lengths[try_frame][-1]
                    this_start = False
            else:
                if sequence_string[go:go + 3] in here_start_codons:
                    orf_lengths[try_frame].append(3)
                    this_start = True
                else:
                    pass
    return sorted(orf_lengths.values(), key=lambda x: -sum(x))[0]


def simulate_fq_simple(
        from_fasta_file, out_dir, out_name=None, is_circular=False, sim_read_len=100, sim_read_jump_size=None,
        generate_paired=False, paired_insert_size=300, generate_spot_num=None, generate_depth=None,
        resume=True):
    """
    :param from_fasta_file:
    :param out_dir:
    :param out_name:
    :param is_circular:
    :param sim_read_len:
    :param sim_read_jump_size: int; mutually exclusive with generate_spot_num, generate_depth; randomly off
    :param generate_paired:
    :param paired_insert_size:
    :param randomly:
    :param generate_spot_num: int; mutually exclusive with sim_read_jump_size, generate_depth; randomly on
    :param generate_depth: int; mutually exclusive with sim_read_jump_size, generate_spot_num; randomly on
    :param resume: continue
    :return:
    """
    if bool(sim_read_jump_size) + bool(generate_spot_num) + bool(generate_depth) == 0:
        raise Exception("One of sim_read_jump_size, generate_spot_num, generate_depth must be given!")
    elif bool(sim_read_jump_size) + bool(generate_spot_num) + bool(generate_depth) > 1:
        raise Exception("Parameters sim_read_jump_size, generate_spot_num, generate_depth are mutually exclusive!")
    if out_name:
        if generate_paired:
            to_fq_files = [os.path.join(out_dir, out_name + "_1.fq"), os.path.join(out_dir, out_name + "_2.fq")]
        else:
            to_fq_files = [os.path.join(out_dir, out_name)]
    else:
        if generate_paired:
            to_fq_files = [os.path.join(out_dir, os.path.basename(from_fasta_file[:-6] + "_1.fq")),
                           os.path.join(out_dir, os.path.basename(from_fasta_file[:-6] + "_2.fq"))]
        else:
            to_fq_files = [os.path.join(out_dir, os.path.basename(from_fasta_file[:-5] + "fq"))]
    if not (resume and os.path.exists(to_fq_files[-1])):
        count_read = 1
        if generate_paired:
            if sim_read_jump_size:
                with open(to_fq_files[0] + ".Temp", "w") as output_handler_1:
                    with open(to_fq_files[1] + ".Temp", "w") as output_handler_2:
                        for from_record in SequenceList(from_fasta_file):
                            from_sequence = from_record.seq
                            if is_circular:
                                from_sequence += from_sequence[:paired_insert_size - 1]
                            from_complement = complementary_seq(from_sequence)
                            len_seq = len(from_sequence)
                            for go_base in range(0, len(from_sequence) - paired_insert_size + 1, sim_read_jump_size):
                                output_handler_1.write("".join(["@", str(count_read), "\n",
                                                                from_sequence[go_base: go_base + sim_read_len], "\n",
                                                                "+", str(count_read), "\n",
                                                                "G" * sim_read_len, "\n"]))
                                go_base = len_seq - (go_base + paired_insert_size)
                                output_handler_2.write("".join(["@", str(count_read), "\n",
                                                                from_complement[go_base: go_base + sim_read_len], "\n",
                                                                "+", str(count_read), "\n",
                                                                "G" * sim_read_len, "\n"]))
                                count_read += 1
            elif generate_spot_num or generate_depth:
                records = SequenceList(from_fasta_file)
                total_len = sum([len(this_r.seq) for this_r in records])
                if generate_depth:
                    generate_spot_num = math.ceil(generate_depth * total_len / (sim_read_len * 2))
                start_ids = []
                cat_all_seqs = []
                accumulated_len = 0
                for from_record in records:
                    from_sequence = from_record.seq
                    if is_circular:
                        from_sequence += from_sequence[:paired_insert_size - 1]
                    start_ids.extend(
                        list(range(accumulated_len, accumulated_len + len(from_sequence) - paired_insert_size + 1)))
                    accumulated_len += len(from_sequence)
                    cat_all_seqs.append(from_sequence)
                cat_all_seqs = "".join(cat_all_seqs)
                cat_all_seqs_rev = complementary_seq(cat_all_seqs)
                chosen_start_ids = [random.choice(start_ids) for foo in range(generate_spot_num)]
                with open(to_fq_files[0] + ".Temp", "w") as output_handler_1:
                    with open(to_fq_files[1] + ".Temp", "w") as output_handler_2:
                        for go_base in chosen_start_ids:
                            output_handler_1.write("".join(["@", str(count_read), "\n",
                                                            cat_all_seqs[go_base: go_base + sim_read_len], "\n",
                                                            "+", str(count_read), "\n",
                                                            "G" * sim_read_len, "\n"]))
                            go_base = accumulated_len - (go_base + paired_insert_size)
                            output_handler_2.write("".join(["@", str(count_read), "\n",
                                                            cat_all_seqs_rev[go_base: go_base + sim_read_len], "\n",
                                                            "+", str(count_read), "\n",
                                                            "G" * sim_read_len, "\n"]))
                            count_read += 1

        else:
            if sim_read_jump_size:
                with open(to_fq_files[0] + ".Temp", "w") as output_handler:
                    for from_record in SequenceList(from_fasta_file):
                        from_sequence = from_record.seq
                        if is_circular:
                            from_sequence += from_sequence[:sim_read_len - 1]
                        for go_base in range(0, len(from_sequence) - sim_read_len + 1, sim_read_jump_size):
                            output_handler.write("".join(["@", str(count_read), "\n",
                                                          from_sequence[go_base: go_base + sim_read_len], "\n",
                                                          "+", str(count_read), "\n",
                                                          "G" * sim_read_len, "\n"]))
                            count_read += 1
            elif generate_spot_num or generate_depth:
                records = SequenceList(from_fasta_file)
                total_len = sum([len(this_r.seq) for this_r in records])
                if generate_depth:
                    generate_spot_num = math.ceil(generate_depth * total_len / (sim_read_len * 2))
                start_ids = []
                cat_all_seqs = []
                accumulated_len = 0
                for from_record in records:
                    from_sequence = from_record.seq
                    if is_circular:
                        from_sequence += from_sequence[:sim_read_len - 1]
                    start_ids.extend(
                        list(range(accumulated_len, accumulated_len + len(from_sequence) - sim_read_len + 1)))
                    accumulated_len += len(from_sequence)
                    cat_all_seqs.append(from_sequence)
                cat_all_seqs = "".join(cat_all_seqs)
                chosen_start_ids = [random.choice(start_ids) for foo in range(generate_spot_num)]
                with open(to_fq_files[0] + ".Temp", "w") as output_handler:
                    for go_base in chosen_start_ids:
                        output_handler.write("".join(["@", str(count_read), "\n",
                                                      cat_all_seqs[go_base: go_base + sim_read_len], "\n",
                                                      "+", str(count_read), "\n",
                                                      "G" * sim_read_len, "\n"]))
                        count_read += 1
        for to_fq_f in to_fq_files:
            os.rename(to_fq_f + ".Temp", to_fq_f)


def get_read_len_mean_max_count(fq_or_fq_files, maximum_n_reads, sampling_percent=1.):
    if type(fq_or_fq_files) is str:
        fq_or_fq_files = [fq_or_fq_files]
    read_lengths = []
    all_counts = []
    if sampling_percent == 1:
        for fq_f in fq_or_fq_files:
            count_r = 0
            for seq in fq_simple_generator(fq_f):
                count_r += 1
                read_lengths.append(len(seq.strip("N")))
                if count_r >= maximum_n_reads:
                    break
            all_counts.append(count_r)
    else:
        sampling_percent = int(1 / sampling_percent)
        for fq_f in fq_or_fq_files:
            count_r = 0
            for seq in fq_simple_generator(fq_f):
                count_r += 1
                if count_r % sampling_percent == 0:
                    read_lengths.append(len(seq.strip("N")))
                if count_r >= maximum_n_reads:
                    break
            all_counts.append(count_r)
    return sum(read_lengths) / len(read_lengths), max(read_lengths), all_counts


def get_read_quality_info(fq_files, maximum_n_reads, min_quality_score, log_handler,
                          maximum_ignore_percent=0.05, sampling_percent=1.):
    if sampling_percent < 1.:
        sampling_percent = int(1 / sampling_percent)
        all_quality_chars_list = []
        for fq_f in fq_files:
            count_r = 0
            this_fq_generator = fq_simple_generator(fq_f, go_to_line=3)
            for quality_str in this_fq_generator:
                if count_r % sampling_percent == 0:
                    all_quality_chars_list.append(quality_str)
                count_r += 1
                if count_r >= maximum_n_reads:
                    break
    else:
        #  sampling_top_n_reads
        all_quality_chars_list = []
        for fq_f in fq_files:
            count_r = 0
            this_fq_generator = fq_simple_generator(fq_f, go_to_line=3)
            for quality_str in this_fq_generator:
                all_quality_chars_list.append(quality_str)
                count_r += 1
                if count_r >= maximum_n_reads:
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
    log_handler.info("Identified quality encoding format = " + the_form)
    if max_quality > the_c_max:
        log_handler.warning("Max quality score " + repr(chr(max_quality)) +
                            "(" + str(max_quality) + ":" + str(max_quality - the_c_min + the_s_min) +
                            ") in your fastq file exceeds the usual range " + str((the_c_min, the_c_max)))
    if min_quality < the_c_min:
        log_handler.warning("Min quality score " + repr(chr(min_quality)) +
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
        log_handler.info("Resetting '--min-quality-score " + str(low_quality_score + the_s_min - the_c_min) + "'")
    if trimmed_quality_chars:
        log_handler.info("Trimming bases with qualities (" + "%.2f" % (ignore_percent * 100) + "%): " +
                         str(ord(min("".join(trimmed_quality_chars)))) + ".." + str(
            ord(max("".join(trimmed_quality_chars)))) +
                         "  " + "".join(trimmed_quality_chars))
    trimmed_quality_chars = "".join(trimmed_quality_chars)

    # calculate mean error rate
    all_quality_char_dict = {in_quality_char: all_quality_chars.count(in_quality_char)
                             for in_quality_char in set(all_quality_chars)}
    error_prob_func = chose_error_prob_func[the_form]
    mean_error_rate = sum([error_prob_func(in_quality_char) * all_quality_char_dict[in_quality_char]
                           for in_quality_char in all_quality_char_dict]) / len_quality_chars_total
    log_handler.info("Mean error rate = " + str(round(mean_error_rate, 4)))

    return "[" + trimmed_quality_chars + "]", mean_error_rate  # , post_trimming_mean


def get_paired_and_unpaired_reads(input_fq_1, input_fq_2, output_p_1, output_p_2, output_u_1, output_u_2):
    file_h_1_lines = open(input_fq_1, 'r').readlines()
    file_h_2_handler = open(input_fq_2, 'r')
    names = {}
    # common_parts = [file_h_1_lines[0].split()[0].split('#')[0].split(".")]
    # len_parts = len(common_parts)
    split_by_dot = False
    split_by_slash = False
    simple_digit_head = False
    if file_h_1_lines and file_h_1_lines[0]:
        if file_h_1_lines[0][-3] == "/" and file_h_1_lines[0][-2].isdigit():
            split_by_slash = True
            for i in range(0, len(file_h_1_lines), 4):
                if file_h_1_lines[i].startswith("@"):
                    this_n = file_h_1_lines[i].split("/")[0]
                    names[this_n] = i
        elif file_h_1_lines[0][1:-1].isdigit():
            simple_digit_head = True
            for i in range(0, len(file_h_1_lines), 4):
                if file_h_1_lines[i].startswith("@"):
                    this_n = file_h_1_lines[i][1:-1]
                    names[this_n] = i
        else:
            first_n = file_h_1_lines[0].split()[0].split('#')[0].split(".")
            split_by_dot = len(first_n) > 2
            if split_by_dot:
                for i in range(0, len(file_h_1_lines), 4):
                    if file_h_1_lines[i].startswith("@"):
                        this_n = ".".join(file_h_1_lines[i].split()[0].split('#')[0].split(".")[:2])
                        names[this_n] = i
            else:
                for i in range(0, len(file_h_1_lines), 4):
                    if file_h_1_lines[i].startswith("@"):
                        this_n = file_h_1_lines[i].split()[0].split('#')[0]
                        names[this_n] = i
    out_paired_h_1 = open(output_p_1 + '.temp', 'w')
    out_paired_h_2 = open(output_p_2 + '.temp', 'w')
    out_unpaired_h_1 = open(output_u_1 + '.temp', 'w')
    out_unpaired_h_2 = open(output_u_2 + '.temp', 'w')
    this_line = file_h_2_handler.readline()
    while this_line:
        if this_line.startswith("@"):
            if split_by_slash:
                this_name = this_line.split("/")[0]
            elif split_by_dot:
                this_name = ".".join(this_line.split()[0].split('#')[0].split(".")[:2])
            elif simple_digit_head:
                this_name = this_line[1:-1]
            else:
                this_name = this_line.split()[0].split('#')[0]
            if this_name in names:
                here_id = names[this_name]
                out_paired_h_1.writelines(file_h_1_lines[here_id:here_id + 4])
                out_paired_h_2.write(this_line)
                for k in range(3):
                    out_paired_h_2.write(file_h_2_handler.readline())
                this_line = file_h_2_handler.readline()
                del names[this_name]
            else:
                out_unpaired_h_2.write(this_line)
                for k in range(3):
                    out_unpaired_h_2.write(file_h_2_handler.readline())
                this_line = file_h_2_handler.readline()
        else:
            this_line = file_h_2_handler.readline()
    left_ids = set(names.values())
    for i in range(0, len(file_h_1_lines), 4):
        if i in left_ids:
            out_unpaired_h_1.writelines(file_h_1_lines[i:i + 4])
    out_paired_h_1.close()
    out_paired_h_2.close()
    out_unpaired_h_1.close()
    out_unpaired_h_2.close()
    os.rename(output_p_1 + '.temp', output_p_1)
    os.rename(output_p_2 + '.temp', output_p_2)
    os.rename(output_u_1 + '.temp', output_u_1)
    os.rename(output_u_2 + '.temp', output_u_2)


# Q = -10*log10(P)
def score_2_error_prob_sanger_quality(quality_score):
    return 10 ** (-quality_score / 10.)


# Q = -10*log10(P/(1-P))
def score_2_error_prob_solexa_quality(quality_score):
    this_sanger_prob = score_2_error_prob_sanger_quality(quality_score)
    return this_sanger_prob / (1 + this_sanger_prob)


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
