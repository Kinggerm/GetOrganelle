import os
import sys
from itertools import combinations, product
from sympy import Symbol, solve, lambdify
from scipy import optimize

path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *
from Library.statistical_func import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]

major_version, minor_version = sys.version_info[:2]
if major_version == 2 and minor_version >= 7:
    python_version = "2.7+"
    RecursionError = RuntimeError
elif major_version == 3 and minor_version >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)


class Assembly:
    def __init__(self, graph_file, min_cov=0., max_cov=inf):
        self.vertex_info = {}
        self.__kmer = 127
        if graph_file.endswith(".gfa"):
            self.parse_gfa(graph_file, min_cov=min_cov, max_cov=max_cov)
        else:
            self.parse_fastg(graph_file, min_cov=min_cov, max_cov=max_cov)
        self.vertex_clusters = []
        self.update_vertex_clusters()
        self.tagged_vertices = {}
        self.vertex_to_copy = {}
        self.vertex_to_float_copy = {}
        self.copy_to_vertex = {}
        self.__inverted_repeat_vertex = {}

    def __repr__(self):
        res = []
        for v in sorted(self.vertex_info):
            res.append(">" + v + "__" + str(self.vertex_info[v]["len"]) + "__" + str(self.vertex_info[v]["cov"]))
            for e in (False, True):
                if len(self.vertex_info[v]["connections"][e]):
                    res.append("(" + ["head", "tail"][e] + ":")
                    res.append(",".join([next_v + "_" + ["head", "tail"][next_e]
                                         for next_v, next_e in self.vertex_info[v]["connections"][e]]))
                    res.append(")")
            res.append("\n")
        return "".join(res)

    def parse_gfa(self, gfa_file, min_cov=0., max_cov=inf):
        with open(gfa_file) as gfa_open:
            kmer_values = set()
            for line in gfa_open:
                if line.startswith("S\t"):
                    flag, vertex_name, sequence, seq_len, seq_num = line.strip().split("\t")
                    seq_len = int(seq_len.split(":")[-1])
                    seq_num = int(seq_num.split(":")[-1])
                    seq_cov = seq_num/float(seq_len)
                    if min_cov <= seq_cov <= max_cov:
                        self.vertex_info[vertex_name] = {"len": seq_len, "cov": seq_cov,
                                                         "connections": {True: set(), False: set()},
                                                         "seq": {True: sequence, False: complementary_seq(sequence)}}
                        if vertex_name.isdigit():
                            self.vertex_info[vertex_name]["long"] = \
                                "EDGE_" + vertex_name + "_length_" + str(seq_len) + "_cov_" + str(round(seq_cov, 5))
                elif line.startswith("L\t"):
                    flag, vertex_1, end_1, vertex_2, end_2, kmer_val = line.strip().split("\t")
                    # "head"~False, "tail"~True
                    end_1 = {"+": True, "-": False}[end_1]
                    end_2 = {"+": False, "-": True}[end_2]
                    kmer_values.add(kmer_val)
                    self.vertex_info[vertex_1]["connections"][end_1].add((vertex_2, end_2))
                    self.vertex_info[vertex_2]["connections"][end_2].add((vertex_1, end_1))
            if len(kmer_values) != 1:
                raise Exception("Multiple overlap values: " + ",".join(sorted(kmer_values)))
            else:
                self.__kmer = int(kmer_values.pop()[:-1])

    def parse_fastg(self, fastg_file, min_cov=0., max_cov=inf):
        fastg_matrix = SequenceList(fastg_file)
        # initialize names; only accept vertex that are formally stored, skip those that are only mentioned after ":"
        for i, seq in enumerate(fastg_matrix):
            if ":" in seq.label:
                this_vertex_str, next_vertices_str = seq.label.strip(";").split(":")
            else:
                this_vertex_str, next_vertices_str = seq.label.strip(";"), ""
            v_tag, vertex_name, l_tag, vertex_len, c_tag, vertex_cov = this_vertex_str.strip("'").split("_")
            # skip vertices with cov out of bounds
            vertex_cov = float(vertex_cov)
            if not (min_cov <= vertex_cov <= max_cov):
                continue
            if vertex_name not in self.vertex_info:
                self.vertex_info[vertex_name] = {"len": int(vertex_len),
                                                 "cov": vertex_cov,
                                                 "long": this_vertex_str.strip("'")}
            if "connections" not in self.vertex_info[vertex_name]:
                self.vertex_info[vertex_name]["connections"] = {True: set(), False: set()}  # "head"~False, "tail"~True

        # adding other info based on existed names
        for i, seq in enumerate(fastg_matrix):
            if ":" in seq.label:
                this_vertex_str, next_vertices_str = seq.label.strip(";").split(":")
            else:
                this_vertex_str, next_vertices_str = seq.label.strip(";"), ""
            v_tag, vertex_name, l_tag, vertex_len, c_tag, vertex_cov = this_vertex_str.strip("'").split("_")
            # skip vertices that not in self.vertex_info: 1. with cov out of bounds
            if vertex_name in self.vertex_info:
                # connections
                this_end = not this_vertex_str.endswith("'")
                if next_vertices_str:
                    for next_vertex_str in next_vertices_str.split(","):
                        next_name = next_vertex_str.strip("'").split("_")[1]
                        if next_name in self.vertex_info:
                            next_end = next_vertex_str.endswith("'")
                            # Adding connection information (edge) to both of the related vertices
                            # even it is only mentioned once in some SPAdes output files
                            self.vertex_info[vertex_name]["connections"][this_end].add((next_name, next_end))
                            self.vertex_info[next_name]["connections"][next_end].add((vertex_name, this_end))
                # sequence
                if "seq" not in self.vertex_info[vertex_name]:
                    self.vertex_info[vertex_name]["seq"] = {}
                    if this_end:
                        self.vertex_info[vertex_name]["seq"][True] = seq.seq
                        self.vertex_info[vertex_name]["seq"][False] = complementary_seq(seq.seq)
                    else:
                        self.vertex_info[vertex_name]["seq"][True] = complementary_seq(seq.seq)
                        self.vertex_info[vertex_name]["seq"][False] = seq.seq

        """detect kmer"""
        ## find initial kmer candidate values
        initial_kmer = set()
        for vertex_name in self.vertex_info:
            if sum([len(self.vertex_info[vertex_name]["connections"][this_e]) for this_e in (True, False)]) != 0:
                for this_e in (True, False):
                    for next_name, next_end in self.vertex_info[vertex_name]["connections"][this_e]:
                        for test_k in range(21, 128, 2):
                            this_seq = self.vertex_info[vertex_name]["seq"][this_e][-test_k:]
                            next_seq = self.vertex_info[next_name]["seq"][not next_end][:test_k]
                            if this_seq == next_seq:
                                initial_kmer.add(test_k)
                        break
                    if initial_kmer:
                        break
            if initial_kmer:
                break
        ## check all edges
        testing_vertices = set(self.vertex_info)
        while initial_kmer and testing_vertices:
            vertex_name = testing_vertices.pop()
            for this_end in (True, False):
                for next_name, next_end in self.vertex_info[vertex_name]["connections"][this_end]:
                    for test_k in list(initial_kmer):
                        this_seq = self.vertex_info[vertex_name]["seq"][this_end][-test_k:]
                        next_seq = self.vertex_info[next_name]["seq"][not next_end][:test_k]
                        if this_seq != next_seq:
                            initial_kmer.discard(test_k)
        if len(initial_kmer) >= 1:
            self.__kmer = max(initial_kmer)
        else:
            raise Exception("No kmer detected!")

    def write_to_fastg(self, out_file, check_postfix=True):
        if check_postfix and not out_file.endswith(".fastg"):
            out_file += ".fastg"
        out_matrix = SequenceList()
        for vertex_name in self.vertex_info:
            try:
                this_name = self.vertex_info[vertex_name]["long"]
                for this_end in (False, True):
                    seq_name = [this_name, ("", "'")[not this_end]]
                    if self.vertex_info[vertex_name]["connections"][this_end]:
                        seq_name.append(":")
                        connect_str = ",".join([self.vertex_info[n_v]["long"] + ("", "'")[n_e]
                                                for n_v, n_e in self.vertex_info[vertex_name]["connections"][this_end]])
                        seq_name.append(connect_str)
                    seq_name.append(";")
                    out_matrix.append(Sequence("".join(seq_name), self.vertex_info[vertex_name]["seq"][this_end]))
            except KeyError as e:
                if str(e) == "'long'":
                    raise Exception("Merged graph cannot be written as fastg format file, please try gfa format!")
        out_matrix.interleaved = 70
        out_matrix.write_fasta(out_file)

    def write_to_fasta(self, out_file, check_postfix=True):
        if check_postfix and not out_file.endswith(".fasta"):
            out_file += ".fasta"
        out_matrix = SequenceList()
        for vertex_name in self.vertex_info:
            out_matrix.append(Sequence(vertex_name, self.vertex_info[vertex_name]["seq"][True]))
        out_matrix.interleaved = 70
        out_matrix.write_fasta(out_file)

    def write_to_gfa(self, out_file, check_postfix=True):
        if check_postfix and not out_file.endswith(".gfa"):
            out_file += ".gfa"
        out_file_handler = open(out_file, "w")
        for vertex_name in self.vertex_info:
            out_file_handler.write("\t".join([
                "S", vertex_name, self.vertex_info[vertex_name]["seq"][True],
                "LN:i:"+str(self.vertex_info[vertex_name]["len"]),
                "RC:i:"+str(int(self.vertex_info[vertex_name]["len"]*self.vertex_info[vertex_name]["cov"]))
                ]) + "\n")
        recorded_connections = set()
        for vertex_name in self.vertex_info:
            for this_end in (False, True):
                for next_v, next_e in self.vertex_info[vertex_name]["connections"][this_end]:
                    this_con = tuple(sorted([(vertex_name, this_end), (next_v, next_e)]))
                    if this_con not in recorded_connections:
                        recorded_connections.add(this_con)
                        out_file_handler.write("\t".join([
                            "L", vertex_name, ("-", "+")[this_end], next_v, ("-", "+")[not next_e],
                            str(self.__kmer) + "M"
                        ]) + "\n")

    def write_out_tags(self, modes, out_file):
        tagged_vertices = set()
        for mode in modes:
            tagged_vertices |= self.tagged_vertices[mode]
        tagged_vertices = sorted(tagged_vertices)
        lines = [["EDGE", "database", "database_weight", "loci"]]
        for this_vertex in tagged_vertices:
            if "tags" in self.vertex_info[this_vertex]:
                all_tags = self.vertex_info[this_vertex]["tags"]
                all_tag_list = sorted(all_tags)
                all_weights = self.vertex_info[this_vertex]["weight"]
                lines.append([this_vertex,
                              ";".join(all_tag_list),
                              ";".join([tag_n + "(" + str(all_weights[tag_n]) + ")" for tag_n in all_tag_list]),
                              ";".join([",".join(sorted(all_tags[tag_n])) for tag_n in all_tag_list])])
            else:
                here_tags = {tag_n for tag_n in modes if this_vertex in self.tagged_vertices[tag_n]}
                lines.append([this_vertex,
                              ";".join(sorted(here_tags)),
                              "", ""])
        open(out_file, "w").writelines(["\t".join(line) + "\n" for line in lines])

    def kmer(self):
        return int(self.__kmer)

    def update_vertex_clusters(self):
        self.vertex_clusters = []
        vertices = set(self.vertex_info)
        for this_vertex in vertices:
            connecting_those = set()
            for connected_set in self.vertex_info[this_vertex]["connections"].values():
                for next_v, next_d in connected_set:
                    for go_to_set, cluster in enumerate(self.vertex_clusters):
                        if next_v in cluster:
                            connecting_those.add(go_to_set)
            if not connecting_those:
                self.vertex_clusters.append({this_vertex})
            elif len(connecting_those) == 1:
                self.vertex_clusters[connecting_those.pop()].add(this_vertex)
            else:
                sorted_those = sorted(connecting_those, reverse=True)
                self.vertex_clusters[sorted_those[-1]].add(this_vertex)
                for go_to_set in sorted_those[:-1]:
                    for that_vertex in self.vertex_clusters[go_to_set]:
                        self.vertex_clusters[sorted_those[-1]].add(that_vertex)
                    del self.vertex_clusters[go_to_set]

    def remove_vertex(self, vertices, update_cluster=True):
        # sys.stdout.write(">>>>" + str(vertices) + " removed!!!!!\n")
        for vertex_name in vertices:
            # sys.stdout.write(">>>>" + vertex_name + " removed!!!!!")
            for this_end, connected_set in self.vertex_info[vertex_name]["connections"].items():
                for next_v, next_e in connected_set:
                    self.vertex_info[next_v]["connections"][next_e].remove((vertex_name, this_end))
            del self.vertex_info[vertex_name]
            for tag in self.tagged_vertices:
                if vertex_name in self.tagged_vertices[tag]:
                    self.tagged_vertices[tag].remove(vertex_name)
            if vertex_name in self.vertex_to_copy:
                this_copy = self.vertex_to_copy[vertex_name]
                self.copy_to_vertex[this_copy].remove(vertex_name)
                if not self.copy_to_vertex[this_copy]:
                    del self.copy_to_vertex[this_copy]
                del self.vertex_to_copy[vertex_name]
                del self.vertex_to_float_copy[vertex_name]
        if update_cluster:
            self.update_vertex_clusters()
        self.__inverted_repeat_vertex = {}

    def detect_parallel_vertices(self, limited_vertices=None):
        if not limited_vertices:
            limiting = False
            limited_vertices = sorted(self.vertex_info)
        else:
            limiting = True
            limited_vertices = sorted(limited_vertices)
        all_both_ends = {}
        for limited_vertex in limited_vertices:
            this_cons = self.vertex_info[limited_vertex]["connections"]
            connect_1 = this_cons[True]
            connect_2 = this_cons[False]
            if connect_1 and connect_2:
                this_ends_raw = [tuple(sorted(connect_1)), tuple(sorted(connect_2))]
                this_ends = sorted(this_ends_raw)
                direction_remained = this_ends_raw == this_ends
                this_ends = tuple(this_ends)
                if this_ends not in all_both_ends:
                    all_both_ends[this_ends] = set()
                all_both_ends[this_ends].add((limited_vertex, direction_remained))
        if limiting:
            limited_vertex_set = set(limited_vertices)
            for each_vertex in self.vertex_info:
                if each_vertex not in limited_vertex_set:
                    this_cons = self.vertex_info[each_vertex]["connections"]
                    connect_1 = this_cons[True]
                    connect_2 = this_cons[False]
                    if connect_1 and connect_2:
                        this_ends_raw = [tuple(sorted(connect_1)), tuple(sorted(connect_2))]
                        this_ends = sorted(this_ends_raw)
                        direction_remained = this_ends_raw == this_ends
                        this_ends = tuple(this_ends)
                        if this_ends in all_both_ends:
                            all_both_ends[this_ends].add((each_vertex, direction_remained))
        return [vertices for vertices in all_both_ends.values() if len(vertices) > 1]

    def is_sequential_repeat(self, vertex, switch_nearby_vertices=True):
        connection_set_t = self.vertex_info[vertex]["connections"][True]
        connection_set_f = self.vertex_info[vertex]["connections"][False]
        inner_circle = []
        
        def search_leakage(start_v, start_e, terminating_end_set):
            in_pipe_leak = False
            circle_in_between = []
            in_vertex_ends = set()
            in_vertex_ends.add((start_v, start_e))
            in_searching_con = [(start_v, not start_e)]
            while in_searching_con:
                in_search_v, in_search_e = in_searching_con.pop(0)
                if (in_search_v, in_search_e) in terminating_end_set:
                    # start from the same (next_t_v, next_t_e), merging to two different ends of connection_set_f
                    if circle_in_between:
                        in_pipe_leak = True
                        break
                    else:
                        circle_in_between.append(((start_v, start_e), (in_search_v, in_search_e)))
                elif (in_search_v, in_search_e) in connection_set_t:
                    in_pipe_leak = True
                    break
                else:
                    for n_in_search_v, n_in_search_e in self.vertex_info[in_search_v]["connections"][in_search_e]:
                        if (n_in_search_v, n_in_search_e) in in_vertex_ends:
                            pass
                        else:
                            in_vertex_ends.add((n_in_search_v, n_in_search_e))
                            in_searching_con.append((n_in_search_v, not n_in_search_e))
            if not in_pipe_leak:
                return circle_in_between
            else:
                return []
        
        # branching ends
        if len(connection_set_t) == len(connection_set_f) == 2:
            for next_t_v, next_t_e in list(connection_set_t):
                this_inner_circle = search_leakage(next_t_v, next_t_e, connection_set_f)
                if this_inner_circle:
                    # check leakage in reverse direction
                    reverse_v, reverse_e = this_inner_circle[0][1]
                    not_leak = search_leakage(reverse_v, reverse_e, connection_set_t)
                    if not_leak:
                        inner_circle.extend(this_inner_circle)
            if switch_nearby_vertices:
                # switch nearby vertices
                # keep those prone to be located in the trunk road of the inverted repeat
                if len(inner_circle) == 1:
                    inner_circle.append([])
                    for next_v, next_e in list(connection_set_t) + list(connection_set_f):
                        if (next_v, next_e) not in inner_circle[0]:
                            inner_circle[1].append((next_v, next_e))
                    inner_circle[1] = tuple(inner_circle[1])
                if len(inner_circle) > 1:
                    del inner_circle[0]
            return inner_circle
        else:
            return inner_circle

    def merge_all_possible_vertices(self, limited_vertices=None, copy_tags=True):
        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)
        while limited_vertices:
            this_vertex = limited_vertices.pop()
            for this_end in (True, False):
                connected_set = self.vertex_info[this_vertex]["connections"][this_end]
                if len(connected_set) == 1:
                    next_vertex, next_end = list(connected_set)[0]
                    if len(self.vertex_info[next_vertex]["connections"][next_end]) == 1 and this_vertex != next_vertex:
                        # reverse the names
                        if this_end:
                            if next_end:
                                new_vertex = this_vertex + "_" + "_".join(next_vertex.split("_")[::-1])
                            else:
                                new_vertex = this_vertex + "_" + next_vertex
                        else:
                            if next_end:
                                new_vertex = next_vertex + "_" + this_vertex
                            else:
                                new_vertex = "_".join(next_vertex.split("_")[::-1]) + "_" + this_vertex

                        limited_vertices.remove(next_vertex)
                        limited_vertices.append(new_vertex)
                        # initialization
                        self.vertex_info[new_vertex] = deepcopy(self.vertex_info[this_vertex])
                        if "long" in self.vertex_info[new_vertex]:
                            del self.vertex_info[new_vertex]["long"]
                        # connections
                        self.vertex_info[new_vertex]["connections"][this_end] \
                            = deepcopy(self.vertex_info[next_vertex]["connections"][not next_end])
                        if (this_vertex, not this_end) in self.vertex_info[new_vertex]["connections"][this_end]:
                            self.vertex_info[new_vertex]["connections"][this_end].remove((this_vertex, not this_end))
                            self.vertex_info[new_vertex]["connections"][this_end].add((new_vertex, not this_end))
                        for new_end in (True, False):
                            for n_n_v, n_n_e in self.vertex_info[new_vertex]["connections"][new_end]:
                                self.vertex_info[n_n_v]["connections"][n_n_e].add((new_vertex, new_end))
                        # len & cov
                        this_len = self.vertex_info[this_vertex]["len"]
                        next_len = self.vertex_info[next_vertex]["len"]
                        this_cov = self.vertex_info[this_vertex]["cov"]
                        next_cov = self.vertex_info[next_vertex]["cov"]
                        self.vertex_info[new_vertex]["len"] = this_len + next_len - self.__kmer
                        self.vertex_info[new_vertex]["cov"] = \
                            ((this_len - self.__kmer + 1) * this_cov + (next_len - self.__kmer + 1) * next_cov) \
                            / ((this_len - self.__kmer + 1) + (next_len - self.__kmer + 1))
                        self.vertex_info[new_vertex]["seq"][this_end] \
                            += self.vertex_info[next_vertex]["seq"][not next_end][self.__kmer:]
                        self.vertex_info[new_vertex]["seq"][not this_end] \
                            = self.vertex_info[next_vertex]["seq"][next_end][:-self.__kmer] \
                            + self.vertex_info[this_vertex]["seq"][not this_end]
                        # tags
                        if copy_tags:
                            if "tags" in self.vertex_info[next_vertex]:
                                if "tags" not in self.vertex_info[new_vertex]:
                                    self.vertex_info[new_vertex]["tags"] = deepcopy(self.vertex_info[next_vertex]["tags"])
                                else:
                                    for mode in self.vertex_info[next_vertex]["tags"]:
                                        if mode not in self.vertex_info[new_vertex]["tags"]:
                                            self.vertex_info[new_vertex]["tags"][mode] \
                                                = deepcopy(self.vertex_info[next_vertex]["tags"][mode])
                                        else:
                                            self.vertex_info[new_vertex]["tags"][mode] \
                                                |= self.vertex_info[next_vertex]["tags"][mode]
                            if "weight" in self.vertex_info[next_vertex]:
                                if "weight" not in self.vertex_info[new_vertex]:
                                    self.vertex_info[new_vertex]["weight"] \
                                        = deepcopy(self.vertex_info[next_vertex]["weight"])
                                else:
                                    for mode in self.vertex_info[next_vertex]["weight"]:
                                        if mode not in self.vertex_info[new_vertex]["weight"]:
                                            self.vertex_info[new_vertex]["weight"][mode] \
                                                = self.vertex_info[next_vertex]["weight"][mode]
                                        else:
                                            self.vertex_info[new_vertex]["weight"][mode] \
                                                += self.vertex_info[next_vertex]["weight"][mode]
                            for mode in self.tagged_vertices:
                                if this_vertex in self.tagged_vertices[mode]:
                                    self.tagged_vertices[mode].add(new_vertex)
                                    self.tagged_vertices[mode].remove(this_vertex)
                                if next_vertex in self.tagged_vertices[mode]:
                                    self.tagged_vertices[mode].add(new_vertex)
                                    self.tagged_vertices[mode].remove(next_vertex)
                        self.remove_vertex([this_vertex, next_vertex], update_cluster=False)
                        break
        self.update_vertex_clusters()

    def estimate_copy_and_depth_by_cov(self, limited_vertices=None, given_average_cov=None, mode="cp",
                                       log_handler=None, verbose=True, debug=False):
        if mode == "cp":
            max_majority_cov = 2
        elif mode == "mt":
            max_majority_cov = 4
        elif mode == "nr":
            max_majority_cov = 1
        elif mode == "all":
            max_majority_cov = 100
        else:
            max_majority_cov = 100
        if not given_average_cov:
            if not limited_vertices:
                limited_vertices = sorted(self.vertex_info)
            else:
                limited_vertices = sorted(limited_vertices)
            previous_val = {0.}
            new_val = -1.
            while new_val not in previous_val:
                previous_val.add(new_val)
                # estimate baseline depth
                total_product = 0.
                total_len = 0
                for vertex_name in limited_vertices:
                    this_len = (self.vertex_info[vertex_name]["len"] - self.__kmer + 1) \
                               * self.vertex_to_copy.get(vertex_name, 1)
                    this_cov = self.vertex_info[vertex_name]["cov"] / self.vertex_to_copy.get(vertex_name, 1)
                    total_len += this_len
                    total_product += this_len * this_cov
                new_val = total_product/total_len
                # adjust this_copy according to new baseline depth
                for vertex_name in self.vertex_info:
                    if vertex_name in self.vertex_to_copy:
                        old_copy = self.vertex_to_copy[vertex_name]
                        self.copy_to_vertex[old_copy].remove(vertex_name)
                        if not self.copy_to_vertex[old_copy]:
                            del self.copy_to_vertex[old_copy]
                    this_float_copy = self.vertex_info[vertex_name]["cov"] / new_val
                    this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_cov)
                    self.vertex_to_float_copy[vertex_name] = this_float_copy
                    self.vertex_to_copy[vertex_name] = this_copy
                    if this_copy not in self.copy_to_vertex:
                        self.copy_to_vertex[this_copy] = set()
                    self.copy_to_vertex[this_copy].add(vertex_name)
            if debug or verbose:
                if log_handler:
                    log_handler.info("updating average target kmer-coverage: " + str(round(new_val, 2)))
                else:
                    sys.stdout.write("updating average target kmer-coverage: " + str(round(new_val, 2)) + "\n")
            return new_val
        else:
            # adjust this_copy according to user-defined depth
            for vertex_name in self.vertex_info:
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    if not self.copy_to_vertex[old_copy]:
                        del self.copy_to_vertex[old_copy]
                this_float_copy = self.vertex_info[vertex_name]["cov"] / given_average_cov
                this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_cov)
                self.vertex_to_float_copy[vertex_name] = this_float_copy
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)
            return given_average_cov

    def estimate_copy_and_depth_precisely(self, verbose=True, maximum_copy_num=10, return_new_graphs=True,
                                          log_handler=None, debug=False):

        def get_formula(from_vertex, from_end, to_vertex, to_end):
            result_form = vertex_to_symbols[from_vertex]
            for next_v, next_e in self.vertex_info[from_vertex]["connections"][from_end]:
                if (next_v, next_e) != (to_vertex, to_end):  # and next_v in vertices_set:
                    recorded_ends.add((next_v, next_e))
                    result_form -= get_formula(next_v, next_e, from_vertex, from_end)
            return result_form

        # for compatibility between scipy and sympy
        def least_square_function_v(x):
            return least_square_function(*tuple(x))

        """ create constraints by creating inequations: the copy of every contig has to be >= 1 """
        def constraint_min_function(x):
            replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
            expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
            min_copy = np.array([1.001] * len(all_v_symbols) + [2.001] * len(extra_symbol_to_str))
            return expression_array - min_copy

        def constraint_min_function_br(x):
            replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
            expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
            min_copy = np.array([1.0] * len(all_v_symbols) + [2.0] * len(extra_symbol_to_str))
            return expression_array - min_copy

        def constraint_int_function(x):
            replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
            expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
            # diff = np.array([0] * len(all_symbols))
            return sum([abs(every_copy - int(every_copy)) for every_copy in expression_array])

        def minimize_brute_force(func, range_list, constraint_list, round_digit=4, display_p=True,
                                 in_log_handler=log_handler):
            # time0 = time.time()
            best_fun_val = inf
            best_para_val = []
            count_round = 0
            count_valid = 0
            for value_set in product(*[list(this_range) for this_range in range_list]):
                count_round += 1
                is_valid_set = True
                for cons in constraint_list:
                    if cons["type"] == "ineq":
                        try:
                            if (cons["fun"](value_set) < 0).any():
                                is_valid_set = False
                                break
                        except TypeError:
                            is_valid_set = False
                            break
                    elif cons["type"] == "eq":
                        try:
                            if cons["fun"](value_set) != 0:
                                is_valid_set = False
                                break
                        except TypeError:
                            is_valid_set = False
                            break
                if not is_valid_set:
                    continue
                count_valid += 1
                this_fun_val = round(func(value_set), round_digit)
                if this_fun_val < best_fun_val:
                    best_para_val = [value_set]
                    best_fun_val = this_fun_val
                elif this_fun_val == best_fun_val:
                    best_para_val.append(value_set)
                else:
                    pass
            if log_handler:
                if debug or display_p:
                    log_handler.info("Brute rounds: " + str(count_valid) + "/" + str(count_round))
                    log_handler.info("Brute best function value: " + str(best_fun_val))
                if debug:
                    log_handler.info("Best solution: " + str(best_para_val))
            else:
                if debug or display_p:
                    sys.stdout.write("Brute rounds: " + str(count_valid) + "/" + str(count_round) + "\n")
                    sys.stdout.write("Brute best function value: " + str(best_fun_val) + "\n")
                if debug:
                    sys.stdout.write("Best solution: " + str(best_para_val) + "\n")
            return best_para_val

        vertices_list = sorted(self.vertex_info)

        """ create constraints by creating multivariate equations """
        vertex_to_symbols = {vertex_name: Symbol("V"+vertex_name, integer=True)  # positive=True)
                             for vertex_name in vertices_list}
        symbols_to_vertex = {vertex_to_symbols[vertex_name]: vertex_name for vertex_name in vertices_list}
        all_v_symbols = list(symbols_to_vertex)
        formulae = []
        recorded_ends = set()
        for vertex_name in vertices_list:
            for this_end in (True, False):
                if (vertex_name, this_end) not in recorded_ends:
                    recorded_ends.add((vertex_name, this_end))
                    connection_set = self.vertex_info[vertex_name]["connections"][this_end]
                    if connection_set:  # len([n_v for n_v, n_e in connection_set if n_v in vertices_set]):
                        this_formula = vertex_to_symbols[vertex_name]
                        for n_v, n_e in connection_set:
                            # if n_v in vertices_set:
                            recorded_ends.add((n_v, n_e))
                            try:
                                this_formula -= get_formula(n_v, n_e, vertex_name, this_end)
                            except RecursionError:
                                direct = ["_tail", "_head"]
                                if log_handler:
                                    log_handler.error("Formulating for: " +
                                                      n_v + direct[n_e] + vertex_name + direct[this_end] + "!")
                                else:
                                    sys.stdout.write("Formulating for: " +
                                                      n_v + direct[n_e] + vertex_name + direct[this_end] + "!\n")
                                raise RecursionError
                        formulae.append(this_formula)

        # add following extra limitation
        # set cov_bubbles = x*near_by_cov, x is an integer
        extra_str_to_symbol = {}
        extra_symbol_to_str = {}
        for vertex_name in vertices_list:
            inner_circle_pairs = self.is_sequential_repeat(vertex_name)
            for ((from_v, from_e), (to_v, to_e)) in inner_circle_pairs:
                new_str = "E" + str(len(extra_str_to_symbol))
                extra_str_to_symbol[new_str] = Symbol(new_str, integer=True)
                extra_symbol_to_str[extra_str_to_symbol[new_str]] = new_str
                formulae.append(vertex_to_symbols[vertex_name] -
                                vertex_to_symbols[from_v] * extra_str_to_symbol[new_str])
        if debug:
            if log_handler:
                log_handler.info("formulae: " + str(formulae) + "\n")
            else:
                sys.stdout.write("formulae: " + str(formulae) + "\n")

        # solve the equations
        all_symbols = all_v_symbols + list(extra_symbol_to_str)
        copy_solution = solve(formulae, all_v_symbols)
        copy_solution = copy_solution if copy_solution else {}
        if type(copy_solution) == list:  # delete 0 containing set
            go_solution = 0
            while go_solution < len(copy_solution):
                if 0 in set(copy_solution[go_solution].values()):
                    del copy_solution[go_solution]
                else:
                    go_solution += 1
        if not copy_solution:
            raise Exception("Incomplete/Complicated target graph (1)!")
        elif type(copy_solution) == list:
            if len(copy_solution) > 2:
                raise Exception("Incomplete/Complicated target graph (2)!")
            else:
                copy_solution = copy_solution[0]

        free_copy_variables = list()
        for symbol_used in all_symbols:
            if symbol_used not in copy_solution:
                free_copy_variables.append(symbol_used)
                copy_solution[symbol_used] = symbol_used

        """ minimizing equation-based copy values and their deviations from coverage-based copy values """
        least_square_expr = 0
        for symbol_used in all_v_symbols:
            # least_square_expr += copy_solution[symbol_used]
            this_vertex = symbols_to_vertex[symbol_used]
            this_copy = self.vertex_to_float_copy[this_vertex]
            least_square_expr += (copy_solution[symbol_used] - this_copy) ** 2  # * self.vertex_info[this_vertex]["len"]
        least_square_function = lambdify(args=free_copy_variables, expr=least_square_expr)

        # for safe running
        if len(free_copy_variables) > 10:
            raise Exception("Free variable > 10 is not accepted yet!")

        if maximum_copy_num ** len(free_copy_variables) < 5E6:
            # sometimes, SLSQP ignores bounds and constraints
            copy_results = minimize_brute_force(func=least_square_function_v,
                                                range_list=[range(1, maximum_copy_num+1)]*len(free_copy_variables),
                                                constraint_list=({'type': 'ineq', 'fun': constraint_min_function_br},
                                                                 {'type': 'eq', 'fun': constraint_int_function}),
                                                display_p=verbose)
        else:
            constraints = ({'type': 'ineq', 'fun': constraint_min_function},
                           {'type': 'eq', 'fun': constraint_int_function})
            copy_results = set()
            best_fun = inf
            opt = {'disp': verbose, "maxiter": 100}
            for initial_copy in range(maximum_copy_num*2 + 1):
                if initial_copy < maximum_copy_num:
                    initials = np.array([initial_copy + 1] * len(free_copy_variables))
                elif initial_copy < maximum_copy_num*2:
                    initials = np.array([random.randint(1, maximum_copy_num)] * len(free_copy_variables))
                else:
                    initials = np.array([self.vertex_to_copy.get(symbols_to_vertex.get(symb, False), 2)
                                         for symb in free_copy_variables])
                bounds = [(1, maximum_copy_num) for foo in range(len(free_copy_variables))]
                try:
                    copy_result = optimize.minimize(fun=least_square_function_v, x0=initials, jac=False,
                                                    method='SLSQP', bounds=bounds, constraints=constraints, options=opt)
                except Exception:
                    continue
                if copy_result.fun < best_fun:
                    best_fun = round(copy_result.fun, 2)
                    copy_results = {tuple(copy_result.x)}
                elif copy_result.fun == best_fun:
                    copy_results.add(tuple(copy_result.x))
                else:
                    pass
            if debug or verbose:
                if log_handler:
                    log_handler.info("Best function value: " + str(best_fun))
                else:
                    sys.stdout.write("Best function value: " + str(best_fun) + "\n")
        if len(copy_results) == 1:
            copy_results = list(copy_results)
        elif len(copy_results) > 1:
            # sort results
            copy_results = sorted(copy_results, key=lambda
                x: sum([(x[go_sym] - self.vertex_to_float_copy[symbols_to_vertex[symbol_used]]) ** 2
                        for go_sym, symbol_used in enumerate(all_v_symbols)]))
        else:
            raise Exception("Incomplete/Complicated target graph (3)!")

        if return_new_graphs:
            """ produce all possible vertex copy combinations """
            final_results = []
            for go_res, copy_result in enumerate(copy_results):
                final_results.append({"graph": deepcopy(self)})
                free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
                                            for i, this_copy in enumerate(copy_result)}

                """ record new copy values """
                for this_symbol in all_v_symbols:
                    vertex_name = symbols_to_vertex[this_symbol]
                    if vertex_name in final_results[go_res]["graph"].vertex_to_copy:
                        old_copy = final_results[go_res]["graph"].vertex_to_copy[vertex_name]
                        final_results[go_res]["graph"].copy_to_vertex[old_copy].remove(vertex_name)
                        if not final_results[go_res]["graph"].copy_to_vertex[old_copy]:
                            del final_results[go_res]["graph"].copy_to_vertex[old_copy]
                    this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
                    assert this_copy > 0, vertex_name
                    final_results[go_res]["graph"].vertex_to_copy[vertex_name] = this_copy
                    if this_copy not in final_results[go_res]["graph"].copy_to_vertex:
                        final_results[go_res]["graph"].copy_to_vertex[this_copy] = set()
                    final_results[go_res]["graph"].copy_to_vertex[this_copy].add(vertex_name)

                """ re-estimate baseline depth """
                total_product = 0.
                total_len = 0
                for vertex_name in vertices_list:
                    this_len = (self.vertex_info[vertex_name]["len"] - self.__kmer + 1) \
                               * final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
                    this_cov = self.vertex_info[vertex_name]["cov"] \
                               / final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
                    total_len += this_len
                    total_product += this_len * this_cov
                final_results[go_res]["cov"] = total_product / total_len
                if log_handler:
                    for vertex_name in vertices_list:
                        log_handler.info("Vertex_" + vertex_name + " #copy: " +
                                         str(final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)))
                    log_handler.info("Average target kmer-coverage(" + str(go_res + 1) + "): " +
                                     str(round(final_results[go_res]["cov"], 2)))
                else:
                    for vertex_name in vertices_list:
                        sys.stdout.write("Vertex_" + vertex_name + " #copy: " +
                                         str(final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)) + "\n")
                    sys.stdout.write("Average target kmer-coverage(" + str(go_res + 1) + "): " +
                                     str(round(final_results[go_res]["cov"], 2)) + "\n")
            return final_results

        else:
            """ produce the first-ranked copy combination """
            free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
                                        for i, this_copy in enumerate(copy_results[0])}

            """ record new copy values """
            for this_symbol in all_v_symbols:
                vertex_name = symbols_to_vertex[this_symbol]
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    if not self.copy_to_vertex[old_copy]:
                        del self.copy_to_vertex[old_copy]
                this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
                assert this_copy > 0
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)

            if debug or verbose:
                """ re-estimate baseline depth """
                total_product = 0.
                total_len = 0
                for vertex_name in vertices_list:
                    this_len = (self.vertex_info[vertex_name]["len"] - self.__kmer + 1) \
                               * self.vertex_to_copy.get(vertex_name, 1)
                    this_cov = self.vertex_info[vertex_name]["cov"] / self.vertex_to_copy.get(vertex_name, 1)
                    total_len += this_len
                    total_product += this_len * this_cov
                new_val = total_product / total_len
                if log_handler:
                    log_handler.info("Average target kmer-coverage: " + str(round(new_val, 2)))
                else:
                    sys.stdout.write("Average target kmer-coverage: " + str(round(new_val, 2)) + "\n")

    def tag_in_between(self, mode):
        # add those in between the tagged vertices to tagged_vertices, which offered the only connection
        updated = True
        candidate_vertices = list(self.vertex_info)
        while updated:
            updated = False
            go_to_v = 0
            while go_to_v < len(candidate_vertices):
                can_v = candidate_vertices[go_to_v]
                if can_v in self.tagged_vertices[mode]:
                    del candidate_vertices[go_to_v]
                    continue
                else:
                    if sum([bool(c_c) for c_c in self.vertex_info[can_v]["connections"].values()]) != 2:
                        del candidate_vertices[go_to_v]
                        continue
                    count_nearby_tagged = []
                    for can_end, can_connect in self.vertex_info[can_v]["connections"].items():
                        for next_v, next_e in can_connect:
                            if next_v in self.tagged_vertices[mode] and \
                                    len(self.vertex_info[next_v]["connections"][next_e]) == 1:
                                count_nearby_tagged.append((next_v, next_e))
                                break
                    if len(count_nearby_tagged) == 2:
                        del candidate_vertices[go_to_v]
                        # add in between
                        self.tagged_vertices[mode].add(can_v)
                        if "weight" not in self.vertex_info[can_v]:
                            self.vertex_info[can_v]["weight"] = {}
                        if mode not in self.vertex_info[can_v]["weight"]:
                            self.vertex_info[can_v]["weight"][mode] = 0
                        self.vertex_info[can_v]["weight"][mode] += 1 * self.vertex_info[can_v]["cov"]
                        # add extra circle
                        near_by_pairs = self.is_sequential_repeat(can_v, switch_nearby_vertices=False)
                        if near_by_pairs:
                            checking_new = []
                            coverage_folds = []
                            for near_by_p in near_by_pairs:
                                for (near_v, near_e) in near_by_p:
                                    if (near_v, near_e) not in count_nearby_tagged:
                                        checking_new.append(near_v)
                                        coverage_folds.append(
                                            round(self.vertex_info[can_v]["cov"] /
                                                  self.vertex_info[near_v]["cov"], 0))
                            for near_v, near_e in count_nearby_tagged:
                                coverage_folds.append(
                                    round(self.vertex_info[can_v]["cov"] /
                                          self.vertex_info[near_v]["cov"], 0))
                            if max(coverage_folds) >= 2:
                                for extra_v_to_add in set(checking_new):
                                    self.tagged_vertices[mode].add(extra_v_to_add)
                                    try:
                                        candidate_vertices.remove(extra_v_to_add)
                                    except ValueError:
                                        pass
                                    if "weight" not in self.vertex_info[extra_v_to_add]:
                                        self.vertex_info[extra_v_to_add]["weight"] = {mode: 0}
                                    self.vertex_info[extra_v_to_add]["weight"][mode] \
                                        += 1 * self.vertex_info[extra_v_to_add]["cov"]
                        updated = True
                        break
                    else:
                        go_to_v += 1

    def parse_tab_file(self, tab_file, mode, type_factor, log_handler=None):
        # parse_csv, every locus only occur in one vertex (removing locations with smaller weight)
        tag_loci = {}
        tab_matrix = [line.strip("\n").split("\t") for line in open(tab_file)][1:]
        for node_record in tab_matrix:
            vertex_name = node_record[0]
            if vertex_name in self.vertex_info:
                matched = node_record[5].split(">>")
                for locus in matched:
                    if "(" in locus:
                        locus_spl = locus.split("(")
                        locus_type = locus_spl[-1].split(",")[1][:-1]
                        if locus_type not in tag_loci:
                            tag_loci[locus_type] = {}
                        locus_name = "(".join(locus_spl[:-1])
                        locus_start, locus_end = locus_spl[-1].split(",")[0].split("-")
                        locus_start, locus_end = int(locus_start), int(locus_end)
                        locus_len = locus_end - locus_start + 1
                        # skip those tags concerning only the overlapping sites
                        if (locus_start == 1 or locus_end == self.vertex_info[vertex_name]["len"]) \
                                and locus_len == self.__kmer:
                            continue
                        if locus_name in tag_loci[locus_type]:
                            new_weight = locus_len * self.vertex_info[vertex_name]["cov"]
                            if new_weight > tag_loci[locus_type][locus_name]["weight"]:
                                tag_loci[locus_type][locus_name] = {"vertex": vertex_name, "len": locus_len,
                                                                    "weight": new_weight}
                        else:
                            tag_loci[locus_type][locus_name] = {"vertex": vertex_name, "len": locus_len,
                                                                "weight": locus_len * self.vertex_info[vertex_name]["cov"]}

        for locus_type in tag_loci:
            self.tagged_vertices[locus_type] = set()
            for locus_name in tag_loci[locus_type]:
                vertex_name = tag_loci[locus_type][locus_name]["vertex"]
                loci_weight = tag_loci[locus_type][locus_name]["weight"]
                # tags
                if "tags" not in self.vertex_info[vertex_name]:
                    self.vertex_info[vertex_name]["tags"] = {}
                if locus_type in self.vertex_info[vertex_name]["tags"]:
                    self.vertex_info[vertex_name]["tags"][locus_type].add(locus_name)
                else:
                    self.vertex_info[vertex_name]["tags"][locus_type] = {locus_name}
                # weight
                if "weight" not in self.vertex_info[vertex_name]:
                    self.vertex_info[vertex_name]["weight"] = {}
                if locus_type in self.vertex_info[vertex_name]["weight"]:
                    self.vertex_info[vertex_name]["weight"][locus_type] += loci_weight
                else:
                    self.vertex_info[vertex_name]["weight"][locus_type] = loci_weight
                self.tagged_vertices[locus_type].add(vertex_name)

        for vertex_name in self.vertex_info:
            if "weight" in self.vertex_info[vertex_name]:
                if len(self.vertex_info[vertex_name]["weight"]) > 1:
                    all_weights = sorted([(loc_type, self.vertex_info[vertex_name]["weight"][loc_type])
                                          for loc_type in self.vertex_info[vertex_name]["weight"]], key=lambda x: -x[1])
                    best_t, best_w = all_weights[0]
                    for next_t, next_w in all_weights[1:]:
                        if next_w * type_factor < best_w:
                            self.tagged_vertices[next_t].remove(vertex_name)

        if len(self.tagged_vertices[mode]) == 0:
            raise Exception("No available target information found in " + tab_file)

    def filter_by_coverage(self, drop_num=1, mode="cp", log_hard_cov_threshold=10.,
                           weight_factor=100., min_sigma_factor=0.1,
                           verbose=False, log_handler=None, debug=False):
        log_hard_cov_threshold = abs(log(log_hard_cov_threshold))
        vertices = sorted(self.vertex_info)
        v_coverages = {this_v: self.vertex_info[this_v]["cov"] / self.vertex_to_copy.get(this_v, 1)
                       for this_v in vertices}
        max_tagged_cov = max([v_coverages[tagged_v] for tagged_v in self.tagged_vertices[mode]])
        # removing coverage with 10 times lower/greater than tagged_cov
        removing_low_cov = [candidate_v
                            for candidate_v in vertices
                            if abs(log(self.vertex_info[candidate_v]["cov"]/max_tagged_cov)) > log_hard_cov_threshold]
        if removing_low_cov:
            if log_handler and (debug or verbose):
                log_handler.info("removing extremely outlying coverage contigs: " + str(removing_low_cov))
            elif verbose or debug:
                sys.stdout.write("removing extremely outlying coverage contigs: " + str(removing_low_cov) + "\n")
            self.remove_vertex(removing_low_cov)
        self.merge_all_possible_vertices()
        vertices = sorted(self.vertex_info)
        v_coverages = {this_v: self.vertex_info[this_v]["cov"] / self.vertex_to_copy.get(this_v, 1)
                       for this_v in vertices}

        coverages = np.array([v_coverages[this_v] for this_v in vertices])
        cover_weights = np.array([(self.vertex_info[this_v]["len"] - self.__kmer) * self.vertex_to_copy.get(this_v, 1)
                                  for this_v in vertices])
        set_cluster = {v_coverages[tagged_v]: 0 for tagged_v in self.tagged_vertices[mode]}

        # gmm_scheme = gmm_with_em_aic(coverages, maximum_cluster=6, cluster_limited=set_cluster,
        #                              min_sigma_factor=min_sigma_factor)
        if log_handler and (debug or verbose):
            log_handler.info("Vertices: " + str(vertices))
            log_handler.info("Coverages: " + str([float("%.1f" % cov_x) for cov_x in coverages]))
        elif verbose or debug:
            sys.stdout.write("Vertices: " + str(vertices) + "\n")
            sys.stdout.write("Coverages: " + str([float("%.1f" % cov_x) for cov_x in coverages]) + "\n")
        gmm_scheme = weighted_gmm_with_em_aic(coverages, data_weights=cover_weights,
                                              maximum_cluster=6, cluster_limited=set_cluster,
                                              min_sigma_factor=min_sigma_factor)
        cluster_num = gmm_scheme["cluster_num"]
        parameters = gmm_scheme["parameters"]
        labels = gmm_scheme["labels"]
        if log_handler and (debug or verbose):
            log_handler.info("Labels: " + str(labels))
        elif verbose or debug:
            sys.stdout.write("Labels: " + str(labels) + "\n")
        selected_label_type = list(
            set([lb for go, lb in enumerate(labels) if vertices[go] in self.tagged_vertices[mode]]))
        if len(selected_label_type) > 1:
            label_weights = {}
            # for lb in selected_label_type:
            #     this_add_up = 0
            #     for go in np.where(labels == lb)[0]:
            #         this_add_up += self.vertex_info[vertices[go]].get("weight", {}).get(mode, 0)
            #     label_weights[lb] = this_add_up
            label_weights = {lb: sum([self.vertex_info[vertices[go]].get("weight", {}).get(mode, 0)
                                      for go in np.where(labels == lb)[0]])
                             for lb in selected_label_type}
            selected_label_type.sort(key=lambda x: -label_weights[x])
            remained_label_type = {selected_label_type[0]}
            for candidate_lb_type in selected_label_type[1:]:
                if label_weights[candidate_lb_type] * weight_factor >= selected_label_type[0]:
                    remained_label_type.add(candidate_lb_type)
                else:
                    break
            extra_kept = set()
            for candidate_lb_type in selected_label_type:
                if candidate_lb_type not in remained_label_type:
                    can_mu = parameters[candidate_lb_type]["mu"]
                    for remained_l in remained_label_type:
                        if abs(can_mu - parameters[remained_l]["mu"]) < 2 * parameters[remained_l]["sigma"]:
                            extra_kept.add(candidate_lb_type)
                            break
            remained_label_type |= extra_kept
        else:
            remained_label_type = {selected_label_type[0]}
        if debug or verbose:
            if log_handler:
                log_handler.info("\t".join(["Mu" + str(go) + ":" + str(parameters[lab_tp]["mu"]) +
                                            " Sigma" + str(go) + ":" + str(parameters[lab_tp]["sigma"])
                                            for go, lab_tp in enumerate(remained_label_type)]))
            else:
                sys.stdout.write("\t".join(["Mu" + str(go) + ":" + str(parameters[lab_tp]["mu"]) +
                                            " Sigma" + str(go) + ":" + str(parameters[lab_tp]["sigma"])
                                            for go, lab_tp in enumerate(remained_label_type)]) + "\n")
        candidate_dropping_label_type = {l_t: inf for l_t in set(range(cluster_num)) - remained_label_type}
        for lab_tp in candidate_dropping_label_type:
            check_mu = parameters[lab_tp]["mu"]
            check_sigma = parameters[lab_tp]["sigma"]
            for remained_l in remained_label_type:
                rem_mu = parameters[remained_l]["mu"]
                rem_sigma = parameters[remained_l]["sigma"]
                this_dist = abs(rem_mu - check_mu) - 2 * (check_sigma + rem_sigma)
                candidate_dropping_label_type[lab_tp] = min(candidate_dropping_label_type[lab_tp], this_dist)
        dropping_type = sorted(candidate_dropping_label_type, key=lambda x: -candidate_dropping_label_type[x])
        dropping_type = dropping_type[:drop_num]
        if debug or verbose:
            if log_handler:
                for lab_tp in dropping_type:
                    if candidate_dropping_label_type[lab_tp] < 0:
                        log_handler.warning("Distinguishable vertices "
                                            + str([vertices[go] for go in np.where(labels == lab_tp)[0]])
                                            + " removed!")
            else:
                for lab_tp in dropping_type:
                    if candidate_dropping_label_type[lab_tp] < 0:
                        sys.stdout.write("Warning: distinguishable vertices "
                                         + str([vertices[go] for go in np.where(labels == lab_tp)[0]])
                                         + " removed!\n")
        vertices_to_del = {vertices[go] for go, lb in enumerate(labels) if lb in set(dropping_type)}
        if vertices_to_del:
            changed = True
            if verbose or debug:
                if log_handler:
                    log_handler.info("removing outlying coverage contigs: " + str(vertices_to_del))
                else:
                    sys.stdout.write("removing outlying coverage contigs: " + str(vertices_to_del) + "\n")
            self.remove_vertex(vertices_to_del)
        else:
            changed = False
        return changed, [(parameters[lab_tp]["mu"], parameters[lab_tp]["sigma"]) for lab_tp in remained_label_type]

    def generate_consensus_vertex(self, vertices, directions, copy_tags=True, check_parallel_vertices=True,
                                  log_handler=None):
        if check_parallel_vertices:
            connection_type = None
            seq_len = None
            assert len(vertices) == len(set(vertices)) == len(directions)
            for go_v, this_v in enumerate(vertices):
                if seq_len:
                    assert seq_len == len(self.vertex_info[this_v]["seq"][True])
                else:
                    seq_len = len(self.vertex_info[this_v]["seq"][True])
                this_cons = self.vertex_info[this_v]["connections"]
                this_ends = tuple([tuple(sorted(this_cons[[directions[go_v]]])),
                                   tuple(sorted(this_cons[not [directions[go_v]]]))])
                if connection_type:
                    assert connection_type == this_ends
                else:
                    connection_type = this_ends

        if len(vertices) > 1:
            new_vertex = "(" + "|".join(vertices) + ")"
            self.vertex_info[new_vertex] = deepcopy(self.vertex_info[vertices[0]])
            self.vertex_info[new_vertex]["cov"] = sum([self.vertex_info[v]["cov"] for v in vertices])
            if "long" in self.vertex_info[new_vertex]:
                del self.vertex_info[new_vertex]["long"]

            for new_end in (True, False):
                for n_n_v, n_n_e in self.vertex_info[new_vertex]["connections"][new_end]:
                    self.vertex_info[n_n_v]["connections"][n_n_e].add((new_vertex, new_end))

            consensus_s = generate_consensus(*[self.vertex_info[v]["seq"][directions[go]] for go, v in enumerate(vertices)])
            self.vertex_info[new_vertex]["seq"][directions[0]] = consensus_s
            self.vertex_info[new_vertex]["seq"][not directions[0]] = complementary_seq(consensus_s)

            # tags
            if copy_tags:
                for other_vertex in vertices[1:]:
                    if "tags" in self.vertex_info[other_vertex]:
                        if "tags" not in self.vertex_info[new_vertex]:
                            self.vertex_info[new_vertex]["tags"] = deepcopy(self.vertex_info[other_vertex]["tags"])
                        else:
                            for mode in self.vertex_info[other_vertex]["tags"]:
                                if mode not in self.vertex_info[new_vertex]["tags"]:
                                    self.vertex_info[new_vertex]["tags"][mode] \
                                        = deepcopy(self.vertex_info[other_vertex]["tags"][mode])
                                else:
                                    self.vertex_info[new_vertex]["tags"][mode] \
                                        |= self.vertex_info[other_vertex]["tags"][mode]
                    if "weight" in self.vertex_info[other_vertex]:
                        if "weight" not in self.vertex_info[new_vertex]:
                            self.vertex_info[new_vertex]["weight"] \
                                = deepcopy(self.vertex_info[other_vertex]["weight"])
                        else:
                            for mode in self.vertex_info[other_vertex]["weight"]:
                                if mode not in self.vertex_info[new_vertex]["weight"]:
                                    self.vertex_info[new_vertex]["weight"][mode] \
                                        = self.vertex_info[other_vertex]["weight"][mode]
                                else:
                                    self.vertex_info[new_vertex]["weight"][mode] \
                                        += self.vertex_info[other_vertex]["weight"][mode]
                    for mode in self.tagged_vertices:
                        if other_vertex in self.tagged_vertices[mode]:
                            self.tagged_vertices[mode].add(new_vertex)
                            self.tagged_vertices[mode].remove(other_vertex)
            self.remove_vertex(vertices)
            if log_handler:
                log_handler.info("Consensus made: " + new_vertex)
            else:
                log_handler.info("Consensus made: " + new_vertex + "\n")

    def processing_polymorphism(self, limited_vertices=None,
                                contamination_depth=5., contamination_similarity=0.95,
                                degenerate=False, degenerate_depth=1.5, degenerate_similarity=0.95,
                                verbose=False, debug=False, log_handler=None):
        parallel_vertices_list = self.detect_parallel_vertices(limited_vertices=limited_vertices)
        if debug:
            if log_handler:
                log_handler.info("detected parallel vertices " + str(parallel_vertices_list))
            else:
                sys.stdout.write("detected parallel vertices " + str(parallel_vertices_list) + "\n")

        degenerate_depth = abs(log(degenerate_depth))
        contamination_depth = abs(log(contamination_depth))
        contamination_dif = 1 - contamination_similarity
        degenerate_dif = 1 - degenerate_similarity

        removing_irrelevant_v = set()
        removing_contaminating_v = set()
        count_contamination_or_degenerate = 0
        for prl_vertices in parallel_vertices_list:
            this_contamination_or_polymorphic = False
            prl_vertices = sorted(prl_vertices, key=lambda x: -self.vertex_info[x[0]]["cov"])
            max_cov_vertex, direction_remained = prl_vertices.pop(0)
            max_cov_seq = self.vertex_info[max_cov_vertex]["seq"][direction_remained]
            max_cov = self.vertex_info[max_cov_vertex]["cov"]
            polymorphic_vertices_with_directions = {(max_cov_vertex, direction_remained)}
            # hard to clearly identify the biological factors
            for this_v, this_direction in prl_vertices:
                this_cov = self.vertex_info[this_v]["cov"]
                if abs(log(this_cov/max_cov)) > contamination_depth:
                    this_seq = self.vertex_info[this_v]["seq"][this_direction]
                    if abs(len(this_seq) - len(max_cov_seq)) / float(len(this_seq)) <= contamination_dif:
                        if len(max(max_cov_seq, this_seq)) > 1E6:
                            removing_irrelevant_v.add(this_v)
                            continue
                        base_dif, proper_end = \
                            find_string_difference(max_cov_seq, this_seq, max(2, int(len(this_seq) * 0.005)))
                        if float(base_dif) / len(this_seq) < contamination_dif:
                            removing_contaminating_v.add(this_v)
                            this_contamination_or_polymorphic = True
                        else:
                            removing_irrelevant_v.add(this_v)
                    else:
                        removing_irrelevant_v.add(this_v)
                elif degenerate and abs(log(this_cov/max_cov)) < degenerate_depth:
                    this_seq = self.vertex_info[this_v]["seq"][this_direction]
                    if abs(len(this_seq) - len(max_cov_seq)) / float(len(this_seq)) <= degenerate_dif:
                        if len(max(max_cov_seq, this_seq)) > 1E6:
                            continue
                        base_dif, proper_end = \
                            find_string_difference(max_cov_seq, this_seq, max(2, int(len(this_seq) * 0.005)))
                        if float(base_dif) / len(this_seq) < degenerate_dif:
                            this_contamination_or_polymorphic = True
                            if len(this_seq) == len(max_cov_seq):
                                polymorphic_vertices_with_directions.add((this_v, this_direction))
                            else:
                                raise Exception("Cannot degenerate inequal-length polymorphic contigs: EDGE_" +
                                                max_cov_vertex + " and EDGE_" + this_v + "!")
                            # else:
                            #     if log_handler:
                            #         log_handler.warning("Polymorphism: EDGE_" + max_cov_vertex + " and EDGE_" +
                            #                             this_v + "!")
                            #     else:
                            #         sys.stdout.write("Warning: Polymorphism: EDGE_" + max_cov_vertex + " and EDGE_" +
                            #                          this_v + "!\n")
            #
            if len(polymorphic_vertices_with_directions) > 1:
                self.generate_consensus_vertex([v_d[0] for v_d in polymorphic_vertices_with_directions],
                                               [v_d[1] for v_d in polymorphic_vertices_with_directions],
                                               check_parallel_vertices=False, log_handler=log_handler)
            count_contamination_or_degenerate += this_contamination_or_polymorphic
        if removing_contaminating_v:
            contaminating_cov = np.array([self.vertex_info[con_v]["cov"] for con_v in removing_contaminating_v])
            contaminating_weight = np.array([len(self.vertex_info[con_v]["seq"][True]) - self.__kmer
                                             for con_v in removing_contaminating_v])
            self.remove_vertex(removing_contaminating_v)
            cont_mean, cont_std = weighted_mean_and_std(contaminating_cov, contaminating_weight)
            cut_off = cont_mean + cont_std
            # remove (low cov and terminal vertex)
            removing_below_cut_off = []
            for del_v in self.vertex_info:
                if self.vertex_info[del_v]["cov"] < cut_off:
                    if sum([bool(cnn) for cnn in self.vertex_info[del_v]["connections"].values()]) < 2:
                        removing_below_cut_off.append(del_v)
            self.remove_vertex(removing_below_cut_off)
            if verbose or debug:
                if log_handler:
                    log_handler.info("removing contaminating vertices: " + " ".join(list(removing_contaminating_v)))
                    log_handler.info("removing contaminating-like vertices: " + " ".join(list(removing_below_cut_off)))
                else:
                    sys.stdout.write("removing contaminating vertices: " + " ".join(list(removing_contaminating_v)) + "\n")
                    sys.stdout.write("removing contaminating-like vertices: " + " ".join(list(removing_below_cut_off)) + "\n")
        if removing_irrelevant_v:
            self.remove_vertex(removing_irrelevant_v)
            if verbose or debug:
                if log_handler:
                    log_handler.info("removing parallel vertices: " + " ".join(list(removing_irrelevant_v)))
                else:
                    sys.stdout.write("removing parallel vertices: " + " ".join(list(removing_irrelevant_v)) + "\n")
        if count_contamination_or_degenerate >= 4:
            if log_handler:
                log_handler.warning("The graph might suffer from contamination or polymorphism!")
            else:
                sys.stdout.write("Warning: The graph might suffer from contamination or polymorphism!")

    def find_target_graph(self, tab_file, mode="cp", type_factor=3, weight_factor=100.0,
                          max_copy=8, min_sigma_factor=0.1,
                          log_hard_cov_threshold=10., contamination_depth=5., contamination_similarity=0.95,
                          degenerate=True, degenerate_depth=1.5, degenerate_similarity=0.95,
                          broken_graph_allowed=False, temp_graph=None, verbose=True,
                          log_handler=None, debug=False):

        self.parse_tab_file(tab_file, mode=mode, type_factor=type_factor, log_handler=log_handler)
        new_assembly = deepcopy(self)
        try:
            new_assembly.merge_all_possible_vertices()
            new_assembly.tag_in_between(mode=mode)
            new_assembly.processing_polymorphism(contamination_depth=contamination_depth,
                                                 contamination_similarity=contamination_similarity,
                                                 degenerate=False, verbose=verbose, debug=debug, log_handler=log_handler)
            # new_assembly.estimate_copy_and_depth_by_cov(new_assembly.tagged_vertices[mode],
            #                                             log_handler=log_handler, verbose=verbose, mode=mode, debug=debug)
            changed = True
            while changed:
                changed = False
                cluster_trimmed = True
                while cluster_trimmed:
                    # remove low coverages
                    first_round = True
                    delete_those_vertices = set()
                    parameters = []
                    while first_round or delete_those_vertices:
                        changed, parameters = new_assembly.filter_by_coverage(mode=mode, weight_factor=weight_factor,
                                                                              log_hard_cov_threshold=log_hard_cov_threshold,
                                                                              min_sigma_factor=min_sigma_factor,
                                                                              log_handler=log_handler,
                                                                              verbose=verbose, debug=debug)
                        new_assembly.estimate_copy_and_depth_by_cov(new_assembly.tagged_vertices[mode], debug=debug,
                                                                    log_handler=log_handler, verbose=verbose, mode=mode)
                        first_round = False

                    cluster_trimmed = False
                    # chose the target cluster (best rank)
                    if len(new_assembly.vertex_clusters) == 0:
                        raise Exception("No available target components detected!")
                    elif len(new_assembly.vertex_clusters) == 1:
                        pass
                    else:
                        cluster_weights = [sum([new_assembly.vertex_info[x_v]["weight"][mode]
                                                for x_v in x
                                                if "weight" in new_assembly.vertex_info[x_v]
                                                and mode in new_assembly.vertex_info[x_v]["weight"]])
                                           for x in new_assembly.vertex_clusters]
                        best = max(cluster_weights)
                        best_id = cluster_weights.index(best)
                        temp_cluster_weights = deepcopy(cluster_weights)
                        del temp_cluster_weights[best_id]
                        second = max(temp_cluster_weights)
                        if best < second * weight_factor:
                            if temp_graph:
                                new_assembly.write_to_gfa(temp_graph)
                                new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                            raise Exception("Multiple isolated target components detected! Broken or contamination?")
                        for j, w in enumerate(cluster_weights):
                            if w == second:
                                for del_v in new_assembly.vertex_clusters[j]:
                                    if del_v in new_assembly.tagged_vertices[mode]:
                                        new_cov = new_assembly.vertex_info[del_v]["cov"]
                                        for mu, sigma in parameters:
                                            if abs(new_cov - mu) < sigma:
                                                if temp_graph:
                                                    new_assembly.write_to_gfa(temp_graph)
                                                    new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                                                raise Exception("Complicated graph: please check around EDGE_" + del_v + "!"
                                                                "\ntags: " + str(new_assembly.vertex_info[del_v]["tags"][mode]))

                        # remove other clusters
                        vertices_to_del = set()
                        for go_cl, v_2_del in enumerate(new_assembly.vertex_clusters):
                            if go_cl != best_id:
                                vertices_to_del |= v_2_del
                        if vertices_to_del:
                            if verbose or debug:
                                if log_handler:
                                    log_handler.info("removing other clusters: " + str(vertices_to_del))
                                else:
                                    sys.stdout.write("removing other clusters: " + str(vertices_to_del) + "\n")
                            new_assembly.remove_vertex(vertices_to_del)
                            cluster_trimmed = True
                            changed = True

                # merge vertices
                new_assembly.merge_all_possible_vertices()
                new_assembly.tag_in_between(mode=mode)

                # no terminal vertices allowed
                first_round = True
                delete_those_vertices = set()
                while first_round or delete_those_vertices:
                    first_round = False
                    delete_those_vertices = set()
                    for vertex_name in new_assembly.vertex_info:
                        # both ends must have edge(s)
                        if sum([bool(len(cn)) for cn in new_assembly.vertex_info[vertex_name]["connections"].values()]) != 2:
                            if vertex_name in new_assembly.tagged_vertices[mode]:
                                if temp_graph:
                                    new_assembly.write_to_gfa(temp_graph)
                                    new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                                raise Exception("Incomplete/Complicated graph: please check around EDGE_" + vertex_name + "!")
                            else:
                                delete_those_vertices.add(vertex_name)
                    if delete_those_vertices:
                        if verbose or debug:
                            if log_handler:
                                log_handler.info("removing terminal contigs: " + str(delete_those_vertices))
                            else:
                                sys.stdout.write("removing terminal contigs: " + str(delete_those_vertices) + "\n")
                        new_assembly.remove_vertex(delete_those_vertices)
                        changed = True

                # merge vertices
                new_assembly.merge_all_possible_vertices()
                new_assembly.processing_polymorphism(contamination_depth=contamination_depth,
                                                     contamination_similarity=contamination_similarity,
                                                     degenerate=False, degenerate_depth=degenerate_depth,
                                                     degenerate_similarity=degenerate_similarity,
                                                     verbose=verbose, debug=debug, log_handler=log_handler)
                new_assembly.tag_in_between(mode=mode)

            if temp_graph:
                new_assembly.write_to_gfa(temp_graph)
                new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
            new_assembly.processing_polymorphism(contamination_depth=contamination_depth,
                                                 contamination_similarity=contamination_similarity,
                                                 degenerate=degenerate, degenerate_depth=degenerate_depth,
                                                 degenerate_similarity=degenerate_similarity,
                                                 verbose=verbose, debug=debug, log_handler=log_handler)
            new_assembly.merge_all_possible_vertices()
            if temp_graph:
                new_assembly.write_to_gfa(temp_graph)
                new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")

            # create idealized vertices and edges
            new_average_cov = new_assembly.estimate_copy_and_depth_by_cov(log_handler=log_handler, verbose=verbose,
                                                                          mode="all", debug=debug)
            try:
                if verbose:
                    if log_handler:
                        log_handler.info("Estimating copy and depth precisely ...")
                    else:
                        sys.stdout.write("Estimating copy and depth precisely ...\n")
                copy_results = new_assembly.estimate_copy_and_depth_precisely(maximum_copy_num=max_copy,
                                                                              log_handler=log_handler, verbose=verbose,
                                                                              debug=debug)
            except RecursionError:
                if temp_graph:
                    new_assembly.write_to_gfa(temp_graph)
                    new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                raise Exception("Complicated target graph! Detecting path(s) failed!\n")
            else:
                return copy_results
        except KeyboardInterrupt as e:
            if temp_graph:
                new_assembly.write_to_gfa(temp_graph)
                new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
            raise KeyboardInterrupt

    def get_all_circular_paths(self, mode="cp", library_info=None, log_handler=None):

        def circular_directed_graph_solver(ongoing_path, next_connections, vertices_left):
            # print("-----------------------------")
            # print("ongoing_path", ongoing_path)
            # print("next_connect", next_connections)
            # print("vertices_lef", vertices_left)
            for next_vertex, next_end in next_connections:
                # print("next_vertex", next_vertex)
                if next_vertex in vertices_left:
                    new_path = deepcopy(ongoing_path)
                    new_left = deepcopy(vertices_left)
                    new_path.append((next_vertex, not next_end))
                    new_left[next_vertex] -= 1
                    if not new_left[next_vertex]:
                        del new_left[next_vertex]
                    new_connections = self.vertex_info[next_vertex]["connections"][not next_end]
                    if not new_left:
                        if (start_vertex, False) in new_connections:
                            paths.add(tuple(new_path))
                            return
                        else:
                            return
                    else:
                        circular_directed_graph_solver(new_path, new_connections, new_left)

        paths = set()
        # start from a single copy vertex
        # print(self.copy_to_vertex)
        if 1 not in self.copy_to_vertex:
            raise Exception("No single copy region?! Detecting path(s) failed!\n")
        start_vertex = sorted(self.copy_to_vertex[1], key=lambda x: -self.vertex_info[x]["len"])[0]
        # each contig stored format:
        first_path = [(start_vertex, True)]
        first_connections = self.vertex_info[start_vertex]["connections"][True]
        vertex_to_copy = deepcopy(self.vertex_to_copy)
        del vertex_to_copy[start_vertex]
        circular_directed_graph_solver(first_path, first_connections, vertex_to_copy)

        if not paths:
            raise Exception("Detecting path(s) from remaining graph failed!\n")
        else:

            # sorting path by average distance among multi-copy loci
            # the highest would be more symmetrical IR, which turns out to be more reasonable
            sorted_paths = []
            total_len = len(list(paths)[0])
            record_pattern = False
            for this_path in paths:
                acc_dist = 0
                for copy_num in self.copy_to_vertex:
                    if copy_num > 2:
                        record_pattern = True
                        for vertex_name in self.copy_to_vertex[copy_num]:
                            loc_ids = [go_to_id for go_to_id, (v, e) in enumerate(this_path) if v == vertex_name]
                            for id_a, id_b in combinations(loc_ids, 2):
                                acc_dist += min((id_a - id_b) % total_len, (id_b - id_a) % total_len)
                sorted_paths.append((this_path, acc_dist))
            if record_pattern:
                sorted_paths.sort(key=lambda x: -x[1])
                pattern_dict = {acc_distance: ad_id + 1
                                for ad_id, acc_distance in enumerate(sorted(set([x[1] for x in sorted_paths]), reverse=True))}
                if len(pattern_dict) > 1:
                    if log_handler:
                        log_handler.warning("Multiple repeat patterns appeared in your data, "
                                            "a more balanced pattern (always the repeat_pattern1) would be suggested "
                                            "for plastomes with inverted repeats!")
                    else:
                        sys.stdout.write("Warning: Multiple repeat patterns appeared in your data, "
                                         "a more balanced pattern (always the repeat_pattern1) would be suggested "
                                         "for plastomes with inverted repeats!\n")
                sorted_paths = [(this_path, ".repeat_pattern" + str(pattern_dict[acc_distance]))
                                for this_path, acc_distance in sorted_paths]
            else:
                sorted_paths = [(this_path, "") for this_path in sorted(paths)]

            if mode == "cp":
                if len(sorted_paths) > 2 and not (100000 < len(self.export_path(sorted_paths[0][0]).seq) < 200000):
                    if log_handler:
                        log_handler.warning("Multiple paths with abnormal plastome length produced!")
                        log_handler.warning("Please check the assembly graph and selected graph to confirm.")
                    else:
                        sys.stdout.write("Warning: Multiple paths with abnormal plastome length produced!\n")
                        sys.stdout.write("Please check the assembly graph and selected graph to confirm.\n")
                elif len(sorted_paths) > 2:
                    if log_handler:
                        log_handler.warning("Multiple paths produced!")
                        log_handler.warning("Please check the existence of those isomers "
                                            "by using reads mapping (library information) or longer reads.")
                    else:
                        sys.stdout.write("Warning: Multiple paths produced!\n")
                        sys.stdout.write("Please check the existence of those isomers by "
                                         "using reads mapping (library information) or longer reads.\n")
                elif len(sorted_paths) > 1:
                    if log_handler:
                        log_handler.warning("More than one path produced ...")
                        log_handler.warning("Please check the final result to confirm whether they are "
                                            "simply flip-flop configurations!")
                    else:
                        sys.stdout.write("More than one path produced ...\n")
                        sys.stdout.write("Please check the final result to confirm whether they are "
                                         "simply flip-flop configurations!\n")
            return sorted_paths

    def export_path(self, in_path):
        seq_names = []
        seq_segments = []
        for this_vertex, this_end in in_path:
            seq_segments.append(self.vertex_info[this_vertex]["seq"][this_end][self.__kmer:])
            seq_names.append(this_vertex + ("-", "+")[this_end])
        # if not circular
        if (in_path[0][0], not in_path[0][1]) not in self.vertex_info[in_path[-1][0]]["connections"][in_path[-1][1]]:
            seq_segments[0] = self.vertex_info[in_path[0][0]]["seq"][in_path[0][1]][:self.__kmer] + seq_segments[0]
        else:
            seq_names[-1] += "(circular)"
        return Sequence(",".join(seq_names), "".join(seq_segments))
