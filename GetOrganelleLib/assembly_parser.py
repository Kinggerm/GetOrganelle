import os
import sys
import time
import random
from copy import deepcopy
from itertools import combinations, product
try:
    from sympy import Symbol, solve, lambdify
    from scipy import optimize
except ImportError:
    def Symbol(foo, integer):
        raise ImportError("No module named sympy")
    def solve(foo1, foo2):
        raise ImportError("No module named sympy")
    def lambdify(args=None, expr=None):
        raise ImportError("No module named sympy")
    class optimize:
        def __init__(self):
            pass
        def minimize(self, fun=None, x0=None, jac=None, method=None, bounds=None, constraints=None, options=None):
            raise ImportError("No module named scipy")

path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.statistical_func import *
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


class Vertex(object):
    def __init__(self, v_name, length=None, coverage=None, forward_seq=None, reverse_seq=None,
                 tail_connections=None, head_connections=None, fastg_form_long_name=None):
        """
        :param v_name:
        :param length:
        :param coverage:
        :param forward_seq:
        :param reverse_seq:
        :param tail_connections: set()
        :param head_connections: set()
        :param fastg_form_long_name:
        self.seq={True: FORWARD_SEQ, False: REVERSE_SEQ}
        self.connections={True: tail_connection_set, False: head_connection_set}
        """
        self.name = v_name
        self.len = length
        self.cov = coverage
        """ True: forward, False: reverse """
        if forward_seq and reverse_seq:
            assert forward_seq == complementary_seq(reverse_seq), "forward_seq != complementary_seq(reverse_seq)"
            self.seq = {True: forward_seq, False: reverse_seq}
        elif forward_seq:
            self.seq = {True: forward_seq, False: complementary_seq(forward_seq)}
        elif reverse_seq:
            self.seq = {True: complementary_seq(reverse_seq), False: reverse_seq}
        else:
            self.seq = {True: None, False: None}
        """ True: tail, False: head """
        self.connections = {True: set(), False: set()}
        if tail_connections:
            self.connections[True] = tail_connections
        if head_connections:
            self.connections[False] = head_connections
        self.fastg_form_name = fastg_form_long_name
        self.other_attr = {}

    def fill_fastg_form_name(self, check_valid=False):
        if check_valid:
            if not str(self.name).isdigit():
                raise ValueError("Invalid vertex name for fastg format!")
            if not isinstance(self.len, int):
                raise ValueError("Invalid vertex length for fastg format!")
            if not (isinstance(self.cov, int) or isinstance(self.cov, float)):
                raise ValueError("Invalid vertex coverage for fastg format!")
        self.fastg_form_name = \
            "EDGE_" + str(self.name) + "_length_" + str(self.len) + "_cov_" + str(round(self.cov, 5))

    def is_self_loop(self):
        return (self.name, False) in self.connections[True]


class VertexInfo(dict):
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            if not isinstance(val, Vertex):
                raise ValueError("Value must be a Vertex type! Current: " + str(type(val)))
        dict.__init__(kwargs)

    def __setitem__(self, key, val):
        if not isinstance(val, Vertex):
            raise ValueError("Value must be a Vertex type! Current: " + str(type(val)))
        dict.__setitem__(self, key, val)


class Assembly(object):
    def __init__(self, graph_file=None, min_cov=0., max_cov=inf):
        """
        :param graph_file:
        :param min_cov:
        :param max_cov:
        """
        self.vertex_info = VertexInfo()
        self.__kmer = 127
        if graph_file:
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
            res.append(">" + v + "__" + str(self.vertex_info[v].len) + "__" + str(self.vertex_info[v].cov))
            for e in (False, True):
                if len(self.vertex_info[v].connections[e]):
                    res.append("(" + ["head", "tail"][e] + ":")
                    res.append(",".join([next_v + "_" + ["head", "tail"][next_e]
                                         for next_v, next_e in self.vertex_info[v].connections[e]]))
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
                        self.vertex_info[vertex_name] = Vertex(vertex_name, seq_len, seq_cov, sequence)
                        # self.vertex_info[vertex_name] = {"len": seq_len, "cov": seq_cov,
                        #                                  "connections": {True: set(), False: set()},
                        #                                  "seq": {True: sequence, False: complementary_seq(sequence)}}
                        if vertex_name.isdigit():
                            self.vertex_info[vertex_name].fill_fastg_form_name()
                            # self.vertex_info[vertex_name]["long"] = \
                            #     "EDGE_" + vertex_name + "_length_" + str(seq_len) + "_cov_" + str(round(seq_cov, 5))
                elif line.startswith("L\t"):
                    flag, vertex_1, end_1, vertex_2, end_2, kmer_val = line.strip().split("\t")
                    # "head"~False, "tail"~True
                    end_1 = {"+": True, "-": False}[end_1]
                    end_2 = {"+": False, "-": True}[end_2]
                    kmer_values.add(kmer_val)
                    self.vertex_info[vertex_1].connections[end_1].add((vertex_2, end_2))
                    self.vertex_info[vertex_2].connections[end_2].add((vertex_1, end_1))
            if len(kmer_values) != 1:
                raise ProcessingGraphFailed("Multiple overlap values: " + ",".join(sorted(kmer_values)))
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
                self.vertex_info[vertex_name] = Vertex(vertex_name, int(vertex_len), vertex_cov,
                                                       fastg_form_long_name=this_vertex_str.strip("'"))
                # self.vertex_info[vertex_name] = {"len": int(vertex_len),
                #                                  "cov": vertex_cov,
                #                                  "long": this_vertex_str.strip("'")}
            # if "connections" not in self.vertex_info[vertex_name]:
            #     self.vertex_info[vertex_name]["connections"] = {True: set(), False: set()}

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
                            self.vertex_info[vertex_name].connections[this_end].add((next_name, next_end))
                            self.vertex_info[next_name].connections[next_end].add((vertex_name, this_end))
                # sequence
                if not self.vertex_info[vertex_name].seq[True]:
                    # self.vertex_info[vertex_name]["seq"] = {}
                    if this_end:
                        self.vertex_info[vertex_name].seq[True] = seq.seq
                        self.vertex_info[vertex_name].seq[False] = complementary_seq(seq.seq)
                    else:
                        self.vertex_info[vertex_name].seq[True] = complementary_seq(seq.seq)
                        self.vertex_info[vertex_name].seq[False] = seq.seq

        """detect general kmer"""
        ## find initial kmer candidate values
        initial_kmer = set()
        no_connection_at_all = True
        for vertex_name in self.vertex_info:
            if sum([len(self.vertex_info[vertex_name].connections[this_e]) for this_e in (True, False)]) != 0:
                no_connection_at_all = False
                for this_e in (True, False):
                    for next_name, next_end in self.vertex_info[vertex_name].connections[this_e]:
                        for test_k in range(21, 128, 2):
                            this_seq = self.vertex_info[vertex_name].seq[this_e][-test_k:]
                            next_seq = self.vertex_info[next_name].seq[not next_end][:test_k]
                            if this_seq == next_seq:
                                initial_kmer.add(test_k)
                        break
                    if initial_kmer:
                        break
            if initial_kmer:
                break
        if no_connection_at_all:
            self.__kmer = 0
        else:
            ## check all edges
            testing_vertices = set(self.vertex_info)
            while initial_kmer and testing_vertices:
                vertex_name = testing_vertices.pop()
                for this_end in (True, False):
                    for next_name, next_end in self.vertex_info[vertex_name].connections[this_end]:
                        for test_k in list(initial_kmer):
                            this_seq = self.vertex_info[vertex_name].seq[this_end][-test_k:]
                            next_seq = self.vertex_info[next_name].seq[not next_end][:test_k]
                            if this_seq != next_seq:
                                initial_kmer.discard(test_k)
            if len(initial_kmer) >= 1:
                self.__kmer = max(initial_kmer)
            else:
                raise ProcessingGraphFailed("No kmer detected!")

    def new_graph_with_vertex_reseeded(self, start_from=1):
        those_vertices = sorted(self.vertex_info)
        new_graph = Assembly()
        name_trans = {those_vertices[go - start_from]: str(go)
                      for go in range(start_from, start_from + len(those_vertices))}
        for old_name in those_vertices:
            new_name = name_trans[old_name]
            this_v_info = deepcopy(self.vertex_info[old_name])
            this_v_info.connections = {True: set(), False: set()}
            for this_end in self.vertex_info[old_name].connections:
                for next_name, next_end in self.vertex_info[old_name].connections[this_end]:
                    this_v_info.connections[this_end].add((name_trans[next_name], next_end))
            this_v_info.fill_fastg_form_name()
            new_graph.vertex_info[new_name] = this_v_info
        return new_graph, name_trans

    def write_to_fastg(self, out_file, check_postfix=True,
                       rename_if_needed=False, out_renaming_table=None, echo_rename_warning=False, log_handler=None):
        if check_postfix and not out_file.endswith(".fastg"):
            out_file += ".fastg"
        try:
            out_matrix = SequenceList()
            for vertex_name in self.vertex_info:
                this_name = self.vertex_info[vertex_name].fastg_form_name
                for this_end in (False, True):
                    seq_name = [this_name, ("", "'")[not this_end]]
                    if self.vertex_info[vertex_name].connections[this_end]:
                        seq_name.append(":")
                        connect_str = ",".join([self.vertex_info[n_v].fastg_form_name + ("", "'")[n_e]
                                                for n_v, n_e in self.vertex_info[vertex_name].connections[this_end]])
                        seq_name.append(connect_str)
                    seq_name.append(";")
                    out_matrix.append(Sequence("".join(seq_name), self.vertex_info[vertex_name].seq[this_end]))
            out_matrix.interleaved = 70
            out_matrix.write_fasta(out_file)
        except TypeError:
            if rename_if_needed:
                if echo_rename_warning:
                    if log_handler:
                        log_handler.info("Graph converted to new fastg with original Vertex names lost.")
                    else:
                        sys.stdout.write("Graph converted to new fastg with original Vertex names lost.")
                new_graph, name_trans = self.new_graph_with_vertex_reseeded()
                new_graph.write_to_fastg(out_file, check_postfix=False)
                if out_renaming_table:
                    with open(out_renaming_table + ".Temp", "w") as out_table:
                        for old_name in sorted(name_trans):
                            out_table.write(old_name + "\t" + name_trans[old_name] + "\n")
                    os.rename(out_renaming_table + ".Temp", out_renaming_table)
                    if echo_rename_warning:
                        if log_handler:
                            log_handler.info("Table (original Vertex names -> new Vertex names) written to " +
                                             out_renaming_table + ".")
                        else:
                            sys.stdout.write("Table (original Vertex names -> new Vertex names) written to " +
                                             out_renaming_table + ".")
            else:
                raise ProcessingGraphFailed(
                    "Merged graph cannot be written as fastg format file, please try gfa format!")

    def write_to_fasta(self, out_file, check_postfix=True):
        if check_postfix and not out_file.endswith(".fasta"):
            out_file += ".fasta"
        out_matrix = SequenceList()
        for vertex_name in self.vertex_info:
            out_matrix.append(Sequence(vertex_name, self.vertex_info[vertex_name].seq[True]))
        out_matrix.interleaved = 70
        out_matrix.write_fasta(out_file)

    def write_to_gfa(self, out_file, check_postfix=True):
        if check_postfix and not out_file.endswith(".gfa"):
            out_file += ".gfa"
        out_file_handler = open(out_file, "w")
        for vertex_name in self.vertex_info:
            out_file_handler.write("\t".join([
                "S", vertex_name, self.vertex_info[vertex_name].seq[True],
                "LN:i:"+str(self.vertex_info[vertex_name].len),
                "RC:i:"+str(int(self.vertex_info[vertex_name].len*self.vertex_info[vertex_name].cov))
                ]) + "\n")
        recorded_connections = set()
        for vertex_name in self.vertex_info:
            for this_end in (False, True):
                for next_v, next_e in self.vertex_info[vertex_name].connections[this_end]:
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
            if "tags" in self.vertex_info[this_vertex].other_attr:
                all_tags = self.vertex_info[this_vertex].other_attr["tags"]
                all_tag_list = sorted(all_tags)
                all_weights = self.vertex_info[this_vertex].other_attr["weight"]
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
            for connected_set in self.vertex_info[this_vertex].connections.values():
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
        for vertex_name in vertices:
            for this_end, connected_set in self.vertex_info[vertex_name].connections.items():
                for next_v, next_e in connected_set:
                    self.vertex_info[next_v].connections[next_e].remove((vertex_name, this_end))
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
            this_cons = self.vertex_info[limited_vertex].connections
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
                    this_cons = self.vertex_info[each_vertex].connections
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

    def is_sequential_repeat(self, search_vertex_name, return_pair_in_the_trunk_path=True):
        if search_vertex_name not in self.vertex_info:
            raise ProcessingGraphFailed("Vertex name " + search_vertex_name + " not found!")
        connection_set_t = self.vertex_info[search_vertex_name].connections[True]
        connection_set_f = self.vertex_info[search_vertex_name].connections[False]
        all_pairs_of_inner_circles = []
        
        def path_without_leakage(start_v, start_e, terminating_end_set):
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
                    for n_in_search_v, n_in_search_e in self.vertex_info[in_search_v].connections[in_search_e]:
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
            for next_t_v, next_t_e in sorted(connection_set_t):
                this_inner_circle = path_without_leakage(next_t_v, next_t_e, connection_set_f)
                if this_inner_circle:
                    # check leakage in reverse direction
                    reverse_v, reverse_e = this_inner_circle[0][1]
                    not_leak = path_without_leakage(reverse_v, reverse_e, connection_set_t)
                    if not_leak:
                        all_pairs_of_inner_circles.extend(this_inner_circle)
            # sort pairs by average depths(?)
            all_pairs_of_inner_circles.sort(
                key=lambda x: (self.vertex_info[x[0][0]].cov + self.vertex_info[x[1][0]].cov))
            if all_pairs_of_inner_circles and return_pair_in_the_trunk_path:
                # switch nearby vertices
                # keep those prone to be located in the "trunk road" of the repeat
                single_pair_in_main_path = []
                if len(all_pairs_of_inner_circles) == 1:
                    for next_v, next_e in sorted(connection_set_t) + sorted(connection_set_f):
                        if (next_v, next_e) not in all_pairs_of_inner_circles[0]:
                            single_pair_in_main_path.append((next_v, next_e))
                    single_pair_in_main_path = tuple(single_pair_in_main_path)
                else:
                    # two circles share this sequential repeat,
                    # return the one with a smaller average depth(?)
                    single_pair_in_main_path = tuple(all_pairs_of_inner_circles[0])
                return single_pair_in_main_path
            return all_pairs_of_inner_circles
        else:
            return all_pairs_of_inner_circles

    def merge_all_possible_vertices(self, limited_vertices=None, copy_tags=True):
        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)
        merged = False
        while limited_vertices:
            this_vertex = limited_vertices.pop()
            for this_end in (True, False):
                connected_set = self.vertex_info[this_vertex].connections[this_end]
                if len(connected_set) == 1:
                    next_vertex, next_end = list(connected_set)[0]
                    if len(self.vertex_info[next_vertex].connections[next_end]) == 1 and this_vertex != next_vertex:
                        # reverse the names
                        merged = True
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
                        self.vertex_info[new_vertex].fastg_form_name = None
                        # if "long" in self.vertex_info[new_vertex]:
                        #     del self.vertex_info[new_vertex]["long"]
                        # connections
                        self.vertex_info[new_vertex].connections[this_end] \
                            = deepcopy(self.vertex_info[next_vertex].connections[not next_end])
                        if (this_vertex, not this_end) in self.vertex_info[new_vertex].connections[this_end]:
                            self.vertex_info[new_vertex].connections[this_end].remove((this_vertex, not this_end))
                            self.vertex_info[new_vertex].connections[this_end].add((new_vertex, not this_end))
                        for new_end in (True, False):
                            for n_n_v, n_n_e in self.vertex_info[new_vertex].connections[new_end]:
                                self.vertex_info[n_n_v].connections[n_n_e].add((new_vertex, new_end))
                        # len & cov
                        this_len = self.vertex_info[this_vertex].len
                        next_len = self.vertex_info[next_vertex].len
                        this_cov = self.vertex_info[this_vertex].cov
                        next_cov = self.vertex_info[next_vertex].cov
                        self.vertex_info[new_vertex].len = this_len + next_len - self.__kmer
                        self.vertex_info[new_vertex].cov = \
                            ((this_len - self.__kmer + 1) * this_cov + (next_len - self.__kmer + 1) * next_cov) \
                            / ((this_len - self.__kmer + 1) + (next_len - self.__kmer + 1))
                        self.vertex_info[new_vertex].seq[this_end] \
                            += self.vertex_info[next_vertex].seq[not next_end][self.__kmer:]
                        self.vertex_info[new_vertex].seq[not this_end] \
                            = self.vertex_info[next_vertex].seq[next_end][:-self.__kmer] \
                            + self.vertex_info[this_vertex].seq[not this_end]
                        # tags
                        if copy_tags:
                            if "tags" in self.vertex_info[next_vertex].other_attr:
                                if "tags" not in self.vertex_info[new_vertex].other_attr:
                                    self.vertex_info[new_vertex].other_attr["tags"] = deepcopy(self.vertex_info[next_vertex].other_attr["tags"])
                                else:
                                    for mode in self.vertex_info[next_vertex].other_attr["tags"]:
                                        if mode not in self.vertex_info[new_vertex].other_attr["tags"]:
                                            self.vertex_info[new_vertex].other_attr["tags"][mode] \
                                                = deepcopy(self.vertex_info[next_vertex].other_attr["tags"][mode])
                                        else:
                                            self.vertex_info[new_vertex].other_attr["tags"][mode] \
                                                |= self.vertex_info[next_vertex].other_attr["tags"][mode]
                            if "weight" in self.vertex_info[next_vertex].other_attr:
                                if "weight" not in self.vertex_info[new_vertex].other_attr:
                                    self.vertex_info[new_vertex].other_attr["weight"] \
                                        = deepcopy(self.vertex_info[next_vertex].other_attr["weight"])
                                else:
                                    for mode in self.vertex_info[next_vertex].other_attr["weight"]:
                                        if mode not in self.vertex_info[new_vertex].other_attr["weight"]:
                                            self.vertex_info[new_vertex].other_attr["weight"][mode] \
                                                = self.vertex_info[next_vertex].other_attr["weight"][mode]
                                        else:
                                            self.vertex_info[new_vertex].other_attr["weight"][mode] \
                                                += self.vertex_info[next_vertex].other_attr["weight"][mode]
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
        return merged

    def estimate_copy_and_depth_by_cov(self, limited_vertices=None, given_average_cov=None, mode="embplant_pt",
                                       re_initialize=False, log_handler=None, verbose=True, debug=False):
        if mode == "embplant_pt":
            max_majority_cov = 2
        elif mode == "other_pt":
            max_majority_cov = 10
        elif mode == "embplant_mt":
            max_majority_cov = 4
        elif mode == "embplant_nr":
            max_majority_cov = 2
        elif mode == "animal_mt":
            max_majority_cov = 4
        elif mode == "fungus_mt":
            max_majority_cov = 8
        elif mode == "all":
            max_majority_cov = 100
        else:
            max_majority_cov = 100

        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)

        if re_initialize:
            for vertex_name in limited_vertices:
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    self.vertex_to_copy[vertex_name] = 1
                    self.vertex_to_float_copy[vertex_name] = 1.
                    if 1 not in self.copy_to_vertex:
                        self.copy_to_vertex[1] = set()
                    self.copy_to_vertex[1].add(vertex_name)

        if not given_average_cov:
            previous_val = {0.}
            new_val = -1.
            while round(new_val, 5) not in previous_val:
                previous_val.add(round(new_val, 5))
                # estimate baseline depth
                total_product = 0.
                total_len = 0
                for vertex_name in limited_vertices:
                    this_len = (self.vertex_info[vertex_name].len - self.__kmer + 1) \
                               * self.vertex_to_copy.get(vertex_name, 1)
                    this_cov = self.vertex_info[vertex_name].cov / self.vertex_to_copy.get(vertex_name, 1)
                    total_len += this_len
                    total_product += this_len * this_cov
                new_val = total_product/total_len
                # print("new val: ", new_val)
                # adjust this_copy according to new baseline depth
                for vertex_name in self.vertex_info:
                    if vertex_name in self.vertex_to_copy:
                        old_copy = self.vertex_to_copy[vertex_name]
                        self.copy_to_vertex[old_copy].remove(vertex_name)
                        if not self.copy_to_vertex[old_copy]:
                            del self.copy_to_vertex[old_copy]
                    this_float_copy = self.vertex_info[vertex_name].cov / new_val
                    this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_cov)
                    self.vertex_to_float_copy[vertex_name] = this_float_copy
                    self.vertex_to_copy[vertex_name] = this_copy
                    if this_copy not in self.copy_to_vertex:
                        self.copy_to_vertex[this_copy] = set()
                    self.copy_to_vertex[this_copy].add(vertex_name)
            if debug or verbose:
                if log_handler:
                    log_handler.info("updating average " + mode + " kmer-coverage: " + str(round(new_val, 2)))
                else:
                    sys.stdout.write("updating average " + mode + " kmer-coverage: " + str(round(new_val, 2)) + "\n")
            # print("return ", new_val)
            return new_val
        else:
            # adjust this_copy according to user-defined depth
            for vertex_name in self.vertex_info:
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    if not self.copy_to_vertex[old_copy]:
                        del self.copy_to_vertex[old_copy]
                this_float_copy = self.vertex_info[vertex_name].cov / given_average_cov
                this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_cov)
                self.vertex_to_float_copy[vertex_name] = this_float_copy
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)
            return given_average_cov

    def estimate_copy_and_depth_precisely(self, maximum_copy_num=10, broken_graph_allowed=False,
                                          return_new_graphs=True, verbose=True, log_handler=None, debug=False,
                                          target_name_for_log="target"):

        def get_formula(from_vertex, from_end, back_to_vertex, back_to_end):
            result_form = vertex_to_symbols[from_vertex]
            # if itself (from_vertex == back_to_vertex) form a loop, skipped
            if from_vertex != back_to_vertex:
                for next_v, next_e in self.vertex_info[from_vertex].connections[from_end]:
                    # if itself (next_v == from_vertex) form a loop, add a pseudo vertex
                    if (next_v, next_e) == (from_vertex, not from_end):
                        pseudo_self_circle_str = "P" + from_vertex
                        if pseudo_self_circle_str not in extra_str_to_symbol:
                            extra_str_to_symbol[pseudo_self_circle_str] = Symbol(pseudo_self_circle_str, integer=True)
                            extra_symbol_to_str[extra_str_to_symbol[pseudo_self_circle_str]] = pseudo_self_circle_str
                        result_form -= (extra_str_to_symbol[pseudo_self_circle_str] - 1)
                    elif (next_v, next_e) != (back_to_vertex, back_to_end):
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

        def constraint_min_function_for_customized_brute(x):
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
            if in_log_handler:
                if debug or display_p:
                    in_log_handler.info("Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round))
                    in_log_handler.info("Brute best function value: " + str(best_fun_val))
                if debug:
                    in_log_handler.info("Best solution: " + str(best_para_val))
            else:
                if debug or display_p:
                    sys.stdout.write("Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round) + "\n")
                    sys.stdout.write("Brute best function value: " + str(best_fun_val) + "\n")
                if debug:
                    sys.stdout.write("Best solution: " + str(best_para_val) + "\n")
            return best_para_val

        vertices_list = sorted(self.vertex_info)
        if len(vertices_list) == 1:
            cov_ = self.vertex_info[vertices_list[0]].cov
            if return_new_graphs:
                return [{"graph": deepcopy(self), "cov": cov_}]
            else:
                if log_handler:
                    log_handler.info("Average " + target_name_for_log +" kmer-coverage = " + str(round(cov_, 2)))
                else:
                    sys.stdout.write("Average " + target_name_for_log +" kmer-coverage = " + str(round(cov_, 2)) + "\n")
                return

        """ create constraints by creating multivariate equations """
        vertex_to_symbols = {vertex_name: Symbol("V"+vertex_name, integer=True)  # positive=True)
                             for vertex_name in vertices_list}
        symbols_to_vertex = {vertex_to_symbols[vertex_name]: vertex_name for vertex_name in vertices_list}
        extra_str_to_symbol = {}
        extra_symbol_to_str = {}
        formulae = []
        recorded_ends = set()
        for vertex_name in vertices_list:
            for this_end in (True, False):
                if (vertex_name, this_end) not in recorded_ends:
                    recorded_ends.add((vertex_name, this_end))
                    connection_set = self.vertex_info[vertex_name].connections[this_end]
                    if connection_set:  # len([n_v for n_v, n_e in connection_set if n_v in vertices_set]):
                        this_formula = vertex_to_symbols[vertex_name]
                        for n_v, n_e in connection_set:
                            # if n_v in vertices_set:
                            recorded_ends.add((n_v, n_e))
                            direct = ["_tail", "_head"]
                            try:
                                this_formula -= get_formula(n_v, n_e, vertex_name, this_end)
                                if verbose:
                                    if log_handler:
                                        log_handler.info("formulating for: " + n_v + direct[n_e] + "->" +
                                                         vertex_name + direct[this_end] + ": " + str(this_formula))
                                    else:
                                        sys.stdout.write("formulating for: " + n_v + direct[n_e] + "->" +
                                                         vertex_name + direct[this_end] + ": " + str(this_formula)+"\n")
                            except RecursionError:

                                if log_handler:
                                    log_handler.warning("formulating for: " + n_v + direct[n_e] + "->" +
                                                        vertex_name + direct[this_end] + " failed!")
                                else:
                                    sys.stdout.write("formulating for: " + n_v + direct[n_e] + "->" +
                                                     vertex_name + direct[this_end] + " failed!\n")
                                raise ProcessingGraphFailed("RecursionError!")
                        formulae.append(this_formula)
                    elif broken_graph_allowed:
                        # Extra limitation to force terminal vertex to have only one copy, to avoid over-estimation
                        # Under-estimation would not be a problem here,
                        # because the True-multiple-copy vertex would simply have no other connections,
                        # or failed in the following estimation if it does
                        formulae.append(vertex_to_symbols[vertex_name] - 1)

        # add following extra limitation
        # set cov_sequential_repeat = x*near_by_cov, x is an integer
        for vertex_name in vertices_list:
            single_pair_in_the_trunk_path = self.is_sequential_repeat(vertex_name)
            if single_pair_in_the_trunk_path:
                (from_v, from_e), (to_v, to_e) = single_pair_in_the_trunk_path
                # from_v and to_v are already in the "trunk path", if they are the same,
                # the graph is like two circles sharing the same sequential repeat, no need to add this limitation
                if from_v != to_v:
                    new_str = "E" + str(len(extra_str_to_symbol))
                    extra_str_to_symbol[new_str] = Symbol(new_str, integer=True)
                    extra_symbol_to_str[extra_str_to_symbol[new_str]] = new_str
                    formulae.append(vertex_to_symbols[vertex_name] -
                                    vertex_to_symbols[from_v] * extra_str_to_symbol[new_str])
        if verbose or debug:
            if log_handler:
                log_handler.info("formulae: " + str(formulae))
            else:
                sys.stdout.write("formulae: " + str(formulae) + "\n")

        # solve the equations
        all_v_symbols = list(symbols_to_vertex)
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
            raise ProcessingGraphFailed("Incomplete/Complicated " + target_name_for_log + " graph (1)!")
        elif type(copy_solution) == list:
            if len(copy_solution) > 2:
                raise ProcessingGraphFailed("Incomplete/Complicated " + target_name_for_log + " graph (2)!")
            else:
                copy_solution = copy_solution[0]

        free_copy_variables = list()
        for symbol_used in all_symbols:
            if symbol_used not in copy_solution:
                free_copy_variables.append(symbol_used)
                copy_solution[symbol_used] = symbol_used
        if verbose:
            if log_handler:
                log_handler.info("Copy equations: " + str(copy_solution))
            else:
                sys.stdout.write("Copy equations: " + str(copy_solution) + "\n")

        # """ minimizing equation-based copy values and their deviations from coverage-based copy values """
        """ minimizing equation-based copy's deviations from coverage-based copy values """
        least_square_expr = 0
        for symbol_used in all_v_symbols:
            # least_square_expr += copy_solution[symbol_used]
            this_vertex = symbols_to_vertex[symbol_used]
            this_copy = self.vertex_to_float_copy[this_vertex]
            least_square_expr += (copy_solution[symbol_used] - this_copy) ** 2  # * self.vertex_info[this_vertex]["len"]
        least_square_function = lambdify(args=free_copy_variables, expr=least_square_expr)

        # for safe running
        if len(free_copy_variables) > 10:
            raise ProcessingGraphFailed("Free variable > 10 is not accepted yet!")

        if maximum_copy_num ** len(free_copy_variables) < 5E6:
            # sometimes, SLSQP ignores bounds and constraints
            copy_results = minimize_brute_force(
                func=least_square_function_v, range_list=[range(1, maximum_copy_num+1)]*len(free_copy_variables),
                constraint_list=({'type': 'ineq', 'fun': constraint_min_function_for_customized_brute},
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
            raise ProcessingGraphFailed("Incomplete/Complicated " + target_name_for_log +" graph (3)!")

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
                    if this_copy <= 0:
                        raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
                    final_results[go_res]["graph"].vertex_to_copy[vertex_name] = this_copy
                    if this_copy not in final_results[go_res]["graph"].copy_to_vertex:
                        final_results[go_res]["graph"].copy_to_vertex[this_copy] = set()
                    final_results[go_res]["graph"].copy_to_vertex[this_copy].add(vertex_name)

                """ re-estimate baseline depth """
                total_product = 0.
                total_len = 0
                for vertex_name in vertices_list:
                    this_len = (self.vertex_info[vertex_name].len - self.__kmer + 1) \
                               * final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
                    this_cov = self.vertex_info[vertex_name].cov \
                               / final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
                    total_len += this_len
                    total_product += this_len * this_cov
                final_results[go_res]["cov"] = total_product / total_len
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
                if this_copy <= 0:
                    raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)

            if debug or verbose:
                """ re-estimate baseline depth """
                total_product = 0.
                total_len = 0
                for vertex_name in vertices_list:
                    this_len = (self.vertex_info[vertex_name].len - self.__kmer + 1) \
                               * self.vertex_to_copy.get(vertex_name, 1)
                    this_cov = self.vertex_info[vertex_name].cov / self.vertex_to_copy.get(vertex_name, 1)
                    total_len += this_len
                    total_product += this_len * this_cov
                new_val = total_product / total_len
                if log_handler:
                    log_handler.info("Average " + target_name_for_log +" kmer-coverage = " + str(round(new_val, 2)))
                else:
                    sys.stdout.write("Average " + target_name_for_log +" kmer-coverage = " + str(round(new_val, 2)) + "\n")

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
                    if sum([bool(c_c) for c_c in self.vertex_info[can_v].connections.values()]) != 2:
                        del candidate_vertices[go_to_v]
                        continue
                    count_nearby_tagged = []
                    for can_end, can_connect in self.vertex_info[can_v].connections.items():
                        for next_v, next_e in can_connect:
                            # candidate_v is the only output vertex to next_v
                            if next_v in self.tagged_vertices[mode] and \
                                    len(self.vertex_info[next_v].connections[next_e]) == 1:
                                count_nearby_tagged.append((next_v, next_e))
                                break
                    if len(count_nearby_tagged) == 2:
                        del candidate_vertices[go_to_v]
                        # add in between
                        self.tagged_vertices[mode].add(can_v)
                        if "weight" not in self.vertex_info[can_v].other_attr:
                            self.vertex_info[can_v].other_attr["weight"] = {}
                        if mode not in self.vertex_info[can_v].other_attr["weight"]:
                            self.vertex_info[can_v].other_attr["weight"][mode] = 0
                        self.vertex_info[can_v].other_attr["weight"][mode] += 1 * self.vertex_info[can_v].cov
                        # add extra circle
                        near_by_pairs = self.is_sequential_repeat(can_v, return_pair_in_the_trunk_path=False)
                        if near_by_pairs:
                            checking_new = []
                            coverage_folds = []
                            for near_by_p in near_by_pairs:
                                for (near_v, near_e) in near_by_p:
                                    if (near_v, near_e) not in count_nearby_tagged:
                                        checking_new.append(near_v)
                                        coverage_folds.append(
                                            round(self.vertex_info[can_v].cov /
                                                  self.vertex_info[near_v].cov, 0))
                            for near_v, near_e in count_nearby_tagged:
                                coverage_folds.append(
                                    round(self.vertex_info[can_v].cov /
                                          self.vertex_info[near_v].cov, 0))
                            if max(coverage_folds) >= 2:
                                for extra_v_to_add in set(checking_new):
                                    self.tagged_vertices[mode].add(extra_v_to_add)
                                    try:
                                        candidate_vertices.remove(extra_v_to_add)
                                    except ValueError:
                                        pass
                                    if "weight" not in self.vertex_info[extra_v_to_add].other_attr:
                                        self.vertex_info[extra_v_to_add].other_attr["weight"] = {mode: 0}
                                    self.vertex_info[extra_v_to_add].other_attr["weight"][mode] \
                                        += 1 * self.vertex_info[extra_v_to_add].cov
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
                        if (locus_start == 1 or locus_end == self.vertex_info[vertex_name].len) \
                                and locus_len == self.__kmer:
                            continue
                        if locus_name in tag_loci[locus_type]:
                            new_weight = locus_len * self.vertex_info[vertex_name].cov
                            if new_weight > tag_loci[locus_type][locus_name]["weight"]:
                                tag_loci[locus_type][locus_name] = {"vertex": vertex_name, "len": locus_len,
                                                                    "weight": new_weight}
                        else:
                            tag_loci[locus_type][locus_name] = {"vertex": vertex_name, "len": locus_len,
                                                                "weight": locus_len * self.vertex_info[vertex_name].cov}

        for locus_type in tag_loci:
            self.tagged_vertices[locus_type] = set()
            for locus_name in tag_loci[locus_type]:
                vertex_name = tag_loci[locus_type][locus_name]["vertex"]
                loci_weight = tag_loci[locus_type][locus_name]["weight"]
                # tags
                if "tags" not in self.vertex_info[vertex_name].other_attr:
                    self.vertex_info[vertex_name].other_attr["tags"] = {}
                if locus_type in self.vertex_info[vertex_name].other_attr["tags"]:
                    self.vertex_info[vertex_name].other_attr["tags"][locus_type].add(locus_name)
                else:
                    self.vertex_info[vertex_name].other_attr["tags"][locus_type] = {locus_name}
                # weight
                if "weight" not in self.vertex_info[vertex_name].other_attr:
                    self.vertex_info[vertex_name].other_attr["weight"] = {}
                if locus_type in self.vertex_info[vertex_name].other_attr["weight"]:
                    self.vertex_info[vertex_name].other_attr["weight"][locus_type] += loci_weight
                else:
                    self.vertex_info[vertex_name].other_attr["weight"][locus_type] = loci_weight
                self.tagged_vertices[locus_type].add(vertex_name)

        for vertex_name in self.vertex_info:
            if "weight" in self.vertex_info[vertex_name].other_attr:
                if len(self.vertex_info[vertex_name].other_attr["weight"]) > 1:
                    all_weights = sorted([(loc_type, self.vertex_info[vertex_name].other_attr["weight"][loc_type])
                                          for loc_type in self.vertex_info[vertex_name].other_attr["weight"]], key=lambda x: -x[1])
                    best_t, best_w = all_weights[0]
                    for next_t, next_w in all_weights[1:]:
                        if next_w * type_factor < best_w:
                            self.tagged_vertices[next_t].remove(vertex_name)

        if len(self.tagged_vertices[mode]) == 0:
            raise Exception("No available " + mode + " information found in " + tab_file)

    def filter_by_coverage(self, drop_num=1, mode="embplant_pt", log_hard_cov_threshold=10.,
                           weight_factor=100., min_sigma_factor=0.1, min_cluster=1,
                           verbose=False, log_handler=None, debug=False):
        changed = False
        log_hard_cov_threshold = abs(log(log_hard_cov_threshold))
        vertices = sorted(self.vertex_info)
        v_coverages = {this_v: self.vertex_info[this_v].cov / self.vertex_to_copy.get(this_v, 1)
                       for this_v in vertices}
        max_tagged_cov = max([v_coverages[tagged_v] for tagged_v in self.tagged_vertices[mode]])
        # removing coverage with 10 times lower/greater than tagged_cov
        removing_low_cov = [candidate_v
                            for candidate_v in vertices
                            if abs(log(self.vertex_info[candidate_v].cov/max_tagged_cov)) > log_hard_cov_threshold]
        if removing_low_cov:
            if log_handler and (debug or verbose):
                log_handler.info("removing extremely outlying coverage contigs: " + str(removing_low_cov))
            elif verbose or debug:
                sys.stdout.write("removing extremely outlying coverage contigs: " + str(removing_low_cov) + "\n")
            self.remove_vertex(removing_low_cov)
            changed = True
        merged = self.merge_all_possible_vertices()
        if merged:
            changed = True
        vertices = sorted(self.vertex_info)
        v_coverages = {this_v: self.vertex_info[this_v].cov / self.vertex_to_copy.get(this_v, 1)
                       for this_v in vertices}

        coverages = np.array([v_coverages[this_v] for this_v in vertices])
        cover_weights = np.array([(self.vertex_info[this_v].len - self.__kmer) * self.vertex_to_copy.get(this_v, 1)
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
                                              minimum_cluster=min_cluster, maximum_cluster=6,
                                              cluster_limited=set_cluster, min_sigma_factor=min_sigma_factor)
        cluster_num = gmm_scheme["cluster_num"]
        parameters = gmm_scheme["parameters"]
        # for debug
        # print('testing', end="\n")
        # for temp in parameters:
        #     print("  ", temp, end="\n")
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
            label_weights = {lb: sum([self.vertex_info[vertices[go]].other_attr.get("weight", {}).get(mode, 0)
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
        return changed, [(parameters[lab_tp]["mu"], parameters[lab_tp]["sigma"]) for lab_tp in remained_label_type]

    def generate_consensus_vertex(self, vertices, directions, copy_tags=True, check_parallel_vertices=True,
                                  log_handler=None):
        if check_parallel_vertices:
            connection_type = None
            seq_len = None
            if not len(vertices) == len(set(vertices)) == len(directions):
                raise ProcessingGraphFailed("Cannot generate consensus (1)!")
            for go_v, this_v in enumerate(vertices):
                if seq_len:
                    if seq_len != len(self.vertex_info[this_v].seq[True]):
                        raise ProcessingGraphFailed("Cannot generate consensus (2)!")
                else:
                    seq_len = len(self.vertex_info[this_v].seq[True])
                this_cons = self.vertex_info[this_v].connections
                this_ends = tuple([tuple(sorted(this_cons[[directions[go_v]]])),
                                   tuple(sorted(this_cons[not [directions[go_v]]]))])
                if connection_type:
                    if connection_type != this_ends:
                        raise ProcessingGraphFailed("Cannot generate consensus (3)!")
                else:
                    connection_type = this_ends

        if len(vertices) > 1:
            new_vertex = "(" + "|".join(vertices) + ")"
            self.vertex_info[new_vertex] = deepcopy(self.vertex_info[vertices[0]])
            self.vertex_info[new_vertex].cov = sum([self.vertex_info[v].cov for v in vertices])
            self.vertex_info[new_vertex].fastg_form_name = None
            # if "long" in self.vertex_info[new_vertex]:
            #     del self.vertex_info[new_vertex]["long"]

            for new_end in (True, False):
                for n_n_v, n_n_e in self.vertex_info[new_vertex].connections[new_end]:
                    self.vertex_info[n_n_v].connections[n_n_e].add((new_vertex, new_end))

            consensus_s = generate_consensus(*[self.vertex_info[v].seq[directions[go]] for go, v in enumerate(vertices)])
            self.vertex_info[new_vertex].seq[directions[0]] = consensus_s
            self.vertex_info[new_vertex].seq[not directions[0]] = complementary_seq(consensus_s)

            # tags
            if copy_tags:
                for other_vertex in vertices[1:]:
                    if "tags" in self.vertex_info[other_vertex].other_attr:
                        if "tags" not in self.vertex_info[new_vertex].other_attr:
                            self.vertex_info[new_vertex].other_attr["tags"] = \
                                deepcopy(self.vertex_info[other_vertex].other_attr["tags"])
                        else:
                            for mode in self.vertex_info[other_vertex].other_attr["tags"]:
                                if mode not in self.vertex_info[new_vertex].other_attr["tags"]:
                                    self.vertex_info[new_vertex].other_attr["tags"][mode] \
                                        = deepcopy(self.vertex_info[other_vertex].other_attr["tags"][mode])
                                else:
                                    self.vertex_info[new_vertex].other_attr["tags"][mode] \
                                        |= self.vertex_info[other_vertex].other_attr["tags"][mode]
                    if "weight" in self.vertex_info[other_vertex].other_attr:
                        if "weight" not in self.vertex_info[new_vertex].other_attr:
                            self.vertex_info[new_vertex].other_attr["weight"] \
                                = deepcopy(self.vertex_info[other_vertex].other_attr["weight"])
                        else:
                            for mode in self.vertex_info[other_vertex].other_attr["weight"]:
                                if mode not in self.vertex_info[new_vertex].other_attr["weight"]:
                                    self.vertex_info[new_vertex].other_attr["weight"][mode] \
                                        = self.vertex_info[other_vertex].other_attr["weight"][mode]
                                else:
                                    self.vertex_info[new_vertex].other_attr["weight"][mode] \
                                        += self.vertex_info[other_vertex].other_attr["weight"][mode]
                    for mode in self.tagged_vertices:
                        if other_vertex in self.tagged_vertices[mode]:
                            self.tagged_vertices[mode].add(new_vertex)
                            self.tagged_vertices[mode].remove(other_vertex)
            self.remove_vertex(vertices)
            if log_handler:
                log_handler.info("Consensus made: " + new_vertex)
            else:
                log_handler.info("Consensus made: " + new_vertex + "\n")

    def processing_polymorphism(self, mode, limited_vertices=None,
                                contamination_depth=3., contamination_similarity=0.95,
                                degenerate=False, degenerate_depth=1.5, degenerate_similarity=0.98, warning_count=4,
                                only_keep_max_cov=False, verbose=False, debug=False, log_handler=None):
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
        count_using_only_max = 0
        sub_sampling = 10000
        half_kmer = int(self.__kmer/2)
        for prl_vertices in parallel_vertices_list:
            this_contamination_or_polymorphic = False
            this_using_only_max = False
            prl_vertices = sorted(prl_vertices, key=lambda x: -self.vertex_info[x[0]].cov)
            max_cov_vertex, direction_remained = prl_vertices.pop(0)
            max_cov_seq = self.vertex_info[max_cov_vertex].seq[direction_remained]
            max_cov = self.vertex_info[max_cov_vertex].cov
            polymorphic_vertices_with_directions = {(max_cov_vertex, direction_remained)}
            # hard to clearly identify the biological factors
            for this_v, this_direction in prl_vertices:
                this_seq = self.vertex_info[this_v].seq[this_direction]
                this_cov = self.vertex_info[this_v].cov
                if abs(log(this_cov/max_cov)) > contamination_depth:
                    if abs(len(this_seq) - len(max_cov_seq)) / float(len(this_seq)) <= contamination_dif:
                        # too long to calculate, too long to be polymorphic
                        if max(len(max_cov_seq), len(this_seq)) > 1E5:
                            removing_irrelevant_v.add(this_v)
                        else:
                            seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                            seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                            base_dif, proper_end = find_string_difference(
                                seq_1, seq_2, max(2, int(len(seq_1) * 0.005)))
                            if float(base_dif) * 2 / (len(seq_1) + len(seq_2)) < contamination_dif:
                                removing_contaminating_v.add(this_v)
                                this_contamination_or_polymorphic = True
                            else:
                                removing_irrelevant_v.add(this_v)
                    else:
                        removing_irrelevant_v.add(this_v)
                elif degenerate and abs(log(this_cov/max_cov)) < degenerate_depth:
                    if abs(len(this_seq) - len(max_cov_seq)) / float(len(this_seq)) <= degenerate_dif:
                        # too long to calculate, too long to be polymorphic
                        if max(len(max_cov_seq), len(this_seq)) > 1E5:
                            continue
                        else:
                            seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                            seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                            base_dif, proper_end = find_string_difference(
                                seq_1, seq_2, max(2, int(len(seq_1) * 0.005)))
                            if float(base_dif) * 2 / (len(seq_1) + len(seq_2)) < degenerate_dif:
                                this_contamination_or_polymorphic = True
                                if len(this_seq) == len(max_cov_seq):
                                    polymorphic_vertices_with_directions.add((this_v, this_direction))
                                elif only_keep_max_cov:
                                    removing_irrelevant_v.add(this_v)
                                    this_using_only_max = True
                                else:
                                    raise ProcessingGraphFailed(
                                        "Cannot degenerate inequal-length polymorphic contigs: EDGE_" +
                                        max_cov_vertex + " and EDGE_" + this_v + "!")
                            elif only_keep_max_cov and float(base_dif)*2/(len(seq_1) + len(seq_2)) < contamination_dif:
                                removing_irrelevant_v.add(this_v)
                                this_contamination_or_polymorphic = True
                                this_using_only_max = True
                    elif only_keep_max_cov:
                        seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                        seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                        base_dif, proper_end = find_string_difference(
                            seq_1, seq_2, max(2, int(len(seq_1) * 0.005)))
                        if float(base_dif) * 2 / (len(seq_1) + len(seq_2)) < contamination_dif:
                            removing_irrelevant_v.add(this_v)
                            this_contamination_or_polymorphic = True
                            this_using_only_max = True
                elif only_keep_max_cov:
                    seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                    seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                    base_dif, proper_end = find_string_difference(
                        seq_1, seq_2, max(2, int(len(seq_1) * 0.005)))
                    if float(base_dif) * 2 / (len(seq_1) + len(seq_2)) < contamination_dif:
                        removing_irrelevant_v.add(this_v)
                        this_contamination_or_polymorphic = True
                        this_using_only_max = True
            if len(polymorphic_vertices_with_directions) > 1:
                self.generate_consensus_vertex([v_d[0] for v_d in polymorphic_vertices_with_directions],
                                               [v_d[1] for v_d in polymorphic_vertices_with_directions],
                                               check_parallel_vertices=False, log_handler=log_handler)
            count_contamination_or_degenerate += this_contamination_or_polymorphic
            count_using_only_max += this_using_only_max
        if removing_contaminating_v:
            contaminating_cov = np.array([self.vertex_info[con_v].cov for con_v in removing_contaminating_v])
            contaminating_weight = np.array([len(self.vertex_info[con_v].seq[True]) - self.__kmer
                                             for con_v in removing_contaminating_v])
            self.remove_vertex(removing_contaminating_v)
            cont_mean, cont_std = weighted_mean_and_std(contaminating_cov, contaminating_weight)
            cut_off_min, cut_off_max = cont_mean - cont_std, cont_mean + cont_std
            # remove (low cov and terminal vertex)
            removing_below_cut_off = []
            for del_v in self.vertex_info:
                if cut_off_min < self.vertex_info[del_v].cov < cut_off_max:
                    if sum([bool(cnn) for cnn in self.vertex_info[del_v].connections.values()]) < 2 \
                            and del_v not in self.tagged_vertices[mode]:
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
        if count_contamination_or_degenerate >= warning_count:
            if log_handler:
                log_handler.warning("The graph might suffer from contamination or polymorphism!")
                if count_using_only_max:
                    log_handler.warning("Only the contig with the max cov was kept for each of those " +
                                        str(count_using_only_max) + " polymorphic loci.")
            else:
                sys.stdout.write("Warning: The graph might suffer from contamination or polymorphism!")
                if count_using_only_max:
                    sys.stdout.write("Warning: Only the contig with the max cov was kept for each of those " +
                                     str(count_using_only_max) + " polymorphic loci.\n")

    def find_target_graph(self, tab_file, mode="embplant_pt", type_factor=3, weight_factor=100.0,
                          max_copy=8, min_sigma_factor=0.1, expected_max_size=inf, expected_min_size=0,
                          log_hard_cov_threshold=10., contamination_depth=3., contamination_similarity=0.95,
                          degenerate=True, degenerate_depth=1.5, degenerate_similarity=0.98, only_keep_max_cov=True,
                          broken_graph_allowed=False, temp_graph=None, verbose=True,
                          read_len_for_log=None, kmer_for_log=None,
                          log_handler=None, debug=False):

        def log_target_res(final_res_combinations_inside):
            echo_graph_id = int(bool(len(final_res_combinations_inside) - 1))
            for go_res, final_res_one in enumerate(final_res_combinations_inside):
                this_graph = final_res_combinations_inside[go_res]["graph"]
                this_k_cov = round(final_res_combinations_inside[go_res]["cov"], 3)
                if read_len_for_log and kmer_for_log:
                    this_b_cov = round(this_k_cov * read_len_for_log / (read_len_for_log - kmer_for_log + 1), 3)
                else:
                    this_b_cov = None
                if log_handler:
                    if echo_graph_id:
                        log_handler.info("Graph " + str(go_res + 1))
                    for vertex_set in sorted(this_graph.vertex_clusters):
                        copies_in_a_set = {this_graph.vertex_to_copy[v_name] for v_name in vertex_set}
                        if copies_in_a_set != {1}:
                            for vertex_name in sorted(vertex_set):
                                log_handler.info("Vertex_" + vertex_name + " #copy = " +
                                                 str(this_graph.vertex_to_copy.get(vertex_name, 1)))

                    log_handler.info("Average " + mode + " kmer-coverage" +
                                     ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + str(this_k_cov))
                    if this_b_cov:
                        log_handler.info("Average " + mode + " base-coverage" +
                                         ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + str(this_b_cov))
                else:
                    if echo_graph_id:
                        sys.stdout.write("Graph " + str(go_res + 1) + "\n")
                    for vertex_set in sorted(this_graph.vertex_clusters):
                        copies_in_a_set = {this_graph.vertex_to_copy[v_name] for v_name in vertex_set}
                        if copies_in_a_set != {1}:
                            for vertex_name in sorted(vertex_set):
                                sys.stdout.write("Vertex_" + vertex_name + " #copy = " +
                                                 str(this_graph.vertex_to_copy.get(vertex_name, 1)) + "\n")
                    sys.stdout.write("Average " + mode + " kmer-coverage" +
                                     ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + str(this_k_cov) + "\n")
                    if this_b_cov:
                        sys.stdout.write("Average " + mode + " base-coverage" +
                                         ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + str(this_b_cov) + "\n")

        if broken_graph_allowed:
            weight_factor = 10000.

        self.parse_tab_file(tab_file, mode=mode, type_factor=type_factor, log_handler=log_handler)
        new_assembly = deepcopy(self)
        is_reasonable_res = False
        data_contains_outlier = False
        try:
            while not is_reasonable_res:
                is_reasonable_res = True
                # if verbose or debug:
                #     if log_handler:
                #         log_handler.info("tagged vertices: " + str(sorted(new_assembly.tagged_vertices[mode])))
                #         log_handler.info("tagged coverage: " +
                #                          str(["%.1f"%new_assembly.vertex_info[log_v].cov
                #                               for log_v in sorted(new_assembly.tagged_vertices[mode])]))
                #     else:
                #         sys.stdout.write("tagged vertices: " + str(sorted(new_assembly.tagged_vertices[mode])) + "\n")
                #         log_handler.info("tagged coverage: " +
                #                          str(["%.1f"%new_assembly.vertex_info[log_v].cov
                #                               for log_v in sorted(new_assembly.tagged_vertices[mode])]) + "\n")
                new_assembly.merge_all_possible_vertices()
                new_assembly.tag_in_between(mode=mode)
                new_assembly.processing_polymorphism(mode=mode, contamination_depth=contamination_depth,
                                                     contamination_similarity=contamination_similarity,
                                                     degenerate=False, verbose=verbose, debug=debug,
                                                     log_handler=log_handler)
                changed = True
                count_large_round = 0
                while changed:
                    count_large_round += 1
                    if verbose or debug:
                        if log_handler:
                            log_handler.info("===================== " + str(count_large_round) + " =====================")
                        else:
                            sys.stdout.write("===================== " + str(count_large_round) + " =====================\n")
                    changed = False
                    cluster_trimmed = True
                    while cluster_trimmed:
                        # remove low coverages
                        first_round = True
                        delete_those_vertices = set()
                        parameters = []
                        this_del = False
                        new_assembly.estimate_copy_and_depth_by_cov(new_assembly.tagged_vertices[mode], debug=debug,
                                                                    log_handler=log_handler, verbose=verbose, mode=mode)
                        while first_round or delete_those_vertices or this_del:
                            if data_contains_outlier:
                                this_del, parameters =\
                                    new_assembly.filter_by_coverage(mode=mode, weight_factor=weight_factor,
                                                                    log_hard_cov_threshold=log_hard_cov_threshold,
                                                                    min_sigma_factor=min_sigma_factor,
                                                                    min_cluster=2, log_handler=log_handler,
                                                                    verbose=verbose, debug=debug)
                                data_contains_outlier = False
                            else:
                                this_del, parameters = \
                                    new_assembly.filter_by_coverage(mode=mode, weight_factor=weight_factor,
                                                                    log_hard_cov_threshold=log_hard_cov_threshold,
                                                                    min_sigma_factor=min_sigma_factor,
                                                                    log_handler=log_handler, verbose=verbose,
                                                                    debug=debug)
                            # if verbose or debug:
                            #     if log_handler:
                            #         log_handler.info("tagged vertices: " + str(sorted(new_assembly.tagged_vertices[mode])))
                            #         log_handler.info("tagged coverage: " +
                            #                          str(["%.1f"%new_assembly.vertex_info[log_v].cov
                            #                               for log_v in sorted(new_assembly.tagged_vertices[mode])]))
                            #     else:
                            #         sys.stdout.write("tagged vertices: " + str(sorted(new_assembly.tagged_vertices[mode])) + "\n")
                            #         log_handler.info("tagged coverage: " +
                            #                          str(["%.1f"%new_assembly.vertex_info[log_v].cov
                            #                               for log_v in sorted(new_assembly.tagged_vertices[mode])]) + "\n")
                            new_assembly.estimate_copy_and_depth_by_cov(new_assembly.tagged_vertices[mode], debug=debug,
                                                                        log_handler=log_handler, verbose=verbose, mode=mode)
                            first_round = False

                        cluster_trimmed = False

                        if len(new_assembly.vertex_clusters) == 0:
                            raise ProcessingGraphFailed("No available " + mode + " components detected!")
                        elif len(new_assembly.vertex_clusters) == 1:
                            pass
                        else:
                            cluster_weights = [sum([new_assembly.vertex_info[x_v].other_attr["weight"][mode]
                                                    for x_v in x
                                                    if "weight" in new_assembly.vertex_info[x_v].other_attr
                                                    and mode in new_assembly.vertex_info[x_v].other_attr["weight"]])
                                               for x in new_assembly.vertex_clusters]
                            best = max(cluster_weights)
                            best_id = cluster_weights.index(best)
                            if broken_graph_allowed:
                                id_remained = {best_id}
                                for j, w in enumerate(cluster_weights):
                                    if w * weight_factor > best:
                                        id_remained.add(j)
                                    else:
                                        for del_v in new_assembly.vertex_clusters[j]:
                                            if del_v in new_assembly.tagged_vertices[mode]:
                                                new_cov = new_assembly.vertex_info[del_v].cov
                                                for mu, sigma in parameters:
                                                    if abs(new_cov - mu) < sigma:
                                                        id_remained.add(j)
                                                        break
                                            if j in id_remained:
                                                break
                            else:
                                # chose the target cluster (best rank)
                                id_remained = {best_id}
                                temp_cluster_weights = deepcopy(cluster_weights)
                                del temp_cluster_weights[best_id]
                                second = max(temp_cluster_weights)
                                if best < second * weight_factor:
                                    if temp_graph:
                                        new_assembly.write_to_gfa(temp_graph)
                                        new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                                    raise ProcessingGraphFailed("Multiple isolated " + mode + " components detected! "
                                                                "Broken or contamination?")
                                for j, w in enumerate(cluster_weights):
                                    if w == second:
                                        for del_v in new_assembly.vertex_clusters[j]:
                                            if del_v in new_assembly.tagged_vertices[mode]:
                                                new_cov = new_assembly.vertex_info[del_v].cov
                                                # for debug
                                                # print(new_cov)
                                                # print(parameters)
                                                for mu, sigma in parameters:
                                                    if abs(new_cov - mu) < sigma:
                                                        if temp_graph:
                                                            new_assembly.write_to_gfa(temp_graph)
                                                            new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                                                        raise ProcessingGraphFailed(
                                                            "Complicated graph: please check around EDGE_" + del_v + "!"
                                                            "# tags: " +
                                                            str(new_assembly.vertex_info[del_v].other_attr["tags"][mode]))

                            # remove other clusters
                            vertices_to_del = set()
                            for go_cl, v_2_del in enumerate(new_assembly.vertex_clusters):
                                if go_cl not in id_remained:
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

                    # no tip contigs allowed
                    if broken_graph_allowed:
                        pass
                    else:
                        first_round = True
                        delete_those_vertices = set()
                        while first_round or delete_those_vertices:
                            first_round = False
                            delete_those_vertices = set()
                            for vertex_name in new_assembly.vertex_info:
                                # both ends must have edge(s)
                                if sum([bool(len(cn))
                                        for cn in new_assembly.vertex_info[vertex_name].connections.values()]) != 2:
                                    if vertex_name in new_assembly.tagged_vertices[mode]:
                                        if temp_graph:
                                            new_assembly.write_to_gfa(temp_graph)
                                            new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                                        raise ProcessingGraphFailed(
                                            "Incomplete/Complicated graph: please check around EDGE_"+vertex_name + "!")
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

                    # # merge vertices
                    # new_assembly.merge_all_possible_vertices()
                    # new_assembly.tag_in_between(mode=mode)
                    # break self-connection if necessary
                    # for vertex_name in new_assembly.vertex_info:
                    #     if (vertex_name, True) in
                    # -> not finished!!

                    # merge vertices
                    new_assembly.merge_all_possible_vertices()
                    new_assembly.processing_polymorphism(mode=mode, contamination_depth=contamination_depth,
                                                         contamination_similarity=contamination_similarity,
                                                         degenerate=False, degenerate_depth=degenerate_depth,
                                                         degenerate_similarity=degenerate_similarity,
                                                         verbose=verbose, debug=debug, log_handler=log_handler)
                    new_assembly.tag_in_between(mode=mode)

                if temp_graph:
                    new_assembly.write_to_gfa(temp_graph)
                    new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                new_assembly.processing_polymorphism(mode=mode, contamination_depth=contamination_depth,
                                                     contamination_similarity=contamination_similarity,
                                                     degenerate=degenerate, degenerate_depth=degenerate_depth,
                                                     degenerate_similarity=degenerate_similarity,
                                                     warning_count=1, only_keep_max_cov=only_keep_max_cov,
                                                     verbose=verbose, debug=debug, log_handler=log_handler)
                new_assembly.merge_all_possible_vertices()
                if temp_graph:
                    new_assembly.write_to_gfa(temp_graph)
                    new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")

                # create idealized vertices and edges
                try:
                    new_average_cov = new_assembly.estimate_copy_and_depth_by_cov(log_handler=log_handler,
                                                                                  verbose=verbose,
                                                                                  mode="all", debug=debug)
                    if verbose:
                        if log_handler:
                            log_handler.info("Estimating copy and depth precisely ...")
                        else:
                            sys.stdout.write("Estimating copy and depth precisely ...\n")
                    final_res_combinations = new_assembly.estimate_copy_and_depth_precisely(
                        maximum_copy_num=max_copy, broken_graph_allowed=broken_graph_allowed, log_handler=log_handler,
                        verbose=verbose, debug=debug)
                    if verbose:
                        if log_handler:
                            log_handler.info(str(len(final_res_combinations)) + " candidate graph(s) generated.")
                        else:
                            sys.stdout.write(str(len(final_res_combinations)) + " candidate graph(s) generated.\n")
                    absurd_copy_nums = True
                    no_single_copy = True
                    while absurd_copy_nums:
                        go_graph = 0
                        while go_graph < len(final_res_combinations):
                            this_assembly_g = final_res_combinations[go_graph]["graph"]
                            this_parallel_v_sets = [v_set for v_set in this_assembly_g.detect_parallel_vertices()]
                            this_parallel_names = set([v_n for v_set in this_parallel_v_sets for v_n, v_e in v_set])
                            if 1 not in this_assembly_g.copy_to_vertex:
                                if verbose or debug:
                                    if log_handler:
                                        for vertex_name in sorted(this_assembly_g.vertex_info):
                                            log_handler.info(
                                                "Vertex_" + vertex_name + " #copy = " +
                                                str(this_assembly_g.vertex_to_copy.get(vertex_name, 1)))
                                        log_handler.info("Removing this graph without single copy contigs.")
                                    else:
                                        for vertex_name in sorted(this_assembly_g.vertex_info):
                                            sys.stdout.write(
                                                "Vertex_" + vertex_name + " #copy = " +
                                                str(this_assembly_g.vertex_to_copy.get(vertex_name, 1)) + "\n")
                                        sys.stdout.write("Removing this graph without single copy contigs.\n")
                                del final_res_combinations[go_graph]
                            else:
                                no_single_copy = False
                                this_absurd = True
                                for single_copy_v in this_assembly_g.copy_to_vertex[1]:
                                    if single_copy_v not in this_parallel_names:
                                        this_absurd = False
                                draft_size_estimates = 0
                                for inside_v in this_assembly_g.vertex_info:
                                    draft_size_estimates += \
                                        (this_assembly_g.vertex_info[inside_v].len - this_assembly_g.kmer()) * \
                                        this_assembly_g.vertex_to_copy[inside_v]
                                if not this_absurd or expected_min_size < draft_size_estimates < expected_max_size:
                                    absurd_copy_nums = False
                                    go_graph += 1
                                else:
                                    if verbose or debug:
                                        if log_handler:
                                            log_handler.info(
                                                "Removing graph with draft size: " + str(draft_size_estimates))
                                        else:
                                            sys.stdout.write(
                                                "Removing graph with draft size: " + str(draft_size_estimates) + "\n")
                                    # add all combinations
                                    for index_set in generate_index_combinations([len(v_set)
                                                                                  for v_set in this_parallel_v_sets]):
                                        new_possible_graph = deepcopy(this_assembly_g)
                                        dropping_names = []
                                        for go_set, this_v_set in enumerate(this_parallel_v_sets):
                                            keep_this = index_set[go_set]
                                            for go_ve, (this_name, this_end) in enumerate(this_v_set):
                                                if go_ve != keep_this:
                                                    dropping_names.append(this_name)
                                        # if log_handler:
                                        #     log_handler.info("Dropping vertices " + " ".join(dropping_names))
                                        # else:
                                        #     log_handler.info("Dropping vertices " + "".join(dropping_names) + "\n")
                                        new_possible_graph.remove_vertex(dropping_names)
                                        new_possible_graph.merge_all_possible_vertices()
                                        new_possible_graph.estimate_copy_and_depth_by_cov(
                                            log_handler=log_handler, verbose=verbose, mode="all", debug=debug)
                                        final_res_combinations.extend(
                                            new_possible_graph.estimate_copy_and_depth_precisely(
                                                maximum_copy_num=max_copy, broken_graph_allowed=broken_graph_allowed,
                                                log_handler=log_handler, verbose=verbose, debug=debug))

                                    del final_res_combinations[go_graph]
                        if not final_res_combinations and absurd_copy_nums:
                            # if absurd_copy_nums:
                            #     raise ProcessingGraphFailed("Complicated graph! Detecting path(s) failed!")
                            # else:
                            raise ProcessingGraphFailed("Complicated graph! Detecting path(s) failed!")
                    if no_single_copy:
                        raise ProcessingGraphFailed("No single copy region?! Detecting path(s) failed!")
                except ImportError as e:
                    raise e
                except (RecursionError, Exception) as e:
                    if broken_graph_allowed:
                        unlabelled_contigs = [check_v for check_v in list(new_assembly.vertex_info)
                                              if check_v not in new_assembly.tagged_vertices[mode]]
                        if unlabelled_contigs:
                            if verbose or debug:
                                if log_handler:
                                    log_handler.info("removing unlabelled contigs: " + str(unlabelled_contigs))
                                else:
                                    sys.stdout.write("removing unlabelled contigs: " + str(unlabelled_contigs) + "\n")
                            new_assembly.remove_vertex(unlabelled_contigs)
                            new_assembly.merge_all_possible_vertices()
                        else:
                            # delete all previous connections if all present contigs are labelled
                            for del_v_connection in new_assembly.vertex_info:
                                new_assembly.vertex_info[del_v_connection].connections = {True: set(),
                                                                                          False: set()}
                            new_assembly.update_vertex_clusters()
                        new_average_cov = new_assembly.estimate_copy_and_depth_by_cov(
                            re_initialize=True, log_handler=log_handler, verbose=verbose, mode="all", debug=debug)
                        outer_continue = False
                        for remove_all_connections in (False, True):
                            if remove_all_connections:  # delete all previous connections
                                for del_v_connection in new_assembly.vertex_info:
                                    new_assembly.vertex_info[del_v_connection].connections = {True: set(),
                                                                                              False: set()}
                            new_assembly.update_vertex_clusters()
                            try:
                                here_max_copy = 1 if remove_all_connections else max_copy
                                final_res_combinations = new_assembly.estimate_copy_and_depth_precisely(
                                    maximum_copy_num=here_max_copy, broken_graph_allowed=True, log_handler=log_handler,
                                    verbose=verbose, debug=debug)
                            except ImportError as e:
                                raise e
                            except Exception as e:
                                if verbose or debug:
                                    if log_handler:
                                        log_handler.info(str(e))
                                    else:
                                        sys.stdout.write(str(e) + "\n")
                                continue
                            test_first_g = final_res_combinations[0]["graph"]
                            if 1 in test_first_g.copy_to_vertex:
                                single_copy_percent = sum([test_first_g.vertex_info[s_v].len
                                                           for s_v in test_first_g.copy_to_vertex[1]]) \
                                                      / float(sum([test_first_g.vertex_info[a_v].len
                                                                   for a_v in test_first_g.vertex_info]))
                                if single_copy_percent < 0.5:
                                    if verbose:
                                        if log_handler:
                                            log_handler.warning(
                                                "Result with single copy vertex percentage < 50% is "
                                                "unacceptable, continue dropping suspicious vertices ...")
                                        else:
                                            sys.stdout.write(
                                                "Warning: Result with single copy vertex percentage < 50% is "
                                                "unacceptable, continue dropping suspicious vertices ...")
                                    data_contains_outlier = True
                                    is_reasonable_res = False
                                    outer_continue = True
                                    break
                                else:
                                    log_target_res(final_res_combinations)
                                    return final_res_combinations
                            else:
                                if verbose:
                                    if log_handler:
                                        log_handler.warning("Result with single copy vertex percentage < 50% is "
                                                            "unacceptable, continue dropping suspicious vertices ...")
                                    else:
                                        sys.stdout.write("Warning: Result with single copy vertex percentage < 50% is "
                                                         "unacceptable, continue dropping suspicious vertices ...")
                                data_contains_outlier = True
                                is_reasonable_res = False
                                outer_continue = True
                                break
                        if outer_continue:
                            continue
                    elif temp_graph:
                        new_assembly.write_to_gfa(temp_graph)
                        new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
                        raise ProcessingGraphFailed("Complicated " + mode + " graph! Detecting path(s) failed!")
                    else:
                        raise e
                else:
                    test_first_g = final_res_combinations[0]["graph"]
                    if 1 in test_first_g.copy_to_vertex:
                        single_copy_percent = sum([test_first_g.vertex_info[s_v].len
                                                   for s_v in test_first_g.copy_to_vertex[1]]) \
                                              / float(sum([test_first_g.vertex_info[a_v].len
                                                           for a_v in test_first_g.vertex_info]))
                        if single_copy_percent < 0.5:
                            if verbose:
                                if log_handler:
                                    log_handler.warning("Result with single copy vertex percentage < 50% is "
                                                        "unacceptable, continue dropping suspicious vertices ...")
                                else:
                                    sys.stdout.write("Warning: Result with single copy vertex percentage < 50% is "
                                                     "unacceptable, continue dropping suspicious vertices ...")
                            data_contains_outlier = True
                            is_reasonable_res = False
                            continue
                        else:
                            log_target_res(final_res_combinations)
                            return final_res_combinations
                    else:
                        if verbose:
                            if log_handler:
                                log_handler.warning("Result with single copy vertex percentage < 50% is "
                                                    "unacceptable, continue dropping suspicious vertices ...")
                            else:
                                sys.stdout.write("Warning: Result with single copy vertex percentage < 50% is "
                                                 "unacceptable, continue dropping suspicious vertices ...")
                        data_contains_outlier = True
                        is_reasonable_res = False
                        continue
        except KeyboardInterrupt as e:
            if temp_graph:
                new_assembly.write_to_gfa(temp_graph)
                new_assembly.write_out_tags([mode], temp_graph[:-5] + "csv")
            raise KeyboardInterrupt

    def get_all_circular_paths(self, mode="embplant_pt", library_info=None, log_handler=None):

        def circular_directed_graph_solver(ongoing_path, next_connections, vertices_left, check_all_kinds):
            # print("-----------------------------")
            # print("ongoing_path", ongoing_path)
            # print("next_connect", next_connections)
            # print("vertices_lef", vertices_left)
            if not vertices_left:
                new_path = deepcopy(ongoing_path)
                if check_all_kinds:
                    rev_path = [(this_v, not this_e) for this_v, this_e in new_path[::-1]]
                    this_path_derived = [new_path, rev_path]
                    for change_start in range(1, len(new_path)):
                        this_path_derived.append(new_path[change_start:] + new_path[:change_start])
                        this_path_derived.append(rev_path[change_start:] + rev_path[:change_start])
                    standardized_path = tuple(sorted(this_path_derived)[0])
                    paths.add(tuple(standardized_path))
                else:
                    paths.add(tuple(new_path))
                return

            for next_vertex, next_end in next_connections:
                # print("next_vertex", next_vertex)
                if next_vertex in vertices_left:
                    new_path = deepcopy(ongoing_path)
                    new_left = deepcopy(vertices_left)
                    new_path.append((next_vertex, not next_end))
                    new_left[next_vertex] -= 1
                    if not new_left[next_vertex]:
                        del new_left[next_vertex]
                    new_connections = self.vertex_info[next_vertex].connections[not next_end]
                    if not new_left:
                        if (start_vertex, False) in new_connections:
                            if check_all_kinds:
                                rev_path = [(this_v, not this_e) for this_v, this_e in new_path[::-1]]
                                this_path_derived = [new_path, rev_path]
                                for change_start in range(1, len(new_path)):
                                    this_path_derived.append(new_path[change_start:] + new_path[:change_start])
                                    this_path_derived.append(rev_path[change_start:] + rev_path[:change_start])
                                standardized_path = tuple(sorted(this_path_derived)[0])
                                paths.add(tuple(standardized_path))
                            else:
                                paths.add(tuple(new_path))
                            return
                        else:
                            return
                    else:
                        circular_directed_graph_solver(new_path, new_connections, new_left, check_all_kinds)

        paths = set()

        # print(self.copy_to_vertex)

        if 1 not in self.copy_to_vertex:
            do_check_all_start_kinds = True
            start_vertex = sorted(self.vertex_info, key=lambda x: -self.vertex_info[x].len)[0]
        else:
            # start from a single copy vertex, no need to check all kinds of start vertex
            do_check_all_start_kinds = False
            start_vertex = sorted(self.copy_to_vertex[1], key=lambda x: -self.vertex_info[x].len)[0]
        # each contig stored format:
        first_path = [(start_vertex, True)]
        first_connections = self.vertex_info[start_vertex].connections[True]
        vertex_to_copy = deepcopy(self.vertex_to_copy)
        vertex_to_copy[start_vertex] -= 1
        if vertex_to_copy[start_vertex] <= 0:
            del vertex_to_copy[start_vertex]
        circular_directed_graph_solver(first_path, first_connections, vertex_to_copy, do_check_all_start_kinds)

        if not paths:
            raise ProcessingGraphFailed("Detecting path(s) from remaining graph failed!")
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
                                            "for plastomes with the canonical IR!")
                    else:
                        sys.stdout.write("Warning: Multiple repeat patterns appeared in your data, "
                                         "a more balanced pattern (always the repeat_pattern1) would be suggested "
                                         "for plastomes with the canonical IR!\n")
                    sorted_paths = [(this_path, ".repeat_pattern" + str(pattern_dict[acc_distance]))
                                    for this_path, acc_distance in sorted_paths]
                else:
                    sorted_paths = [(this_path, "") for this_path in sorted(paths)]
            else:
                sorted_paths = [(this_path, "") for this_path in sorted(paths)]

            if mode == "embplant_pt":
                if len(sorted_paths) > 2 and not (100000 < len(self.export_path(sorted_paths[0][0]).seq) < 200000):
                    if log_handler:
                        log_handler.warning("Multiple circular genome structures with abnormal length produced!")
                        log_handler.warning("Please check the assembly graph and selected graph to confirm.")
                    else:
                        sys.stdout.write("Warning: Multiple circular genome structures with abnormal length produced!\n")
                        sys.stdout.write("Please check the assembly graph and selected graph to confirm.\n")
                elif len(sorted_paths) > 2:
                    if log_handler:
                        log_handler.warning("Multiple circular genome structures produced!")
                        log_handler.warning("Please check the existence of those isomers "
                                            "by using reads mapping (library information) or longer reads.")
                    else:
                        sys.stdout.write("Warning: Multiple circular genome structures produced!\n")
                        sys.stdout.write("Please check the existence of those isomers by "
                                         "using reads mapping (library information) or longer reads.\n")
                elif len(sorted_paths) > 1:
                    if log_handler:
                        log_handler.warning("More than one circular genome structure produced ...")
                        log_handler.warning("Please check the final result to confirm whether they are "
                                            "simply flip-flop configurations!")
                    else:
                        sys.stdout.write("More than one circular genome structure produced ...\n")
                        sys.stdout.write("Please check the final result to confirm whether they are "
                                         "simply flip-flop configurations!\n")
            return sorted_paths

    def get_all_paths(self, mode="embplant_pt", log_handler=None):

        def standardize_paths(raw_paths):
            here_standardized_path = []
            for part_path in raw_paths:
                rev_part = [(this_v, not this_e) for this_v, this_e in part_path[::-1]]
                if (part_path[0][0], not part_path[0][1]) \
                        in self.vertex_info[part_path[-1][0]].connections[part_path[-1][1]]:
                    # circular
                    this_part_derived = [part_path, rev_part]
                    for change_start in range(1, len(part_path)):
                        this_part_derived.append(part_path[change_start:] + part_path[:change_start])
                        this_part_derived.append(rev_part[change_start:] + rev_part[:change_start])
                    standard_part = tuple(sorted(this_part_derived, key=lambda x: smart_trans_for_sort(x))[0])
                else:
                    standard_part = tuple(sorted([part_path, rev_part], key=lambda x: smart_trans_for_sort(x))[0])
                here_standardized_path.append(standard_part)
            return tuple(sorted(here_standardized_path, key=lambda x: smart_trans_for_sort(x)))

        def directed_graph_solver(ongoing_paths, next_connections, vertices_left, in_all_start_ve):
            # print("-----------------------------")
            # print("ongoing_path", ongoing_path)
            # print("next_connect", next_connections)
            # print("vertices_lef", vertices_left)
            # print("vertices_lef", len(vertices_left))
            if not vertices_left:
                new_paths = deepcopy(ongoing_paths)
                path_paris.append([new_paths, standardize_paths(new_paths)])
                return

            find_next = False
            for next_vertex, next_end in next_connections:
                # print("next_vertex", next_vertex, next_end)
                if next_vertex in vertices_left:
                    find_next = True
                    new_paths = deepcopy(ongoing_paths)
                    new_left = deepcopy(vertices_left)
                    new_paths[-1].append((next_vertex, not next_end))
                    new_left[next_vertex] -= 1
                    if not new_left[next_vertex]:
                        del new_left[next_vertex]
                    new_connections = self.vertex_info[next_vertex].connections[not next_end]
                    if not new_left:
                        path_paris.append([new_paths, standardize_paths(new_paths)])
                        return
                    else:
                        directed_graph_solver(new_paths, new_connections, new_left, in_all_start_ve)
            if not find_next:
                new_all_start_ve = deepcopy(in_all_start_ve)
                while new_all_start_ve:
                    new_start_vertex, new_start_end = new_all_start_ve.pop(0)
                    if new_start_vertex in vertices_left:
                        new_paths = deepcopy(ongoing_paths)
                        new_left = deepcopy(vertices_left)
                        new_paths.append([(new_start_vertex, new_start_end)])
                        new_left[new_start_vertex] -= 1
                        if not new_left[new_start_vertex]:
                            del new_left[new_start_vertex]
                        new_connections = self.vertex_info[new_start_vertex].connections[new_start_end]
                        if not new_left:
                            path_paris.append([new_paths, standardize_paths(new_paths)])
                            return
                        else:
                            directed_graph_solver(new_paths, new_connections, new_left, new_all_start_ve)
                            break
                if not new_all_start_ve:
                    return

        path_paris = list()
        # start from a terminal vertex in an open graph/subgraph
        #         or a single copy vertex in a closed graph/subgraph
        all_start_v_e = []
        start_vertices = set()
        for go_set, v_set in enumerate(self.vertex_clusters):
            is_closed = True
            for test_vertex_n in sorted(v_set):
                for test_end in (False, True):
                    if not self.vertex_info[test_vertex_n].connections[test_end]:
                        is_closed = False
                        if test_vertex_n not in start_vertices:
                            all_start_v_e.append((test_vertex_n, not test_end))
                            start_vertices.add(test_vertex_n)
            if is_closed:
                if 1 in self.copy_to_vertex[1] and bool(self.copy_to_vertex[1] & v_set):
                    single_copy_v = sorted(self.copy_to_vertex[1] & v_set, key=lambda x: -self.vertex_info[x].len)[0]
                    all_start_v_e.append((single_copy_v, True))
                else:
                    longest_v = sorted(v_set, key=lambda x: -self.vertex_info[x].len)[0]
                    all_start_v_e.append((longest_v, True))
        all_start_v_e.sort(key=lambda x: (smart_trans_for_sort(x[0]), x[1]))
        # start from a self-loop vertex in an open/closed graph/subgraph
        for go_set, v_set in enumerate(self.vertex_clusters):
            for test_vertex_n in sorted(v_set):
                if self.vertex_info[test_vertex_n].is_self_loop():
                    all_start_v_e.append((test_vertex_n, True))
                    all_start_v_e.append((test_vertex_n, False))

        start_v_e = all_start_v_e.pop(0)
        first_path = [[start_v_e]]
        first_connections = self.vertex_info[start_v_e[0]].connections[start_v_e[1]]
        vertex_to_copy = deepcopy(self.vertex_to_copy)
        vertex_to_copy[start_v_e[0]] -= 1
        if not vertex_to_copy[start_v_e[0]]:
            del vertex_to_copy[start_v_e[0]]
        directed_graph_solver(first_path, first_connections, vertex_to_copy, all_start_v_e)

        standardized_path_unique_set = set([this_path_pair[1] for this_path_pair in path_paris])
        paths = []
        for raw_path, standardized_path in path_paris:
            if standardized_path in standardized_path_unique_set:
                paths.append(raw_path)
                standardized_path_unique_set.remove(standardized_path)

        if not paths:
            raise ProcessingGraphFailed("Detecting path(s) from remaining graph failed!")
        else:
            sorted_paths = []
            # total_len = len(list(set(paths))[0])
            record_pattern = False
            for this_path in paths:
                acc_dist = 0
                for copy_num in self.copy_to_vertex:
                    if copy_num > 2:
                        for vertex_name in self.copy_to_vertex[copy_num]:
                            for this_p_part in this_path:
                                loc_ids = [go_to_id for go_to_id, (v, e) in enumerate(this_p_part) if v == vertex_name]
                                if len(loc_ids) > 1:
                                    record_pattern = True
                                    if (this_p_part[0][0], not this_p_part[0][1]) \
                                            in self.vertex_info[this_p_part[-1][0]].connections[this_p_part[-1][1]]:
                                        # circular
                                        part_len = len(this_p_part)
                                        for id_a, id_b in combinations(loc_ids, 2):
                                            acc_dist += min((id_a - id_b) % part_len, (id_b - id_a) % part_len)
                                    else:
                                        for id_a, id_b in combinations(loc_ids, 2):
                                            acc_dist += id_b - id_a
                sorted_paths.append((this_path, acc_dist))
            if record_pattern:
                sorted_paths.sort(key=lambda x: -x[1])
                pattern_dict = {acc_distance: ad_id + 1
                                for ad_id, acc_distance
                                in enumerate(sorted(set([x[1] for x in sorted_paths]), reverse=True))}
                if len(pattern_dict) > 1:
                    if log_handler:
                        if mode == "embplant_pt":
                            log_handler.warning("Multiple repeat patterns appeared in your data, "
                                                "a more balanced pattern (always the repeat_pattern1) would be "
                                                "suggested for plastomes with inverted repeats!")
                        else:
                            log_handler.warning("Multiple repeat patterns appeared in your data.")
                    else:
                        if mode == "embplant_pt":
                            sys.stdout.write("Warning: Multiple repeat patterns appeared in your data, "
                                             "a more balanced pattern (always the repeat_pattern1) would be suggested "
                                             "for plastomes with inverted repeats!\n")
                        else:
                            sys.stdout.write("Warning: Multiple repeat patterns appeared in your data.\n")
                    sorted_paths = [(this_path, ".repeat_pattern" + str(pattern_dict[acc_distance]))
                                    for this_path, acc_distance in sorted_paths]
                else:
                    sorted_paths = [(this_path, "") for this_path in sorted(paths)]
            else:
                sorted_paths = [(this_path, "") for this_path in sorted(paths)]

            if mode == "embplant_pt":
                if len(sorted_paths) > 2 and \
                        not (100000 < sum([len(self.export_path(part_p).seq) for part_p in sorted_paths[0][0]]) < 200000):
                    if log_handler:
                        log_handler.warning("Multiple structures (gene order) with abnormal plastome length produced!")
                        log_handler.warning("Please check the assembly graph and selected graph to confirm.")
                    else:
                        sys.stdout.write(
                            "Warning: Multiple structures (gene order) with abnormal plastome length produced!\n")
                        sys.stdout.write("Please check the assembly graph and selected graph to confirm.\n")
                elif len(sorted_paths) > 2:
                    if log_handler:
                        log_handler.warning("Multiple structures (gene order) produced!")
                        log_handler.warning("Please check the existence of those isomers "
                                            "by using reads mapping (library information) or longer reads.")
                    else:
                        sys.stdout.write("Warning: Multiple structures (gene order) produced!\n")
                        sys.stdout.write("Please check the existence of those isomers by "
                                         "using reads mapping (library information) or longer reads.\n")
                elif len(sorted_paths) > 1:
                    if log_handler:
                        log_handler.warning("More than one structure (gene order) produced ...")
                        log_handler.warning("Please check the final result to confirm whether they are "
                                            "simply flip-flop configurations!")
                    else:
                        sys.stdout.write("More than one structure (gene order) produced ...\n")
                        sys.stdout.write("Please check the final result to confirm whether they are "
                                         "simply flip-flop configurations!\n")
            return sorted_paths

    def export_path(self, in_path):
        seq_names = []
        seq_segments = []
        for this_vertex, this_end in in_path:
            seq_segments.append(self.vertex_info[this_vertex].seq[this_end][self.__kmer:])
            seq_names.append(this_vertex + ("-", "+")[this_end])
        # if not circular
        if (in_path[0][0], not in_path[0][1]) not in self.vertex_info[in_path[-1][0]].connections[in_path[-1][1]]:
            seq_segments[0] = self.vertex_info[in_path[0][0]].seq[in_path[0][1]][:self.__kmer] + seq_segments[0]
        else:
            seq_names[-1] += "(circular)"
        return Sequence(",".join(seq_names), "".join(seq_segments))


class ProcessingGraphFailed(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def generate_index_combinations(index_list):
    if not index_list:
        yield []
    else:
        for go_id in range(index_list[0]):
            for next_ids in generate_index_combinations(index_list[1:]):
                yield [go_id] + next_ids


def smart_trans_for_sort(candidate_item):
    if type(candidate_item) in (tuple, list):
        return [smart_trans_for_sort(this_sub) for this_sub in candidate_item]
    elif type(candidate_item) == bool:
        return not candidate_item
    else:
        all_e = candidate_item.split("_")
        for go_e, this_ele in enumerate(all_e):
            try:
                all_e[go_e] = int(this_ele)
            except ValueError:
                try:
                    all_e[go_e] = float(this_ele)
                except ValueError:
                    pass
        return all_e


def get_graph_coverage_dict_simple(fasta_matrix, is_fastg):
    if is_fastg:
        coverages = {}
        for fastg_name in fasta_matrix[0]:
            this_coverage = float(fastg_name.split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0])
            coverages[fastg_name.split('_')[1]] = this_coverage
        return coverages
    else:
        return {}


def average_weighted_np_free(vals, weights):
    return sum([val * weights[go_v] for go_v, val in enumerate(vals)])/float(sum(weights))


def weighted_mean_and_std_np_free(values, weights):
    mean = average_weighted_np_free(values, weights=weights)
    std = average_weighted_np_free([(val-mean)**2 for val in values], weights=weights)**0.5
    return mean, std


def get_graph_coverages_range_simple(fasta_matrix, drop_low_percent=0.10, drop_high_percent=0.40):
    coverages = []
    lengths = []
    for fastg_name in fasta_matrix[0]:
        this_coverage = float(fastg_name.split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0])
        this_length = int(fastg_name.split('length_')[1].split('_cov_')[0])
        coverages.append(this_coverage)
        lengths.append(this_length)
    weights = [inside_cov * lengths[go_v] for go_v, inside_cov in enumerate(coverages)]
    sum_weights = sum(weights)
    assert drop_low_percent + drop_high_percent < 1
    # get coverage low threshold
    forward_sorted_cov_id = sorted(list(range(len(coverages))), key=lambda x: coverages[x])
    low_threshold = 0
    low_accumulated = 0
    low_weights_threshold = drop_low_percent * sum_weights
    while forward_sorted_cov_id and low_accumulated < low_weights_threshold:
        this_id = forward_sorted_cov_id.pop(0)
        low_threshold = coverages[this_id]
        low_accumulated += low_threshold * lengths[this_id]
    # get coverage high threshold
    reverse_sorted_cov_id = sorted(list(range(len(coverages))), key=lambda x: -coverages[x])
    high_threshold = inf
    high_accumulated = 0
    high_weights_threshold = drop_high_percent * sum_weights
    while reverse_sorted_cov_id and high_accumulated < high_weights_threshold:
        this_id = reverse_sorted_cov_id.pop(0)
        high_threshold = coverages[this_id]
        high_accumulated += high_threshold * lengths[this_id]
    # remove pairs outside [low_threshold, high_threshold]
    check_v = 0
    while check_v < len(coverages):
        if low_threshold <= coverages[check_v] <= high_threshold:
            check_v += 1
        else:
            del coverages[check_v]
            del lengths[check_v]
    cov_mean, cov_std = weighted_mean_and_std_np_free(coverages, lengths)
    return max(cov_mean - cov_std, min(coverages)), cov_mean, min(cov_mean + cov_std, max(coverages))


def filter_fastg_by_depth_simple(fas_file, max_depth, min_depth, log_handler=None, out_dir=None):
    bounded = (max_depth and float(max_depth) and float(max_depth) != inf) or (min_depth and float(min_depth))
    if fas_file.endswith('.fastg') and bounded:
        max_depth = float(max_depth)
        min_depth = float(min_depth)
        time0 = time.time()
        fastg_matrix = read_fasta(fas_file)
        new_fastg_matrix = [[], [], fastg_matrix[2]]
        for i in range(len(fastg_matrix[0])):
            if max_depth > float(fastg_matrix[0][i].split('cov_')[1].split(':')[0].split(';')[0].split('\'')[0]) \
                    >= min_depth:
                new_fastg_matrix[0].append(fastg_matrix[0][i])
                new_fastg_matrix[1].append(fastg_matrix[1][i])
        out_fasta = '.'.join(fas_file.split('.')[:-1]) + '.depth' + str(min_depth) + "-" + str(max_depth) + \
                    '.' + fas_file.split('.')[-1]
        if out_dir:
            out_fasta = os.path.join(out_dir, os.path.basename(out_fasta))
        write_fasta(out_file=out_fasta, matrix=new_fastg_matrix, overwrite=True)
        if log_handler:
            log_handler.info('filtering by depth cost: ' + str(round(time.time() - time0, 2)))
        else:
            sys.stdout.write('\nfiltering by depth cost: '+str(round(time.time()-time0, 2)) + "\n")
        return out_fasta
    else:
        return fas_file

