import time
import os
import sys
import random
from copy import deepcopy
from itertools import combinations
try:
    from sympy import Symbol, solve, lambdify
    from scipy import optimize, stats
    from numpy import array, dot
except ImportError:
    not_optimized = True
else:
    not_optimized = False
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *
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
    def __init__(self, fastg_file):
        self.vertex_info = {}
        self.__kmer = 127
        self.parse_fastg(fastg_file)
        self.vertex_clusters = []
        self.update_vertex_clusters()
        self.tagged_vertices = {}
        self.tag_loci = {}
        self.vertex_to_copy = {}
        self.copy_to_vertex = {}

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

    def parse_fastg(self, fastg_file):
        fastg_matrix = SequenceList(fastg_file)
        for i, seq in enumerate(fastg_matrix):
            if ":" in seq.label:
                this_vertex_str, next_vertices_str = seq.label.strip(";").split(":")
            else:
                this_vertex_str, next_vertices_str = seq.label.strip(";"), ""
            v_tag, vertex_name, l_tag, vertex_len, c_tag, vertex_cov = this_vertex_str.strip("'").split("_")
            if vertex_name not in self.vertex_info:
                self.vertex_info[vertex_name] = {"len": int(vertex_len),
                                                 "cov": float(vertex_cov),
                                                 "long": this_vertex_str.strip("'")}
            # connections
            if "connections" not in self.vertex_info[vertex_name]:
                self.vertex_info[vertex_name]["connections"] = {True: set(), False: set()}  # "tail"~True, "head"~False
            this_end = not this_vertex_str.endswith("'")
            self.vertex_info[vertex_name]["connections"][this_end] = set()
            if next_vertices_str:
                for next_vertex_str in next_vertices_str.split(","):
                    next_name = next_vertex_str.strip("'").split("_")[1]
                    next_end = next_vertex_str.endswith("'")
                    self.vertex_info[vertex_name]["connections"][this_end].add((next_name, next_end))
            # sequence
            if "seq" not in self.vertex_info[vertex_name]:
                self.vertex_info[vertex_name]["seq"] = {}
                if this_end:
                    self.vertex_info[vertex_name]["seq"][True] = seq.seq
                    self.vertex_info[vertex_name]["seq"][False] = complementary_seq(seq.seq)
                else:
                    self.vertex_info[vertex_name]["seq"][True] = complementary_seq(seq.seq)
                    self.vertex_info[vertex_name]["seq"][False] = seq.seq
        """delete vertices that not recorded"""
        for vertex_name in self.vertex_info:
            for this_end in (True, False):
                for next_name, next_end in list(self.vertex_info[vertex_name]["connections"][this_end]):
                    if next_name not in self.vertex_info:
                        self.vertex_info[vertex_name]["connections"][this_end].remove((next_name, next_end))
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

    def write_to_fastg(self, out_file):
        if not out_file.endswith(".fastg"):
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
        out_matrix.write_fasta(out_file)

    def write_to_gfa(self, out_file):
        if not out_file.endswith(".gfa"):
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
            # skip self.tag_loci
            if vertex_name in self.vertex_to_copy:
                this_copy = self.vertex_to_copy[vertex_name]
                self.copy_to_vertex[this_copy].remove(vertex_name)
                if not self.copy_to_vertex[this_copy]:
                    del self.copy_to_vertex[this_copy]
                del self.vertex_to_copy[vertex_name]
        if update_cluster:
            self.update_vertex_clusters()

    def merge_all_possible_vertices(self, limited_vertices=None, copy_tags=False):
        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)
        while limited_vertices:
            this_vertex = limited_vertices.pop()
            for this_end in (True, False):
                # print("this_vertex", this_vertex, this_end)
                connected_set = self.vertex_info[this_vertex]["connections"][this_end]
                if len(connected_set) == 1:
                    next_vertex, next_end = list(connected_set)[0]
                    # print("next_vertex", next_vertex, next_end)
                    # print("next_con_len", len(self.vertex_info[next_vertex]["connections"][next_end]))
                    if len(self.vertex_info[next_vertex]["connections"][next_end]) == 1 and this_vertex != next_vertex:
                        new_vertex = this_vertex + "_" + next_vertex  # + ("+", "-")[next_end]
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
                        self.remove_vertex([this_vertex, next_vertex], update_cluster=False)
                        break
        self.update_vertex_clusters()

    def estimate_copy_and_depth_by_cov(self, limited_vertices=None, given_average_cov=None, display=True):
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
                    this_copy = max(1, int(round(self.vertex_info[vertex_name]["cov"] / new_val, 0)))
                    self.vertex_to_copy[vertex_name] = this_copy
                    if this_copy not in self.copy_to_vertex:
                        self.copy_to_vertex[this_copy] = set()
                    self.copy_to_vertex[this_copy].add(vertex_name)
            if display:
                sys.stdout.write("updating average target cov: " + str(round(new_val, 4)) + "\n")
            return new_val
        else:
            # adjust this_copy according to user-defined depth
            for vertex_name in self.vertex_info:
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    if not self.copy_to_vertex[old_copy]:
                        del self.copy_to_vertex[old_copy]
                this_copy = max(1, int(round(self.vertex_info[vertex_name]["cov"] / given_average_cov, 0)))
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)
            return given_average_cov

    def estimate_copy_and_depth_precisely(self, limited_vertices=None, mode=None, display=True, maximum_copy_num=10):
        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)
        limited_vertices_set = set(limited_vertices)

        previous_val = {0.}
        new_val = -1.
        while new_val not in previous_val:
            previous_val.add(new_val)
            """ create constraints by creating multivariate equations """
            vertex_to_symbols = {vertex_name: Symbol("V"+vertex_name, integer=True)  # positive=True)
                                 for vertex_name in limited_vertices}
            symbols_to_vertex = {vertex_to_symbols[vertex_name]: vertex_name for vertex_name in limited_vertices}
            all_symbols = list(symbols_to_vertex)
            formulae = []
            def get_formula(from_vertex, from_end, to_vertex, to_end):
                result_form = vertex_to_symbols[from_vertex]
                for next_v, next_e in self.vertex_info[from_vertex]["connections"][from_end]:
                    if (next_v, next_e) != (to_vertex, to_end) and next_v in limited_vertices_set:
                        result_form -= get_formula(next_v, next_e, from_vertex, from_end)
                return result_form
            for vertex_name in limited_vertices:
                for this_end in (True, False):
                    connection_set = self.vertex_info[vertex_name]["connections"][this_end]
                    if len([n_v for n_v, n_e in connection_set if n_v in limited_vertices_set]):
                        this_formula = vertex_to_symbols[vertex_name]
                        for n_v, n_e in connection_set:
                            if n_v in limited_vertices_set:
                                this_formula -= get_formula(n_v, n_e, vertex_name, this_end)
                        formulae.append(this_formula)

            copy_solution = solve(formulae, all_symbols)
            copy_solution = copy_solution if copy_solution else {}
            for symbol_here in copy_solution:
                if copy_solution[symbol_here] == 0:
                    raise Exception("Incomplete/Complicated target graph!")
            free_copy_variables = list()
            for symbol_used in all_symbols:
                if symbol_used not in copy_solution:
                    free_copy_variables.append(symbol_used)
                    copy_solution[symbol_used] = symbol_used
            # the copy of every contig has to be >= 1
            coefficients = array(
                [[copy_solution[this_symbol].coeff(symbol_used) for symbol_used in free_copy_variables]
                 for this_symbol in all_symbols])
            min_copy = array([1.00001] * len(all_symbols))
            # print(self.vertex_to_copy)
            # if mode == "cp":
            #     for s_id, symbol_lim in enumerate(all_symbols):

            constraints = {'type': 'ineq', 'fun': lambda x: dot(coefficients, x) - min_copy}

            """ minimizing deviations from draft copy values """
            least_square_expr = sum([(solution_f - self.vertex_to_copy[symbols_to_vertex[symbol_used]]) ** 2
                                     for symbol_used, solution_f in copy_solution.items()])
            least_square_function = lambdify(args=free_copy_variables, expr=least_square_expr)
            # for compatibility between scipy and sympy
            def least_square_function_v(x):
                return least_square_function(*tuple(x))
            opt = {'disp': display, "maxiter": 1000}
            initials = array([random.randint(1, maximum_copy_num)] * len(free_copy_variables))
            bounds = [(1, maximum_copy_num)] * len(free_copy_variables)
            result = optimize.minimize(fun=least_square_function_v, x0=initials,
                                       method='SLSQP', bounds=bounds, constraints=constraints, options=opt)
            free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
                                        for i, this_copy in enumerate(result.x)}
            
            """ record new copy values """
            for this_symbol in all_symbols:
                vertex_name = symbols_to_vertex[this_symbol]
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    if not self.copy_to_vertex[old_copy]:
                        del self.copy_to_vertex[old_copy]
                this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
                # print(vertex_name, this_copy)
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)

            """ re-estimate baseline depth """
            total_product = 0.
            total_len = 0
            for vertex_name in limited_vertices:
                this_len = (self.vertex_info[vertex_name]["len"] - self.__kmer + 1) \
                           * self.vertex_to_copy.get(vertex_name, 1)
                this_cov = self.vertex_info[vertex_name]["cov"] / self.vertex_to_copy.get(vertex_name, 1)
                total_len += this_len
                total_product += this_len * this_cov
            new_val = total_product / total_len
        if display:
            sys.stdout.write("updating average target cov: " + str(round(new_val, 4)) + "\n")
        return new_val

    def find_target_graph(self, tab_file, mode="cp", weight_factor=100.0, depth_factor=None, max_copy=8,
                          temp_graph=None, display=True, debug=False):
        # parse_csv, every locus only occur in one vertex (removing locations with smaller weight)
        self.tag_loci = {}
        tab_matrix = [line.strip("\n").split("\t") for line in open(tab_file)][1:]
        for node_record in tab_matrix:
            vertex_name = node_record[0]
            if vertex_name in self.vertex_info:
                matched = node_record[5].split(">>")
                for locus in matched:
                    if "(" in locus:
                        locus_spl = locus.split("(")
                        locus_type = locus_spl[-1].split(",")[1][:-1]
                        if locus_type == mode:
                            locus_name = "(".join(locus_spl[:-1])
                            locus_start, locus_end = locus_spl[-1].split(",")[0].split("-")
                            locus_len = int(locus_end) - int(locus_start) + 1
                            if locus_name in self.tag_loci:
                                new_weight = locus_len * self.vertex_info[vertex_name]["cov"]
                                if new_weight > self.tag_loci[locus_name]["weight"]:
                                    self.tag_loci[locus_name] = {"vertex": vertex_name,
                                                                 "len": locus_len,
                                                                 "weight": new_weight}
                            else:
                                self.tag_loci[locus_name] = {"vertex": vertex_name,
                                                             "len": locus_len,
                                                             "weight": locus_len * self.vertex_info[vertex_name]["cov"]}
        # weight_tag = mode + "-weight"
        self.tagged_vertices[mode] = set()
        for locus_name in self.tag_loci:
            vertex_name = self.tag_loci[locus_name]["vertex"]
            loci_weight = self.tag_loci[locus_name]["weight"]
            # tags
            if "tags" not in self.vertex_info[vertex_name]:
                self.vertex_info[vertex_name]["tags"] = {}
            if mode in self.vertex_info[vertex_name]["tags"]:
                self.vertex_info[vertex_name]["tags"][mode].add(locus_name)
            else:
                self.vertex_info[vertex_name]["tags"][mode] = {locus_name}
            # weight
            if "weight" not in self.vertex_info[vertex_name]:
                self.vertex_info[vertex_name]["weight"] = {}
            if mode in self.vertex_info[vertex_name]["weight"]:
                self.vertex_info[vertex_name]["weight"][mode] += loci_weight
            else:
                self.vertex_info[vertex_name]["weight"][mode] = loci_weight
            self.tagged_vertices[mode].add(vertex_name)
        if len(self.tagged_vertices[mode]) == 0:
            raise Exception("No available target information found in " + tab_file)

        def update_depth_factor(limited_vs, average_cov=None):
            if len(self.tagged_vertices[mode]) == 1 or average_cov is None:
                df = 10
                if display:
                    sys.stdout.write("using initial depth-factor: " + str(df) + "\n")
                return df
            else:
                df = 3 * stats.tstd([self.vertex_info[this_v]["cov"]/self.vertex_to_copy.get(this_v, 1)
                                     for this_v in limited_vs])/average_cov
                if display:
                    sys.stdout.write("estimated depth-factor: " + str(df) + "\n")
                return df
        if depth_factor is None:
            auto_depth = True
        else:
            auto_depth = False

        new_assembly = deepcopy(self)
        changed = True
        initial_depth = True
        average_target_cov = 0.
        while changed:
            changed = False
            cluster_trimmed = True
            while cluster_trimmed:
                # remove low coverages
                first_round = True
                delete_those_vertices = set()
                while first_round or delete_those_vertices:
                    if auto_depth:
                        if initial_depth:
                            depth_factor = update_depth_factor(self.tagged_vertices[mode])
                            initial_depth = False
                        else:
                            # depth_factor = update_depth_factor(self.vertex_info, average_target_cov)
                            depth_factor = update_depth_factor(self.tagged_vertices[mode], average_target_cov)
                    delete_those_vertices = set()
                    average_target_cov = new_assembly.estimate_copy_and_depth_by_cov(new_assembly.tagged_vertices[mode],
                                                                                     display=display)
                    for vertex_name in new_assembly.vertex_info:
                        if new_assembly.vertex_info[vertex_name]["cov"] * depth_factor < average_target_cov \
                                or new_assembly.vertex_info[vertex_name]["cov"]/(max_copy*depth_factor) > average_target_cov:
                            delete_those_vertices.add(vertex_name)
                    if delete_those_vertices:
                        if display and debug:
                            sys.stdout.write("removing low coverage contigs: " + str(delete_those_vertices) + "\n")
                        changed = True
                        new_assembly.remove_vertex(delete_those_vertices)
                    first_round = False

                cluster_trimmed = False
                # chose the target cluster (best rank)
                if len(new_assembly.vertex_clusters) == 0:
                    raise Exception("No available target graph detected!")
                elif len(new_assembly.vertex_clusters) == 1:
                    pass
                else:
                    cluster_weights = [sum([new_assembly.vertex_info[x_v]["weight"][mode]
                                            for x_v in x if "weight" in new_assembly.vertex_info[x_v]])
                                       for x in new_assembly.vertex_clusters]
                    # print("cluster_weights", cluster_weights)
                    best = max(cluster_weights)
                    best_id = cluster_weights.index(best)
                    temp_cluster_weights = deepcopy(cluster_weights)
                    del temp_cluster_weights[best_id]
                    second = max(temp_cluster_weights)
                    if best < second * weight_factor:
                        if temp_graph:
                            new_assembly.write_to_gfa(temp_graph)
                        raise Exception("Complicated target graph!")
                    for j, w in enumerate(cluster_weights):
                        if w == second:
                            for del_v in new_assembly.vertex_clusters[j]:
                                if del_v in new_assembly.tagged_vertices[mode]:
                                    if new_assembly.vertex_info[del_v]["cov"] * depth_factor >= average_target_cov:
                                        if temp_graph:
                                            new_assembly.write_to_gfa(temp_graph)
                                        raise Exception("Complicated target graph: please check EDGE_" + del_v + "!")

                    # remove other clusters
                    vertices_to_del = set()
                    for go_cl, v_2_del in enumerate(new_assembly.vertex_clusters):
                        if go_cl != best_id:
                            vertices_to_del |= v_2_del
                    if vertices_to_del:
                        if display and debug:
                            sys.stdout.write("removing other clusters: " + str(vertices_to_del) + "\n")
                        new_assembly.remove_vertex(vertices_to_del)
                        cluster_trimmed = True
                        changed = True

            # merge vertices
            new_assembly.merge_all_possible_vertices()

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
                            raise Exception("Incomplete/Complicated target graph!")
                        else:
                            delete_those_vertices.add(vertex_name)
                if delete_those_vertices:
                    if display and debug:
                        sys.stdout.write("removing terminal contigs: " + str(delete_those_vertices) + "\n")
                    new_assembly.remove_vertex(delete_those_vertices)
                    changed = True

            # merge vertices
            new_assembly.merge_all_possible_vertices()

        if temp_graph:
            new_assembly.write_to_gfa(temp_graph)

        # create idealized vertices and edges
        new_average_cov = new_assembly.estimate_copy_and_depth_by_cov(display=display)
        if not_optimized:
            if display:
                sys.stdout.write("Warning: numpy/scipy/sympy not installed, using coverage information only!\n")
        else:
            try:
                if display:
                    sys.stdout.write("Estimating copy and depth precisely ...\n")
                new_average_cov = new_assembly.estimate_copy_and_depth_precisely(maximum_copy_num=max_copy,
                                                                                 display=display)
            except RecursionError:
                if temp_graph:
                    new_assembly.write_to_gfa(temp_graph)
                raise Exception("Complicated target graph! Detecting path(s) failed!\n")
        return {"graph": new_assembly, "cov": new_average_cov}

    def get_all_circular_paths(self):

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
                    new_path.append((next_vertex, next_end))
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
        if 1 not in self.copy_to_vertex:
            raise Exception("No single copy region?! Detecting path(s) failed!\n")
        start_vertex = sorted(self.copy_to_vertex[1], key=lambda x: -self.vertex_info[x]["len"])[0]
        first_path = [(start_vertex, False)]
        first_connections = self.vertex_info[start_vertex]["connections"][not False]
        vertex_to_copy = deepcopy(self.vertex_to_copy)
        del vertex_to_copy[start_vertex]
        circular_directed_graph_solver(first_path, first_connections, vertex_to_copy)

        if not paths:
            raise Exception("Complicated target graph! Detecting path(s) failed!\n")
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
                pattern_dict = {acc_distance: ad_id for
                                ad_id, acc_distance in enumerate(sorted(set([x[1] for x in sorted_paths])))}
                sorted_paths = [(this_path, ".nesting_repeat_pattern_" + str(pattern_dict[acc_distance]))
                                for this_path, acc_distance in sorted_paths]
            else:
                sorted_paths = [(this_path, "") for this_path in sorted(paths)]
            return sorted_paths

    def export_path(self, in_path):
        seq_names = []
        seq_segments = []
        for this_vertex, this_end in in_path:
            seq_segments.append(self.vertex_info[this_vertex]["seq"][not this_end][self.__kmer:])
            seq_names.append(this_vertex + ("+", "-")[not this_end])
        # if not circular
        if in_path[0] not in self.vertex_info[in_path[-1][0]]["connections"][not in_path[-1][1]]:
            seq_segments[0] = self.vertex_info[in_path[0][0]]["seq"][not in_path[0][1]][:self.__kmer] + seq_segments[0]
        else:
            seq_names[-1] += "(circular)"
        return Sequence(",".join(seq_names), "".join(seq_segments))

