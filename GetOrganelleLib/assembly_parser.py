import os
import sys
from itertools import combinations, product
from hashlib import sha256
from collections import OrderedDict

try:
    from sympy import Symbol, solve, lambdify
    from sympy import log as symlog
    # from scipy import optimize
except ImportError:
    def Symbol(foo, integer):
        raise ImportError("Failed in 'from sympy import Symbol, solve, lambdify, log'!")


    def solve(foo1, foo2):
        raise ImportError("Failed in 'from sympy import Symbol, solve, lambdify, log'!")


    def lambdify(args=None, expr=None):
        raise ImportError("Failed in 'from sympy import Symbol, solve, lambdify, log'!")


    def symlog(foo):
        raise ImportError("Failed in 'from sympy import Symbol, solve, lambdify, log'!")

try:
    from gekko import GEKKO
except ImportError:
    def GEKKO(remote):
        raise ImportError("Failed in 'from gekko import GEKKO'!")


    # class optimize:
    #     def __init__(self):
    #         pass
    #
    #     def minimize(self, fun=None, x0=None, jac=None, method=None, bounds=None, constraints=None, options=None):
    #         raise ImportError("Failed in 'from scipy import optimize'!")

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
# for test
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.statistical_func import *
from GetOrganelleLib.pipe_control_func import log_target_res

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
import random
from copy import deepcopy

MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    python_version = "2.7+"
    RecursionError = RuntimeError
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)
ECHO_DIRECTION = ["_tail", "_head"]
VERTEX_DIRECTION_STR_TO_BOOL = {"+": True, "-": False}
VERTEX_DIRECTION_BOOL_TO_STR = {True: "+", False: "-"}


class ProcessingGraphFailed(Exception):
    def __init__(self, value=""):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Vertex(object):
    def __init__(self, v_name, length=None, coverage=None, forward_seq=None, reverse_seq=None,
                 tail_connections=None, head_connections=None, fastg_form_long_name=None):
        """
        :param v_name: str
        :param length: int
        :param coverage: float
        :param forward_seq: str
        :param reverse_seq: str
        :param tail_connections: OrderedDict()
        :param head_connections: OrderedDict()
        :param fastg_form_long_name: str
        self.seq={True: FORWARD_SEQ, False: REVERSE_SEQ}
        self.connections={True: tail_connection, False: head_connection}
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
        elif forward_seq is None:
            self.seq = {True: None, False: None}
        else:
            self.seq = {True: None, False: None}

        # True: tail, False: head
        self.connections = {True: OrderedDict(), False: OrderedDict()}
        assert tail_connections is None or isinstance(tail_connections, OrderedDict), \
            "tail_connections must be an OrderedDict()"
        assert head_connections is None or isinstance(head_connections, OrderedDict), \
            "head_connections must be an OrderedDict()"
        if tail_connections:
            self.connections[True] = tail_connections
        if head_connections:
            self.connections[False] = head_connections
        self.fastg_form_name = fastg_form_long_name
        self.merging_history = MergingHistory([(v_name, True)])
        self.other_attr = {}

    def __repr__(self):
        return self.name

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

    def is_terminal(self):
        return not (self.connections[True] and self.connections[False])

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
        val.name = key
        dict.__setitem__(self, key, val)


class MergingHistory(object):
    def __init__(self, history_or_path=None):
        self.__list = []
        if history_or_path:
            for each_item in history_or_path:
                is_vertex = isinstance(each_item, tuple) and len(each_item) == 2 and isinstance(each_item[1], bool)
                is_hist = isinstance(each_item, MergingHistory)
                assert is_vertex or is_hist
                if is_vertex:
                    self.__list.append(each_item)
                else:
                    self.__list.extend(each_item.list())

    def add(self, new_history_or_vertex, add_new_to_front=False, reverse_the_new=False):
        is_vertex = isinstance(new_history_or_vertex, tuple) and len(new_history_or_vertex) == 2
        is_hist = isinstance(new_history_or_vertex, MergingHistory)
        assert is_vertex or is_hist
        if add_new_to_front:
            if reverse_the_new:
                if is_vertex:
                    self.__list.insert(0, (new_history_or_vertex[0], not new_history_or_vertex[1]))
                else:
                    self.__list = list(-new_history_or_vertex) + self.__list
            else:
                if is_vertex:
                    self.__list.insert(0, new_history_or_vertex)
                else:
                    self.__list = list(new_history_or_vertex) + self.__list
        else:
            if reverse_the_new:
                if is_vertex:
                    self.__list.append((new_history_or_vertex[0], not new_history_or_vertex[1]))
                else:
                    self.__list.extend(list(-new_history_or_vertex))
            else:
                if is_vertex:
                    self.__list.append(new_history_or_vertex)
                else:
                    self.__list.extend(list(new_history_or_vertex))

    def __neg__(self):
        return MergingHistory([(each_vertex[0], not each_vertex[1]) for each_vertex in self.__list[::-1]])

    def __iter__(self):
        for item in self.__list:
            yield item

    def __str__(self):
        return "_".join([str(each_vertex[0]) if isinstance(each_vertex, tuple) else str(each_vertex)
                         for each_vertex in self.__list])

    def reverse(self):
        self.__list = [(each_vertex[0], not each_vertex[1]) for each_vertex in self.__list[::-1]]

    def path_list(self):
        return list(self.__list)
        # return [each_vertex.path_list() if isinstance(each_vertex, MergingHistory) else each_vertex
        #         for each_vertex in self.__list]

    def vertex_set(self):
        v_set = set()
        for each_item in self.__list:
            if isinstance(each_item[0], ConsensusHistory):
                v_set.update(each_item[0].vertex_set())
            else:
                v_set.add(each_item[0])
        return v_set


class ConsensusHistory(object):
    def __init__(self, merging_items):
        """
        :param merging_items: [(name1, vertex1), (name2, vertex2)]
        """
        self.__list = []
        for each_item in merging_items:
            assert isinstance(each_item, tuple) and len(each_item) == 2
            self.__list.append(each_item)

    def add_consensus(self, vertex):
        """
        :param vertex: (name, vertex)
        :return:
        """
        assert isinstance(vertex, tuple) and len(vertex) == 2
        self.__list.append(vertex)

    def __str__(self):
        sorted_list = sorted(self.__list, key=lambda x: (str(x[0]), x[1]))
        # format of consensus name, similar to newick
        return "(" + "|".join([str(each_hist[0]) + VERTEX_DIRECTION_BOOL_TO_STR[each_hist[1]]
                               for each_hist in sorted_list]) + ")"

    def vertex_set(self):
        v_set = set()
        for each_item in self.__list:
            if isinstance(each_item[0], MergingHistory):
                v_set.update(each_item[0].vertex_set())
            else:
                v_set.add(each_item[0])
        return v_set


class SimpleAssembly(object):
    def __init__(self, graph_file=None, min_cov=0., max_cov=inf):
        """
        :param graph_file:
        :param min_cov:
        :param max_cov:
        """
        self.vertex_info = VertexInfo()
        # self.__overlap = overlap
        self.__uni_overlap = None
        if graph_file:
            if graph_file.endswith(".gfa"):
                self.parse_gfa(graph_file, min_cov=min_cov, max_cov=max_cov)
            else:
                self.parse_fastg(graph_file, min_cov=min_cov, max_cov=max_cov)

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

    def __bool__(self):
        return bool(self.vertex_info)

    def __iter__(self):
        for vertex in sorted(self.vertex_info):
            yield self.vertex_info[vertex]

    def parse_gfa(self, gfa_file, default_cov=1., min_cov=0., max_cov=inf):
        with open(gfa_file) as gfa_open:
            kmer_values = set()
            line = gfa_open.readline()
            gfa_version_number = "1.0"
            if line.startswith("H\t"):
                for element in line.strip().split("\t")[1:]:
                    element = element.split(":")
                    element_tag, element_type, element_description = element[0], element[1], ":".join(element[2:])
                    if element_tag == "VN":
                        gfa_version_number = element_description
            gfa_open.seek(0)
            if gfa_version_number == "1.0":
                for line in gfa_open:
                    if line.startswith("S\t"):
                        elements = line.strip().split("\t")
                        elements.pop(0)  # record_type
                        vertex_name = elements.pop(0)  # segment name
                        sequence = elements.pop(0)
                        seq_len_tag = None
                        kmer_count = None
                        seq_depth_tag = None
                        sh_256_val = None
                        other_attributes = {}
                        for element in elements:
                            element = element.split(":")  # element_tag, element_type, element_description
                            # skip RC/FC
                            if element[0].upper() == "LN":
                                seq_len_tag = int(element[-1])
                            elif element[0].upper() == "KC":
                                kmer_count = int(element[-1])
                            elif element[0].upper() == "RC":  # took read counts as kmer counts
                                kmer_count = int(element[-1])
                            elif element[0].upper() == "DP":
                                seq_depth_tag = float(element[-1])
                            elif element[0].upper() == "RD":  # took read depth as seq_depth_tag counts
                                seq_depth_tag = int(element[-1])
                            elif element[0].upper() == "SH":
                                sh_256_val = ":".join(element[2:])
                            elif element[0].upper() == "UR":
                                seq_file_path = element[-1]
                                if os.path.isfile(seq_file_path):
                                    if sequence == "*":
                                        sequence = "".join([sub_seq.strip() for sub_seq in open(seq_file_path)])
                                    else:
                                        tag_seq = "".join([sub_seq.strip() for sub_seq in open(seq_file_path)])
                                        if tag_seq != sequence:
                                            raise ProcessingGraphFailed(
                                                vertex_name + " sequences from different sources!")
                                else:
                                    raise ProcessingGraphFailed(
                                        seq_file_path + " for " + vertex_name + " does not exist!")
                            else:
                                other_attributes[element[0].upper()] = element[-1]
                        seq_len = len(sequence)
                        if seq_len_tag is not None and seq_len != seq_len_tag:
                            raise ProcessingGraphFailed(vertex_name + " has unmatched sequence length as noted!")
                        if sh_256_val is not None and sh_256_val != sha256(sequence):
                            raise ProcessingGraphFailed(vertex_name + " has unmatched sha256 value as noted!")
                        if kmer_count is not None or seq_depth_tag is not None:
                            if kmer_count is not None:
                                seq_depth = kmer_count / float(seq_len)
                            elif seq_depth_tag is not None:
                                seq_depth = seq_depth_tag
                            if min_cov <= seq_depth <= max_cov:
                                vert = Vertex(vertex_name, seq_len, seq_depth, sequence)
                                vert.other_attr = other_attributes
                                self.vertex_info[vertex_name] = vert
                                if vertex_name.isdigit():
                                    self.vertex_info[vertex_name].fill_fastg_form_name()
                        else:
                            vert = Vertex(vertex_name, seq_len, default_cov, sequence)
                            vert.other_attr = other_attributes
                            self.vertex_info[vertex_name] = vert
                gfa_open.seek(0)
                for line in gfa_open:
                    if line.startswith("L\t"):
                        elements = line.strip().split("\t")
                        elements.pop(0)  # flag
                        vertex_1 = elements.pop(0)
                        end_1 = elements.pop(0)
                        vertex_2 = elements.pop(0)
                        end_2 = elements.pop(0)
                        alignment_cigar = elements.pop(0)
                        # "head"~False, "tail"~True
                        if vertex_1 in self.vertex_info and vertex_2 in self.vertex_info:
                            # different notation as to our Assembly system
                            end_1 = {"+": True, "-": False}[end_1]
                            end_2 = {"+": False, "-": True}[end_2]
                            if not (alignment_cigar[:-1].isdigit() and alignment_cigar[-1] == "M"):
                                raise ProcessingGraphFailed(
                                    "Unmatched overlaps are currently not supported: " + alignment_cigar)
                            else:
                                overlap = int(alignment_cigar[:-1])
                            self.vertex_info[vertex_1].connections[end_1][(vertex_2, end_2)] = overlap
                            self.vertex_info[vertex_2].connections[end_2][(vertex_1, end_1)] = overlap
                            kmer_values.add(overlap)
            elif gfa_version_number == "2.0":
                for line in gfa_open:
                    if line.startswith("S\t"):
                        elements = line.strip().split("\t")
                        elements.pop(0)  # record_type
                        vertex_name = elements.pop(0)  # segment name
                        int(elements.pop(0))  # seq_len_tag
                        sequence = elements.pop(0)
                        seq_len_tag = None
                        kmer_count = None
                        seq_depth_tag = None
                        sh_256_val = None
                        other_attributes = {}
                        for element in elements:
                            element = element.split(":")  # element_tag, element_type, element_description
                            # skip RC/FC
                            if element[0].upper() == "KC":
                                kmer_count = int(element[-1])
                            elif element[0].upper() == "RC":  # took read counts as kmer counts
                                kmer_count = int(element[-1])
                            elif element[0].upper() == "DP":
                                seq_depth_tag = float(element[-1])
                            elif element[0].upper() == "RD":  # took read depth as seq_depth_tag counts
                                seq_depth_tag = int(element[-1])
                            elif element[0].upper() == "SH":
                                sh_256_val = ":".join(element[2:])
                            elif element[0].upper() == "UR":
                                seq_file_path = element[-1]
                                if os.path.isfile(seq_file_path):
                                    if sequence == "*":
                                        sequence = "".join([sub_seq.strip() for sub_seq in open(seq_file_path)])
                                    else:
                                        tag_seq = "".join([sub_seq.strip() for sub_seq in open(seq_file_path)])
                                        if tag_seq != sequence:
                                            raise ProcessingGraphFailed(
                                                vertex_name + " sequences from different sources!")
                                else:
                                    raise ProcessingGraphFailed(
                                        seq_file_path + " for " + vertex_name + " does not exist!")
                            else:
                                other_attributes[element[0].upper()] = element[-1]
                        seq_len = len(sequence)
                        if seq_len_tag is not None and seq_len != seq_len_tag:
                            raise ProcessingGraphFailed(vertex_name + " has unmatched sequence length as noted!")
                        if sh_256_val is not None and sh_256_val != sha256(sequence):
                            raise ProcessingGraphFailed(vertex_name + " has unmatched sha256 value as noted!")
                        if kmer_count is not None or seq_depth_tag is not None:
                            if kmer_count is not None:
                                seq_depth = kmer_count / float(seq_len)
                            elif seq_depth_tag is not None:
                                seq_depth = seq_depth_tag
                            if min_cov <= seq_depth <= max_cov:
                                self.vertex_info[vertex_name] = Vertex(vertex_name, seq_len, seq_depth, sequence)
                                self.vertex_info[vertex_name].other_attr = other_attributes
                                if vertex_name.isdigit():
                                    self.vertex_info[vertex_name].fill_fastg_form_name()
                        else:
                            self.vertex_info[vertex_name] = Vertex(vertex_name, seq_len, default_cov, sequence)
                            self.vertex_info[vertex_name].other_attr = other_attributes
                gfa_open.seek(0)
                for line in gfa_open:
                    if line.startswith("E\t"):  # gfa2 uses E
                        elements = line.strip().split("\t")
                        elements.pop(0)  # flag
                        vertex_1 = elements.pop(0)
                        end_1 = elements.pop(0)
                        vertex_2 = elements.pop(0)
                        end_2 = elements.pop(0)
                        alignment_cigar = elements.pop(0)
                        # "head"~False, "tail"~True
                        if vertex_1 in self.vertex_info and vertex_2 in self.vertex_info:
                            end_1 = {"+": True, "-": False}[end_1]
                            end_2 = {"+": False, "-": True}[end_2]
                            # kmer_values.add(alignment_cigar)
                            if not (alignment_cigar[:-1].isdigit() and alignment_cigar[-1] == "M"):
                                raise ProcessingGraphFailed(
                                    "Unmatched overlaps are currently not supported: " + alignment_cigar)
                            else:
                                overlap = int(alignment_cigar[:-1])
                            self.vertex_info[vertex_1].connections[end_1][(vertex_2, end_2)] = overlap
                            self.vertex_info[vertex_2].connections[end_2][(vertex_1, end_1)] = overlap
                            kmer_values.add(overlap)
            else:
                raise ProcessingGraphFailed("Unrecognized GFA version number: " + gfa_version_number)
            if len(kmer_values) == 1:
                self.__uni_overlap = kmer_values.pop()
            # if len(kmer_values) == 0:
            #     self.__overlap = None
            # elif len(kmer_values) > 1:
            #     raise ProcessingGraphFailed("Multiple overlap values: " + ",".join(sorted(kmer_values)))
            # else:
            #     self.__overlap = int(kmer_values.pop()[:-1])

    def parse_fastg(self, fastg_file, min_cov=0., max_cov=inf):
        fastg_matrix = SequenceList(fastg_file)
        # initialize names; only accept vertex that are formally stored, skip those that are only mentioned after ":"
        for i, seq in enumerate(fastg_matrix):
            if ":" in seq.label:
                this_vertex_str, next_vertices_str = seq.label.strip(";").split(":")
            else:
                this_vertex_str, next_vertices_str = seq.label.strip(";"), ""
            v_tag, vertex_name, l_tag, vertex_len, c_tag, vertex_cov = this_vertex_str.strip("'").split("_")
            # skip vertices_set with cov out of bounds
            vertex_cov = float(vertex_cov)
            if not (min_cov <= vertex_cov <= max_cov):
                continue
            if vertex_name not in self.vertex_info:
                self.vertex_info[vertex_name] = Vertex(vertex_name, int(vertex_len), vertex_cov,
                                                       fastg_form_long_name=this_vertex_str.strip("'"))
        # adding other info based on existed names
        for i, seq in enumerate(fastg_matrix):
            if ":" in seq.label:
                this_vertex_str, next_vertices_str = seq.label.strip(";").split(":")
            else:
                this_vertex_str, next_vertices_str = seq.label.strip(";"), ""
            v_tag, vertex_name, l_tag, vertex_len, c_tag, vertex_cov = this_vertex_str.strip("'").split("_")
            # skip vertices_set that not in self.vertex_info: 1. with cov out of bounds
            if vertex_name in self.vertex_info:
                # connections
                this_end = not this_vertex_str.endswith("'")
                if next_vertices_str:
                    for next_vertex_str in next_vertices_str.split(","):
                        next_name = next_vertex_str.strip("'").split("_")[1]
                        if next_name in self.vertex_info:
                            next_end = next_vertex_str.endswith("'")
                            # Adding connection information (edge) to both of the related vertices_set
                            # even it is only mentioned once in some SPAdes output files
                            self.vertex_info[vertex_name].connections[this_end][(next_name, next_end)] = 0  # None?
                            self.vertex_info[next_name].connections[next_end][(vertex_name, this_end)] = 0
                # sequence
                if not self.vertex_info[vertex_name].seq[True]:
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
            general_overlap = 0
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
                general_overlap = max(initial_kmer)
            else:
                general_overlap = 0
                # raise ProcessingGraphFailed("No kmer detected!")
        # assign general kmer to all edges
        for vertex_name in self.vertex_info:
            for this_end in (True, False):
                for next_tuple in self.vertex_info[vertex_name].connections[this_end]:
                    self.vertex_info[vertex_name].connections[this_end][next_tuple] = general_overlap
        self.__uni_overlap = general_overlap

    def uni_overlap(self):
        if self.__uni_overlap is None:
            return None
        else:
            return int(self.__uni_overlap)

    def write_to_fasta(self, out_file, interleaved=None, check_postfix=True):
        if check_postfix and not out_file.endswith(".fasta"):
            out_file += ".fasta"
        out_matrix = SequenceList()
        for vertex_name in self.vertex_info:
            out_matrix.append(Sequence(vertex_name, self.vertex_info[vertex_name].seq[True]))
        out_matrix.interleaved = 70
        out_matrix.write_fasta(out_file, interleaved=interleaved)

    def write_to_gfa(self, out_file, check_postfix=True, other_attr=None):
        """
        :param out_file: str
        :param check_postfix: bool
        :param other_attr: dict, e.g. {"CL":"z", "C2":"z"}
        """
        if check_postfix and not out_file.endswith(".gfa"):
            out_file += ".gfa"
        if not other_attr:
            other_attr = {}
        out_file_handler = open(out_file, "w")
        for vertex_name in self.vertex_info:
            out_file_handler.write("\t".join(
                [
                    "S", vertex_name, self.vertex_info[vertex_name].seq[True],
                    "LN:i:" + str(self.vertex_info[vertex_name].len),
                    "RC:i:" + str(int(self.vertex_info[vertex_name].len * self.vertex_info[vertex_name].cov))] +
                [
                    "%s:%s:%s" % (attr_name, attr_type, self.vertex_info[vertex_name].other_attr.get(attr_name, ""))
                    for attr_name, attr_type in other_attr.items()
                    if self.vertex_info[vertex_name].other_attr.get(attr_name, False)
                ]) + "\n")
        recorded_connections = set()
        for vertex_name in self.vertex_info:
            for this_end in (False, True):
                for (next_v, next_e), this_overlap in self.vertex_info[vertex_name].connections[this_end].items():
                    this_con = tuple(sorted([(vertex_name, this_end), (next_v, next_e)]))
                    if this_con not in recorded_connections:
                        recorded_connections.add(this_con)
                        out_file_handler.write("\t".join([
                            "L", vertex_name, ("-", "+")[this_end], next_v, ("-", "+")[not next_e],
                            str(this_overlap) + "M"
                        ]) + "\n")


class Assembly(SimpleAssembly):
    def __init__(self, graph_file=None, min_cov=0., max_cov=inf, uni_overlap=None):
        """
        :param graph_file:
        :param min_cov:
        :param max_cov:
        """
        super(Assembly, self).__init__(graph_file=graph_file, min_cov=min_cov, max_cov=max_cov)
        if uni_overlap:
            self.__uni_overlap = uni_overlap
        else:
            self.__uni_overlap = super(Assembly, self).uni_overlap()
        self.vertex_clusters = []
        self.update_vertex_clusters()
        self.tagged_vertices = {}
        self.vertex_to_copy = {}
        self.vertex_to_float_copy = {}
        self.copy_to_vertex = {}
        self.__inverted_repeat_vertex = {}

        # Note: traversome move palindromic_repeats here rather than inside get_all_circular_paths/get_all_paths

    # def uni_overlap(self):
    #     if self.__uni_overlap is None:
    #         return None
    #     else:
    #         return int(self.__uni_overlap)

    def new_graph_with_vertex_reseeded(self, start_from=1):
        those_vertices = sorted(self.vertex_info)
        new_graph = Assembly(uni_overlap=self.__uni_overlap)
        name_trans = {those_vertices[go - start_from]: str(go)
                      for go in range(start_from, start_from + len(those_vertices))}
        for old_name in those_vertices:
            new_name = name_trans[old_name]
            this_v_info = deepcopy(self.vertex_info[old_name])
            this_v_info.name = new_name
            this_v_info.connections = {True: OrderedDict(), False: OrderedDict()}
            for this_end in self.vertex_info[old_name].connections:
                for next_name, next_end in self.vertex_info[old_name].connections[this_end]:
                    this_v_info.connections[this_end][(name_trans[next_name], next_end)] = \
                        self.vertex_info[old_name].connections[this_end][(next_name, next_end)]
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
                        sys.stdout.write("Graph converted to new fastg with original Vertex names lost.\n")
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
                                             out_renaming_table + ".\n")
            else:
                raise ProcessingGraphFailed(
                    "Merged graph cannot be written as fastg format file, please try gfa format!")

    def write_out_tags(self, db_names, out_file):
        tagged_vertices = set()
        for db_n in db_names:
            tagged_vertices |= self.tagged_vertices[db_n]
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
                here_tags = {tag_n for tag_n in db_names if this_vertex in self.tagged_vertices[tag_n]}
                lines.append([this_vertex,
                              ";".join(sorted(here_tags)),
                              "", ""])
        open(out_file, "w").writelines(["\t".join(line) + "\n" for line in lines])

    def update_orf_total_len(self, limited_vertices=None):
        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)
        for vertex_name in limited_vertices:
            self.vertex_info[vertex_name].other_attr["orf"] = {}
            for direction in (True, False):
                this_orf_lens = get_orf_lengths(self.vertex_info[vertex_name].seq[direction])
                self.vertex_info[vertex_name].other_attr["orf"][direction] = {"lengths": this_orf_lens,
                                                                              "sum_len": sum(this_orf_lens)}

    # def update_vertex_clusters(self):
    #     self.vertex_clusters = []
    #     vertices_set = sorted(self.vertex_info)
    #     for this_vertex in vertices_set:
    #         connecting_those = set()
    #         for connected_set in self.vertex_info[this_vertex].connections.values():
    #             for next_v, next_d in connected_set:
    #                 for go_to_set, cluster in enumerate(self.vertex_clusters):
    #                     if next_v in cluster:
    #                         connecting_those.add(go_to_set)
    #         if not connecting_those:
    #             self.vertex_clusters.append({this_vertex})
    #         elif len(connecting_those) == 1:
    #             self.vertex_clusters[connecting_those.pop()].add(this_vertex)
    #         else:
    #             sorted_those = sorted(connecting_those, reverse=True)
    #             self.vertex_clusters[sorted_those[-1]].add(this_vertex)
    #             for go_to_set in sorted_those[:-1]:
    #                 for that_vertex in self.vertex_clusters[go_to_set]:
    #                     self.vertex_clusters[sorted_those[-1]].add(that_vertex)
    #                 del self.vertex_clusters[go_to_set]

    def update_vertex_clusters(self):
        """
        faster than original v1.7.7.0 update_vertex_clusters algorithm. 2022-12-18.
        """
        self.vertex_clusters = []
        candidate_vs = set(self.vertex_info)
        while candidate_vs:
            new_root = candidate_vs.pop()
            self.vertex_clusters.append({new_root})
            waiting_vs = set([next_v
                              for this_e in (True, False)
                              for next_v, next_e in self.vertex_info[new_root].connections[this_e]
                              if next_v in candidate_vs])
            while candidate_vs and waiting_vs:
                next_v = waiting_vs.pop()
                self.vertex_clusters[-1].add(next_v)
                candidate_vs.discard(next_v)
                for next_e in (True, False):
                    for n_next_v, n_next_e in self.vertex_info[next_v].connections[next_e]:
                        if n_next_v in candidate_vs:
                            waiting_vs.add(n_next_v)
        # for reproducible, not necessary for some cases
        self.vertex_clusters.sort(key=lambda x: max(x))

    def get_clusters(self, limited_vertices=None):
        if limited_vertices is None:
            candidate_vs = set(self.vertex_info)
        else:
            candidate_vs = set(limited_vertices)
        vertex_clusters = []
        while candidate_vs:
            new_root = candidate_vs.pop()
            vertex_clusters.append({new_root})
            waiting_vs = set([next_v
                              for this_e in (True, False)
                              for next_v, next_e in self.vertex_info[new_root].connections[this_e]
                              if next_v in candidate_vs])
            while candidate_vs and waiting_vs:
                next_v = waiting_vs.pop()
                vertex_clusters[-1].add(next_v)
                candidate_vs.discard(next_v)
                for next_e in (True, False):
                    for n_next_v, n_next_e in self.vertex_info[next_v].connections[next_e]:
                        if n_next_v in candidate_vs:
                            waiting_vs.add(n_next_v)
        # for reproducible, not necessary for some cases
        return sorted(vertex_clusters, key=lambda x: max(x))

    def check_connected(self, vertices_set):
        """
        a fast algorithm modified from update_vertex_clusters
        :param vertices_set:
        :return:
        """
        candidate_vs = set(vertices_set)
        while candidate_vs:
            new_root = candidate_vs.pop()
            waiting_vs = set([next_v
                              for this_e in (True, False)
                              for next_v, next_e in self.vertex_info[new_root].connections[this_e]
                              if next_v in candidate_vs])
            while candidate_vs and waiting_vs:
                next_v = waiting_vs.pop()
                candidate_vs.discard(next_v)
                for next_e in (True, False):
                    for n_next_v, n_next_e in self.vertex_info[next_v].connections[next_e]:
                        if n_next_v in candidate_vs:
                            waiting_vs.add(n_next_v)
            if candidate_vs and not waiting_vs:
                return False
        return True

    def remove_vertex(self, vertices, update_cluster=True):
        for vertex_name in vertices:
            for this_end, connected_dict in list(self.vertex_info[vertex_name].connections.items()):
                for next_v, next_e in list(connected_dict):
                    del self.vertex_info[next_v].connections[next_e][(vertex_name, this_end)]
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
            # if vertex_name in self.merging_history:
            #     del self.merging_history[vertex_name]
        if update_cluster:
            self.update_vertex_clusters()
        self.__inverted_repeat_vertex = {}

    def rename_vertex(self, old_vertex, new_vertex, update_cluster=True):
        assert old_vertex != new_vertex
        assert new_vertex not in self.vertex_info, new_vertex + " exists!"
        self.vertex_info[new_vertex] = deepcopy(self.vertex_info[old_vertex])
        self.vertex_info[new_vertex].name = new_vertex
        # merging history removed
        self.vertex_info[new_vertex].merging_history = MergingHistory([(new_vertex, True)])
        for this_end in (True, False):
            for next_v, next_e in list(self.vertex_info[new_vertex].connections[this_end]):
                self.vertex_info[next_v].connections[next_e][(new_vertex, this_end)] = \
                    self.vertex_info[next_v].connections[next_e][(old_vertex, this_end)]
                del self.vertex_info[next_v].connections[next_e][(old_vertex, this_end)]
        for tag in self.tagged_vertices:
            if old_vertex in self.tagged_vertices[tag]:
                self.tagged_vertices[tag].add(new_vertex)
                self.tagged_vertices[tag].remove(old_vertex)
        if old_vertex in self.vertex_to_copy:
            this_copy = self.vertex_to_copy[old_vertex]
            self.copy_to_vertex[this_copy].remove(old_vertex)
            self.copy_to_vertex[this_copy].add(new_vertex)
            self.vertex_to_copy[new_vertex] = self.vertex_to_copy[old_vertex]
            del self.vertex_to_copy[old_vertex]
            self.vertex_to_float_copy[new_vertex] = self.vertex_to_float_copy[old_vertex]
            del self.vertex_to_float_copy[old_vertex]
        if self.vertex_info[old_vertex].fastg_form_name:
            split_long_name = self.vertex_info[old_vertex].fastg_form_name.split("_")
            self.vertex_info[new_vertex].fastg_form_name = \
                "_".join([split_long_name[0], new_vertex] + split_long_name[2:])
        del self.vertex_info[old_vertex]
        if update_cluster:
            for go_c, v_cluster in enumerate(self.vertex_clusters):
                if old_vertex in v_cluster:
                    self.vertex_clusters[go_c].remove(old_vertex)
                    self.vertex_clusters[go_c].add(new_vertex)
        # if old_vertex in self.merging_history:
        #     self.merging_history[new_vertex] = self.merging_history[old_vertex]
        #     del self.merging_history[old_vertex]

    def detect_parallel_vertices(self, limited_vertices=None, detect_neighbors=True):
        if not limited_vertices:
            limiting = False
            limited_vertices = sorted(self.vertex_info)
        else:
            limiting = True
            limited_vertices = sorted(limited_vertices)
        all_both_ends = {}
        for vertex_name in limited_vertices:
            this_cons = self.vertex_info[vertex_name].connections
            connect_1 = this_cons[True]
            connect_2 = this_cons[False]
            if connect_1 and connect_2:
                this_ends_raw = [tuple(sorted(connect_1)), tuple(sorted(connect_2))]
                this_ends = sorted(this_ends_raw)
                direction_remained = this_ends_raw == this_ends
                this_ends = tuple(this_ends)
                if this_ends not in all_both_ends:
                    all_both_ends[this_ends] = set()
                all_both_ends[this_ends].add((vertex_name, direction_remained))
        if limiting and detect_neighbors:
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

    def find_pair_closing_the_path(self, start_v, start_e, terminating_end_set, starting_end_set):
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
            elif (in_search_v, in_search_e) in starting_end_set:
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

    def is_sequential_repeat(self, search_vertex_name, return_pair_in_the_trunk_path=True):
        if search_vertex_name not in self.vertex_info:
            raise ProcessingGraphFailed("Vertex name " + search_vertex_name + " not found!")
        connection_set_t = self.vertex_info[search_vertex_name].connections[True]
        connection_set_f = self.vertex_info[search_vertex_name].connections[False]
        all_pairs_of_inner_circles = []

        # branching ends
        if len(connection_set_t) == len(connection_set_f) == 2:
            for next_t_v, next_t_e in list(connection_set_t):
                this_inner_circle = self.find_pair_closing_the_path(next_t_v, next_t_e, connection_set_f,
                                                                    connection_set_t)
                if this_inner_circle:
                    # check leakage in reverse direction
                    reverse_v, reverse_e = this_inner_circle[0][1]
                    not_leak = self.find_pair_closing_the_path(reverse_v, reverse_e, connection_set_t, connection_set_f)
                    if not_leak:
                        all_pairs_of_inner_circles.extend(this_inner_circle)
            # sort pairs by average depths(?)
            all_pairs_of_inner_circles.sort(
                key=lambda x: (self.vertex_info[x[0][0]].cov + self.vertex_info[x[1][0]].cov))
            if all_pairs_of_inner_circles and return_pair_in_the_trunk_path:
                # switch nearby vertices_set
                # keep those prone to be located in the "trunk road" of the repeat
                single_pair_in_main_path = []
                if len(all_pairs_of_inner_circles) == 1:
                    for next_v, next_e in list(connection_set_t) + list(connection_set_f):
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
        # follow variables were not updated because of max_majority_copy should not be considered here
        # self.copy_to_vertex
        # self.vertex_to_copy
        # self.vertex_to_float_copy

        if not limited_vertices:
            limited_vertices = sorted(self.vertex_info)
        else:
            limited_vertices = sorted(limited_vertices)
        merged = False
        # overlap = self.__overlap if self.__overlap else 0
        while limited_vertices:
            this_vertex = limited_vertices.pop()
            for this_end in (True, False):
                connected_dict = self.vertex_info[this_vertex].connections[this_end]
                if len(connected_dict) == 1:
                    (next_vertex, next_end), this_overlap = list(connected_dict.items())[0]
                    if len(self.vertex_info[next_vertex].connections[next_end]) == 1 and this_vertex != next_vertex:
                        # reverse the names
                        merged = True
                        self.vertex_info[this_vertex].merging_history.add((next_vertex, not this_end == next_end),
                                                                          add_new_to_front=not this_end,
                                                                          reverse_the_new=False)
                        new_vertex = str(self.vertex_info[this_vertex].merging_history)

                        limited_vertices.remove(next_vertex)
                        limited_vertices.append(new_vertex)
                        # initialization
                        self.vertex_info[new_vertex] = deepcopy(self.vertex_info[this_vertex])
                        self.vertex_info[new_vertex].name = new_vertex
                        self.vertex_info[new_vertex].fastg_form_name = None
                        # modify connections
                        self.vertex_info[new_vertex].connections[this_end] \
                            = deepcopy(self.vertex_info[next_vertex].connections[not next_end])
                        if (this_vertex, not this_end) in self.vertex_info[new_vertex].connections[this_end]:
                            # forms a circle
                            overlap_x = self.vertex_info[new_vertex].connections[this_end][(this_vertex, not this_end)]
                            del self.vertex_info[new_vertex].connections[this_end][(this_vertex, not this_end)]
                            self.vertex_info[new_vertex].connections[this_end][(new_vertex, not this_end)] = overlap_x
                        for new_end in (True, False):
                            for n_n_v, n_n_e in self.vertex_info[new_vertex].connections[new_end]:
                                self.vertex_info[n_n_v].connections[n_n_e][(new_vertex, new_end)] = \
                                    self.vertex_info[new_vertex].connections[new_end][(n_n_v, n_n_e)]
                        # len & cov
                        this_len = self.vertex_info[this_vertex].len
                        next_len = self.vertex_info[next_vertex].len
                        this_cov = self.vertex_info[this_vertex].cov
                        next_cov = self.vertex_info[next_vertex].cov
                        self.vertex_info[new_vertex].len = this_len + next_len - this_overlap
                        self.vertex_info[new_vertex].cov = \
                            ((this_len - this_overlap + 1) * this_cov + (next_len - this_overlap + 1) * next_cov) \
                            / ((this_len - this_overlap + 1) + (next_len - this_overlap + 1))
                        self.vertex_info[new_vertex].seq[this_end] \
                            += self.vertex_info[next_vertex].seq[not next_end][this_overlap:]
                        self.vertex_info[new_vertex].seq[not this_end] \
                            = self.vertex_info[next_vertex].seq[next_end][:next_len - this_overlap] \
                              + self.vertex_info[this_vertex].seq[not this_end]

                        # follow variables were not updated because of max_majority_copy should not be considered here
                        # self.copy_to_vertex
                        # self.vertex_to_copy
                        # self.vertex_to_float_copy
                        # average_cov = this_cov / self.vertex_to_float_copy[this_vertex]
                        # this_float_copy = self.vertex_info[new_vertex].cov / average_cov
                        # this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_copy)
                        # self.vertex_to_float_copy[new_vertex] = this_float_copy
                        # self.vertex_to_copy[new_vertex] = this_copy
                        # if this_copy not in self.copy_to_vertex:
                        #     self.copy_to_vertex[this_copy] = set()
                        # self.copy_to_vertex[this_copy].add(new_vertex)

                        # tags
                        if copy_tags:
                            if "tags" in self.vertex_info[next_vertex].other_attr:
                                if "tags" not in self.vertex_info[new_vertex].other_attr:
                                    self.vertex_info[new_vertex].other_attr["tags"] = \
                                        deepcopy(self.vertex_info[next_vertex].other_attr["tags"])
                                else:
                                    for db_n in self.vertex_info[next_vertex].other_attr["tags"]:
                                        if db_n not in self.vertex_info[new_vertex].other_attr["tags"]:
                                            self.vertex_info[new_vertex].other_attr["tags"][db_n] \
                                                = deepcopy(self.vertex_info[next_vertex].other_attr["tags"][db_n])
                                        else:
                                            self.vertex_info[new_vertex].other_attr["tags"][db_n] \
                                                |= self.vertex_info[next_vertex].other_attr["tags"][db_n]
                            if "weight" in self.vertex_info[next_vertex].other_attr:
                                if "weight" not in self.vertex_info[new_vertex].other_attr:
                                    self.vertex_info[new_vertex].other_attr["weight"] \
                                        = deepcopy(self.vertex_info[next_vertex].other_attr["weight"])
                                else:
                                    for db_n in self.vertex_info[next_vertex].other_attr["weight"]:
                                        if db_n not in self.vertex_info[new_vertex].other_attr["weight"]:
                                            self.vertex_info[new_vertex].other_attr["weight"][db_n] \
                                                = self.vertex_info[next_vertex].other_attr["weight"][db_n]
                                        else:
                                            self.vertex_info[new_vertex].other_attr["weight"][db_n] \
                                                += self.vertex_info[next_vertex].other_attr["weight"][db_n]
                            for db_n in self.tagged_vertices:
                                if this_vertex in self.tagged_vertices[db_n]:
                                    self.tagged_vertices[db_n].add(new_vertex)
                                    self.tagged_vertices[db_n].remove(this_vertex)
                                if next_vertex in self.tagged_vertices[db_n]:
                                    self.tagged_vertices[db_n].add(new_vertex)
                                    self.tagged_vertices[db_n].remove(next_vertex)
                        self.remove_vertex([this_vertex, next_vertex], update_cluster=False)
                        break
        self.update_vertex_clusters()
        return merged

    def estimate_copy_and_depth_by_cov(self, limited_vertices=None, given_average_cov=None, mode="embplant_pt",
                                       min_sigma=0.1, re_initialize=False, log_handler=None, verbose=True, debug=False):
        """
        :param limited_vertices:
        :param given_average_cov:
        :param mode:
        :param min_sigma: when only one sample
        :param re_initialize:
        :param log_handler:
        :param verbose:
        :param debug:
        :return:
        """
        # overlap = self.__overlap if self.__overlap else 0
        # those are all empirical values
        # TODO: GetOrganelle need a better algorithm for target filtering
        # if mode == "embplant_pt":
        #     max_majority_copy = 2
        # elif mode == "other_pt":
        #     max_majority_copy = 10
        # elif mode == "embplant_mt":
        #     max_majority_copy = 4
        # elif mode == "embplant_nr":
        #     max_majority_copy = 2
        # elif mode == "animal_mt":
        #     # the difference between mt and nucl are usually not that large,
        #     # making the upper boundary useless for excluding non-target but excluding target by mistake
        #     max_majority_copy = 1000
        # elif mode == "fungus_mt":
        #     max_majority_copy = 8
        # elif mode == "fungus_nr":
        #     max_majority_copy = 4
        # elif mode == "all":
        #     max_majority_copy = 100
        # else:
        #     max_majority_copy = 100

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

        cov_ls = []
        len_ls = []
        if not given_average_cov:
            previous_val = {0.}
            new_val = -1.
            new_std = -1.
            # arbitrary setting, without influence of limited_vertices
            min_average_depth = 0.9 * min([self.vertex_info[vertex_n].cov for vertex_n in self.vertex_info])
            while round(new_val, 5) not in previous_val:
                previous_val.add(round(new_val, 5))
                # estimate baseline depth
                total_product = 0.
                total_len = 0
                for vertex_name in limited_vertices:
                    # do we need to exclude the overlap?
                    this_len = self.vertex_info[vertex_name].len * self.vertex_to_copy.get(vertex_name, 1)
                    len_ls.append(this_len)
                    this_cov = self.vertex_info[vertex_name].cov / self.vertex_to_copy.get(vertex_name, 1)
                    cov_ls.append(this_cov)
                    total_len += this_len
                    total_product += this_len * this_cov
                # new_val = total_product / total_len
                new_val = max(total_product / total_len, min_average_depth)
                new_std = (sum([_w * (new_val - _c) ** 2 for _c, _w in zip(cov_ls, len_ls)]) / sum(len_ls)) ** 0.5
                # (sum(len_ls) * (len(cov_ls) - 1) / len(cov_ls))) ** 0.5
                # print("new val: ", new_val)
                # adjust this_copy according to new baseline depth
                for vertex_name in self.vertex_info:
                    if vertex_name in self.vertex_to_copy:
                        old_copy = self.vertex_to_copy[vertex_name]
                        self.copy_to_vertex[old_copy].remove(vertex_name)
                        if not self.copy_to_vertex[old_copy]:
                            del self.copy_to_vertex[old_copy]
                    this_float_copy = self.vertex_info[vertex_name].cov / new_val
                    # this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_copy)
                    this_copy = max(1, int(round(this_float_copy, 0)))
                    self.vertex_to_float_copy[vertex_name] = this_float_copy
                    self.vertex_to_copy[vertex_name] = this_copy
                    if this_copy not in self.copy_to_vertex:
                        self.copy_to_vertex[this_copy] = set()
                    self.copy_to_vertex[this_copy].add(vertex_name)
            if debug or verbose:
                cov_str = " kmer-coverage: " if bool(self.uni_overlap) else " coverage: "
                if log_handler:
                    log_handler.info("updating average " + mode + cov_str + str(round(new_val, 2)))
                else:
                    sys.stdout.write("updating average " + mode + cov_str + str(round(new_val, 2)) + "\n")
            # print("return ", new_val)
            return new_val, new_val * min_sigma if len(limited_vertices) == 1 and new_std == 0. else new_std
        else:
            # adjust this_copy according to user-defined depth
            for vertex_name in limited_vertices:
                if vertex_name in self.vertex_to_copy:
                    old_copy = self.vertex_to_copy[vertex_name]
                    self.copy_to_vertex[old_copy].remove(vertex_name)
                    if not self.copy_to_vertex[old_copy]:
                        del self.copy_to_vertex[old_copy]
                this_float_copy = self.vertex_info[vertex_name].cov / given_average_cov
                # this_copy = min(max(1, int(round(this_float_copy, 0))), max_majority_copy)
                this_copy = max(1, int(round(this_float_copy, 0)))
                self.vertex_to_float_copy[vertex_name] = this_float_copy
                self.vertex_to_copy[vertex_name] = this_copy
                if this_copy not in self.copy_to_vertex:
                    self.copy_to_vertex[this_copy] = set()
                self.copy_to_vertex[this_copy].add(vertex_name)
                cov_ls.append(self.vertex_info[vertex_name].cov / this_copy)
                len_ls.append(self.vertex_info[vertex_name].len * this_copy)
            new_std = (sum([_w * (given_average_cov - _c) ** 2 for _c, _w in zip(cov_ls, len_ls)]) / sum(len_ls)) ** 0.5
            return given_average_cov, \
                   given_average_cov * min_sigma if len(limited_vertices) == 1 and new_std == 0. else new_std

    # def estimate_copy_and_depth_precisely_using_multinomial(
    #         self, expected_average_cov,  # broken_graph_allowed=False,
    #         verbose=False, log_handler=None, debug=False,
    #         target_name_for_log="target", n_iterations=None):
    #     """
    #     Currently problematic because of
    #     frequently reporting
    #     Exception: @error: Solution Not Found
    #
    #     :param expected_average_cov:
    #     :param verbose:
    #     :param log_handler:
    #     :param debug:
    #     :param target_name_for_log:
    #     :return:
    #     """
    #     # TODO: to test whether it's better than least-square
    #
    #     def get_formula(from_vertex, from_end, back_to_vertex, back_to_end, here_record_ends):
    #         result_form = vertex_to_symbols[from_vertex]
    #         here_record_ends.add((from_vertex, from_end))
    #         # if back_to_vertex ~ from_vertex (from_vertex == back_to_vertex) form a loop, skipped
    #         if from_vertex != back_to_vertex:
    #             for next_v, next_e in self.vertex_info[from_vertex].connections[from_end]:
    #                 # if next_v ~ from_vertex (next_v == from_vertex) form a loop, add a pseudo vertex
    #                 if (next_v, next_e) == (from_vertex, not from_end):
    #                     # skip every self-loop 2020-06-23
    #                     # pseudo_self_circle_str = "P" + from_vertex
    #                     # if pseudo_self_circle_str not in extra_str_to_symbol_m2:
    #                     #     extra_str_to_symbol_m2[pseudo_self_circle_str] = Symbol(pseudo_self_circle_str, integer=True)
    #                     #     extra_symbol_to_str_m2[extra_str_to_symbol_m2[pseudo_self_circle_str]] = pseudo_self_circle_str
    #                     # result_form -= (extra_str_to_symbol_m2[pseudo_self_circle_str] - 1)
    #                     pass
    #                 # elif (next_v, next_e) != (back_to_vertex, back_to_end):
    #                 elif (next_v, next_e) not in here_record_ends:
    #                     result_form -= get_formula(next_v, next_e, from_vertex, from_end, here_record_ends)
    #         return result_form
    #
    #     # # for compatibility between scipy and sympy
    #     # def least_square_function_v(x):
    #     #     return least_square_function(*tuple(x))
    #     #
    #     # """ create constraints by creating inequations: the copy of every contig has to be >= 1 """
    #     #
    #     # def constraint_min_function(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     min_copy = np.array([1.001] * len(all_v_symbols) +
    #     #                         [1.001] * len(extra_symbol_to_str_m1) +
    #     #                         [2.001] * len(extra_symbol_to_str_m2))
    #     #     # effect: expression_array >= int(min_copy)
    #     #     return expression_array - min_copy
    #     #
    #     # def constraint_min_function_for_customized_brute(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     min_copy = np.array([1.0] * len(all_v_symbols) +
    #     #                         [1.0] * len(extra_symbol_to_str_m1) +
    #     #                         [2.0] * len(extra_symbol_to_str_m2))
    #     #     # effect: expression_array >= min_copy
    #     #     return expression_array - min_copy
    #
    #     def constraint_min_function_for_gekko(g_vars):
    #         subs_tuples = [(symb_used_, Symbol("g_vars[" + str(go_sym) + "]"))
    #                        for go_sym, symb_used_ in enumerate(free_copy_variables)]
    #         expression_array = [copy_solution[this_sym].subs(subs_tuples) for this_sym in all_symbols]
    #         min_copy = [1] * len(all_v_symbols) + \
    #                    [1] * len(extra_symbol_to_str_m1) + \
    #                    [2] * len(extra_symbol_to_str_m2)
    #         # effect: expression_array >= min_copy
    #         expression = []
    #         if verbose or debug:
    #             for e, c in zip(expression_array, min_copy):
    #                 expression.append(eval(str(e) + ">=" + str(c)))
    #                 log_handler.info("  constraint: " + str(e) + ">=" + str(c))
    #         else:
    #             for e, c in zip(expression_array, min_copy):
    #                 expression.append(eval(str(e) + ">=" + str(c)))
    #         expression = [expr for expr in expression if not isinstance(expr, bool)]
    #         return expression
    #
    #     # def constraint_max_function(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     max_copy = np.array([expected_average_cov] * len(all_v_symbols) +
    #     #                         [expected_average_cov] * len(extra_symbol_to_str_m1) +
    #     #                         [expected_average_cov * 2] * len(extra_symbol_to_str_m2))
    #     #     # effect: expression_array <= max_copy
    #     #     return max_copy - expression_array
    #     #
    #     # def constraint_int_function(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     # diff = np.array([0] * len(all_symbols))
    #     #     return sum([abs(every_copy - int(every_copy)) for every_copy in expression_array])
    #     #
    #     # def minimize_brute_force(func, range_list, constraint_list, round_digit=4, display_p=True,
    #     #                          in_log_handler=log_handler):
    #     #     # time0 = time.time()
    #     #     best_fun_val = inf
    #     #     best_para_val = []
    #     #     count_round = 0
    #     #     count_valid = 0
    #     #     for value_set in product(*[list(this_range) for this_range in range_list]):
    #     #         count_round += 1
    #     #         is_valid_set = True
    #     #         for cons in constraint_list:
    #     #             if cons["type"] == "ineq":
    #     #                 try:
    #     #                     if (cons["fun"](value_set) < 0).any():
    #     #                         is_valid_set = False
    #     #                         # if in_log_handler and (debug or display_p):
    #     #                         #     in_log_handler.info("value_set={} ; illegal ineq constraints".format(value_set))
    #     #                         break
    #     #                 except TypeError:
    #     #                     # if in_log_handler and (debug or display_p):
    #     #                     #     in_log_handler.info("value_set={} ; illegal ineq constraints".format(value_set))
    #     #                     is_valid_set = False
    #     #                     break
    #     #             elif cons["type"] == "eq":
    #     #                 try:
    #     #                     if cons["fun"](value_set) != 0:
    #     #                         is_valid_set = False
    #     #                         # if in_log_handler and (debug or display_p):
    #     #                         #     in_log_handler.info("value_set={} ; illegal eq constraints".format(value_set))
    #     #                         break
    #     #                 except TypeError:
    #     #                     # if in_log_handler and (debug or display_p):
    #     #                     #     in_log_handler.info("value_set={} ; illegal eq constraints".format(value_set))
    #     #                     is_valid_set = False
    #     #                     break
    #     #         if not is_valid_set:
    #     #             continue
    #     #         count_valid += 1
    #     #         this_fun_val = func(value_set)
    #     #         if in_log_handler:
    #     #             if debug or display_p:
    #     #                 in_log_handler.info("value_set={} ; fun_val={}".format(value_set, this_fun_val))
    #     #         this_fun_val = round(this_fun_val, round_digit)
    #     #         if this_fun_val < best_fun_val:
    #     #             best_para_val = [value_set]
    #     #             best_fun_val = this_fun_val
    #     #         elif this_fun_val == best_fun_val:
    #     #             best_para_val.append(value_set)
    #     #         else:
    #     #             pass
    #     #     if in_log_handler:
    #     #         if debug or display_p:
    #     #             in_log_handler.info("Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round))
    #     #             in_log_handler.info("Brute best function value: " + str(best_fun_val))
    #     #         if debug:
    #     #             in_log_handler.info("Best solution: " + str(best_para_val))
    #     #     else:
    #     #         if debug or display_p:
    #     #             sys.stdout.write(
    #     #                 "Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round) + "\n")
    #     #             sys.stdout.write("Brute best function value: " + str(best_fun_val) + "\n")
    #     #         if debug:
    #     #             sys.stdout.write("Best solution: " + str(best_para_val) + "\n")
    #     #     return best_para_val
    #
    #     vertices_list = sorted(self.vertex_info)
    #     if len(vertices_list) == 1:
    #         cov_ = self.vertex_info[vertices_list[0]].cov
    #         # 2022-12-15, remove return_new_graph
    #         # if return_new_graphs:
    #         return [{"graph": deepcopy(self), "cov": cov_}]
    #         # else:
    #         #     if log_handler:
    #         #         log_handler.info("Average " + target_name_for_log + " kmer-coverage = " + str(round(cov_, 2)))
    #         #     else:
    #         #         sys.stdout.write(
    #         #             "Average " + target_name_for_log + " kmer-coverage = " + str(round(cov_, 2)) + "\n")
    #         #     return
    #
    #     # reduce expected_average_cov to reduce computational burden
    #     all_coverages = [self.vertex_info[v_name].cov for v_name in vertices_list]
    #     # max_contig_multiplicity = \
    #     #     min(max_contig_multiplicity, int(2 * math.ceil(max(all_coverages) / min(all_coverages))))
    #     # if verbose:
    #     #     if log_handler:
    #     #         log_handler.info("Maximum multiplicity: " + str(max_contig_multiplicity))
    #     #     else:
    #     #         sys.stdout.write("Maximum multiplicity: " + str(max_contig_multiplicity) + "\n")
    #
    #     """ create constraints by creating multivariate equations """
    #     vertex_to_symbols = {vertex_name: Symbol("V" + vertex_name, integer=True)  # positive=True)
    #                          for vertex_name in vertices_list}
    #     symbols_to_vertex = {vertex_to_symbols[vertex_name]: vertex_name for vertex_name in vertices_list}
    #     extra_str_to_symbol_m1 = {}
    #     extra_str_to_symbol_m2 = {}
    #     extra_symbol_to_str_m1 = {}
    #     extra_symbol_to_str_m2 = {}
    #     extra_symbol_initial_values = {}
    #     formulae = []
    #     recorded_ends = set()
    #     for vertex_name in vertices_list:
    #         for this_end in (True, False):
    #             if (vertex_name, this_end) not in recorded_ends:
    #                 recorded_ends.add((vertex_name, this_end))
    #                 if self.vertex_info[vertex_name].connections[this_end]:
    #                     this_formula = vertex_to_symbols[vertex_name]
    #                     formulized = False
    #                     for n_v, n_e in self.vertex_info[vertex_name].connections[this_end]:
    #                         if (n_v, n_e) not in recorded_ends:
    #                             # if n_v in vertices_set:
    #                             # recorded_ends.add((n_v, n_e))
    #                             try:
    #                                 this_formula -= get_formula(n_v, n_e, vertex_name, this_end, recorded_ends)
    #                                 formulized = True
    #                                 # if verbose:
    #                                 #     if log_handler:
    #                                 #         log_handler.info("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                 #                          vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 #                          str(this_formula))
    #                                 #     else:
    #                                 #         sys.stdout.write("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                 #                          vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 #                          str(this_formula)+"\n")
    #                             except RecursionError:
    #                                 if log_handler:
    #                                     log_handler.warning("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                                         vertex_name + ECHO_DIRECTION[this_end] + " failed!")
    #                                 else:
    #                                     sys.stdout.write("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                                      vertex_name + ECHO_DIRECTION[this_end] + " failed!\n")
    #                                 raise ProcessingGraphFailed("RecursionError!")
    #                     if verbose:
    #                         if log_handler:
    #                             log_handler.info(
    #                                 "formulating for: " + vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 str(this_formula))
    #                         else:
    #                             sys.stdout.write(
    #                                 "formulating for: " + vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 str(this_formula) + "\n")
    #                     if formulized:
    #                         formulae.append(this_formula)
    #                 # 2022-12-13 remove this restriction
    #                 #            because we have a reduce_list_with_gcd for all graph component
    #                 # elif broken_graph_allowed:
    #                 #     # Extra limitation to force terminal vertex to have only one copy, to avoid over-estimation
    #                 #     # Under-estimation would not be a problem here,
    #                 #     # because the True-multiple-copy vertex would simply have no other connections,
    #                 #     # or failed in the following estimation if it does
    #                 #     formulae.append(vertex_to_symbols[vertex_name] - 1)
    #
    #     # add self-loop formulae
    #     self_loop_v = set()
    #     for vertex_name in vertices_list:
    #         if self.vertex_info[vertex_name].is_self_loop():
    #             self_loop_v.add(vertex_name)
    #             if log_handler:
    #                 log_handler.warning("Self-loop contig detected: Vertex_" + vertex_name)
    #             pseudo_self_loop_str = "P" + vertex_name
    #             if pseudo_self_loop_str not in extra_str_to_symbol_m1:
    #                 extra_str_to_symbol_m1[pseudo_self_loop_str] = Symbol(pseudo_self_loop_str, integer=True)
    #                 extra_symbol_to_str_m1[extra_str_to_symbol_m1[pseudo_self_loop_str]] = pseudo_self_loop_str
    #             this_formula = vertex_to_symbols[vertex_name] - extra_str_to_symbol_m1[pseudo_self_loop_str]
    #             extra_symbol_initial_values[extra_str_to_symbol_m1[pseudo_self_loop_str]] = \
    #                 self.vertex_to_copy[vertex_name]
    #             formulae.append(this_formula)
    #             if verbose:
    #                 if log_handler:
    #                     log_handler.info(
    #                         "formulating for: " + vertex_name + ECHO_DIRECTION[True] + ": " + str(this_formula))
    #                 else:
    #                     sys.stdout.write(
    #                         "formulating for: " + vertex_name + ECHO_DIRECTION[True] + ": " + str(this_formula) + "\n")
    #
    #     # add following extra limitation
    #     # set cov_sequential_repeat = x*near_by_cov, x is an integer
    #     for vertex_name in vertices_list:
    #         single_pair_in_the_trunk_path = self.is_sequential_repeat(vertex_name)
    #         if single_pair_in_the_trunk_path:
    #             (from_v, from_e), (to_v, to_e) = single_pair_in_the_trunk_path
    #             # from_v and to_v are already in the "trunk path", if they are the same,
    #             # the graph is like two circles sharing the same sequential repeat, no need to add this limitation
    #             if from_v != to_v:
    #                 new_str = "E" + str(len(extra_str_to_symbol_m1) + len(extra_str_to_symbol_m2))
    #                 if vertex_name in self_loop_v:
    #                     # self-loop vertex is allowed to have the multiplicity of 1
    #                     extra_str_to_symbol_m1[new_str] = Symbol(new_str, integer=True)
    #                     extra_symbol_to_str_m1[extra_str_to_symbol_m1[new_str]] = new_str
    #                     this_formula = vertex_to_symbols[vertex_name] - \
    #                                    vertex_to_symbols[from_v] * extra_str_to_symbol_m1[new_str]
    #                     extra_symbol_initial_values[extra_str_to_symbol_m1[new_str]] = \
    #                         round(self.vertex_to_float_copy[vertex_name] / self.vertex_to_float_copy[from_v])
    #                 else:
    #                     extra_str_to_symbol_m2[new_str] = Symbol(new_str, integer=True)
    #                     extra_symbol_to_str_m2[extra_str_to_symbol_m2[new_str]] = new_str
    #                     this_formula = vertex_to_symbols[vertex_name] - \
    #                                    vertex_to_symbols[from_v] * extra_str_to_symbol_m2[new_str]
    #                     extra_symbol_initial_values[extra_str_to_symbol_m2[new_str]] = \
    #                         round(self.vertex_to_float_copy[vertex_name] / self.vertex_to_float_copy[from_v])
    #                 formulae.append(this_formula)
    #                 if verbose:
    #                     if log_handler:
    #                         log_handler.info("formulating for: " + vertex_name + ": " + str(this_formula))
    #                     else:
    #                         sys.stdout.write("formulating for: " + vertex_name + ": " + str(this_formula) + "\n")
    #
    #     all_v_symbols = list(symbols_to_vertex)
    #     all_symbols = all_v_symbols + list(extra_symbol_to_str_m1) + list(extra_symbol_to_str_m2)
    #     if verbose or debug:
    #         if log_handler:
    #             log_handler.info("formulae: " + str(formulae))
    #         else:
    #             sys.stdout.write("formulae: " + str(formulae) + "\n")
    #     # solve the equations
    #     copy_solution = solve(formulae, all_v_symbols)
    #
    #     copy_solution = copy_solution if copy_solution else {}
    #     if type(copy_solution) == list:  # delete 0 containing set, even for self-loop vertex
    #         go_solution = 0
    #         while go_solution < len(copy_solution):
    #             if 0 in set(copy_solution[go_solution].values()):
    #                 del copy_solution[go_solution]
    #             else:
    #                 go_solution += 1
    #     if not copy_solution:
    #         raise ProcessingGraphFailed("Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (1)!")
    #     elif type(copy_solution) == list:
    #         if len(copy_solution) > 2:
    #             raise ProcessingGraphFailed("Incomplete/Complicated " + target_name_for_log + " graph (2)!")
    #         else:
    #             copy_solution = copy_solution[0]
    #
    #     free_copy_variables = list()
    #     for symbol_used in all_symbols:
    #         if symbol_used not in copy_solution:
    #             free_copy_variables.append(symbol_used)
    #             copy_solution[symbol_used] = symbol_used
    #     if verbose:
    #         if log_handler:
    #             log_handler.info("copy equations: " + str(copy_solution))
    #             log_handler.info("free variables: " + str(free_copy_variables))
    #         else:
    #             sys.stdout.write("copy equations: " + str(copy_solution) + "\n")
    #             sys.stdout.write("free variables: " + str(free_copy_variables) + "\n")
    #
    #     # """ minimizing equation-based copy's deviations from coverage-based copy values """
    #     # least_square_expr = 0
    #     # for symbol_used in all_v_symbols:
    #     #     # least_square_expr += copy_solution[symbol_used]
    #     #     this_vertex = symbols_to_vertex[symbol_used]
    #     #     this_copy = self.vertex_to_float_copy[this_vertex]
    #     #     least_square_expr += (copy_solution[symbol_used] - this_copy) ** 2  # * self.vertex_info[this_vertex]["len"]
    #     # least_square_function = lambdify(args=free_copy_variables, expr=least_square_expr)
    #
    #     if free_copy_variables:
    #         """ Maximize the likelihood of the multinomial distribution of kmers (kmer_cov * contig_len)"""
    #         m = GEKKO(remote=False)
    #         g_vars = m.Array(m.Var,
    #                          len(free_copy_variables),
    #                          lb=1,
    #                          ub=int(4 * math.ceil(max(all_coverages) / min(all_coverages))),
    #                          integer=True)
    #         # initialize free variables
    #         for go_sym, symbol_used in enumerate(free_copy_variables):
    #             if symbol_used in symbols_to_vertex and symbols_to_vertex[symbol_used] in self.vertex_to_copy:
    #                 g_vars[go_sym].value = self.vertex_to_copy[symbols_to_vertex[symbol_used]]
    #             elif symbol_used in extra_symbol_initial_values:
    #                 g_vars[go_sym].value = extra_symbol_initial_values[symbol_used]
    #         replacements = [(symbol_used, Symbol("g_vars[" + str(go_sym) + "]"))
    #                         for go_sym, symbol_used in enumerate(free_copy_variables)]
    #         # account for the influence of the overlap
    #         total_len = 0
    #         multinomial_loglike_list = []
    #         v_to_len = {}
    #         v_to_copy = {}
    #         v_to_real_len = {}
    #         all_obs = []
    #         if self.__uni_overlap:
    #             for symbol_used in all_v_symbols:
    #                 this_vertex = symbols_to_vertex[symbol_used]
    #                 v_to_real_len[this_vertex] = self.vertex_info[this_vertex].len - self.__uni_overlap
    #                 v_to_copy[this_vertex] = eval(str(copy_solution[symbol_used].subs(replacements)))
    #                 v_to_len[this_vertex] = v_to_copy[this_vertex] * v_to_real_len[this_vertex]
    #                 total_len += v_to_len[this_vertex]
    #             for symbol_used in all_v_symbols:
    #                 this_vertex = symbols_to_vertex[symbol_used]
    #                 prob = v_to_len[this_vertex] / total_len
    #                 obs = self.vertex_info[this_vertex].cov * v_to_real_len[this_vertex]
    #                 multinomial_loglike_list.append(m.log(prob) * obs)
    #                 all_obs.append(obs)
    #                 if verbose:
    #                     if log_handler:
    #                         log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
    #                     else:
    #                         sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
    #         else:
    #             for symbol_used in all_v_symbols:
    #                 this_vertex = symbols_to_vertex[symbol_used]
    #                 overlaps = [_ovl
    #                             for _strand in (True, False)
    #                             for _next, _ovl in self.vertex_info[this_vertex].connections[_strand].items()]
    #                 approximate_overlap = average_np_free(overlaps)
    #                 v_to_real_len[this_vertex] = self.vertex_info[this_vertex].len - approximate_overlap
    #                 v_to_copy[this_vertex] = eval(str(copy_solution[symbol_used].subs(replacements)))
    #                 v_to_len[this_vertex] = v_to_copy[this_vertex] * v_to_real_len[this_vertex]
    #                 total_len += v_to_len[this_vertex]
    #             for symbol_used in all_v_symbols:
    #                 this_vertex = symbols_to_vertex[symbol_used]
    #                 prob = v_to_len[this_vertex] / total_len
    #                 obs = self.vertex_info[this_vertex].cov * v_to_real_len[this_vertex]
    #                 multinomial_loglike_list.append(m.log(prob) * obs)
    #                 all_obs.append(obs)
    #                 if verbose:
    #                     if log_handler:
    #                         log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
    #                     else:
    #                         sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
    #         """extra arbitrary restriction to avoid over inflation of copies"""
    #         sum_obs = sum(all_obs)
    #         multinomial_loglike_list.append(-abs(sum_obs / expected_average_cov - total_len))
    #         # """extra restriction to constraint the integer solution for dependant variables"""
    #         # for symbol_used in all_v_symbols:
    #         #     this_vertex = symbols_to_vertex[symbol_used]
    #         #     multinomial_loglike_list.append(sum_obs * (v_to_copy[this_vertex] - int(v_to_copy[this_vertex])) ** 2)
    #         """generate the expression"""
    #         # multinomial_loglike_expr = m.sum(multinomial_loglike_list) will lead to No solution error
    #         multinomial_loglike_expr = sum(multinomial_loglike_list)
    #         exp_str_len = len(str(multinomial_loglike_expr))
    #         if exp_str_len > 15000:  # not allowed by Gekko:APM
    #             num_blocks = math.ceil(exp_str_len / 10000.)
    #             block_size = math.ceil(len(multinomial_loglike_list) / float(num_blocks))
    #             block_list = []
    #             for g_b in range(num_blocks):
    #                 block_list.append(sum(multinomial_loglike_list[g_b * block_size: (g_b + 1)* block_size]))
    #             multinomial_loglike_expr = m.sum(block_list)
    #         # for symbol_used in all_v_symbols:
    #         #     this_vertex = symbols_to_vertex[symbol_used]
    #         #     total_len += eval(str(copy_solution[symbol_used].subs(replacements))) * self.vertex_info[this_vertex].len
    #         # multinomial_like_expr = 0
    #         # for symbol_used in all_v_symbols:
    #         #     this_vertex = symbols_to_vertex[symbol_used]
    #         #     prob = eval(str(copy_solution[symbol_used].subs(replacements))) \
    #         #            * self.vertex_info[this_vertex].len / total_len
    #         #     obs = self.vertex_info[this_vertex].cov * self.vertex_info[this_vertex].len
    #         #     multinomial_like_expr += m.log(prob) * obs
    #         #     if verbose:
    #         #         if log_handler:
    #         #             log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
    #         #         else:
    #         #             sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
    #
    #         # multinomial_like_function = lambdify(args=free_copy_variables, expr=multinomial_like_expr)
    #         m.Equations(constraint_min_function_for_gekko(g_vars))
    #         m.Maximize(multinomial_loglike_expr)
    #         # m.Minimize(least_square_function_v(g_vars))
    #         # 1 for APOPT, 2 for BPOPT, 3 for IPOPT, 0 for all available solvers
    #         # here only 1 and 3 are available
    #         m.options.SOLVER = 1
    #         # setting empirical options
    #         # 5000 costs ~ 150 sec
    #         if n_iterations is None:
    #             n_high_copy = sum([math.log2(self.vertex_to_float_copy[_v])
    #                                for _v in self.vertex_info if self.vertex_to_float_copy[_v] > 2])
    #             n_iterations = 500 + int(len(self.vertex_info) * n_high_copy)
    #         if verbose or debug:
    #             log_handler.info("setting n_iterations=" + str(n_iterations))
    #         m.solver_options = ['minlp_maximum_iterations ' + str(n_iterations),
    #                             # minlp iterations with integer solution
    #                             'minlp_max_iter_with_int_sol ' + str(n_iterations),
    #                             # treat minlp as nlp
    #                             'minlp_as_nlp 0',
    #                             # nlp sub-problem max iterations
    #                             'nlp_maximum_iterations ' + str(n_iterations),
    #                             # 1 = depth first, 2 = breadth first
    #                             'minlp_branch_method 2',
    #                             # maximum deviation from whole number
    #                             'minlp_integer_tol 1.0e-6',
    #                             # covergence tolerance
    #                             'minlp_gap_tol 1.0e-6']
    #         if debug or verbose:
    #             m.solve()
    #         else:
    #             m.solve(disp=False)
    #         # print([x.value[0] for x in g_vars])
    #         copy_results = list([x.value[0] for x in g_vars])
    #
    #         # # for safe running
    #         # if len(free_copy_variables) > 10:
    #         #     raise ProcessingGraphFailed("Free variable > 10 is not accepted yet!")
    #         #
    #         # if expected_average_cov ** len(free_copy_variables) < 5E6:
    #         #     # sometimes, SLSQP ignores bounds and constraints
    #         #     copy_results = minimize_brute_force(
    #         #         func=least_square_function_v, range_list=[range(1, expected_average_cov + 1)] * len(free_copy_variables),
    #         #         constraint_list=({'type': 'ineq', 'fun': constraint_min_function_for_customized_brute},
    #         #                          {'type': 'eq', 'fun': constraint_int_function},
    #         #                          {'type': 'ineq', 'fun': constraint_max_function}),
    #         #         display_p=verbose)
    #         # else:
    #         #     constraints = ({'type': 'ineq', 'fun': constraint_min_function},
    #         #                    {'type': 'eq', 'fun': constraint_int_function},
    #         #                    {'type': 'ineq', 'fun': constraint_max_function})
    #         #     copy_results = set()
    #         #     best_fun = inf
    #         #     opt = {'disp': verbose, "maxiter": 100}
    #         #     for initial_copy in range(expected_average_cov * 2 + 1):
    #         #         if initial_copy < expected_average_cov:
    #         #             initials = np.array([initial_copy + 1] * len(free_copy_variables))
    #         #         elif initial_copy < expected_average_cov * 2:
    #         #             initials = np.array([random.randint(1, expected_average_cov)] * len(free_copy_variables))
    #         #         else:
    #         #             initials = np.array([self.vertex_to_copy.get(symbols_to_vertex.get(symb, False), 2)
    #         #                                  for symb in free_copy_variables])
    #         #         bounds = [(1, expected_average_cov) for foo in range(len(free_copy_variables))]
    #         #         try:
    #         #             copy_result = optimize.minimize(fun=least_square_function_v, x0=initials, jac=False,
    #         #                                             method='SLSQP', bounds=bounds, constraints=constraints, options=opt)
    #         #         except Exception:
    #         #             continue
    #         #         if copy_result.fun < best_fun:
    #         #             best_fun = round(copy_result.fun, 2)
    #         #             copy_results = {tuple(copy_result.x)}
    #         #         elif copy_result.fun == best_fun:
    #         #             copy_results.add(tuple(copy_result.x))
    #         #         else:
    #         #             pass
    #         #     if debug or verbose:
    #         #         if log_handler:
    #         #             log_handler.info("Best function value: " + str(best_fun))
    #         #         else:
    #         #             sys.stdout.write("Best function value: " + str(best_fun) + "\n")
    #         if verbose or debug:
    #             if log_handler:
    #                 log_handler.info("Copy results: " + str(copy_results))
    #             else:
    #                 sys.stdout.write("Copy results: " + str(copy_results) + "\n")
    #         # if len(copy_results) == 1:
    #         #     copy_results = list(copy_results)
    #         # elif len(copy_results) > 1:
    #         #     # draftly sort results by freedom vertices_set
    #         #     copy_results = sorted(copy_results, key=lambda
    #         #         x: sum([(x[go_sym] - self.vertex_to_float_copy[symbols_to_vertex[symb_used]]) ** 2
    #         #                 for go_sym, symb_used in enumerate(free_copy_variables)
    #         #                 if symb_used in symbols_to_vertex]))
    #         # else:
    #         #     raise ProcessingGraphFailed("Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (3)!")
    #     else:
    #         copy_results = []
    #
    #     # if return_new_graphs:
    #     # """ produce all possible vertex copy combinations """
    #     final_results = []
    #     all_copy_sets = set()
    #     # maybe no more multiple results since 2022-12 gekko update
    #     for go_res, copy_result in enumerate([copy_results]):
    #         free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
    #                                     for i, this_copy in enumerate(copy_result)}
    #
    #         """ simplify copy values """
    #         # 2020-02-22 added to avoid multiplicities res such as: [4, 8, 4]
    #         # 2022-12-15 add cluster info to simplify by graph components when the graph is broken
    #         all_copies = []
    #         v_to_cid = {}
    #         for go_id, this_symbol in enumerate(all_v_symbols):
    #             vertex_name = symbols_to_vertex[this_symbol]
    #             v_to_cid[vertex_name] = go_id
    #             this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
    #             if this_copy <= 0:
    #                 raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
    #             all_copies.append(this_copy)
    #         if len(self.vertex_clusters) == 1:
    #             if len(all_copies) == 0:
    #                 raise ProcessingGraphFailed(
    #                     "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
    #             elif len(all_copies) == 1:
    #                 all_copies = [1]
    #             elif min(all_copies) == 1:
    #                 pass
    #             else:
    #                 new_all_copies = reduce_list_with_gcd(all_copies)
    #                 if verbose and new_all_copies != all_copies:
    #                     if log_handler:
    #                         log_handler.info("Estimated copies: " + str(all_copies))
    #                         log_handler.info("Reduced copies: " + str(new_all_copies))
    #                     else:
    #                         sys.stdout.write("Estimated copies: " + str(all_copies) + "\n")
    #                         sys.stdout.write("Reduced copies: " + str(new_all_copies) + "\n")
    #                 all_copies = new_all_copies
    #         else:
    #             for v_cluster in self.vertex_clusters:
    #                 ids = [v_to_cid[_v] for _v in v_cluster]
    #                 component_copies = [all_copies[_id] for _id in ids]
    #                 if len(component_copies) == 0:
    #                     raise ProcessingGraphFailed(
    #                         "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
    #                 elif len(component_copies) == 1:
    #                     component_copies = [1]
    #                 elif min(component_copies) == 1:
    #                     pass
    #                 else:
    #                     new_comp_copies = reduce_list_with_gcd(component_copies)
    #                     if verbose and new_comp_copies != component_copies:
    #                         if log_handler:
    #                             log_handler.info("Estimated copies: " + str(component_copies))
    #                             log_handler.info("Reduced copies: " + str(new_comp_copies))
    #                         else:
    #                             sys.stdout.write("Estimated copies: " + str(component_copies) + "\n")
    #                             sys.stdout.write("Reduced copies: " + str(new_comp_copies) + "\n")
    #                     component_copies = new_comp_copies
    #                 for sequential_id, _id in enumerate(ids):
    #                     all_copies[_id] = component_copies[sequential_id]
    #
    #         all_copies = tuple(all_copies)
    #         if all_copies not in all_copy_sets:
    #             all_copy_sets.add(all_copies)
    #         else:
    #             continue
    #
    #         """ record new copy values """
    #         final_results.append({"graph": deepcopy(self)})
    #         for go_s, this_symbol in enumerate(all_v_symbols):
    #             vertex_name = symbols_to_vertex[this_symbol]
    #             if vertex_name in final_results[go_res]["graph"].vertex_to_copy:
    #                 old_copy = final_results[go_res]["graph"].vertex_to_copy[vertex_name]
    #                 final_results[go_res]["graph"].copy_to_vertex[old_copy].remove(vertex_name)
    #                 if not final_results[go_res]["graph"].copy_to_vertex[old_copy]:
    #                     del final_results[go_res]["graph"].copy_to_vertex[old_copy]
    #             this_copy = all_copies[go_s]
    #             final_results[go_res]["graph"].vertex_to_copy[vertex_name] = this_copy
    #             if this_copy not in final_results[go_res]["graph"].copy_to_vertex:
    #                 final_results[go_res]["graph"].copy_to_vertex[this_copy] = set()
    #             final_results[go_res]["graph"].copy_to_vertex[this_copy].add(vertex_name)
    #
    #         """ re-estimate baseline depth """
    #         total_product = 0.
    #         total_len = 0
    #         for vertex_name in vertices_list:
    #             this_len = self.vertex_info[vertex_name].len \
    #                        * final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
    #             this_cov = self.vertex_info[vertex_name].cov \
    #                        / final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
    #             total_len += this_len
    #             total_product += this_len * this_cov
    #         final_results[go_res]["cov"] = total_product / total_len
    #     return final_results
    #     # else:
    #     #     """ produce the first-ranked copy combination """
    #     #     free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
    #     #                                 for i, this_copy in enumerate(copy_results)}
    #     #
    #     #     """ simplify copy values """  # 2020-02-22 added to avoid multiplicities res such as: [4, 8, 4]
    #     #     all_copies = []
    #     #     for this_symbol in all_v_symbols:
    #     #         vertex_name = symbols_to_vertex[this_symbol]
    #     #         this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
    #     #         if this_copy <= 0:
    #     #             raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
    #     #         all_copies.append(this_copy)
    #     #     if len(all_copies) == 0:
    #     #         raise ProcessingGraphFailed(
    #     #             "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
    #     #     elif len(all_copies) == 1:
    #     #         all_copies = [1]
    #     #     elif min(all_copies) == 1:
    #     #         pass
    #     #     else:
    #     #         new_all_copies = reduce_list_with_gcd(all_copies)
    #     #         if verbose and new_all_copies != all_copies:
    #     #             if log_handler:
    #     #                 log_handler.info("Estimated copies: " + str(all_copies))
    #     #                 log_handler.info("Reduced copies: " + str(new_all_copies))
    #     #             else:
    #     #                 sys.stdout.write("Estimated copies: " + str(all_copies) + "\n")
    #     #                 sys.stdout.write("Reduced copies: " + str(new_all_copies) + "\n")
    #     #         all_copies = new_all_copies
    #     #
    #     #     """ record new copy values """
    #     #     for go_s, this_symbol in enumerate(all_v_symbols):
    #     #         vertex_name = symbols_to_vertex[this_symbol]
    #     #         if vertex_name in self.vertex_to_copy:
    #     #             old_copy = self.vertex_to_copy[vertex_name]
    #     #             self.copy_to_vertex[old_copy].remove(vertex_name)
    #     #             if not self.copy_to_vertex[old_copy]:
    #     #                 del self.copy_to_vertex[old_copy]
    #     #         this_copy = all_copies[go_s]
    #     #         self.vertex_to_copy[vertex_name] = this_copy
    #     #         if this_copy not in self.copy_to_vertex:
    #     #             self.copy_to_vertex[this_copy] = set()
    #     #         self.copy_to_vertex[this_copy].add(vertex_name)
    #     #
    #     #     if debug or verbose:
    #     #         """ re-estimate baseline depth """
    #     #         total_product = 0.
    #     #         total_len = 0
    #     #         for vertex_name in vertices_list:
    #     #             this_len = self.vertex_info[vertex_name].len \
    #     #                        * self.vertex_to_copy.get(vertex_name, 1)
    #     #             this_cov = self.vertex_info[vertex_name].cov / self.vertex_to_copy.get(vertex_name, 1)
    #     #             total_len += this_len
    #     #             total_product += this_len * this_cov
    #     #         new_val = total_product / total_len
    #     #         if log_handler:
    #     #             log_handler.info("Average " + target_name_for_log + " kmer-coverage = " + str(round(new_val, 2)))
    #     #         else:
    #     #             sys.stdout.write(
    #     #                 "Average " + target_name_for_log + " kmer-coverage = " + str(round(new_val, 2)) + "\n")

    def estimate_copy_and_depth_precisely(self, expected_average_cov=None,  # broken_graph_allowed=False,
                                          verbose=False, log_handler=None, debug=False,
                                          target_name_for_log="target", n_iterations=None):
        """
        :param expected_average_cov: not used in the least-square version
        :param verbose:
        :param log_handler:
        :param debug:
        :param target_name_for_log:
        :param n_iterations:
        :return:
        """
        def get_formula(from_vertex, from_end, back_to_vertex, here_record_ends):
            result_form = v_vars[vertices_ids[from_vertex]]
            here_record_ends.add((from_vertex, from_end))
            # if back_to_vertex ~ from_vertex (from_vertex == back_to_vertex) form a loop, skipped
            if from_vertex != back_to_vertex:
                for next_v, next_e in self.vertex_info[from_vertex].connections[from_end]:
                    # if next_v ~ from_vertex (next_v == from_vertex) form a loop, add a pseudo vertex
                    if (next_v, next_e) == (from_vertex, not from_end):
                        # skip every self-loop 2020-06-23
                        # pseudo_self_circle_str = "P" + from_vertex
                        # if pseudo_self_circle_str not in extra_str_to_symbol_m2:
                        #     extra_str_to_symbol_m2[pseudo_self_circle_str] = Symbol(pseudo_self_circle_str, integer=True)
                        #     extra_symbol_to_str_m2[extra_str_to_symbol_m2[pseudo_self_circle_str]] = pseudo_self_circle_str
                        # result_form -= (extra_str_to_symbol_m2[pseudo_self_circle_str] - 1)
                        pass
                    # elif (next_v, next_e) != (back_to_vertex, back_to_end):
                    elif (next_v, next_e) not in here_record_ends:
                        result_form -= get_formula(next_v, next_e, from_vertex, here_record_ends)
            return result_form

        # # for compatibility between scipy and sympy
        # def least_square_function_v(x):
        #     return least_square_expr(*tuple(x))
        #
        # """ create constraints by creating inequations: the copy of every contig has to be >= 1 """
        #
        # def constraint_min_function(x):
        #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
        #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
        #     min_copy = np.array([1.001] * len(all_v_symbols) +
        #                         [1.001] * len(extra_symbol_to_str_m1) +
        #                         [2.001] * len(extra_symbol_to_str_m2))
        #     # effect: expression_array >= int(min_copy)
        #     return expression_array - min_copy
        #
        # def constraint_min_function_for_customized_brute(x):
        #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
        #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
        #     min_copy = np.array([1.0] * len(all_v_symbols) +
        #                         [1.0] * len(extra_symbol_to_str_m1) +
        #                         [2.0] * len(extra_symbol_to_str_m2))
        #     # effect: expression_array >= min_copy
        #     return expression_array - min_copy
        # def constraint_max_function(x):
        #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
        #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
        #     max_copy = np.array([expected_average_cov] * len(all_v_symbols) +
        #                         [expected_average_cov] * len(extra_symbol_to_str_m1) +
        #                         [expected_average_cov * 2] * len(extra_symbol_to_str_m2))
        #     # effect: expression_array <= max_copy
        #     return max_copy - expression_array
        #
        # def constraint_int_function(x):
        #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
        #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
        #     # diff = np.array([0] * len(all_symbols))
        #     return sum([abs(every_copy - int(every_copy)) for every_copy in expression_array])
        #
        # def minimize_brute_force(func, range_list, constraint_list, round_digit=4, display_p=True,
        #                          in_log_handler=log_handler):
        #     # time0 = time.time()
        #     best_fun_val = inf
        #     best_para_val = []
        #     count_round = 0
        #     count_valid = 0
        #     for value_set in product(*[list(this_range) for this_range in range_list]):
        #         count_round += 1
        #         is_valid_set = True
        #         for cons in constraint_list:
        #             if cons["type"] == "ineq":
        #                 try:
        #                     if (cons["fun"](value_set) < 0).any():
        #                         is_valid_set = False
        #                         # if in_log_handler and (debug or display_p):
        #                         #     in_log_handler.info("value_set={} ; illegal ineq constraints".format(value_set))
        #                         break
        #                 except TypeError:
        #                     # if in_log_handler and (debug or display_p):
        #                     #     in_log_handler.info("value_set={} ; illegal ineq constraints".format(value_set))
        #                     is_valid_set = False
        #                     break
        #             elif cons["type"] == "eq":
        #                 try:
        #                     if cons["fun"](value_set) != 0:
        #                         is_valid_set = False
        #                         # if in_log_handler and (debug or display_p):
        #                         #     in_log_handler.info("value_set={} ; illegal eq constraints".format(value_set))
        #                         break
        #                 except TypeError:
        #                     # if in_log_handler and (debug or display_p):
        #                     #     in_log_handler.info("value_set={} ; illegal eq constraints".format(value_set))
        #                     is_valid_set = False
        #                     break
        #         if not is_valid_set:
        #             continue
        #         count_valid += 1
        #         this_fun_val = func(value_set)
        #         if in_log_handler:
        #             if debug or display_p:
        #                 in_log_handler.info("value_set={} ; fun_val={}".format(value_set, this_fun_val))
        #         this_fun_val = round(this_fun_val, round_digit)
        #         if this_fun_val < best_fun_val:
        #             best_para_val = [value_set]
        #             best_fun_val = this_fun_val
        #         elif this_fun_val == best_fun_val:
        #             best_para_val.append(value_set)
        #         else:
        #             pass
        #     if in_log_handler:
        #         if debug or display_p:
        #             in_log_handler.info("Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round))
        #             in_log_handler.info("Brute best function value: " + str(best_fun_val))
        #         if debug:
        #             in_log_handler.info("Best solution: " + str(best_para_val))
        #     else:
        #         if debug or display_p:
        #             sys.stdout.write(
        #                 "Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round) + "\n")
        #             sys.stdout.write("Brute best function value: " + str(best_fun_val) + "\n")
        #         if debug:
        #             sys.stdout.write("Best solution: " + str(best_para_val) + "\n")
        #     return best_para_val
        if verbose:
            log_handler.info("Estimating copy and depth precisely ...")

        vertices_list = sorted(self.vertex_info)
        vertices_ids = {_v: _i for _i, _v in enumerate(vertices_list)}
        if len(vertices_list) == 1:
            cov_ = self.vertex_info[vertices_list[0]].cov
            # 2022-12-15, remove return_new_graph
            # if return_new_graphs:
            return [{"graph": deepcopy(self), "cov": cov_}]
            # else:
            #     if log_handler:
            #         log_handler.info("Average " + target_name_for_log + " kmer-coverage = " + str(round(cov_, 2)))
            #     else:
            #         sys.stdout.write(
            #             "Average " + target_name_for_log + " kmer-coverage = " + str(round(cov_, 2)) + "\n")
            #     return

        # reduce expected_average_cov to reduce computational burden
        all_coverages = [self.vertex_info[v_name].cov for v_name in vertices_list]
        # max_contig_multiplicity = \
        #     min(max_contig_multiplicity, int(2 * math.ceil(max(all_coverages) / min(all_coverages))))
        # if verbose:
        #     if log_handler:
        #         log_handler.info("Maximum multiplicity: " + str(max_contig_multiplicity))
        #     else:
        #         sys.stdout.write("Maximum multiplicity: " + str(max_contig_multiplicity) + "\n")

        """ use local gekko """
        m = GEKKO(remote=False)
        """ create constraints by creating multivariate equations """
        copy_upper_bound = int(4 * math.ceil(max(all_coverages) / min(all_coverages)))
        v_vars = m.Array(m.Var,
                         len(vertices_list),
                         lb=1,
                         ub=copy_upper_bound,
                         integer=True)
        # initialize free variables
        for go_v, v_var in enumerate(v_vars):
            v_var.value = self.vertex_to_copy[vertices_list[go_v]]
        # for go_sym, symbol_used in enumerate(free_copy_variables):
        #     if symbol_used in symbols_to_vertex and symbols_to_vertex[symbol_used] in self.vertex_to_copy:
        #         g_vars[go_sym].value = self.vertex_to_copy[symbols_to_vertex[symbol_used]]
        #     elif symbol_used in extra_symbol_initial_values:
        #         g_vars[go_sym].value = extra_symbol_initial_values[symbol_used]

        # vertex_to_symbols = {vertex_name: Symbol("V" + vertex_name, integer=True)  # positive=True)
        #                      for vertex_name in vertices_list}
        # symbols_to_vertex = {vertex_to_symbols[vertex_name]: vertex_name for vertex_name in vertices_list}
        formulae = []
        recorded_ends = set()
        for go_v, vertex_name in enumerate(vertices_list):
            for this_end in (True, False):
                if (vertex_name, this_end) not in recorded_ends:
                    recorded_ends.add((vertex_name, this_end))
                    if self.vertex_info[vertex_name].connections[this_end]:
                        this_formula = v_vars[go_v]
                        formulized = False
                        for n_v, n_e in self.vertex_info[vertex_name].connections[this_end]:
                            if (n_v, n_e) not in recorded_ends:
                                try:
                                    this_formula -= get_formula(n_v, n_e, vertex_name, recorded_ends)
                                    formulized = True
                                    # if verbose:
                                    #     if log_handler:
                                    #         log_handler.info("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
                                    #                          vertex_name + ECHO_DIRECTION[this_end] + ": " +
                                    #                          str(this_formula))
                                    #     else:
                                    #         sys.stdout.write("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
                                    #                          vertex_name + ECHO_DIRECTION[this_end] + ": " +
                                    #                          str(this_formula)+"\n")
                                except RecursionError:
                                    if log_handler:
                                        log_handler.warning("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
                                                            vertex_name + ECHO_DIRECTION[this_end] + " failed!")
                                    else:
                                        sys.stdout.write("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
                                                         vertex_name + ECHO_DIRECTION[this_end] + " failed!\n")
                                    raise ProcessingGraphFailed("RecursionError!")
                        if verbose:
                            if log_handler:
                                log_handler.info(
                                    "formulating for: " + vertex_name + ECHO_DIRECTION[this_end] + ": " +
                                    str(this_formula))
                            else:
                                sys.stdout.write(
                                    "formulating for: " + vertex_name + ECHO_DIRECTION[this_end] + ": " +
                                    str(this_formula) + "\n")
                        if formulized:
                            formulae.append(this_formula)
                    # 2022-12-13 remove this restriction
                    #            because we have a reduce_list_with_gcd for all graph component
                    # elif broken_graph_allowed:
                    #     # Extra limitation to force terminal vertex to have only one copy, to avoid over-estimation
                    #     # Under-estimation would not be a problem here,
                    #     # because the True-multiple-copy vertex would simply have no other connections,
                    #     # or failed in the following estimation if it does
                    #     formulae.append(vertex_to_symbols[vertex_name] - 1)

        """ add self-loop formulae """
        self_loop_v = OrderedDict()
        for vertex_name in vertices_list:
            if self.vertex_info[vertex_name].is_self_loop():
                self_loop_v[vertex_name] = self.vertex_to_copy[vertex_name]
                if log_handler:
                    log_handler.warning("Self-loop contig detected: Vertex_" + vertex_name)
        #
        if self_loop_v:
            p_vars = m.Array(m.Var,
                             len(self_loop_v),
                             lb=1,
                             ub=copy_upper_bound,
                             integer=True)
            for go_p, (vertex_name, initial_val) in enumerate(self_loop_v.items()):
                # set initial value
                p_vars[go_p].value = initial_val
                # add formulae
                this_formula = v_vars[vertices_ids[vertex_name]] - p_vars[go_p]
                formulae.append(this_formula)
                if verbose and log_handler:
                    log_handler.info(
                        "formulating for: " + vertex_name + ECHO_DIRECTION[True] + ": " + str(this_formula))
        """ add extra restriction on repeats """
        extra_m1 = []
        extra_m2 = []
        for vertex_name in vertices_list:
            single_pair_in_the_trunk_path = self.is_sequential_repeat(vertex_name)
            if single_pair_in_the_trunk_path:
                (from_v, from_e), (to_v, to_e) = single_pair_in_the_trunk_path
                # from_v and to_v are already in the "trunk path", if they are the same,
                # the graph is like two circles sharing the same sequential repeat, no need to add this limitation
                if from_v != to_v:
                    initial_val = round(self.vertex_to_float_copy[vertex_name] / self.vertex_to_float_copy[from_v])
                    if vertex_name in self_loop_v:
                        # self-loop vertex is allowed to have the multiplicity of 1
                        extra_m1.append([vertex_name, from_v, initial_val])
                    else:
                        extra_m2.append([vertex_name, from_v, max(initial_val, 2)])
        if extra_m1:
            m1_vars = m.Array(m.Var,
                              len(extra_m1),
                              lb=1,
                              ub=copy_upper_bound,
                              integer=True)
            for go_m, (vertex_name, from_v, initial_val) in enumerate(extra_m1):
                m1_vars[go_m].value = initial_val
                this_formula = v_vars[vertices_ids[vertex_name]] - v_vars[vertices_ids[from_v]] * m1_vars[go_m]
                formulae.append(this_formula)
                if verbose and log_handler:
                    log_handler.info("formulating for: " + vertex_name + ": " + str(this_formula))
        if extra_m2:
            m2_vars = m.Array(m.Var,
                              len(extra_m2),
                              lb=2,
                              ub=copy_upper_bound,
                              integer=True)
            for go_m, (vertex_name, from_v, initial_val) in enumerate(extra_m2):
                m2_vars[go_m].value = initial_val
                this_formula = v_vars[vertices_ids[vertex_name]] - v_vars[vertices_ids[from_v]] * m2_vars[go_m]
                formulae.append(this_formula)
                if verbose and log_handler:
                    log_handler.info("formulating for: " + vertex_name + ": " + str(this_formula))

        # # solve the equations
        # copy_solution = solve(formulae, all_v_symbols)
        #
        # copy_solution = copy_solution if copy_solution else {}
        # if type(copy_solution) == list:  # delete 0 containing set, even for self-loop vertex
        #     go_solution = 0
        #     while go_solution < len(copy_solution):
        #         if 0 in set(copy_solution[go_solution].values()):
        #             del copy_solution[go_solution]
        #         else:
        #             go_solution += 1
        # if not copy_solution:
        #     raise ProcessingGraphFailed("Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (1)!")
        # elif type(copy_solution) == list:
        #     if len(copy_solution) > 2:
        #         raise ProcessingGraphFailed("Incomplete/Complicated " + target_name_for_log + " graph (2)!")
        #     else:
        #         copy_solution = copy_solution[0]
        #
        # free_copy_variables = list()
        # for symbol_used in all_symbols:
        #     if symbol_used not in copy_solution:
        #         free_copy_variables.append(symbol_used)
        #         copy_solution[symbol_used] = symbol_used
        # if verbose:
        #     if log_handler:
        #         log_handler.info("copy equations: " + str(copy_solution))
        #         log_handler.info("free variables: " + str(free_copy_variables))
        #     else:
        #         sys.stdout.write("copy equations: " + str(copy_solution) + "\n")
        #         sys.stdout.write("free variables: " + str(free_copy_variables) + "\n")

        # """  """
        # least_square_expr = 0
        # for symbol_used in all_v_symbols:
        #     # least_square_expr += copy_solution[symbol_used]
        #     this_vertex = symbols_to_vertex[symbol_used]
        #     this_copy = self.vertex_to_float_copy[this_vertex]
        #     least_square_expr += (copy_solution[symbol_used] - this_copy) ** 2  # * self.vertex_info[this_vertex]["len"]
        # least_square_expr = lambdify(args=free_copy_variables, expr=least_square_expr)

        # if free_copy_variables:
        """ minimizing equation-based copy's deviations from coverage-based copy values """
        # ignore overlap influence
        least_square_list = []
        for go_v, vertex_name in enumerate(vertices_list):
            estimated_copy = self.vertex_to_float_copy[vertex_name]
            least_square_list.append((v_vars[go_v] - estimated_copy) ** 2)
        least_square_expr = sum(least_square_list)
        if verbose or debug:
            log_handler.info("square function: " + str(repr(least_square_expr)))
        # reform least_square_expr if string length > 15000
        exp_str_len = len(str(least_square_expr))
        if exp_str_len > 15000:  # not allowed by Gekko:APM
            num_blocks = math.ceil(exp_str_len / 10000.)
            block_size = math.ceil(len(least_square_list) / float(num_blocks))
            block_list = []
            for g_b in range(num_blocks):
                block_list.append(sum(least_square_list[g_b * block_size: (g_b + 1) * block_size]))
            least_square_expr = m.sum(block_list)

        # account for the influence of the overlap
        # total_len = 0
        # multinomial_loglike_list = []
        # v_to_len = {}
        # v_to_real_len = {}
        # all_obs = []
        # if self.__uni_overlap:
        #     for symbol_used in all_v_symbols:
        #         this_vertex = symbols_to_vertex[symbol_used]
        #         v_to_real_len[this_vertex] = self.vertex_info[this_vertex].len - self.__uni_overlap
        #         v_to_len[this_vertex] = eval(str(copy_solution[symbol_used].subs(replacements))) \
        #                                 * v_to_real_len[this_vertex]
        #         total_len += v_to_len[this_vertex]
        #     for symbol_used in all_v_symbols:
        #         this_vertex = symbols_to_vertex[symbol_used]
        #         prob = v_to_len[this_vertex] / total_len
        #         obs = self.vertex_info[this_vertex].cov * v_to_real_len[this_vertex]
        #         multinomial_loglike_list.append(m.log(prob) * obs)
        #         all_obs.append(obs)
        #         if verbose:
        #             if log_handler:
        #                 log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
        #             else:
        #                 sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
        # else:
        #     for symbol_used in all_v_symbols:
        #         this_vertex = symbols_to_vertex[symbol_used]
        #         overlaps = [_ovl
        #                     for _strand in (True, False)
        #                     for _next, _ovl in self.vertex_info[this_vertex].connections[_strand].items()]
        #         approximate_overlap = average_np_free(overlaps)
        #         v_to_real_len[this_vertex] = self.vertex_info[this_vertex].len - approximate_overlap
        #         v_to_len[this_vertex] = eval(str(copy_solution[symbol_used].subs(replacements)))\
        #                                 * v_to_real_len[this_vertex]
        #         total_len += v_to_len[this_vertex]
        #     for symbol_used in all_v_symbols:
        #         this_vertex = symbols_to_vertex[symbol_used]
        #         prob = v_to_len[this_vertex] / total_len
        #         obs = self.vertex_info[this_vertex].cov * v_to_real_len[this_vertex]
        #         # multinomial_loglike_list.append(m.log(prob) * obs)
        #         all_obs.append(obs)
        #         if verbose:
        #             if log_handler:
        #                 log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
        #             else:
        #                 sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
        # """extra arbitrary restriction to avoid over inflation of copies"""
        # multinomial_loglike_list.append(-abs(sum(all_obs) / expected_average_cov - total_len))
        # multinomial_loglike_expr = m.sum(multinomial_loglike_list)

        # for symbol_used in all_v_symbols:
        #     this_vertex = symbols_to_vertex[symbol_used]
        #     total_len += eval(str(copy_solution[symbol_used].subs(replacements))) * self.vertex_info[this_vertex].len
        # multinomial_like_expr = 0
        # for symbol_used in all_v_symbols:
        #     this_vertex = symbols_to_vertex[symbol_used]
        #     prob = eval(str(copy_solution[symbol_used].subs(replacements))) \
        #            * self.vertex_info[this_vertex].len / total_len
        #     obs = self.vertex_info[this_vertex].cov * self.vertex_info[this_vertex].len
        #     multinomial_like_expr += m.log(prob) * obs
        #     if verbose:
        #         if log_handler:
        #             log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
        #         else:
        #             sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
        m.Equations([f_ == 0 for f_ in formulae])
        m.Minimize(least_square_expr)
        # 1 for APOPT, 2 for BPOPT, 3 for IPOPT, 0 for all available solvers
        # here only 1 and 3 are available
        m.options.SOLVER = 1
        # setting empirical options
        # 5000 costs ~ 150 sec
        # get the variation within the data
        single_variations = []
        for vertex_name in self.vertex_info:
            f_copy = self.vertex_to_float_copy[vertex_name]
            single_variations.append(abs((math.ceil(f_copy) - f_copy) ** 2 - (math.floor(f_copy) - f_copy) ** 2))
        single_variations.sort()
        # largest_var = single_variations[-1] * 4
        if n_iterations is None:
            n_high_copy = sum([math.log2(self.vertex_to_float_copy[_v])
                               for _v in self.vertex_info if self.vertex_to_float_copy[_v] > 2])
            n_iterations = 500 + int(len(self.vertex_info) * n_high_copy)
        if verbose or debug:
            log_handler.info("setting n_iterations=" + str(n_iterations))
        if verbose or debug:
            log_handler.info("setting minlp_gap_tol=%.0e" % single_variations[0])
        m.solver_options = ['minlp_maximum_iterations ' + str(n_iterations),
                            # minlp iterations with integer solution
                            'minlp_max_iter_with_int_sol ' + str(n_iterations),
                            # treat minlp as nlp
                            'minlp_as_nlp 0',
                            # nlp sub-problem max iterations
                            'nlp_maximum_iterations ' + str(n_iterations),
                            # 1 = depth first, 2 = breadth first
                            'minlp_branch_method 2',
                            # maximum deviation from whole number:
                            # amount that a candidate solution variable can deviate from an integer solution
                            # and still be considered an integer
                            'minlp_integer_tol 1.0e-2',
                            # covergence tolerance
                            'minlp_gap_tol %.0e' % single_variations[0]]
        try:
            if debug or verbose:
                m.solve()
            else:
                m.solve(disp=False)
        except Exception as e:
            # TODO adjust parameters according to apm result, currently I do not know how to load apm result
            if "Solution Not Found" in str(e):
                raise ProcessingGraphFailed("Solution not found by apm for current graph!")
            else:
                raise e
        copy_results = list([x.value[0] for x in v_vars])
        if debug or verbose:
            for go_v, vertex_name in enumerate(vertices_list):
                log_handler.info(vertex_name + ": " + str(copy_results[go_v]))
        # """ produce all possible vertex copy combinations """
        # maybe no more multiple results since 2022-12 gekko update
        final_results = []
        all_copy_sets = set()  # removing duplicates in multiple results
        for go_res, copy_result in enumerate([copy_results]):
            """ simplify copy values """
            # 2020-02-22 added to avoid multiplicities res such as: [4, 8, 4]
            # 2022-12-15 add cluster info to simplify by graph components when the graph is broken
            all_copies = copy_result
            if len(self.vertex_clusters) == 1:
                if len(all_copies) == 0:
                    raise ProcessingGraphFailed(
                        "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
                elif len(all_copies) == 1:
                    all_copies = [1]
                elif min(all_copies) == 1:
                    pass
                else:
                    new_all_copies = reduce_list_with_gcd(all_copies)
                    if verbose and new_all_copies != all_copies:
                        if log_handler:
                            log_handler.info("Estimated copies: " + str(all_copies))
                            log_handler.info("Reduced copies: " + str(new_all_copies))
                        else:
                            sys.stdout.write("Estimated copies: " + str(all_copies) + "\n")
                            sys.stdout.write("Reduced copies: " + str(new_all_copies) + "\n")
                    all_copies = new_all_copies
            else:
                for v_cluster in self.vertex_clusters:
                    ids = [vertices_ids[_v] for _v in v_cluster]
                    component_copies = [all_copies[_id] for _id in ids]
                    if len(component_copies) == 0:
                        raise ProcessingGraphFailed(
                            "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
                    elif len(component_copies) == 1:
                        component_copies = [1]
                    elif min(component_copies) == 1:
                        pass
                    else:
                        new_comp_copies = reduce_list_with_gcd(component_copies)
                        if verbose and new_comp_copies != component_copies:
                            if log_handler:
                                log_handler.info("Estimated copies: " + str(component_copies))
                                log_handler.info("Reduced copies: " + str(new_comp_copies))
                            else:
                                sys.stdout.write("Estimated copies: " + str(component_copies) + "\n")
                                sys.stdout.write("Reduced copies: " + str(new_comp_copies) + "\n")
                        component_copies = new_comp_copies
                    for sequential_id, _id in enumerate(ids):
                        all_copies[_id] = component_copies[sequential_id]

            all_copies = tuple(all_copies)
            if all_copies not in all_copy_sets:
                all_copy_sets.add(all_copies)
            else:
                continue

            """ record new copy values """
            final_results.append({"graph": deepcopy(self)})
            for go_v, vertex_name in enumerate(vertices_list):
                if vertex_name in final_results[go_res]["graph"].vertex_to_copy:
                    old_copy = final_results[go_res]["graph"].vertex_to_copy[vertex_name]
                    final_results[go_res]["graph"].copy_to_vertex[old_copy].remove(vertex_name)
                    if not final_results[go_res]["graph"].copy_to_vertex[old_copy]:
                        del final_results[go_res]["graph"].copy_to_vertex[old_copy]
                estimated_copy = all_copies[go_v]
                final_results[go_res]["graph"].vertex_to_copy[vertex_name] = estimated_copy
                if estimated_copy not in final_results[go_res]["graph"].copy_to_vertex:
                    final_results[go_res]["graph"].copy_to_vertex[estimated_copy] = set()
                final_results[go_res]["graph"].copy_to_vertex[estimated_copy].add(vertex_name)

            """ re-estimate baseline depth """
            total_product = 0.
            total_len = 0
            for vertex_name in vertices_list:
                this_len = self.vertex_info[vertex_name].len \
                           * final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
                this_cov = self.vertex_info[vertex_name].cov \
                           / final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
                total_len += this_len
                total_product += this_len * this_cov
            final_results[go_res]["cov"] = total_product / total_len
        return final_results
        # else:
        #     """ produce the first-ranked copy combination """
        #     free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
        #                                 for i, this_copy in enumerate(copy_results)}
        #
        #     """ simplify copy values """  # 2020-02-22 added to avoid multiplicities res such as: [4, 8, 4]
        #     all_copies = []
        #     for this_symbol in all_v_symbols:
        #         vertex_name = symbols_to_vertex[this_symbol]
        #         this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
        #         if this_copy <= 0:
        #             raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
        #         all_copies.append(this_copy)
        #     if len(all_copies) == 0:
        #         raise ProcessingGraphFailed(
        #             "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
        #     elif len(all_copies) == 1:
        #         all_copies = [1]
        #     elif min(all_copies) == 1:
        #         pass
        #     else:
        #         new_all_copies = reduce_list_with_gcd(all_copies)
        #         if verbose and new_all_copies != all_copies:
        #             if log_handler:
        #                 log_handler.info("Estimated copies: " + str(all_copies))
        #                 log_handler.info("Reduced copies: " + str(new_all_copies))
        #             else:
        #                 sys.stdout.write("Estimated copies: " + str(all_copies) + "\n")
        #                 sys.stdout.write("Reduced copies: " + str(new_all_copies) + "\n")
        #         all_copies = new_all_copies
        #
        #     """ record new copy values """
        #     for go_s, this_symbol in enumerate(all_v_symbols):
        #         vertex_name = symbols_to_vertex[this_symbol]
        #         if vertex_name in self.vertex_to_copy:
        #             old_copy = self.vertex_to_copy[vertex_name]
        #             self.copy_to_vertex[old_copy].remove(vertex_name)
        #             if not self.copy_to_vertex[old_copy]:
        #                 del self.copy_to_vertex[old_copy]
        #         this_copy = all_copies[go_s]
        #         self.vertex_to_copy[vertex_name] = this_copy
        #         if this_copy not in self.copy_to_vertex:
        #             self.copy_to_vertex[this_copy] = set()
        #         self.copy_to_vertex[this_copy].add(vertex_name)
        #
        #     if debug or verbose:
        #         """ re-estimate baseline depth """
        #         total_product = 0.
        #         total_len = 0
        #         for vertex_name in vertices_list:
        #             this_len = self.vertex_info[vertex_name].len \
        #                        * self.vertex_to_copy.get(vertex_name, 1)
        #             this_cov = self.vertex_info[vertex_name].cov / self.vertex_to_copy.get(vertex_name, 1)
        #             total_len += this_len
        #             total_product += this_len * this_cov
        #         new_val = total_product / total_len
        #         if log_handler:
        #             log_handler.info("Average " + target_name_for_log + " kmer-coverage = " + str(round(new_val, 2)))
        #         else:
        #             sys.stdout.write(
        #                 "Average " + target_name_for_log + " kmer-coverage = " + str(round(new_val, 2)) + "\n")

    # def estimate_copy_and_depth_precisely_sympy(self, expected_average_cov=None,  # broken_graph_allowed=False,
    #                                             verbose=False, log_handler=None, debug=False,
    #                                             target_name_for_log="target", n_iterations=None):
    #
    #     def get_formula(from_vertex, from_end, back_to_vertex, here_record_ends):
    #         result_form = vertex_to_symbols[from_vertex]
    #         here_record_ends.add((from_vertex, from_end))
    #         # if back_to_vertex ~ from_vertex (from_vertex == back_to_vertex) form a loop, skipped
    #         if from_vertex != back_to_vertex:
    #             for next_v, next_e in self.vertex_info[from_vertex].connections[from_end]:
    #                 # if next_v ~ from_vertex (next_v == from_vertex) form a loop, add a pseudo vertex
    #                 if (next_v, next_e) == (from_vertex, not from_end):
    #                     # skip every self-loop 2020-06-23
    #                     # pseudo_self_circle_str = "P" + from_vertex
    #                     # if pseudo_self_circle_str not in extra_str_to_symbol_m2:
    #                     #     extra_str_to_symbol_m2[pseudo_self_circle_str] = Symbol(pseudo_self_circle_str, integer=True)
    #                     #     extra_symbol_to_str_m2[extra_str_to_symbol_m2[pseudo_self_circle_str]] = pseudo_self_circle_str
    #                     # result_form -= (extra_str_to_symbol_m2[pseudo_self_circle_str] - 1)
    #                     pass
    #                 # elif (next_v, next_e) != (back_to_vertex, back_to_end):
    #                 elif (next_v, next_e) not in here_record_ends:
    #                     result_form -= get_formula(next_v, next_e, from_vertex, here_record_ends)
    #         return result_form
    #
    #     # # for compatibility between scipy and sympy
    #     # def least_square_function_v(x):
    #     #     return least_square_expr(*tuple(x))
    #     #
    #     # """ create constraints by creating inequations: the copy of every contig has to be >= 1 """
    #     #
    #     # def constraint_min_function(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     min_copy = np.array([1.001] * len(all_v_symbols) +
    #     #                         [1.001] * len(extra_symbol_to_str_m1) +
    #     #                         [2.001] * len(extra_symbol_to_str_m2))
    #     #     # effect: expression_array >= int(min_copy)
    #     #     return expression_array - min_copy
    #     #
    #     # def constraint_min_function_for_customized_brute(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     min_copy = np.array([1.0] * len(all_v_symbols) +
    #     #                         [1.0] * len(extra_symbol_to_str_m1) +
    #     #                         [2.0] * len(extra_symbol_to_str_m2))
    #     #     # effect: expression_array >= min_copy
    #     #     return expression_array - min_copy
    #
    #     def constraint_min_function_for_gekko(g_vars):
    #         subs_tuples = [(symb_used_, Symbol("g_vars[" + str(go_sym) + "]"))
    #                        for go_sym, symb_used_ in enumerate(free_copy_variables)]
    #         expression_array = [copy_solution[this_sym].subs(subs_tuples) for this_sym in all_symbols]
    #         min_copy = [1] * len(all_v_symbols) + \
    #                    [1] * len(extra_symbol_to_str_m1) + \
    #                    [2] * len(extra_symbol_to_str_m2)
    #         # effect: expression_array >= min_copy
    #         expression = []
    #         if verbose or debug:
    #             for e, c in zip(expression_array, min_copy):
    #                 expression.append(eval(str(e) + ">=" + str(c)))
    #                 log_handler.info("  constraint: " + str(e) + ">=" + str(c))
    #         else:
    #             for e, c in zip(expression_array, min_copy):
    #                 expression.append(eval(str(e) + ">=" + str(c)))
    #         expression = [expr for expr in expression if not isinstance(expr, bool)]
    #         return expression
    #
    #     # def constraint_max_function(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     max_copy = np.array([expected_average_cov] * len(all_v_symbols) +
    #     #                         [expected_average_cov] * len(extra_symbol_to_str_m1) +
    #     #                         [expected_average_cov * 2] * len(extra_symbol_to_str_m2))
    #     #     # effect: expression_array <= max_copy
    #     #     return max_copy - expression_array
    #     #
    #     # def constraint_int_function(x):
    #     #     replacements = [(symbol_used, x[go_sym]) for go_sym, symbol_used in enumerate(free_copy_variables)]
    #     #     expression_array = np.array([copy_solution[this_sym].subs(replacements) for this_sym in all_symbols])
    #     #     # diff = np.array([0] * len(all_symbols))
    #     #     return sum([abs(every_copy - int(every_copy)) for every_copy in expression_array])
    #     #
    #     # def minimize_brute_force(func, range_list, constraint_list, round_digit=4, display_p=True,
    #     #                          in_log_handler=log_handler):
    #     #     # time0 = time.time()
    #     #     best_fun_val = inf
    #     #     best_para_val = []
    #     #     count_round = 0
    #     #     count_valid = 0
    #     #     for value_set in product(*[list(this_range) for this_range in range_list]):
    #     #         count_round += 1
    #     #         is_valid_set = True
    #     #         for cons in constraint_list:
    #     #             if cons["type"] == "ineq":
    #     #                 try:
    #     #                     if (cons["fun"](value_set) < 0).any():
    #     #                         is_valid_set = False
    #     #                         # if in_log_handler and (debug or display_p):
    #     #                         #     in_log_handler.info("value_set={} ; illegal ineq constraints".format(value_set))
    #     #                         break
    #     #                 except TypeError:
    #     #                     # if in_log_handler and (debug or display_p):
    #     #                     #     in_log_handler.info("value_set={} ; illegal ineq constraints".format(value_set))
    #     #                     is_valid_set = False
    #     #                     break
    #     #             elif cons["type"] == "eq":
    #     #                 try:
    #     #                     if cons["fun"](value_set) != 0:
    #     #                         is_valid_set = False
    #     #                         # if in_log_handler and (debug or display_p):
    #     #                         #     in_log_handler.info("value_set={} ; illegal eq constraints".format(value_set))
    #     #                         break
    #     #                 except TypeError:
    #     #                     # if in_log_handler and (debug or display_p):
    #     #                     #     in_log_handler.info("value_set={} ; illegal eq constraints".format(value_set))
    #     #                     is_valid_set = False
    #     #                     break
    #     #         if not is_valid_set:
    #     #             continue
    #     #         count_valid += 1
    #     #         this_fun_val = func(value_set)
    #     #         if in_log_handler:
    #     #             if debug or display_p:
    #     #                 in_log_handler.info("value_set={} ; fun_val={}".format(value_set, this_fun_val))
    #     #         this_fun_val = round(this_fun_val, round_digit)
    #     #         if this_fun_val < best_fun_val:
    #     #             best_para_val = [value_set]
    #     #             best_fun_val = this_fun_val
    #     #         elif this_fun_val == best_fun_val:
    #     #             best_para_val.append(value_set)
    #     #         else:
    #     #             pass
    #     #     if in_log_handler:
    #     #         if debug or display_p:
    #     #             in_log_handler.info("Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round))
    #     #             in_log_handler.info("Brute best function value: " + str(best_fun_val))
    #     #         if debug:
    #     #             in_log_handler.info("Best solution: " + str(best_para_val))
    #     #     else:
    #     #         if debug or display_p:
    #     #             sys.stdout.write(
    #     #                 "Brute valid/candidate rounds: " + str(count_valid) + "/" + str(count_round) + "\n")
    #     #             sys.stdout.write("Brute best function value: " + str(best_fun_val) + "\n")
    #     #         if debug:
    #     #             sys.stdout.write("Best solution: " + str(best_para_val) + "\n")
    #     #     return best_para_val
    #     if verbose:
    #         log_handler.info("Estimating copy and depth precisely ...")
    #
    #     vertices_list = sorted(self.vertex_info)
    #     if len(vertices_list) == 1:
    #         cov_ = self.vertex_info[vertices_list[0]].cov
    #         # 2022-12-15, remove return_new_graph
    #         # if return_new_graphs:
    #         return [{"graph": deepcopy(self), "cov": cov_}]
    #         # else:
    #         #     if log_handler:
    #         #         log_handler.info("Average " + target_name_for_log + " kmer-coverage = " + str(round(cov_, 2)))
    #         #     else:
    #         #         sys.stdout.write(
    #         #             "Average " + target_name_for_log + " kmer-coverage = " + str(round(cov_, 2)) + "\n")
    #         #     return
    #
    #     # reduce expected_average_cov to reduce computational burden
    #     all_coverages = [self.vertex_info[v_name].cov for v_name in vertices_list]
    #     # max_contig_multiplicity = \
    #     #     min(max_contig_multiplicity, int(2 * math.ceil(max(all_coverages) / min(all_coverages))))
    #     # if verbose:
    #     #     if log_handler:
    #     #         log_handler.info("Maximum multiplicity: " + str(max_contig_multiplicity))
    #     #     else:
    #     #         sys.stdout.write("Maximum multiplicity: " + str(max_contig_multiplicity) + "\n")
    #
    #     """ create constraints by creating multivariate equations """
    #     vertex_to_symbols = {vertex_name: Symbol("V" + vertex_name, integer=True)  # positive=True)
    #                          for vertex_name in vertices_list}
    #     symbols_to_vertex = {vertex_to_symbols[vertex_name]: vertex_name for vertex_name in vertices_list}
    #     extra_str_to_symbol_m1 = {}
    #     extra_str_to_symbol_m2 = {}
    #     extra_symbol_to_str_m1 = {}
    #     extra_symbol_to_str_m2 = {}
    #     extra_symbol_initial_values = {}
    #     formulae = []
    #     recorded_ends = set()
    #     for vertex_name in vertices_list:
    #         for this_end in (True, False):
    #             if (vertex_name, this_end) not in recorded_ends:
    #                 recorded_ends.add((vertex_name, this_end))
    #                 if self.vertex_info[vertex_name].connections[this_end]:
    #                     this_formula = vertex_to_symbols[vertex_name]
    #                     formulized = False
    #                     for n_v, n_e in self.vertex_info[vertex_name].connections[this_end]:
    #                         if (n_v, n_e) not in recorded_ends:
    #                             # if n_v in vertices_set:
    #                             # recorded_ends.add((n_v, n_e))
    #                             try:
    #                                 this_formula -= get_formula(n_v, n_e, vertex_name, recorded_ends)
    #                                 formulized = True
    #                                 # if verbose:
    #                                 #     if log_handler:
    #                                 #         log_handler.info("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                 #                          vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 #                          str(this_formula))
    #                                 #     else:
    #                                 #         sys.stdout.write("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                 #                          vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 #                          str(this_formula)+"\n")
    #                             except RecursionError:
    #                                 if log_handler:
    #                                     log_handler.warning("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                                         vertex_name + ECHO_DIRECTION[this_end] + " failed!")
    #                                 else:
    #                                     sys.stdout.write("formulating for: " + n_v + ECHO_DIRECTION[n_e] + "->" +
    #                                                      vertex_name + ECHO_DIRECTION[this_end] + " failed!\n")
    #                                 raise ProcessingGraphFailed("RecursionError!")
    #                     if verbose:
    #                         if log_handler:
    #                             log_handler.info(
    #                                 "formulating for: " + vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 str(this_formula))
    #                         else:
    #                             sys.stdout.write(
    #                                 "formulating for: " + vertex_name + ECHO_DIRECTION[this_end] + ": " +
    #                                 str(this_formula) + "\n")
    #                     if formulized:
    #                         formulae.append(this_formula)
    #                 # 2022-12-13 remove this restriction
    #                 #            because we have a reduce_list_with_gcd for all graph component
    #                 # elif broken_graph_allowed:
    #                 #     # Extra limitation to force terminal vertex to have only one copy, to avoid over-estimation
    #                 #     # Under-estimation would not be a problem here,
    #                 #     # because the True-multiple-copy vertex would simply have no other connections,
    #                 #     # or failed in the following estimation if it does
    #                 #     formulae.append(vertex_to_symbols[vertex_name] - 1)
    #
    #     # add self-loop formulae
    #     self_loop_v = set()
    #     for vertex_name in vertices_list:
    #         if self.vertex_info[vertex_name].is_self_loop():
    #             self_loop_v.add(vertex_name)
    #             if log_handler:
    #                 log_handler.warning("Self-loop contig detected: Vertex_" + vertex_name)
    #             pseudo_self_loop_str = "P" + vertex_name
    #             if pseudo_self_loop_str not in extra_str_to_symbol_m1:
    #                 extra_str_to_symbol_m1[pseudo_self_loop_str] = Symbol(pseudo_self_loop_str, integer=True)
    #                 extra_symbol_to_str_m1[extra_str_to_symbol_m1[pseudo_self_loop_str]] = pseudo_self_loop_str
    #             this_formula = vertex_to_symbols[vertex_name] - extra_str_to_symbol_m1[pseudo_self_loop_str]
    #             extra_symbol_initial_values[extra_str_to_symbol_m1[pseudo_self_loop_str]] = \
    #                 self.vertex_to_copy[vertex_name]
    #             formulae.append(this_formula)
    #             if verbose:
    #                 if log_handler:
    #                     log_handler.info(
    #                         "formulating for: " + vertex_name + ECHO_DIRECTION[True] + ": " + str(this_formula))
    #                 else:
    #                     sys.stdout.write(
    #                         "formulating for: " + vertex_name + ECHO_DIRECTION[True] + ": " + str(this_formula) + "\n")
    #
    #     # add following extra limitation
    #     # set cov_sequential_repeat = x*near_by_cov, x is an integer
    #     for vertex_name in vertices_list:
    #         single_pair_in_the_trunk_path = self.is_sequential_repeat(vertex_name)
    #         if single_pair_in_the_trunk_path:
    #             (from_v, from_e), (to_v, to_e) = single_pair_in_the_trunk_path
    #             # from_v and to_v are already in the "trunk path", if they are the same,
    #             # the graph is like two circles sharing the same sequential repeat, no need to add this limitation
    #             if from_v != to_v:
    #                 new_str = "E" + str(len(extra_str_to_symbol_m1) + len(extra_str_to_symbol_m2))
    #                 if vertex_name in self_loop_v:
    #                     # self-loop vertex is allowed to have the multiplicity of 1
    #                     extra_str_to_symbol_m1[new_str] = Symbol(new_str, integer=True)
    #                     extra_symbol_to_str_m1[extra_str_to_symbol_m1[new_str]] = new_str
    #                     this_formula = vertex_to_symbols[vertex_name] - \
    #                                    vertex_to_symbols[from_v] * extra_str_to_symbol_m1[new_str]
    #                     extra_symbol_initial_values[extra_str_to_symbol_m1[new_str]] = \
    #                         round(self.vertex_to_float_copy[vertex_name] / self.vertex_to_float_copy[from_v])
    #                 else:
    #                     extra_str_to_symbol_m2[new_str] = Symbol(new_str, integer=True)
    #                     extra_symbol_to_str_m2[extra_str_to_symbol_m2[new_str]] = new_str
    #                     this_formula = vertex_to_symbols[vertex_name] - \
    #                                    vertex_to_symbols[from_v] * extra_str_to_symbol_m2[new_str]
    #                     extra_symbol_initial_values[extra_str_to_symbol_m2[new_str]] = \
    #                         round(self.vertex_to_float_copy[vertex_name] / self.vertex_to_float_copy[from_v])
    #                 formulae.append(this_formula)
    #                 if verbose:
    #                     if log_handler:
    #                         log_handler.info("formulating for: " + vertex_name + ": " + str(this_formula))
    #                     else:
    #                         sys.stdout.write("formulating for: " + vertex_name + ": " + str(this_formula) + "\n")
    #
    #     all_v_symbols = list(symbols_to_vertex)
    #     all_symbols = all_v_symbols + list(extra_symbol_to_str_m1) + list(extra_symbol_to_str_m2)
    #     if verbose or debug:
    #         if log_handler:
    #             log_handler.info("formulae: " + str(formulae))
    #         else:
    #             sys.stdout.write("formulae: " + str(formulae) + "\n")
    #     # solve the equations
    #     copy_solution = solve(formulae, all_v_symbols)
    #
    #     copy_solution = copy_solution if copy_solution else {}
    #     if type(copy_solution) == list:  # delete 0 containing set, even for self-loop vertex
    #         go_solution = 0
    #         while go_solution < len(copy_solution):
    #             if 0 in set(copy_solution[go_solution].values()):
    #                 del copy_solution[go_solution]
    #             else:
    #                 go_solution += 1
    #     if not copy_solution:
    #         raise ProcessingGraphFailed("Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (1)!")
    #     elif type(copy_solution) == list:
    #         if len(copy_solution) > 2:
    #             raise ProcessingGraphFailed("Incomplete/Complicated " + target_name_for_log + " graph (2)!")
    #         else:
    #             copy_solution = copy_solution[0]
    #
    #     free_copy_variables = list()
    #     for symbol_used in all_symbols:
    #         if symbol_used not in copy_solution:
    #             free_copy_variables.append(symbol_used)
    #             copy_solution[symbol_used] = symbol_used
    #     if verbose:
    #         if log_handler:
    #             log_handler.info("copy equations: " + str(copy_solution))
    #             log_handler.info("free variables: " + str(free_copy_variables))
    #         else:
    #             sys.stdout.write("copy equations: " + str(copy_solution) + "\n")
    #             sys.stdout.write("free variables: " + str(free_copy_variables) + "\n")
    #
    #     # """  """
    #     # least_square_expr = 0
    #     # for symbol_used in all_v_symbols:
    #     #     # least_square_expr += copy_solution[symbol_used]
    #     #     this_vertex = symbols_to_vertex[symbol_used]
    #     #     this_copy = self.vertex_to_float_copy[this_vertex]
    #     #     least_square_expr += (copy_solution[symbol_used] - this_copy) ** 2  # * self.vertex_info[this_vertex]["len"]
    #     # least_square_expr = lambdify(args=free_copy_variables, expr=least_square_expr)
    #
    #     if free_copy_variables:
    #         """ minimizing equation-based copy's deviations from coverage-based copy values """
    #         # ignore overlap influence
    #         m = GEKKO(remote=False)
    #         g_vars = m.Array(m.Var,
    #                          len(free_copy_variables),
    #                          lb=1,
    #                          ub=int(4 * math.ceil(max(all_coverages) / min(all_coverages))),
    #                          integer=True)
    #         # initialize free variables
    #         for go_sym, symbol_used in enumerate(free_copy_variables):
    #             if symbol_used in symbols_to_vertex and symbols_to_vertex[symbol_used] in self.vertex_to_copy:
    #                 g_vars[go_sym].value = self.vertex_to_copy[symbols_to_vertex[symbol_used]]
    #             elif symbol_used in extra_symbol_initial_values:
    #                 g_vars[go_sym].value = extra_symbol_initial_values[symbol_used]
    #         # get the variation within the data
    #         single_variations = []
    #         for vertex_name in self.vertex_info:
    #             f_copy = self.vertex_to_float_copy[vertex_name]
    #             single_variations.append(abs((math.ceil(f_copy) - f_copy) ** 2 - (math.floor(f_copy) - f_copy) ** 2))
    #         single_variations.sort()
    #         # largest_var = single_variations[-1] * 4
    #
    #         replacements = [(symbol_used, Symbol("g_vars[" + str(go_sym) + "]"))
    #                         for go_sym, symbol_used in enumerate(free_copy_variables)]
    #         least_square_list = []
    #         for symbol_used in all_v_symbols:
    #             this_vertex = symbols_to_vertex[symbol_used]
    #             this_copy = self.vertex_to_float_copy[this_vertex]
    #             symbol_copy = eval(str(copy_solution[symbol_used].subs(replacements)))
    #             least_square_list.append((symbol_copy - this_copy) ** 2)
    #             # not working
    #             # # constraint the number to be integer
    #             # least_square_list.append(largest_var * (symbol_copy - int(symbol_copy)) ** 2)
    #         # least_square_expr = sum(least_square_list) will lead to no solution error for many variables
    #         least_square_expr = sum(least_square_list)
    #         if verbose or debug:
    #             log_handler.info("square function: " + str(repr(least_square_expr)))
    #         # reform least_square_expr if string length > 15000
    #         exp_str_len = len(str(least_square_expr))
    #         if exp_str_len > 15000:  # not allowed by Gekko:APM
    #             num_blocks = math.ceil(exp_str_len / 10000.)
    #             block_size = math.ceil(len(least_square_list) / float(num_blocks))
    #             block_list = []
    #             for g_b in range(num_blocks):
    #                 block_list.append(sum(least_square_list[g_b * block_size: (g_b + 1) * block_size]))
    #             least_square_expr = m.sum(block_list)
    #
    #         # account for the influence of the overlap
    #         # total_len = 0
    #         # multinomial_loglike_list = []
    #         # v_to_len = {}
    #         # v_to_real_len = {}
    #         # all_obs = []
    #         # if self.__uni_overlap:
    #         #     for symbol_used in all_v_symbols:
    #         #         this_vertex = symbols_to_vertex[symbol_used]
    #         #         v_to_real_len[this_vertex] = self.vertex_info[this_vertex].len - self.__uni_overlap
    #         #         v_to_len[this_vertex] = eval(str(copy_solution[symbol_used].subs(replacements))) \
    #         #                                 * v_to_real_len[this_vertex]
    #         #         total_len += v_to_len[this_vertex]
    #         #     for symbol_used in all_v_symbols:
    #         #         this_vertex = symbols_to_vertex[symbol_used]
    #         #         prob = v_to_len[this_vertex] / total_len
    #         #         obs = self.vertex_info[this_vertex].cov * v_to_real_len[this_vertex]
    #         #         multinomial_loglike_list.append(m.log(prob) * obs)
    #         #         all_obs.append(obs)
    #         #         if verbose:
    #         #             if log_handler:
    #         #                 log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
    #         #             else:
    #         #                 sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
    #         # else:
    #         #     for symbol_used in all_v_symbols:
    #         #         this_vertex = symbols_to_vertex[symbol_used]
    #         #         overlaps = [_ovl
    #         #                     for _strand in (True, False)
    #         #                     for _next, _ovl in self.vertex_info[this_vertex].connections[_strand].items()]
    #         #         approximate_overlap = average_np_free(overlaps)
    #         #         v_to_real_len[this_vertex] = self.vertex_info[this_vertex].len - approximate_overlap
    #         #         v_to_len[this_vertex] = eval(str(copy_solution[symbol_used].subs(replacements)))\
    #         #                                 * v_to_real_len[this_vertex]
    #         #         total_len += v_to_len[this_vertex]
    #         #     for symbol_used in all_v_symbols:
    #         #         this_vertex = symbols_to_vertex[symbol_used]
    #         #         prob = v_to_len[this_vertex] / total_len
    #         #         obs = self.vertex_info[this_vertex].cov * v_to_real_len[this_vertex]
    #         #         # multinomial_loglike_list.append(m.log(prob) * obs)
    #         #         all_obs.append(obs)
    #         #         if verbose:
    #         #             if log_handler:
    #         #                 log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
    #         #             else:
    #         #                 sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
    #         # """extra arbitrary restriction to avoid over inflation of copies"""
    #         # multinomial_loglike_list.append(-abs(sum(all_obs) / expected_average_cov - total_len))
    #         # multinomial_loglike_expr = m.sum(multinomial_loglike_list)
    #
    #         # for symbol_used in all_v_symbols:
    #         #     this_vertex = symbols_to_vertex[symbol_used]
    #         #     total_len += eval(str(copy_solution[symbol_used].subs(replacements))) * self.vertex_info[this_vertex].len
    #         # multinomial_like_expr = 0
    #         # for symbol_used in all_v_symbols:
    #         #     this_vertex = symbols_to_vertex[symbol_used]
    #         #     prob = eval(str(copy_solution[symbol_used].subs(replacements))) \
    #         #            * self.vertex_info[this_vertex].len / total_len
    #         #     obs = self.vertex_info[this_vertex].cov * self.vertex_info[this_vertex].len
    #         #     multinomial_like_expr += m.log(prob) * obs
    #         #     if verbose:
    #         #         if log_handler:
    #         #             log_handler.info("   >" + this_vertex + "\t" + str(obs))  # + "\t" + str(prob)
    #         #         else:
    #         #             sys.stdout.write("   >" + this_vertex + "\t" + str(obs) + "\n")
    #         m.Equations(constraint_min_function_for_gekko(g_vars))
    #         m.Minimize(least_square_expr)
    #         # 1 for APOPT, 2 for BPOPT, 3 for IPOPT, 0 for all available solvers
    #         # here only 1 and 3 are available
    #         m.options.SOLVER = 1
    #         # setting empirical options
    #         # 5000 costs ~ 150 sec
    #         if n_iterations is None:
    #             n_high_copy = sum([math.log2(self.vertex_to_float_copy[_v])
    #                                for _v in self.vertex_info if self.vertex_to_float_copy[_v] > 2])
    #             n_iterations = 500 + int(len(self.vertex_info) * n_high_copy)
    #         if verbose or debug:
    #             log_handler.info("setting n_iterations=" + str(n_iterations))
    #         if verbose or debug:
    #             log_handler.info("setting minlp_gap_tol=%.0e" % single_variations[0])
    #         m.solver_options = ['minlp_maximum_iterations ' + str(n_iterations),
    #                             # minlp iterations with integer solution
    #                             'minlp_max_iter_with_int_sol ' + str(n_iterations),
    #                             # treat minlp as nlp
    #                             'minlp_as_nlp 0',
    #                             # nlp sub-problem max iterations
    #                             'nlp_maximum_iterations ' + str(n_iterations),
    #                             # 1 = depth first, 2 = breadth first
    #                             'minlp_branch_method 2',
    #                             # maximum deviation from whole number:
    #                             # amount that a candidate solution variable can deviate from an integer solution
    #                             # and still be considered an integer
    #                             'minlp_integer_tol 1.0e-2',
    #                             # covergence tolerance
    #                             'minlp_gap_tol %.0e' % single_variations[0]]
    #         if debug or verbose:
    #             m.solve()
    #         else:
    #             m.solve(disp=False)
    #         # print([x.value[0] for x in g_vars])
    #         copy_results = list([x.value[0] for x in g_vars])
    #
    #         # # for safe running
    #         # if len(free_copy_variables) > 10:
    #         #     raise ProcessingGraphFailed("Free variable > 10 is not accepted yet!")
    #         #
    #         # if expected_average_cov ** len(free_copy_variables) < 5E6:
    #         #     # sometimes, SLSQP ignores bounds and constraints
    #         #     copy_results = minimize_brute_force(
    #         #         func=least_square_function_v, range_list=[range(1, expected_average_cov + 1)] * len(free_copy_variables),
    #         #         constraint_list=({'type': 'ineq', 'fun': constraint_min_function_for_customized_brute},
    #         #                          {'type': 'eq', 'fun': constraint_int_function},
    #         #                          {'type': 'ineq', 'fun': constraint_max_function}),
    #         #         display_p=verbose)
    #         # else:
    #         #     constraints = ({'type': 'ineq', 'fun': constraint_min_function},
    #         #                    {'type': 'eq', 'fun': constraint_int_function},
    #         #                    {'type': 'ineq', 'fun': constraint_max_function})
    #         #     copy_results = set()
    #         #     best_fun = inf
    #         #     opt = {'disp': verbose, "maxiter": 100}
    #         #     for initial_copy in range(expected_average_cov * 2 + 1):
    #         #         if initial_copy < expected_average_cov:
    #         #             initials = np.array([initial_copy + 1] * len(free_copy_variables))
    #         #         elif initial_copy < expected_average_cov * 2:
    #         #             initials = np.array([random.randint(1, expected_average_cov)] * len(free_copy_variables))
    #         #         else:
    #         #             initials = np.array([self.vertex_to_copy.get(symbols_to_vertex.get(symb, False), 2)
    #         #                                  for symb in free_copy_variables])
    #         #         bounds = [(1, expected_average_cov) for foo in range(len(free_copy_variables))]
    #         #         try:
    #         #             copy_result = optimize.minimize(fun=least_square_function_v, x0=initials, jac=False,
    #         #                                             method='SLSQP', bounds=bounds, constraints=constraints, options=opt)
    #         #         except Exception:
    #         #             continue
    #         #         if copy_result.fun < best_fun:
    #         #             best_fun = round(copy_result.fun, 2)
    #         #             copy_results = {tuple(copy_result.x)}
    #         #         elif copy_result.fun == best_fun:
    #         #             copy_results.add(tuple(copy_result.x))
    #         #         else:
    #         #             pass
    #         #     if debug or verbose:
    #         #         if log_handler:
    #         #             log_handler.info("Best function value: " + str(best_fun))
    #         #         else:
    #         #             sys.stdout.write("Best function value: " + str(best_fun) + "\n")
    #         if verbose or debug:
    #             if log_handler:
    #                 log_handler.info("Copy results: " + str(copy_results))
    #             else:
    #                 sys.stdout.write("Copy results: " + str(copy_results) + "\n")
    #         # if len(copy_results) == 1:
    #         #     copy_results = list(copy_results)
    #         # elif len(copy_results) > 1:
    #         #     # draftly sort results by freedom vertices_set
    #         #     copy_results = sorted(copy_results, key=lambda
    #         #         x: sum([(x[go_sym] - self.vertex_to_float_copy[symbols_to_vertex[symb_used]]) ** 2
    #         #                 for go_sym, symb_used in enumerate(free_copy_variables)
    #         #                 if symb_used in symbols_to_vertex]))
    #         # else:
    #         #     raise ProcessingGraphFailed("Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (3)!")
    #     else:
    #         copy_results = []
    #
    #     # if return_new_graphs:
    #     # """ produce all possible vertex copy combinations """
    #     final_results = []
    #     all_copy_sets = set()
    #     # maybe no more multiple results since 2022-12 gekko update
    #     for go_res, copy_result in enumerate([copy_results]):
    #         free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
    #                                     for i, this_copy in enumerate(copy_result)}
    #
    #         """ simplify copy values """
    #         # 2020-02-22 added to avoid multiplicities res such as: [4, 8, 4]
    #         # 2022-12-15 add cluster info to simplify by graph components when the graph is broken
    #         all_copies = []
    #         v_to_cid = {}
    #         for go_id, this_symbol in enumerate(all_v_symbols):
    #             vertex_name = symbols_to_vertex[this_symbol]
    #             v_to_cid[vertex_name] = go_id
    #             this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
    #             if this_copy <= 0:
    #                 raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
    #             all_copies.append(this_copy)
    #         if len(self.vertex_clusters) == 1:
    #             if len(all_copies) == 0:
    #                 raise ProcessingGraphFailed(
    #                     "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
    #             elif len(all_copies) == 1:
    #                 all_copies = [1]
    #             elif min(all_copies) == 1:
    #                 pass
    #             else:
    #                 new_all_copies = reduce_list_with_gcd(all_copies)
    #                 if verbose and new_all_copies != all_copies:
    #                     if log_handler:
    #                         log_handler.info("Estimated copies: " + str(all_copies))
    #                         log_handler.info("Reduced copies: " + str(new_all_copies))
    #                     else:
    #                         sys.stdout.write("Estimated copies: " + str(all_copies) + "\n")
    #                         sys.stdout.write("Reduced copies: " + str(new_all_copies) + "\n")
    #                 all_copies = new_all_copies
    #         else:
    #             for v_cluster in self.vertex_clusters:
    #                 ids = [v_to_cid[_v] for _v in v_cluster]
    #                 component_copies = [all_copies[_id] for _id in ids]
    #                 if len(component_copies) == 0:
    #                     raise ProcessingGraphFailed(
    #                         "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
    #                 elif len(component_copies) == 1:
    #                     component_copies = [1]
    #                 elif min(component_copies) == 1:
    #                     pass
    #                 else:
    #                     new_comp_copies = reduce_list_with_gcd(component_copies)
    #                     if verbose and new_comp_copies != component_copies:
    #                         if log_handler:
    #                             log_handler.info("Estimated copies: " + str(component_copies))
    #                             log_handler.info("Reduced copies: " + str(new_comp_copies))
    #                         else:
    #                             sys.stdout.write("Estimated copies: " + str(component_copies) + "\n")
    #                             sys.stdout.write("Reduced copies: " + str(new_comp_copies) + "\n")
    #                     component_copies = new_comp_copies
    #                 for sequential_id, _id in enumerate(ids):
    #                     all_copies[_id] = component_copies[sequential_id]
    #
    #         all_copies = tuple(all_copies)
    #         if all_copies not in all_copy_sets:
    #             all_copy_sets.add(all_copies)
    #         else:
    #             continue
    #
    #         """ record new copy values """
    #         final_results.append({"graph": deepcopy(self)})
    #         for go_s, this_symbol in enumerate(all_v_symbols):
    #             vertex_name = symbols_to_vertex[this_symbol]
    #             if vertex_name in final_results[go_res]["graph"].vertex_to_copy:
    #                 old_copy = final_results[go_res]["graph"].vertex_to_copy[vertex_name]
    #                 final_results[go_res]["graph"].copy_to_vertex[old_copy].remove(vertex_name)
    #                 if not final_results[go_res]["graph"].copy_to_vertex[old_copy]:
    #                     del final_results[go_res]["graph"].copy_to_vertex[old_copy]
    #             this_copy = all_copies[go_s]
    #             final_results[go_res]["graph"].vertex_to_copy[vertex_name] = this_copy
    #             if this_copy not in final_results[go_res]["graph"].copy_to_vertex:
    #                 final_results[go_res]["graph"].copy_to_vertex[this_copy] = set()
    #             final_results[go_res]["graph"].copy_to_vertex[this_copy].add(vertex_name)
    #
    #         """ re-estimate baseline depth """
    #         total_product = 0.
    #         total_len = 0
    #         for vertex_name in vertices_list:
    #             this_len = self.vertex_info[vertex_name].len \
    #                        * final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
    #             this_cov = self.vertex_info[vertex_name].cov \
    #                        / final_results[go_res]["graph"].vertex_to_copy.get(vertex_name, 1)
    #             total_len += this_len
    #             total_product += this_len * this_cov
    #         final_results[go_res]["cov"] = total_product / total_len
    #     return final_results
    #     # else:
    #     #     """ produce the first-ranked copy combination """
    #     #     free_copy_variables_dict = {free_copy_variables[i]: int(this_copy)
    #     #                                 for i, this_copy in enumerate(copy_results)}
    #     #
    #     #     """ simplify copy values """  # 2020-02-22 added to avoid multiplicities res such as: [4, 8, 4]
    #     #     all_copies = []
    #     #     for this_symbol in all_v_symbols:
    #     #         vertex_name = symbols_to_vertex[this_symbol]
    #     #         this_copy = int(copy_solution[this_symbol].evalf(subs=free_copy_variables_dict, chop=True))
    #     #         if this_copy <= 0:
    #     #             raise ProcessingGraphFailed("Cannot identify copy number of " + vertex_name + "!")
    #     #         all_copies.append(this_copy)
    #     #     if len(all_copies) == 0:
    #     #         raise ProcessingGraphFailed(
    #     #             "Incomplete/Complicated/Unsolvable " + target_name_for_log + " graph (4)!")
    #     #     elif len(all_copies) == 1:
    #     #         all_copies = [1]
    #     #     elif min(all_copies) == 1:
    #     #         pass
    #     #     else:
    #     #         new_all_copies = reduce_list_with_gcd(all_copies)
    #     #         if verbose and new_all_copies != all_copies:
    #     #             if log_handler:
    #     #                 log_handler.info("Estimated copies: " + str(all_copies))
    #     #                 log_handler.info("Reduced copies: " + str(new_all_copies))
    #     #             else:
    #     #                 sys.stdout.write("Estimated copies: " + str(all_copies) + "\n")
    #     #                 sys.stdout.write("Reduced copies: " + str(new_all_copies) + "\n")
    #     #         all_copies = new_all_copies
    #     #
    #     #     """ record new copy values """
    #     #     for go_s, this_symbol in enumerate(all_v_symbols):
    #     #         vertex_name = symbols_to_vertex[this_symbol]
    #     #         if vertex_name in self.vertex_to_copy:
    #     #             old_copy = self.vertex_to_copy[vertex_name]
    #     #             self.copy_to_vertex[old_copy].remove(vertex_name)
    #     #             if not self.copy_to_vertex[old_copy]:
    #     #                 del self.copy_to_vertex[old_copy]
    #     #         this_copy = all_copies[go_s]
    #     #         self.vertex_to_copy[vertex_name] = this_copy
    #     #         if this_copy not in self.copy_to_vertex:
    #     #             self.copy_to_vertex[this_copy] = set()
    #     #         self.copy_to_vertex[this_copy].add(vertex_name)
    #     #
    #     #     if debug or verbose:
    #     #         """ re-estimate baseline depth """
    #     #         total_product = 0.
    #     #         total_len = 0
    #     #         for vertex_name in vertices_list:
    #     #             this_len = self.vertex_info[vertex_name].len \
    #     #                        * self.vertex_to_copy.get(vertex_name, 1)
    #     #             this_cov = self.vertex_info[vertex_name].cov / self.vertex_to_copy.get(vertex_name, 1)
    #     #             total_len += this_len
    #     #             total_product += this_len * this_cov
    #     #         new_val = total_product / total_len
    #     #         if log_handler:
    #     #             log_handler.info("Average " + target_name_for_log + " kmer-coverage = " + str(round(new_val, 2)))
    #     #         else:
    #     #             sys.stdout.write(
    #     #                 "Average " + target_name_for_log + " kmer-coverage = " + str(round(new_val, 2)) + "\n")

    def tag_in_between(self, database_n=None):
        """add those in between the tagged vertices_set to tagged_vertices, which offered the only connection"""
        if database_n is None:
            db_types = sorted(self.tagged_vertices)
        else:
            db_types = [database_n]
        for db_n in db_types:
            updated = True
            candidate_vertices = list(self.vertex_info)
            while updated:
                updated = False
                go_to_v = 0
                while go_to_v < len(candidate_vertices):
                    can_v = candidate_vertices[go_to_v]
                    if can_v in self.tagged_vertices[db_n]:
                        del candidate_vertices[go_to_v]
                        continue
                    else:
                        if sum([bool(c_c) for c_c in self.vertex_info[can_v].connections.values()]) != 2:
                            del candidate_vertices[go_to_v]
                            continue
                        count_nearby_tagged = []
                        for can_end in (True, False):
                            for next_v, next_e in self.vertex_info[can_v].connections[can_end]:
                                # candidate_v is the only output vertex to next_v
                                if next_v in self.tagged_vertices[db_n] and \
                                        len(self.vertex_info[next_v].connections[next_e]) == 1:
                                    count_nearby_tagged.append((next_v, next_e))
                                    break
                        if len(count_nearby_tagged) == 2:
                            del candidate_vertices[go_to_v]
                            # add in between
                            self.tagged_vertices[db_n].add(can_v)
                            if "weight" not in self.vertex_info[can_v].other_attr:
                                self.vertex_info[can_v].other_attr["weight"] = {}
                            if db_n not in self.vertex_info[can_v].other_attr["weight"]:
                                self.vertex_info[can_v].other_attr["weight"][db_n] = 0
                            self.vertex_info[can_v].other_attr["weight"][db_n] += 1 * self.vertex_info[can_v].cov
                            if db_n != "embplant_mt":
                                # Adding extra circle - the contig in-between the sequential repeats
                                # To avoid risk of tagging mt as pt by mistake,
                                # the repeated contig must be at least 2 folds of the nearby tagged contigs
                                near_by_pairs = self.is_sequential_repeat(can_v, return_pair_in_the_trunk_path=False)
                                if near_by_pairs:
                                    checking_new = []
                                    coverage_folds = []
                                    for near_by_p in near_by_pairs:
                                        for (near_v, near_e) in near_by_p:
                                            if (near_v, near_e) not in count_nearby_tagged:
                                                checking_new.append(near_v)
                                                # comment out for improper design: if the untagged is mt
                                                # coverage_folds.append(
                                                #     round(self.vertex_info[can_v].cov /
                                                #           self.vertex_info[near_v].cov, 0))
                                    for near_v, near_e in count_nearby_tagged:
                                        coverage_folds.append(
                                            round(self.vertex_info[can_v].cov /
                                                  self.vertex_info[near_v].cov, 0))
                                    # if coverage folds is
                                    if max(coverage_folds) >= 2:
                                        for extra_v_to_add in set(checking_new):
                                            self.tagged_vertices[db_n].add(extra_v_to_add)
                                            try:
                                                candidate_vertices.remove(extra_v_to_add)
                                            except ValueError:
                                                pass
                                            # when a contig has no weights
                                            if "weight" not in self.vertex_info[extra_v_to_add].other_attr:
                                                self.vertex_info[extra_v_to_add].other_attr["weight"] = {db_n: 0}
                                            # when a contig has weights of other database
                                            if db_n not in self.vertex_info[extra_v_to_add].other_attr["weight"]:
                                                self.vertex_info[extra_v_to_add].other_attr["weight"][db_n] = 0
                                            self.vertex_info[extra_v_to_add].other_attr["weight"][db_n] \
                                                += 1 * self.vertex_info[extra_v_to_add].cov
                            updated = True
                            break
                        else:
                            go_to_v += 1

    def parse_tab_file(self,
                       tab_file,
                       database_name,
                       type_factor,
                       max_gene_gap=300,
                       max_cov_diff=3.,
                       log_handler=None,
                       append_info=False,
                       verbose=False):
        """
        :param tab_file:
        :param database_name:
        :param type_factor:
        :param max_gene_gap:
        :param max_cov_diff:
        :param log_handler:
        :param append_info: not recommended, keep the original information in the vertex info
        :param verbose:
        :return:
        """
        # # parse_csv, every locus only occur in one vertex (removing locations with smaller weight)
        # 2022-12-22 modified for v2: locus can occur in multiple vertices that are linearly continuous

        # 1. parsing to tag_loci
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
                                and self.uni_overlap() and locus_len <= self.uni_overlap():
                            continue
                        # 2022-12-22 added
                        if locus_name not in tag_loci[locus_type]:
                            tag_loci[locus_type][locus_name] = []
                        tag_loci[locus_type][locus_name].append(
                            {"vertex": vertex_name,
                             "len": locus_len,
                             "weight": locus_len * self.vertex_info[vertex_name].cov})
                        # if locus_name in tag_loci[locus_type]:
                        #     new_weight = locus_len * self.vertex_info[vertex_name].cov
                        #     if new_weight > tag_loci[locus_type][locus_name]["weight"]:
                        #         tag_loci[locus_type][locus_name] = {"vertex": vertex_name, "len": locus_len,
                        #                                             "weight": new_weight}
                        # else:
                        #     tag_loci[locus_type][locus_name] = \
                        #         {"vertex": vertex_name,
                        #          "len": locus_len,
                        #          "weight": locus_len * self.vertex_info[vertex_name].cov}

        # 2022-12-22~24 added
        # 2. remove redundant tags can occur in multiple vertices that are not linearly continuous
        #    Under current version, there is no hit-start-end information from csv produced by slim_graph.py
        #    for 1) easy modification and 2) compatible with older versions
        #    So we have to guess the order of vertices in the linear gene
        # TODO: in the future, the best solution is annotating the assembly graph accurately

        sum_tag_loci = {}
        idx_v_cluster = False
        v_to_cluster = {}
        len_cluster = len(self.vertex_clusters)
        # import time
        # time0 = time.time()
        # gmm_time = 0
        for locus_type in tag_loci:
            sum_tag_loci[locus_type] = {}
            for locus_name in tag_loci[locus_type]:
                if len(tag_loci[locus_type][locus_name]) == 1:
                    sum_tag_loci[locus_type][locus_name] = {"vertex": [tag_loci[locus_type][locus_name][0]["vertex"]],
                                                            "weight": [tag_loci[locus_type][locus_name][0]["weight"]]}
                else:
                    if not idx_v_cluster:
                        for go_c, v_clusters in enumerate(self.vertex_clusters):
                            for v_name in v_clusters:
                                v_to_cluster[v_name] = go_c
                    # 2023-01-07 added
                    single_locus_info = tag_loci[locus_type][locus_name]
                    # 2.1 to speed up, remove tags (de-weight) out of the main connected component
                    # if len(single_locus_info) > 10:
                    if verbose and log_handler:
                        log_handler.info("  de-weighting minor-component tags " + locus_type + ":" + locus_name)
                    g_weights = [0.] * len_cluster
                    cluster_to_info_id = {c_id: [] for c_id in range(len_cluster)}
                    for go_r, rec in enumerate(single_locus_info):
                        cluster_id = v_to_cluster[rec["vertex"]]
                        g_weights[cluster_id] += rec["weight"]
                        cluster_to_info_id[cluster_id].append(go_r)
                    max_g_w = max(g_weights)
                    rm_r_ids = []
                    for go_c in range(len_cluster):
                        # arbitrary weight different between connected components
                        if g_weights[go_c] * 20 < max_g_w:
                            rm_r_ids.extend(cluster_to_info_id[go_c])
                    rm_r_ids.sort(reverse=True)
                    if verbose and log_handler:
                        log_handler.info("      " + str(len(rm_r_ids)) + "/" + str(len(single_locus_info)) +
                                         " de-weighted: " + str([single_locus_info[_r]["vertex"] for _r in rm_r_ids]))
                    for go_r in rm_r_ids:
                        del single_locus_info[go_r]

                    # if len(self._get_tagged_merged_paths([_rec["vertex"] for _rec in single_locus_info])) > 1:
                    # 2.2 mark tags of minor coverage as negative
                    if verbose and log_handler:
                        log_handler.info("  negatizing tags based coverage " + locus_type + ":" + locus_name)
                    single_locus_info.sort(key=lambda x: -x["weight"])
                    vertices = [x["vertex"] for x in single_locus_info]
                    # self.get_clusters(limited_vertices=vertices)
                    # maybe increase the vertex weight in the main component
                    coverages = [self.vertex_info[x["vertex"]].cov for x in single_locus_info]
                    # v_weights = [x["weight"] for x in single_locus_info]
                    v_weights = [self.vertex_info[x["vertex"]].len for x in single_locus_info]
                    if verbose and log_handler:
                        log_handler.info("      vertices: " + str(vertices) + "; depths: " + str(coverages) +
                                         "; weights: " + str(v_weights))
                    # timex = time.time()
                    # most time consuming step
                    gmm_scheme = weighted_gmm_with_em_aic(
                        data_array=coverages,
                        data_weights=v_weights,
                        maximum_cluster=5,
                        log_handler=log_handler,
                        verbose_log=verbose)
                    # print(time.time() - timex)
                    # gmm_time += time.time() - timex
                    labels = gmm_scheme["labels"]
                    if log_handler and verbose:
                        log_handler.info("      labels: " + str(list(labels)))
                    l_weights = [0.] * gmm_scheme["cluster_num"]
                    for go_r, lb in enumerate(labels):
                        l_weights[lb] += single_locus_info[go_r]["weight"]
                    selected_lb = l_weights.index(max(l_weights))
                    selected_param = gmm_scheme["parameters"][selected_lb]
                    selected_mu, selected_sigma = selected_param["mu"], selected_param["sigma"]
                    keep_lbs = {go_l
                                for go_l, params in enumerate(gmm_scheme["parameters"])
                                if params["mu"] - selected_mu > -2 * max(selected_sigma, params["sigma"])}
                    # rm_idx = sorted([go_r for go_r, lb in enumerate(labels) if lb not in keep_lbs], reverse=True)
                    ne_idx = [go_r for go_r, lb in enumerate(labels) if lb not in keep_lbs]
                    if verbose and log_handler:
                        log_handler.info("      " + str(len(ne_idx)) + "/" + str(len(single_locus_info)) +
                                         " negatized: " + str([single_locus_info[_r]["vertex"] for _r in ne_idx]))
                    # for go_r in rm_idx:
                    #     del single_locus_info[go_r]
                    for go_r in ne_idx:
                        single_locus_info[go_r]["weight"] = -1 * abs(single_locus_info[go_r]["weight"])

                    # 2.3. remove redundant tags that occur in parallel vertices
                    if verbose and log_handler:
                        log_handler.info("  negatizing parallel tags: " + locus_type + ":" + locus_name)
                    v_to_sl_id = {rec["vertex"]: go_r for go_r, rec in enumerate(single_locus_info)}
                    parallel_vertices_list = self.detect_parallel_vertices(
                        limited_vertices=list(v_to_sl_id),
                        detect_neighbors=False)
                    if parallel_vertices_list:
                        rm_r_ids = set()
                        ne_r_ids = set()
                        for prl_vertices_set in parallel_vertices_list:
                            # sort by weight, then coverage
                            prl_vertices = sorted(
                                prl_vertices_set,
                                key=lambda x: (
                                    -single_locus_info[v_to_sl_id[x[0]]]["weight"],
                                    -self.vertex_info[x[0]].cov))
                            up_v, up_e = prl_vertices[0]
                            up_id = v_to_sl_id[up_v]
                            up_lb = labels[up_id]
                            up_sigma = gmm_scheme["parameters"][up_lb]["sigma"]
                            up_cov = self.vertex_info[up_v].cov
                            for de_name, de_end in prl_vertices[1:]:
                                de_id = v_to_sl_id[de_name]
                                de_lb = labels[de_id]
                                de_sigma = gmm_scheme["parameters"][de_lb]["sigma"]
                                de_cov = self.vertex_info[de_name].cov
                                if abs(de_cov - up_cov) < 2 * max(up_sigma, de_sigma) or \
                                        single_locus_info[de_id]["weight"] / de_cov > \
                                        single_locus_info[up_id]["weight"] / up_cov:
                                    # to be conserved
                                    rm_r_ids.add(de_id)
                                else:
                                    ne_r_ids.add(de_id)
                        if verbose and log_handler:
                            log_handler.info("    (" + str(len(ne_r_ids)) + "+" + str(len(rm_r_ids)) +
                                             ")/" + str(len(single_locus_info)) +
                                             " negatized: " +
                                             str([single_locus_info[_r]["vertex"] for _r in ne_r_ids]) +
                                             " de-weighted: " +
                                             str([single_locus_info[_r]["vertex"] for _r in rm_r_ids]))
                        for go_r in ne_r_ids:
                            single_locus_info[go_r]["weight"] = -1 * abs(single_locus_info[go_r]["weight"])
                        for rm_id in sorted(rm_r_ids, reverse=True):
                            del single_locus_info[rm_id]

                    # 2.4 search for the linear tags maximize the total gene weight
                    if verbose and log_handler:
                        log_handler.info("  linearize " + locus_type + ":" + locus_name)
                    sum_tag_loci[locus_type][locus_name] = \
                        self._find_linear_tags(single_locus_info, max_gene_gap, max_cov_diff, verbose, log_handler)
        # print("gmm cost", gmm_time)
        # print("tagging cost", time.time() - time0)
        # 3. assign information in sum_tag_loci to contigs.other_attr
        # 2022-12-22 modified
        if not append_info:
            # clean previous info
            for vertex_name in self.vertex_info:
                self.vertex_info[vertex_name].other_attr["tags"] = {}
        for locus_type in sum_tag_loci:
            self.tagged_vertices[locus_type] = set()
            self.tagged_vertices[locus_type + "-"] = set()  # negative type
            for locus_name in sum_tag_loci[locus_type]:
                # 2022-12-22 modified
                for vertex_name, loci_weight in zip(sum_tag_loci[locus_type][locus_name]["vertex"],
                                                    sum_tag_loci[locus_type][locus_name]["weight"]):
                    # vertex_name = tag_loci[locus_type][locus_name]["vertex"]
                    # loci_weight = tag_loci[locus_type][locus_name]["weight"]
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
                    # self.tagged_vertices[locus_type].add(vertex_name)

        # 4. clarify locus_type for each contig by comparing weights, and add to self.tagged_vertices
        for vertex_name in self.vertex_info:
            if "weight" in self.vertex_info[vertex_name].other_attr:
                all_weights = [(loc_type, self.vertex_info[vertex_name].other_attr["weight"][loc_type])
                               for loc_type in self.vertex_info[vertex_name].other_attr["weight"]]
                if len(self.vertex_info[vertex_name].other_attr["weight"]) > 1:
                    all_weights.sort(key=lambda x: -x[1])
                best_t, best_w = all_weights[0]
                if best_w > 0:
                    self.tagged_vertices[best_t].add(vertex_name)
                    for next_t, next_w in all_weights[1:]:
                        if next_w * type_factor >= best_w:
                            self.tagged_vertices[next_t].add(vertex_name)
                elif best_w <= 0:
                    for next_t, next_w in all_weights:
                        if next_w < 0:
                            self.tagged_vertices[next_t + "-"].add(vertex_name)

        if database_name not in self.tagged_vertices or len(self.tagged_vertices[database_name]) == 0:
            raise ProcessingGraphFailed("No available " + database_name + " information found in " + tab_file)
        # print("parsing cost", time.time() - time0)

    def _get_tagged_merged_paths(self, tagged_vs):
        vs_to_merge = set(tagged_vs)
        merged_paths = []
        while vs_to_merge:
            check_v = vs_to_merge.pop()
            extend_e = True
            this_path = [(check_v, extend_e)]
            while True:
                this_v, this_e = this_path[-1]
                next_con_tagged = [(_v, _e)
                                   for _v, _e in self.vertex_info[this_v].connections[this_e]
                                   if _v in vs_to_merge]
                # print("this_path", this_path)
                # print("next_con_tagged", next_con_tagged)
                if len(next_con_tagged) == 1:
                    next_v, next_e = next_con_tagged[0]
                    back_con_tagged = [(_v, _e)
                                       for _v, _e in self.vertex_info[next_v].connections[next_e]
                                       if _v == this_v or _v in vs_to_merge]
                    # print("back_con_tagged", back_con_tagged)
                    if back_con_tagged == [(this_v, this_e)]:
                        this_path.append((next_v, not next_e))
                        vs_to_merge.discard(next_v)
                        continue
                if extend_e:
                    this_path = [(_v, not _e) for _v, _e in this_path[::-1]]
                    extend_e = False
                else:
                    break
            merged_paths.append(this_path)
        return merged_paths

    def _find_linear_tags(self, tag_locus_info, max_gene_gap, max_cov_diff, verbose=False, log_handler=None):

        def _try_merge(_current_p_id, _rev_p, _next_p_id, potential_start, middle_gap_p=None):
            # if not the first one to extend/merge
            if count_keep > 0:
                _c_opt = deepcopy(raw_opt)
            else:
                _c_opt = c_opt
            if _rev_p:
                _c_opt["paths"][_current_p_id] = [(_v, not _e)
                                                  for _v, _e in _c_opt["paths"][_current_p_id][::-1]]
            # if check_gene:
            #     print("_c_opt['paths'][_current_p_id]", _c_opt["paths"][_current_p_id])
            #     print('_c_opt["paths"][_next_p_id]', _c_opt["paths"][_next_p_id])
            _nv, _ne = potential_start

            middle_gap_p = [] if not middle_gap_p else middle_gap_p
            if _c_opt["paths"][_next_p_id][0] == potential_start:
                _c_opt["paths"][_current_p_id].extend(middle_gap_p)
                # merge go_p and next_p
                _c_opt["paths"][_current_p_id].extend(_c_opt["paths"][_next_p_id])
                # if check_gene:
                #     print("    merge forward, new paths", _c_opt["paths"][_current_p_id])
            elif _c_opt["paths"][_next_p_id][-1] == (_nv, not _ne):
                _c_opt["paths"][_current_p_id].extend(middle_gap_p)
                # merge go_p and next_p in the reverse
                rev_next_p = [(_v, not _e) for _v, _e in _c_opt["paths"][_next_p_id][::-1]]
                _c_opt["paths"][_current_p_id].extend(rev_next_p)
                # if check_gene:
                #     print("    merge reverse, new paths", _c_opt["paths"][_current_p_id])
            else:
                return False
            del _c_opt["paths"][_next_p_id]
            # update path id
            _c_opt["path_id"] = {_v: _p_id
                                 for _p_id, _p in enumerate(_c_opt["paths"])
                                 for _v, _e in _p}
            # update tuple
            _c_opt["tuple"] = self.standardize_paths(_c_opt["paths"])
            # if not the first one to extend/merge
            if count_keep > 0:
                candidate_options.append(_c_opt)
            return True

        tagged_vs = sorted([_rec["vertex"] for _rec in tag_locus_info if _rec["weight"] > 0])

        # check_gene = False
        # if "323673" in tagged_vs:
        #     check_gene = True

        # merge tagged vertices into paths before linear searching
        merged_paths = self._get_tagged_merged_paths(tagged_vs)
        candidate_options = [{}]
        candidate_options[0]["tuple"] = self.standardize_paths(merged_paths)
        candidate_options[0]["paths"] = [list(_p) for _p in candidate_options[0]["tuple"]]
        candidate_options[0]["path_id"] = {_v: _p_id
                                           for _p_id, _p in enumerate(candidate_options[0]["paths"])
                                           for _v, _e in _p}
        if len(candidate_options[0]["tuple"]) > 1:
            # tagged_set = set(tagged_vs)
            intermediate_combinations = set([])  # to avoid repeated calculation

            # # start_v = sorted(tag_loci[locus_type][locus_name],
            # #                  key=lambda x: (-x["weight"], x["vertex"]))[0]["vertex"]
            # # start_v = tagged_vs[0]
            # # sv_id = tagged_vs.index(start_v)
            # candidate_options = [{"paths": [[(_v, True)] for _v in tagged_vs],
            #                       "path_id": {_v: p_id for p_id, _v in enumerate(tagged_vs)}
            #                       }]
            # # palindromic repeats does not matter, just cause duplicates
            # candidate_options[0]["tuple"] = self.standardize_paths(candidate_options[0]["paths"])

            go_candidate = 0
            while go_candidate < len(candidate_options):
                # if check_gene:
                #     print("go_candidate", go_candidate)
                #     input("")
                if candidate_options[go_candidate]["tuple"] in intermediate_combinations:
                    del candidate_options[go_candidate]  # searched
                else:
                    c_opt = candidate_options[go_candidate]
                    # if check_gene:
                    #     print("c_opt (" + str(len(c_opt["paths"])) + "):", c_opt["paths"])
                    #     input("")
                    intermediate_combinations.add(c_opt["tuple"])
                    extended = False
                    count_keep = 0
                    raw_opt = deepcopy(c_opt)
                    for go_p, this_path in enumerate(list(raw_opt["paths"])):
                        # Problematic
                        # # only start from single-end terminal path, middle path will be extended anyway
                        # next_con_pair = []
                        # next_con_ls_pair = []
                        # next_con_ls_tagged_pair = []
                        # for rev_p in (False, True):
                        #     if rev_p:
                        #         extend_v, extend_e = this_path[0]
                        #         extend_e = not extend_e
                        #     else:
                        #         extend_v, extend_e = this_path[-1]
                        #     next_connections = self.vertex_info[extend_v].connections[extend_e]
                        #     next_con_pair.append(next_connections)
                        #     # constraint the coverage change
                        #     next_connect_ls = [(_n, _e)
                        #                        for _n, _e in next_connections
                        #                        if self.vertex_info[_n].cov / max_cov_diff
                        #                        < self.vertex_info[extend_v].cov
                        #                        < max_cov_diff * self.vertex_info[_n].cov]
                        #     next_con_ls_pair.append(next_connect_ls)
                        #     next_con_ls_tagged_pair.append([_n for _n, _e in next_connect_ls if _n in raw_opt["path_id"]])
                        # if len(this_path) == 1 and bool(next_con_ls_tagged_pair[0]) == bool(next_con_ls_tagged_pair[1]):
                        #     # 1. middle path will be extended by other starts
                        #     # 2. double-ended terminal or self-loop will be isolated anyway
                        #     # so skip those senarios
                        #     continue
                        # rev_p = bool(next_con_ls_tagged_pair[1])

                        for rev_p in (False, True):
                            # if check_gene:
                            #     print("  go_p", go_p, rev_p, this_path)
                            #     # print("  next_con_ls_tagged_pair:", next_con_ls_tagged_pair)
                            #     input("")
                            if rev_p:
                                # palindromic repeats does not matter, just cause duplicates
                                this_path = [(this_v, not this_e) for this_v, this_e in this_path[::-1]]

                            # Problematic
                            # next_connections = next_con_pair[int(rev_p)]
                            # next_connect_ls = next_con_ls_pair[int(rev_p)]

                            extend_v, extend_e = this_path[-1]
                            next_connections = self.vertex_info[extend_v].connections[extend_e]
                            # constraint the coverage change
                            next_connect_ls = [(_n, _e)
                                               for _n, _e in next_connections
                                               if self.vertex_info[_n].cov / max_cov_diff
                                               < self.vertex_info[extend_v].cov
                                               < max_cov_diff * self.vertex_info[_n].cov]
                            # if check_gene:
                            #     print("    next_connections", next_connections)
                            #     input("")
                            if len(next_connect_ls) == 0:
                                continue
                            else:
                                for next_v, next_e in next_connect_ls:
                                    if this_path.count((next_v, not next_e)) \
                                            >= len(self.vertex_info[next_v].connections[next_e]):
                                        # real multiplicity does not matter, just search for the simplest path
                                        # that represent the gene
                                        continue
                                    elif next_v in raw_opt["path_id"]:
                                        next_p_id = raw_opt["path_id"][next_v]
                                        if next_p_id != go_p:  # not self-loop
                                            if _try_merge(go_p, rev_p, next_p_id, (next_v, not next_e)):
                                                # if check_gene:
                                                #     print("    merged with ", next_v, next_p_id, self.vertex_info[next_v].cov)
                                                #     input("")
                                                count_keep += 1
                                                extended = True
                                    else:
                                        # allow gaps
                                        # if check_gene:
                                        #     print("    check gaps")
                                        accumulated_gap = self.vertex_info[next_v].len - \
                                                          next_connections[(next_v, next_e)]
                                        gap_paths = [{"p": [(next_v, not next_e)], "l": accumulated_gap}]
                                        go_g = 0
                                        while go_g < len(gap_paths):
                                            if gap_paths[go_g]["l"] > max_gene_gap:
                                                del gap_paths[go_g]
                                            else:
                                                next_ext_v, next_ext_e = gap_paths[go_g]["p"][-1]
                                                nn_cons = self.vertex_info[next_ext_v].connections[next_ext_e]
                                                # constraint the coverage change
                                                nn_con_ls = [(_n, _e)
                                                             for _n, _e in nn_cons
                                                             if self.vertex_info[_n].cov / max_cov_diff
                                                             < self.vertex_info[next_ext_v].cov
                                                             < max_cov_diff * self.vertex_info[_n].cov]
                                                if len(nn_con_ls) == 0:
                                                    del gap_paths[go_g]
                                                # elif len(nn_con_ls) == 1:
                                                #     nn_v, nn_e = nn_con_ls[0]
                                                #     gap_paths[go_g]["p"].append((nn_v, not nn_e))
                                                #     # either jump to the next gap path option
                                                #     # or add the accumulated gap_length,
                                                #     # which both lead to search termination
                                                #     if nn_v in raw_opt["path_id"]:
                                                #         # nn_p_id = c_opt["path_id"][nn_v]
                                                #         # if nn_p_id == go_p:
                                                #         #     del gap_paths[go_g]
                                                #         go_g += 1
                                                #     else:
                                                #         gap_paths[go_g]["l"] += \
                                                #             self.vertex_info[nn_v].len - nn_cons[(nn_v, nn_e)]
                                                else:
                                                    dup_p = deepcopy(gap_paths[go_g])
                                                    for go_c, (nn_v, nn_e) in enumerate(nn_con_ls):
                                                        if go_c == 0:
                                                            gap_paths[go_g]["p"].append((nn_v, not nn_e))
                                                            if nn_v in raw_opt["path_id"]:
                                                                go_g += 1
                                                            else:
                                                                gap_paths[go_g]["l"] += \
                                                                    self.vertex_info[nn_v].len - nn_cons[(nn_v, nn_e)]
                                                        else:
                                                            if go_c < len(nn_con_ls) - 1:
                                                                this_p = deepcopy(dup_p)
                                                            else:
                                                                this_p = dup_p
                                                            this_p["p"].append((nn_v, not nn_e))
                                                            if nn_v in raw_opt["path_id"]:
                                                                gap_paths.insert(go_g, this_p)
                                                                go_g += 1
                                                            else:
                                                                this_p["l"] += \
                                                                    self.vertex_info[nn_v].len - nn_cons[(nn_v, nn_e)]
                                                                gap_paths.append(this_p)
                                        # if check_gene:
                                        #     print("    gap_paths", gap_paths)
                                        if gap_paths:
                                            if len(gap_paths) == 1:
                                                p_start = nn_v, nn_e = gap_paths[0]["p"][-1]
                                                nn_p_id = raw_opt["path_id"][nn_v]
                                                # if check_gene:
                                                #     print("go_p, rev_p, nn_p_id, p_start, gap_paths[0]['p'][:-1]")
                                                #     print(go_p, rev_p, nn_p_id, p_start, gap_paths[0]["p"][:-1])
                                                if nn_p_id != go_p:
                                                    if _try_merge(go_p, rev_p, nn_p_id, p_start, gap_paths[0]["p"][:-1]):
                                                        # if check_gene:
                                                        #     print("    merged (gap) with ", nn_v, nn_p_id,
                                                        #           self.vertex_info[nn_v].cov)
                                                        #     input("")
                                                        count_keep += 1
                                                        extended = True
                                                    # else:
                                                    #     continue
                                            else:
                                                for go_g, gap_path in enumerate(gap_paths):
                                                    p_start = nn_v, nn_e = gap_path["p"][-1]
                                                    nn_p_id = raw_opt["path_id"][nn_v]
                                                    if nn_p_id == go_p:
                                                        continue
                                                    if _try_merge(go_p, rev_p, nn_p_id, p_start, gap_path["p"][:-1]):
                                                        # if check_gene:
                                                        #     print("    merged (gap) with ", nn_v, nn_p_id,
                                                        #           self.vertex_info[nn_v].cov)
                                                        #     input("")
                                                        count_keep += 1
                                                        extended = True
                            #     if extended:
                            #         break
                            # if extended:
                            #     break
                    if not extended:
                        go_candidate += 1
        tagged_v_w = {_rec["vertex"]: _rec["weight"] for _rec in tag_locus_info}
        # pick the paths with the largest weight
        # sort candidate_options by its maximum path weight (sum of v weights in a path), decreasingly
        candidate_options.sort(key=lambda opt: (max([sum([tagged_v_w.get(_v, 0.)
                                                          for _v, _e, in _p])
                                                     for _p in opt["paths"]]),
                                                opt["tuple"]),
                               reverse=True)
        if verbose:
            for candidate_opt in candidate_options:
                log_handler.info("    paths: " + str(candidate_opt["paths"]))
                log_handler.info("    weights: " + str([sum([tagged_v_w.get(_v, 0.)
                                                              for _v, _e, in _p])
                                                         for _p in candidate_opt["paths"]]))
        best_paths = candidate_options[0]["paths"]
        # pick the path with the largest weight
        best_paths.sort(key=lambda _path: sum([tagged_v_w.get(_v, 0.) for _v, _e, in _path]), reverse=True)
        best_path = best_paths[0]
        if verbose and log_handler:
            log_handler.info("    best_path: " + str(best_path))
        # generate info table
        # labeled_vs = sorted(set([_v for _v, _e in best_path]))
        labeled_vs = set([_v for _v, _e in best_path])
        res_dict = {"vertex": [], "weight": []}
        for this_v in sorted(labeled_vs):
            res_dict["vertex"].append(this_v)
            res_dict["weight"].append(tagged_v_w.get(this_v, 0.))
        for record in tag_locus_info:
            this_v = record["vertex"]
            this_w = record["weight"]
            if this_v not in labeled_vs:
                if this_w < 0:
                    res_dict["vertex"].append(this_v)
                    res_dict["weight"].append(this_w)
                elif self.check_connected(labeled_vs | {this_v}):
                    # connected but not real one should have negative weight
                    res_dict["vertex"].append(this_v)
                    res_dict["weight"].append(-abs(this_w))
        return res_dict

    def filter_by_coverage(self, drop_num=1, database_n="embplant_pt", min_cov_folds=5.,
                           weight_factor=100., min_sigma_factor=0.1, min_cluster=1, terminal_extra_weight=5.,
                           verbose=False, log_handler=None, debug=False):
        changed = False
        if len(self.vertex_info) == 1 and list(self.vertex_info)[0] in self.tagged_vertices[database_n]:
            only_cov = self.vertex_info[list(self.vertex_info)[0]].cov
            return changed, [(only_cov, only_cov * min_sigma_factor)]
        # overlap = self.__overlap if self.__overlap else 0
        # log_min_cov_folds = abs(log(min_cov_folds))
        vertices = sorted(
            self.vertex_info, key=lambda x: (-self.vertex_info[x].other_attr.get("weight", {}).get(database_n, 0), x))
        # 2022-05-06: use the coverage of the contig with the max weight instead of the max coverage
        standard_coverage = self.vertex_info[vertices[0]].cov / self.vertex_to_copy.get(vertices[0], 1)
        if log_handler and (debug or verbose):
            log_handler.info("coverage threshold set: " + str(standard_coverage))
        elif verbose or debug:
            sys.stdout.write("coverage threshold set: " + str(standard_coverage) + "\n")

        # removing coverage with min_cov_folds times lower/greater than tagged_cov
        # if abs(log(self.vertex_info[candidate_v].cov / standard_coverage)) > log_min_cov_folds]
        # 2022-12-05 only remove lower contigs
        removing_low_cov = [candidate_v
                            for candidate_v in vertices
                            if self.vertex_info[candidate_v].cov * min_cov_folds < standard_coverage]
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
            self.estimate_copy_and_depth_by_cov(self.tagged_vertices[database_n],
                                                min_sigma=min_sigma_factor,
                                                debug=debug, log_handler=log_handler,
                                                verbose=verbose, mode=database_n)
        vertices = sorted(self.vertex_info)
        v_coverages = {this_v: self.vertex_info[this_v].cov / self.vertex_to_copy.get(this_v, 1)
                       for this_v in vertices}

        coverages = np.array([v_coverages[this_v] for this_v in vertices])
        cover_weights = np.array([self.vertex_info[this_v].len
                                  # multiply by copy number
                                  * self.vertex_to_copy.get(this_v, 1)
                                  # extra weight to short non-target
                                  * (terminal_extra_weight if self.vertex_info[this_v].is_terminal() else 1)
                                  for this_v in vertices])
        tag_kinds = [tag_kind for tag_kind in self.tagged_vertices if self.tagged_vertices[tag_kind]]
        set_kinds = [tag_kind for tag_kind in tag_kinds if not tag_kind.startswith("-")]
        # introduced 2023-01-11
        ban_kind_set = set([ban_kind for ban_kind in tag_kinds if ban_kind.startswith("-") and ban_kind in set_kinds])
        set_kinds.sort(key=lambda x: x != database_n)
        # force labeled vertex to be in specific cluster, which provide the supervision information for the clustering
        set_cluster = {}
        for go_tag, this_tag in enumerate(set_kinds):
            for v_id, vertex_name in enumerate(vertices):
                if vertex_name in self.tagged_vertices[this_tag]:
                    if v_id not in set_cluster:
                        set_cluster[v_id] = set()
                    set_cluster[v_id].add(go_tag)
        # introduced 2023-01-11
        ban_cluster = {}
        for go_tag, this_tag in enumerate(set_kinds):
            ban_tag = this_tag + "-"
            if ban_tag in ban_kind_set:
                for v_id, vertex_name in enumerate(vertices):
                    if vertex_name in self.tagged_vertices[ban_tag]:
                        if v_id not in ban_cluster:
                            ban_cluster[v_id] = set()
                        ban_cluster[v_id].add(go_tag)
        # min number of clusters
        min_tag_kind = {0}
        for v_id in set_cluster:
            if 0 not in set_cluster[v_id]:
                min_tag_kind |= set_cluster[v_id]
        min_cluster = max(min_cluster, len(min_tag_kind))
        # old way:
        # set_cluster = {v_coverages[tagged_v]: 0 for tagged_v in self.tagged_vertices[mode]}

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
                                              cluster_limited=set_cluster, cluster_bans=ban_cluster,
                                              min_sigma_factor=min_sigma_factor,
                                              log_handler=log_handler, verbose_log=verbose)
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

        # 1
        selected_label_type = list(
            set([lb for go, lb in enumerate(labels) if vertices[go] in self.tagged_vertices[database_n]]))
        if len(selected_label_type) > 1:
            label_weights = {}
            # for lb in selected_label_type:
            #     this_add_up = 0
            #     for go in np.where(labels == lb)[0]:
            #         this_add_up += self.vertex_info[vertices_set[go]].get("weight", {}).get(mode, 0)
            #     label_weights[lb] = this_add_up
            label_weights = {lb: sum([self.vertex_info[vertices[go]].other_attr.get("weight", {}).get(database_n, 0)
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

        # 2
        # exclude_label_type = set()
        # if len(set_kinds) > 1:
        #     for go_l, this_label in enumerate(labels):
        #         for this_tag in set_kinds[1:]:
        #             if vertices_set[go_l] in self.tagged_vertices[this_tag]:
        #                 exclude_label_type.add(this_label)
        #                 break
        # exclude_label_type = sorted(exclude_label_type)
        # if exclude_label_type:
        #     check_ex = 0
        #     while check_ex < len(exclude_label_type):
        #         if exclude_label_type[check_ex] in remained_label_type:
        #             if debug or verbose:
        #                 if log_handler:
        #                     log_handler.info("label " + str(exclude_label_type[check_ex]) + " kept")
        #                 else:
        #                     sys.stdout.write("label " + str(exclude_label_type[check_ex]) + " kept\n")
        #             del exclude_label_type[check_ex]
        #         else:
        #             check_ex += 1

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
        drop_num = max(len(set_kinds) - 1, drop_num)
        dropping_type = dropping_type[:drop_num]
        if debug or verbose:
            if log_handler:
                for lab_tp in dropping_type:
                    if candidate_dropping_label_type[lab_tp] < 0:
                        log_handler.warning("Indistinguishable vertices_set "
                                            + str([vertices[go] for go in np.where(labels == lab_tp)[0]])
                                            + " removed!")
            else:
                for lab_tp in dropping_type:
                    if candidate_dropping_label_type[lab_tp] < 0:
                        sys.stdout.write("Warning: indistinguishable vertices_set "
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

    def exclude_other_hits(self, database_n):
        vertices_to_exclude = []
        for vertex_name in self.vertex_info:
            if "tags" in self.vertex_info[vertex_name].other_attr:
                if database_n in self.vertex_info[vertex_name].other_attr["tags"]:
                    pass
                elif self.vertex_info[vertex_name].other_attr["tags"]:
                    vertices_to_exclude.append(vertex_name)
        if vertices_to_exclude:
            self.remove_vertex(vertices_to_exclude)
            return True
        else:
            return False

    def reduce_to_subgraph(self, bait_vertices, bait_offsets=None,
                           limit_extending_len=None,
                           extending_len_weighted_by_depth=False,
                           verbose=False,
                           log_handler=None):
        """
        :param bait_vertices:
        :param bait_offsets:
        :param limit_extending_len:
        :param extending_len_weighted_by_depth:
        :param verbose:
        :param log_handler:
        :return:
        """
        if bait_offsets is None:
            bait_offsets = {}
        rm_contigs = set()
        rm_sub_ids = []
        # overlap = self.__overlap if self.__overlap else 0
        for go_sub, vertices in enumerate(self.vertex_clusters):
            for vertex in sorted(vertices):
                if vertex in bait_vertices:
                    break
            else:
                rm_sub_ids.append(go_sub)
                rm_contigs.update(vertices)
        # rm vertices_set
        if rm_contigs and verbose and log_handler:
            log_handler.info("removing clusters without baits(" + str(len(rm_contigs)) + "):" + str(rm_contigs))
        self.remove_vertex(rm_contigs, update_cluster=False)
        # rm clusters
        for sub_id in rm_sub_ids[::-1]:
            del self.vertex_clusters[sub_id]
        # searching within a certain length scope
        if limit_extending_len not in (None, inf):
            if extending_len_weighted_by_depth:
                explorers = {(v_n, v_e): (limit_extending_len - bait_offsets.get((v_n, v_e), 0),
                                          self.vertex_info[v_n].cov)
                             for v_n in set(bait_vertices)
                             for v_e in (True, False)}
                best_explored_record = {}
                # explore all minimum distances starting from the bait_vertices
                while True:
                    changed = False
                    for (this_v, this_e), (quota_len, base_cov) in sorted(explorers.items()):
                        # if there's any this_v active: quota_len>0 AND (not_recorded OR recorded_changed))
                        if quota_len > 0 and \
                                (quota_len, base_cov) != best_explored_record.get((this_v, this_e), 0):
                            changed = True
                            best_explored_record[(this_v, this_e)] = (quota_len, base_cov)
                            for (next_v, next_e), this_overlap in self.vertex_info[this_v].connections[this_e].items():
                                # not the starting vertices_set
                                if next_v not in bait_vertices:
                                    new_quota_len = quota_len - (self.vertex_info[next_v].len - this_overlap) * \
                                                    max(1, self.vertex_info[next_v].cov / base_cov)
                                    # if next_v is active: quota_len>0 AND (not_explored OR larger_len))
                                    next_p = (next_v, not next_e)
                                    if new_quota_len > 0 and \
                                            (next_p not in explorers or
                                             # follow the bait contigs with higher coverage:
                                             # replace new_quota_len > explorers[next_p][0]): with
                                             new_quota_len * base_cov > explorers[next_p][0] * explorers[next_p][1]):
                                        explorers[next_p] = (new_quota_len, base_cov)
                    if not changed:
                        break  # if no this_v active, stop the exploring
            else:
                explorers = {(v_n, v_e): limit_extending_len - bait_offsets.get((v_n, v_e), 0)
                             for v_n in set(bait_vertices)
                             for v_e in (True, False)}
                best_explored_record = {}
                # explore all minimum distances starting from the bait_vertices
                while True:
                    changed = False
                    for (this_v, this_e), quota_len in sorted(explorers.items()):
                        # if there's any this_v active: quota_len>0 AND (not_recorded OR recorded_changed))
                        if quota_len > 0 and quota_len != best_explored_record.get((this_v, this_e), None):
                            changed = True
                            best_explored_record[(this_v, this_e)] = quota_len
                            # for this_direction in (True, False):
                            for (next_v, next_e), this_overlap in self.vertex_info[this_v].connections[this_e].items():
                                # not the starting vertices_set
                                if next_v not in bait_vertices:
                                    new_quota_len = quota_len - (self.vertex_info[next_v].len - this_overlap)
                                    # if next_v is active: quota_len>0 AND (not_explored OR larger_len))
                                    next_p = (next_v, not next_e)
                                    if new_quota_len > explorers.get(next_p, 0):
                                        explorers[next_p] = new_quota_len
                    if not changed:
                        break  # if no this_v active, stop the exploring
            accepted = {candidate_v for (candidate_v, candidate_e) in explorers}
            rm_contigs = {candidate_v for candidate_v in self.vertex_info if candidate_v not in accepted}
            self.remove_vertex(rm_contigs, update_cluster=True)

    def generate_consensus_vertex(self, vertices, directions, copy_tags=True, check_parallel_vertices=True,
                                  log_handler=None):
        if self.__uni_overlap is None:
            raise ProcessingGraphFailed("Consensus is only applicable to graph with a universal overlap value!")
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
            self.vertex_info[vertices[0]].merging_history = \
                MergingHistory(
                    [(ConsensusHistory([(v_n, v_e) for v_n, v_e in zip(vertices, directions)]), directions[0])])
            new_vertex = str(self.vertex_info[vertices[0]].merging_history)
            # new_vertex = "(" + "|".join(vertices_set) + ")"
            self.vertex_info[new_vertex] = deepcopy(self.vertex_info[vertices[0]])
            self.vertex_info[new_vertex].name = new_vertex
            self.vertex_info[new_vertex].cov = sum([self.vertex_info[v].cov for v in vertices])
            self.vertex_info[new_vertex].fastg_form_name = None
            # if "long" in self.vertex_info[new_vertex]:
            #     del self.vertex_info[new_vertex]["long"]

            # self.merging_history[new_vertex] = set()
            # for candidate_v in vertices_set:
            #     if candidate_v in self.merging_history:
            #         for sub_v_n in self.merging_history[candidate_v]:
            #             self.merging_history[new_vertex].add(sub_v_n)
            #     else:
            #         self.merging_history[new_vertex].add(candidate_v)
            # for candidate_v in vertices_set:
            #     if candidate_v in self.merging_history:
            #         del self.merging_history[candidate_v]

            for new_end in (True, False):
                for n_n_v, n_n_e in self.vertex_info[new_vertex].connections[new_end]:
                    self.vertex_info[n_n_v].connections[n_n_e][(new_vertex, new_end)] = \
                        self.vertex_info[new_vertex].connections[new_end][(n_n_v, n_n_e)]

            consensus_s = generate_consensus(
                *[self.vertex_info[v].seq[directions[go]] for go, v in enumerate(vertices)])
            self.vertex_info[new_vertex].seq[directions[0]] = consensus_s
            self.vertex_info[new_vertex].seq[not directions[0]] = complementary_seq(consensus_s)
            if copy_tags:
                for db_n in self.tagged_vertices:
                    if vertices[0] in self.tagged_vertices[db_n]:
                        self.tagged_vertices[db_n].add(new_vertex)
                        self.tagged_vertices[db_n].remove(vertices[0])

            # tags
            if copy_tags:
                for other_vertex in vertices[1:]:
                    if "tags" in self.vertex_info[other_vertex].other_attr:
                        if "tags" not in self.vertex_info[new_vertex].other_attr:
                            self.vertex_info[new_vertex].other_attr["tags"] = \
                                deepcopy(self.vertex_info[other_vertex].other_attr["tags"])
                        else:
                            for db_n in self.vertex_info[other_vertex].other_attr["tags"]:
                                if db_n not in self.vertex_info[new_vertex].other_attr["tags"]:
                                    self.vertex_info[new_vertex].other_attr["tags"][db_n] \
                                        = deepcopy(self.vertex_info[other_vertex].other_attr["tags"][db_n])
                                else:
                                    self.vertex_info[new_vertex].other_attr["tags"][db_n] \
                                        |= self.vertex_info[other_vertex].other_attr["tags"][db_n]
                    if "weight" in self.vertex_info[other_vertex].other_attr:
                        if "weight" not in self.vertex_info[new_vertex].other_attr:
                            self.vertex_info[new_vertex].other_attr["weight"] \
                                = deepcopy(self.vertex_info[other_vertex].other_attr["weight"])
                        else:
                            for db_n in self.vertex_info[other_vertex].other_attr["weight"]:
                                if db_n not in self.vertex_info[new_vertex].other_attr["weight"]:
                                    self.vertex_info[new_vertex].other_attr["weight"][db_n] \
                                        = self.vertex_info[other_vertex].other_attr["weight"][db_n]
                                else:
                                    self.vertex_info[new_vertex].other_attr["weight"][db_n] \
                                        += self.vertex_info[other_vertex].other_attr["weight"][db_n]
                    for db_n in self.tagged_vertices:
                        if other_vertex in self.tagged_vertices[db_n]:
                            self.tagged_vertices[db_n].add(new_vertex)
                            self.tagged_vertices[db_n].remove(other_vertex)
            self.remove_vertex(vertices)
            if log_handler:
                log_handler.info("Consensus made: " + new_vertex)
            else:
                log_handler.info("Consensus made: " + new_vertex + "\n")

    def processing_polymorphism(self, database_name, average_depth=None, limited_vertices=None,
                                contamination_depth=3., contamination_similarity=0.95,
                                degenerate=False, degenerate_depth=1.5, degenerate_similarity=0.98, warning_count=4,
                                only_keep_max_cov=False, verbose=False, debug=False, log_handler=None):
        if average_depth is None:
            average_depth, ave_std = self.estimate_copy_and_depth_by_cov(
                self.tagged_vertices[database_name], debug=debug, log_handler=log_handler,
                verbose=verbose, mode=database_name)
        else:
            average_depth = float(average_depth)

        parallel_vertices_list = self.detect_parallel_vertices(limited_vertices=limited_vertices)
        # overlap = self.__overlap if self.__overlap else 0
        if verbose or debug:
            if log_handler:
                log_handler.info("detected parallel vertices_set " + str(parallel_vertices_list))
            else:
                sys.stdout.write("detected parallel vertices_set " + str(parallel_vertices_list) + "\n")

        degenerate_depth = abs(log(degenerate_depth))
        contamination_depth = abs(log(contamination_depth))
        contamination_dif = 1 - contamination_similarity
        degenerate_dif = 1 - degenerate_similarity

        removing_irrelevant_v = set()
        removing_contaminating_v = set()
        count_contamination_or_degenerate = 0
        count_using_only_max = 0
        sub_sampling = 10000
        half_kmer = None if self.uni_overlap() is None else int(self.uni_overlap() / 2)
        for prl_vertices in parallel_vertices_list:
            this_contamination_or_polymorphic = False
            this_using_only_max = False
            # sort by weight, then coverage
            prl_vertices = sorted(
                prl_vertices,
                key=lambda x: (
                    -self.vertex_info[x[0]].other_attr.get("weight", {database_name: 0.}).get(database_name, 0.),
                    -self.vertex_info[x[0]].cov))
            max_cov_vertex, direction_remained = prl_vertices.pop(0)
            max_cov_seq = self.vertex_info[max_cov_vertex].seq[direction_remained]
            max_cov = self.vertex_info[max_cov_vertex].cov
            polymorphic_vertices_with_directions = {(max_cov_vertex, direction_remained)}
            # hard to clearly identify the biological factors
            for this_v, this_direction in prl_vertices:
                this_seq = self.vertex_info[this_v].seq[this_direction]
                this_cov = self.vertex_info[this_v].cov
                if abs(log(this_cov / max_cov)) > contamination_depth:
                    if abs(len(this_seq) - len(max_cov_seq)) / float(len(this_seq)) <= contamination_dif:
                        # too long to calculate, too long to be polymorphic
                        if max(len(max_cov_seq), len(this_seq)) > 1E5 or self.uni_overlap() is None:
                            # removing_irrelevant_v.add(this_v)
                            pass
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
                elif degenerate and abs(log(this_cov / max_cov)) < degenerate_depth:
                    if abs(len(this_seq) - len(max_cov_seq)) * 2 / float(
                            len(this_seq) + len(max_cov_seq)) <= degenerate_dif:
                        # too long to calculate, too long to be polymorphic
                        if max(len(max_cov_seq), len(this_seq)) > 1E5 or self.uni_overlap() is None:
                            continue
                        else:
                            seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                            seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                            if len(seq_1) + len(seq_2) > 6000:
                                base_dif = find_string_difference_regional_kmer_counting(seq_1, seq_2)
                            else:
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
                            elif only_keep_max_cov and float(base_dif) * 2 / (
                                    len(seq_1) + len(seq_2)) < contamination_dif:
                                removing_irrelevant_v.add(this_v)
                                this_contamination_or_polymorphic = True
                                this_using_only_max = True
                    elif only_keep_max_cov:
                        if self.uni_overlap() is None:
                            removing_irrelevant_v.add(this_v)
                            this_contamination_or_polymorphic = True
                            this_using_only_max = True
                        else:
                            seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                            seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                            if len(seq_1) + len(seq_2) > 6000 and \
                                    2 * abs(len(seq_1) - len(seq_2)) / float(len(seq_1) + len(seq_2)) <= degenerate_dif:
                                base_dif = find_string_difference_regional_kmer_counting(seq_1, seq_2)
                            else:
                                base_dif, proper_end = find_string_difference(seq_1, seq_2,
                                                                              max(2, int(len(seq_1) * 0.005)))
                            if float(base_dif) * 2 / (len(seq_1) + len(seq_2)) < contamination_dif:
                                removing_irrelevant_v.add(this_v)
                                this_contamination_or_polymorphic = True
                                this_using_only_max = True
                elif only_keep_max_cov:
                    if self.uni_overlap() is None:
                        removing_irrelevant_v.add(this_v)
                        this_contamination_or_polymorphic = True
                        this_using_only_max = True
                    else:
                        seq_1 = max_cov_seq[half_kmer: min(half_kmer + sub_sampling, len(max_cov_seq) - half_kmer)]
                        seq_2 = this_seq[half_kmer: min(half_kmer + sub_sampling, len(this_seq) - half_kmer)]
                        if len(seq_1) + len(seq_2) > 6000 and \
                                2 * abs(len(seq_1) - len(seq_2)) / float(len(seq_1) + len(seq_2)) <= degenerate_dif:
                            base_dif = find_string_difference_regional_kmer_counting(seq_1, seq_2)
                        else:
                            base_dif, proper_end = find_string_difference(seq_1, seq_2, max(2, int(len(seq_1) * 0.005)))
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
            contaminating_weight = np.array([len(self.vertex_info[con_v].seq[True])
                                             for con_v in removing_contaminating_v])
            for candidate_rm_v in list(removing_contaminating_v):  # fixed in 2022-12-30
                if candidate_rm_v in self.tagged_vertices[database_name]:
                    removing_contaminating_v.remove(candidate_rm_v)
            self.remove_vertex(removing_contaminating_v)
            cont_mean, cont_std = weighted_mean_and_std(contaminating_cov, contaminating_weight)
            cut_off_min, cut_off_max = cont_mean - cont_std, cont_mean + cont_std
            # remove (low cov and terminal vertex)
            removing_below_cut_off = []
            for del_v in self.vertex_info:
                if cut_off_min < self.vertex_info[del_v].cov < cut_off_max:
                    if sum([bool(cnn) for cnn in self.vertex_info[del_v].connections.values()]) < 2 \
                            and del_v not in self.tagged_vertices[database_name]:
                        removing_below_cut_off.append(del_v)
            self.remove_vertex(removing_below_cut_off)
            if verbose or debug:
                if log_handler:
                    log_handler.info("removing contaminating vertices_set: " + " ".join(list(removing_contaminating_v)))
                    log_handler.info("removing contaminating-like vertices_set: " + " ".join(list(removing_below_cut_off)))
                else:
                    sys.stdout.write(
                        "removing contaminating vertices_set: " + " ".join(list(removing_contaminating_v)) + "\n")
                    sys.stdout.write(
                        "removing contaminating-like vertices_set: " + " ".join(list(removing_below_cut_off)) + "\n")
        if removing_irrelevant_v:
            for candidate_rm_v in list(removing_irrelevant_v):
                if candidate_rm_v in self.tagged_vertices[database_name]:
                    removing_irrelevant_v.remove(candidate_rm_v)
            if removing_irrelevant_v:
                self.remove_vertex(removing_irrelevant_v)
                if verbose or debug:
                    # if log_handler:
                    log_handler.info("removing parallel vertices_set: " + " ".join(list(removing_irrelevant_v)))
                    # else:
                    # sys.stdout.write("removing parallel vertices_set: " + " ".join(list(removing_irrelevant_v)) + "\n")
        if count_contamination_or_degenerate >= warning_count:
            # if log_handler:
            log_handler.warning("The graph might suffer from contamination or polymorphism!")
            if count_using_only_max and removing_irrelevant_v:
                log_handler.warning("Only the contig with the max cov was kept for each of those " +
                                    str(count_using_only_max) + " polymorphic loci.")
            # else:
            #     sys.stdout.write("Warning: The graph might suffer from contamination or polymorphism!")
            #     if count_using_only_max:
            #         sys.stdout.write("Warning: Only the contig with the max cov was kept for each of those " +
            #                          str(count_using_only_max) + " polymorphic loci.\n")

    def find_target_graph(self,
                          # tab_file,
                          db_name,
                          mode="embplant_pt",
                          # type_factor=3,
                          weight_factor=100.0,
                          min_sigma_factor=0.1,
                          expected_max_size=inf,
                          expected_min_size=0,
                          # max_contig_multiplicity=8,
                          hard_cov_threshold=5., contamination_depth=3., contamination_similarity=0.95,
                          degenerate=True, degenerate_depth=1.5, degenerate_similarity=0.98, only_keep_max_cov=True,
                          min_single_copy_percent=50, meta=False,
                          broken_graph_allowed=False,
                          selected_graph=None,
                          temp_graph=None, verbose=True,
                          read_len_for_log=None, kmer_for_log=None,
                          log_handler=None, debug=False):
        """
        :param db_name:
        :param mode:
        :param weight_factor:
        :param min_sigma_factor:
        :param expected_max_size:
        :param expected_min_size:
        :param hard_cov_threshold:
        :param contamination_depth:
        :param contamination_similarity:
        :param degenerate:
        :param degenerate_depth:
        :param degenerate_similarity:
        :param only_keep_max_cov:
        :param min_single_copy_percent: [0-100]
        :param broken_graph_allowed:
        :param selected_graph:
        :param temp_graph:
        :param verbose:
        :param read_len_for_log:
        :param kmer_for_log:
        :param log_handler:
        :param debug:
        :return:
        """

        # overlap = self.__overlap if self.__overlap else 0

        # def log_target_res(final_res_combinations_inside):
        #     echo_graph_id = int(bool(len(final_res_combinations_inside) - 1))
        #     for go_res, final_res_one in enumerate(final_res_combinations_inside):
        #         this_graph = final_res_one["graph"]
        #         this_k_cov = round(final_res_one["cov"], 3)
        #         if read_len_for_log and kmer_for_log:
        #             this_b_cov = round(this_k_cov * read_len_for_log / (read_len_for_log - kmer_for_log + 1), 3)
        #         else:
        #             this_b_cov = None
        #         if log_handler:
        #             if echo_graph_id:
        #                 log_handler.info("Graph " + str(go_res + 1))
        #             for vertex_set in sorted(this_graph.vertex_clusters):
        #                 copies_in_a_set = {this_graph.vertex_to_copy[v_name] for v_name in vertex_set}
        #                 if copies_in_a_set != {1}:
        #                     for in_vertex_name in sorted(vertex_set):
        #                         log_handler.info("Vertex_" + in_vertex_name + " #copy = " +
        #                                          str(this_graph.vertex_to_copy.get(in_vertex_name, 1)))
        #             cov_str = " kmer-coverage" if bool(self.uni_overlap()) else " coverage"
        #             log_handler.info("Average " + mode + cov_str +
        #                              ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + "%.1f" % this_k_cov)
        #             if this_b_cov:
        #                 log_handler.info("Average " + mode + " base-coverage" +
        #                                  ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + "%.1f" % this_b_cov)
        #         else:
        #             if echo_graph_id:
        #                 sys.stdout.write("Graph " + str(go_res + 1) + "\n")
        #             for vertex_set in sorted(this_graph.vertex_clusters):
        #                 copies_in_a_set = {this_graph.vertex_to_copy[v_name] for v_name in vertex_set}
        #                 if copies_in_a_set != {1}:
        #                     for in_vertex_name in sorted(vertex_set):
        #                         sys.stdout.write("Vertex_" + in_vertex_name + " #copy = " +
        #                                          str(this_graph.vertex_to_copy.get(in_vertex_name, 1)) + "\n")
        #             cov_str = " kmer-coverage" if bool(self.uni_overlap()) else " coverage"
        #             sys.stdout.write("Average " + mode + cov_str +
        #                              ("(" + str(go_res + 1) + ")") * echo_graph_id + " = " + "%.1f" % this_k_cov + "\n")
        #             if this_b_cov:
        #                 sys.stdout.write("Average " + mode + " base-coverage" + ("(" + str(go_res + 1) + ")") *
        #                                  echo_graph_id + " = " + "%.1f" % this_b_cov + "\n")

        if temp_graph:
            if temp_graph.endswith(".gfa"):
                temp_csv = temp_graph[:-3] + "csv"
            elif temp_graph.endswith(".fastg"):
                temp_csv = temp_graph[:-5] + "csv"
            elif temp_graph.endswith(".fasta"):
                temp_csv = temp_graph[:-5] + "csv"
            else:
                temp_csv = temp_graph + ".csv"
        else:
            temp_csv = None
        count_all_temp = [1]

        def add_temp_id(old_temp_file, extra_str=""):
            if old_temp_file.endswith(".gfa"):
                return old_temp_file[:-4] + extra_str + ".gfa"
            elif old_temp_file.endswith(".csv"):
                return old_temp_file[:-4] + extra_str + ".csv"
            else:
                return old_temp_file + extra_str

        def write_temp_out(_assembly, _database_name, _temp_graph, _temp_csv, step_tag):
            if _temp_graph:
                tmp_graph_1 = add_temp_id(_temp_graph, ".%02d.%s" % (count_all_temp[0], step_tag))
                tmp_csv_1 = add_temp_id(_temp_csv, ".%02d.%s" % (count_all_temp[0], step_tag))
                if verbose:
                    if log_handler:
                        log_handler.info("Writing out temp graph (%s): %s" % (step_tag, tmp_graph_1))
                    else:
                        sys.stdout.write("Writing out temp graph (%s): %s" % (step_tag, tmp_graph_1) + "\n")
                _assembly.write_to_gfa(tmp_graph_1)
                if _database_name in ("embplant_pt", "embplant_mt"):
                    _database_name = ["embplant_pt", "embplant_mt"]
                else:
                    _database_name = [_database_name]
                _assembly.write_out_tags(_database_name, tmp_csv_1)
                count_all_temp[0] += 1

        def write_selected(_assembly, _selected_graph):
            if _selected_graph is None:
                pass
            else:
                # write out selected graph
                log_handler.info("Output selected graph file " + str(_selected_graph))
                _assembly.write_to_gfa(_selected_graph)
                if mode in ("embplant_pt", "embplant_mt"):
                    _this_dbs = ["embplant_pt", "embplant_mt"]
                else:
                    _this_dbs = [mode]
                _assembly.write_out_tags(_this_dbs, _selected_graph[:-3] + "csv")

        def check_remaining_singleton():
            if len(new_assembly.vertex_info) == 1:
                the_only_v = list(new_assembly.vertex_info)[0]
                if the_only_v in new_assembly.tagged_vertices[db_name]:
                    if new_assembly.vertex_info[the_only_v].is_self_loop() or broken_graph_allowed:
                        return True
                    else:
                        if verbose:
                            raise ProcessingGraphFailed("Linear graph: " + the_only_v + "! # tags: " +
                                                        str(new_assembly.vertex_info[the_only_v].other_attr.
                                                            get("tags", {db_name: ""})[db_name]))
                        else:
                            raise ProcessingGraphFailed("Linear graph")

        if broken_graph_allowed and not meta:
            weight_factor = 10000.

        # if meta:
        #     try:
        #         self.parse_tab_file(
        #             tab_file,
        #             database_name=db_name,
        #             type_factor=type_factor,
        #             max_gene_gap=250,
        #             max_cov_diff=hard_cov_threshold,
        #             verbose=verbose,
        #             log_handler=log_handler)
        #     except ProcessingGraphFailed:
        #         return []
        # else:
        #     self.parse_tab_file(
        #         tab_file,
        #         database_name=db_name,
        #         type_factor=type_factor,
        #         max_gene_gap=250,
        #         max_cov_diff=hard_cov_threshold,  #  contamination_depth?
        #         verbose=verbose,
        #         log_handler=log_handler)

        new_assembly = deepcopy(self)
        is_reasonable_res = False
        data_contains_outlier = False
        try:
        # if True:
            while not is_reasonable_res:
                is_reasonable_res = True
                if verbose or debug:
                    if log_handler:
                        log_handler.info("tagged vertices_set: " + str(sorted(new_assembly.tagged_vertices[db_name])))
                        log_handler.info("tagged coverage: " +
                                         str(["%.1f" % new_assembly.vertex_info[log_v].cov
                                              for log_v in sorted(new_assembly.tagged_vertices[db_name])]))
                    else:
                        sys.stdout.write("tagged vertices_set: " + str(sorted(new_assembly.tagged_vertices[db_name]))
                                         + "\n")
                        sys.stdout.write("tagged coverage: " +
                                         str(["%.1f" % new_assembly.vertex_info[log_v].cov
                                              for log_v in sorted(new_assembly.tagged_vertices[db_name])]) + "\n")
                new_assembly.merge_all_possible_vertices()
                new_assembly.tag_in_between()
                write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "a")
                changed = True
                count_large_round = 0
                while changed:
                    count_large_round += 1
                    if verbose or debug:
                        if log_handler:
                            log_handler.info(
                                "===================== " + str(count_large_round) + " =====================")
                        else:
                            sys.stdout.write(
                                "===================== " + str(count_large_round) + " =====================\n")
                    changed = False
                    cluster_trimmed = True
                    while cluster_trimmed:
                        # remove low coverages
                        first_round = True
                        delete_those_vertices = set()
                        # parameters = []
                        this_del = False
                        new_ave_cov, ave_std = new_assembly.estimate_copy_and_depth_by_cov(
                            new_assembly.tagged_vertices[db_name], min_sigma=min_sigma_factor,
                            debug=debug, log_handler=log_handler,
                            verbose=verbose, mode=mode)
                        while first_round or delete_those_vertices or this_del:
                            if data_contains_outlier:
                                this_del, foo_parameters = \
                                    new_assembly.filter_by_coverage(database_n=db_name,
                                                                    weight_factor=weight_factor,
                                                                    min_cov_folds=hard_cov_threshold,
                                                                    min_sigma_factor=min_sigma_factor,
                                                                    min_cluster=2, log_handler=log_handler,
                                                                    verbose=verbose, debug=debug)
                                data_contains_outlier = False
                                if not this_del:
                                    raise ProcessingGraphFailed(
                                        "Unable to generate result with single copy vertex percentage < {}%"
                                            .format(min_single_copy_percent))
                            else:
                                this_del, foo_parameters = \
                                    new_assembly.filter_by_coverage(database_n=db_name,
                                                                    weight_factor=weight_factor,
                                                                    min_cov_folds=hard_cov_threshold,
                                                                    min_sigma_factor=min_sigma_factor,
                                                                    log_handler=log_handler, verbose=verbose,
                                                                    debug=debug)
                            if verbose or debug:
                                if log_handler:
                                    log_handler.info("tagged vertices_set: " +
                                                     str(sorted(new_assembly.tagged_vertices[db_name])))
                                    log_handler.info("tagged coverage: " +
                                                     str(["%.1f" % new_assembly.vertex_info[log_v].cov
                                                          for log_v
                                                          in sorted(new_assembly.tagged_vertices[db_name])]))
                                else:
                                    sys.stdout.write("tagged vertices_set: " +
                                                     str(sorted(new_assembly.tagged_vertices[db_name])) + "\n")
                                    log_handler.info("tagged coverage: " +
                                                     str(["%.1f" % new_assembly.vertex_info[log_v].cov
                                                          for log_v
                                                          in
                                                          sorted(new_assembly.tagged_vertices[db_name])]) + "\n")
                            if this_del and temp_graph:
                                write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "b")
                            new_ave_cov, ave_std = new_assembly.estimate_copy_and_depth_by_cov(
                                new_assembly.tagged_vertices[db_name], min_sigma=min_sigma_factor,
                                debug=debug, log_handler=log_handler,
                                verbose=verbose, mode=mode)
                            first_round = False

                        if new_assembly.exclude_other_hits(database_n=db_name):
                            changed = True

                        cluster_trimmed = False

                        if len(new_assembly.vertex_clusters) == 0:
                            raise ProcessingGraphFailed("No available " + mode + " components detected!")
                        elif len(new_assembly.vertex_clusters) == 1:
                            pass
                        else:
                            cluster_weights = \
                                [sum([new_assembly.vertex_info[x_v].other_attr["weight"][db_name]
                                      for x_v in x
                                      if
                                      "weight" in new_assembly.vertex_info[x_v].other_attr
                                      and
                                      db_name in new_assembly.vertex_info[x_v].other_attr["weight"]])
                                 for x in new_assembly.vertex_clusters]
                            if verbose and log_handler:
                                log_handler.info("cluster_weights: " + str(cluster_weights))
                            best = max(cluster_weights)
                            best_id = cluster_weights.index(best)
                            if broken_graph_allowed:
                                id_remained = {best_id}
                                for j, w in enumerate(cluster_weights):
                                    if w * weight_factor > best:
                                        id_remained.add(j)
                                    else:
                                        for del_v in new_assembly.vertex_clusters[j]:
                                            if del_v in new_assembly.tagged_vertices[db_name]:
                                                new_cov = new_assembly.vertex_info[del_v].cov
                                                # 2023-01-04 modified
                                                if abs(new_cov - new_ave_cov) < 3 * ave_std:
                                                    id_remained.add(j)
                                                # for mu, sigma in parameters:
                                                #     if abs(new_cov - mu) < sigma:
                                                #         id_remained.add(j)
                                                #         break
                                            if j in id_remained:
                                                break
                            else:
                                # chose the target cluster (best rank)
                                id_remained = {best_id}
                                temp_cluster_weights = deepcopy(cluster_weights)
                                del temp_cluster_weights[best_id]
                                second = max(temp_cluster_weights)
                                if best < second * weight_factor:
                                    write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "c")
                                    raise ProcessingGraphFailed(
                                        "Multiple isolated " + mode + " components detected! Broken or contamination?")
                                for j, w in enumerate(cluster_weights):
                                    if w == second:
                                        for del_v in new_assembly.vertex_clusters[j]:
                                            if del_v in new_assembly.tagged_vertices[db_name]:
                                                new_cov = new_assembly.vertex_info[del_v].cov
                                                # for debug
                                                # print(new_cov)
                                                # print(parameters)
                                                # 2023-01-04 modified
                                                if abs(new_cov - new_ave_cov) < 3 * ave_std:
                                                    raise ProcessingGraphFailed(
                                                        "Complicated graph: please check around EDGE_" + del_v + "!"
                                                        "# tags: " +
                                                        str(new_assembly.vertex_info[del_v].other_attr.
                                                            get("tags", {db_name: ""})[db_name]))
                                                else:
                                                    if (verbose or debug) and log_handler:
                                                        log_handler.warning(
                                                            "removing tagged but low-coverage isolated contig: " +
                                                            del_v + ":" +
                                                            str(new_assembly.vertex_info[del_v].other_attr["tags"]))
                                                # for mu, sigma in parameters:
                                                #     if abs(new_cov - mu) < sigma:
                                                #         write_temp_out(new_assembly, db_name,
                                                #                        temp_graph, temp_csv, "d")
                                                #         raise ProcessingGraphFailed(
                                                #             "Complicated graph: please check around EDGE_" + del_v + "!"
                                                #             "# tags: " +
                                                #             str(new_assembly.vertex_info[del_v].other_attr.
                                                #                 get("tags", {db_name: ""})[db_name]))

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

                                write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "e")
                        if check_remaining_singleton():
                            break

                    # merge vertices_set
                    new_assembly.merge_all_possible_vertices()
                    new_assembly.tag_in_between()
                    if check_remaining_singleton():
                        break

                    # no tip contigs allowed
                    if broken_graph_allowed:
                        pass
                    else:
                        total_weight = sum([new_assembly.vertex_info[x_v].other_attr["weight"][db_name]
                                            for x_v in new_assembly.vertex_info
                                            if
                                            "weight" in new_assembly.vertex_info[x_v].other_attr
                                            and
                                            db_name in new_assembly.vertex_info[x_v].other_attr["weight"]])
                        first_round = True
                        delete_those_vertices = set()
                        while first_round or delete_those_vertices:
                            first_round = False
                            delete_those_vertices = set()
                            for _v_n in new_assembly.vertex_info:
                                # both ends must have edge(s)
                                if sum([bool(len(cn))
                                        for cn in new_assembly.vertex_info[_v_n].connections.values()]) != 2:
                                    # To keep a terminal vertex
                                    # 1. tagged
                                    # 2. normal depth (3 sigma)
                                    # 3. enough weight
                                    if verbose and log_handler:
                                        log_handler.info(str(new_ave_cov) + "[" + str(ave_std) + "]~" +
                                                         str(new_assembly.vertex_info[_v_n].cov))
                                        if "weight" in new_assembly.vertex_info[_v_n].other_attr:
                                            log_handler.info("v_weight/total_weight: " +
                                                             str(new_assembly.vertex_info[_v_n].other_attr["weight"].
                                                                 get(db_name, 0.)) + "/" + str(total_weight))

                                    if _v_n in new_assembly.tagged_vertices[db_name]:
                                        if abs(new_ave_cov - new_assembly.vertex_info[_v_n].cov) <= 3 * ave_std and \
                                                "weight" in new_assembly.vertex_info[_v_n].other_attr and \
                                                new_assembly.vertex_info[_v_n].other_attr["weight"].get(db_name, 0.) * \
                                                weight_factor > total_weight:
                                            write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "f")
                                            raise ProcessingGraphFailed(
                                                "Incomplete/Complicated graph: please check around EDGE_" + _v_n + "!")
                                        else:
                                            if (verbose or debug) and log_handler:
                                                log_handler.warning(
                                                    "removing tagged but low-coverage terminal contig: " + _v_n + ":" +
                                                    str(new_assembly.vertex_info[_v_n].other_attr["tags"]))
                                            delete_those_vertices.add(_v_n)
                                    else:
                                        delete_those_vertices.add(_v_n)
                            if delete_those_vertices:
                                if verbose or debug:
                                    if log_handler:
                                        log_handler.info("removing terminal contigs: " + str(delete_those_vertices))
                                    else:
                                        sys.stdout.write(
                                            "removing terminal contigs: " + str(delete_those_vertices) + "\n")
                                new_assembly.remove_vertex(delete_those_vertices)
                                changed = True
                            if check_remaining_singleton():
                                break

                    # merge vertices_set
                    new_assembly.merge_all_possible_vertices()
                    if check_remaining_singleton():
                        break
                    new_ave_cov, ave_std = new_assembly.estimate_copy_and_depth_by_cov(
                        new_assembly.tagged_vertices[db_name], min_sigma=min_sigma_factor,
                        debug=debug, log_handler=log_handler,
                        verbose=verbose, mode=mode)
                    new_assembly.processing_polymorphism(database_name=db_name,
                                                         average_depth=new_ave_cov,
                                                         contamination_depth=contamination_depth,
                                                         contamination_similarity=contamination_similarity,
                                                         degenerate=False, degenerate_depth=degenerate_depth,
                                                         degenerate_similarity=degenerate_similarity,
                                                         verbose=verbose, debug=debug, log_handler=log_handler)
                    new_assembly.tag_in_between()
                    write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "g")

                write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "h")
                if check_remaining_singleton():
                    pass
                else:
                    new_assembly.processing_polymorphism(database_name=db_name,
                                                         contamination_depth=contamination_depth,
                                                         contamination_similarity=contamination_similarity,
                                                         degenerate=degenerate, degenerate_depth=degenerate_depth,
                                                         degenerate_similarity=degenerate_similarity,
                                                         warning_count=1, only_keep_max_cov=only_keep_max_cov,
                                                         verbose=verbose, debug=debug, log_handler=log_handler)
                    new_assembly.merge_all_possible_vertices()
                    write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "i")

                # create idealized vertices_set and edges
                try:
                    new_average_cov, ave_std = new_assembly.estimate_copy_and_depth_by_cov(
                        min_sigma=min_sigma_factor,
                        log_handler=log_handler,
                        verbose=verbose,
                        mode="all",
                        debug=debug)
                    if check_remaining_singleton():
                        final_res_combinations = [{"graph": new_assembly, "cov": new_average_cov}]
                    else:
                        final_res_combinations = new_assembly.estimate_copy_and_depth_precisely(
                            expected_average_cov=new_average_cov,
                            # broken_graph_allowed=broken_graph_allowed,
                            log_handler=log_handler,
                            verbose=verbose, debug=debug)
                        # maybe no more multiple results since 2022-12 gekko update
                        # if verbose:
                        #     if log_handler:
                        #         log_handler.info(str(len(final_res_combinations)) + " candidate graph(s) generated.")
                        #     else:
                        #         sys.stdout.write(str(len(final_res_combinations)) + " candidate graph(s) generated.\n")
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
                                if this_assembly_g.uni_overlap():
                                    for inside_v in this_assembly_g.vertex_info:
                                        draft_size_estimates += \
                                            (this_assembly_g.vertex_info[inside_v].len - this_assembly_g.uni_overlap()) \
                                            * this_assembly_g.vertex_to_copy[inside_v]
                                else:
                                    for inside_v in this_assembly_g.vertex_info:
                                        draft_size_estimates += \
                                            this_assembly_g.vertex_info[inside_v].len * \
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
                                        #     log_handler.info("Dropping vertices_set " + " ".join(dropping_names))
                                        # else:
                                        #     log_handler.info("Dropping vertices_set " + "".join(dropping_names) + "\n")
                                        new_possible_graph.remove_vertex(dropping_names)
                                        new_possible_graph.merge_all_possible_vertices()
                                        new_ave_cov, ave_std = new_possible_graph.estimate_copy_and_depth_by_cov(
                                            min_sigma=min_sigma_factor,
                                            log_handler=log_handler, verbose=verbose, mode="all", debug=debug)

                                        final_res_combinations.extend(
                                            new_possible_graph.estimate_copy_and_depth_precisely(
                                                expected_average_cov=new_ave_cov,
                                                # broken_graph_allowed=broken_graph_allowed,
                                                log_handler=log_handler, verbose=verbose, debug=debug))
                                        write_temp_out(new_possible_graph, db_name, temp_graph, temp_csv, "j")
                                    del final_res_combinations[go_graph]
                        if not final_res_combinations and absurd_copy_nums:
                            # if absurd_copy_nums:
                            #     raise ProcessingGraphFailed("Complicated graph! Detecting path(s) failed!")
                            # else:
                            raise ProcessingGraphFailed("Complicated graph! Detecting path(s) failed!")
                    if no_single_copy:
                        raise ProcessingGraphFailed("No single copy region?! Detecting path(s) failed!")
                except ImportError as e:
                    write_selected(_assembly=new_assembly, _selected_graph=selected_graph)
                    raise e
                # except (RecursionError, Exception) as e:
                # 2022-12-21 remove base class Exception
                except RecursionError as e:  # RecursionError is created by complex graph
                    unlabelled_contigs = [check_v for check_v in list(new_assembly.vertex_info)
                                          if check_v not in new_assembly.tagged_vertices[db_name]]
                    connections_removed = False
                    if unlabelled_contigs:
                        if verbose or debug:
                            if log_handler:
                                log_handler.info("removing unlabelled contigs: " + str(unlabelled_contigs))
                            else:
                                sys.stdout.write("removing unlabelled contigs: " + str(unlabelled_contigs) + "\n")
                        new_assembly.remove_vertex(unlabelled_contigs)
                        new_assembly.merge_all_possible_vertices()
                        write_selected(_assembly=new_assembly, _selected_graph=selected_graph)
                    else:
                        write_selected(_assembly=new_assembly, _selected_graph=selected_graph)
                        if broken_graph_allowed:
                            # delete all previous connections if all present contigs are labelled
                            for del_v_connection in new_assembly.vertex_info:
                                new_assembly.vertex_info[del_v_connection].connections = {True: OrderedDict(),
                                                                                          False: OrderedDict()}
                            # new_assembly.update_vertex_clusters()
                            connections_removed = True
                        else:
                            if verbose and log_handler:
                                log_handler.exception("")
                            raise e

                    new_average_cov, ave_std = new_assembly.estimate_copy_and_depth_by_cov(
                        re_initialize=True, min_sigma=min_sigma_factor,
                        log_handler=log_handler, verbose=verbose, mode="all", debug=debug)
                    outer_continue = False
                    for remove_all_connections in (False, True):
                        # if connections_removed and remove_all_connections:
                        #     is_reasonable_res = False
                        #     outer_continue = True
                        #     break
                        if remove_all_connections and not connections_removed:  # delete all previous connections
                            for del_v_connection in new_assembly.vertex_info:
                                new_assembly.vertex_info[del_v_connection].connections = {True: OrderedDict(),
                                                                                          False: OrderedDict()}
                        # new_assembly.merge_all_possible_vertices()
                        new_assembly.update_vertex_clusters()
                        try:
                            final_res_combinations = new_assembly.estimate_copy_and_depth_precisely(
                                expected_average_cov=new_average_cov,
                                log_handler=log_handler, verbose=verbose, debug=debug)
                        except ImportError as e:
                            raise e
                        except Exception as e:
                            if verbose or debug:
                                if remove_all_connections or connections_removed:
                                    log_handler.info("Unlikely error: " + str(e))
                                else:
                                    log_handler.info(str(e))
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
                                            "unacceptable, continue dropping suspicious vertices_set ...")
                                    else:
                                        sys.stdout.write(
                                            "Warning: Result with single copy vertex percentage < 50% is "
                                            "unacceptable, continue dropping suspicious vertices_set ...")
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
                                                        "unacceptable, continue dropping suspicious vertices_set ...")
                                else:
                                    sys.stdout.write("Warning: Result with single copy vertex percentage < 50% is "
                                                     "unacceptable, continue dropping suspicious vertices_set ...")
                            data_contains_outlier = True
                            is_reasonable_res = False
                            outer_continue = True
                            break
                    if outer_continue:
                        continue
                else:
                    write_selected(_assembly=new_assembly, _selected_graph=selected_graph)
                    test_first_g = final_res_combinations[0]["graph"]
                    if 1 in test_first_g.copy_to_vertex or min_single_copy_percent == 0:
                        single_copy_percent = sum([test_first_g.vertex_info[s_v].len
                                                   for s_v in test_first_g.copy_to_vertex[1]]) \
                                              / float(sum([test_first_g.vertex_info[a_v].len
                                                           for a_v in test_first_g.vertex_info]))
                        if single_copy_percent < min_single_copy_percent / 100.:
                            if verbose:
                                if log_handler:
                                    log_handler.warning("Result with single copy vertex percentage < {}% is "
                                                        "unacceptable, continue dropping suspicious vertices_set ..."
                                                        .format(min_single_copy_percent))
                                else:
                                    sys.stdout.write("Warning: Result with single copy vertex percentage < {}% is "
                                                     "unacceptable, continue dropping suspicious vertices_set ..."
                                                     .format(min_single_copy_percent))
                            data_contains_outlier = True
                            is_reasonable_res = False
                            continue
                        else:
                            log_target_res(final_res_combinations,
                                           log_handler=log_handler,
                                           read_len_for_log=read_len_for_log,
                                           kmer_for_log=kmer_for_log,
                                           universal_overlap=bool(self.uni_overlap()),
                                           mode=mode)
                            return final_res_combinations
                    else:
                        if verbose:
                            if log_handler:
                                log_handler.warning("Result with single copy vertex percentage < {}% is "
                                                    "unacceptable, continue dropping suspicious vertices_set ..."
                                                    .format(min_single_copy_percent))
                            else:
                                sys.stdout.write("Warning: Result with single copy vertex percentage < {}% is "
                                                 "unacceptable, continue dropping suspicious vertices_set ..."
                                                 .format(min_single_copy_percent))
                        data_contains_outlier = True
                        is_reasonable_res = False
                        continue
        except KeyboardInterrupt as e:
            write_temp_out(new_assembly, db_name, temp_graph, temp_csv, "k")
            if log_handler:
                log_handler.exception("")
                raise e
            else:
                raise e

    def peel_subgraph(self, subgraph, mode="", log_handler=None, verbose=False):
        assert isinstance(subgraph, Assembly)
        subgraph_vertices = set()
        for sub_v in subgraph.vertex_info:
            subgraph_vertices.update(subgraph.vertex_info[sub_v].merging_history.vertex_set())
        limited_vertices = set(self.vertex_info) & set(subgraph_vertices)
        if not limited_vertices:
            if log_handler:
                log_handler.warning("No overlapped vertices_set found for peeling!")
            else:
                sys.stdout.write("No overlapped vertices_set found for peeling!\n")
            if verbose:
                if log_handler:
                    log_handler.warning("graph vertices_set: " + str(sorted(self.vertex_info)))
                    log_handler.warning("subgraph vertices_set: " + str(sorted(subgraph.vertex_info)))
                else:
                    sys.stdout.write("graph vertices_set: " + str(sorted(self.vertex_info)))
                    sys.stdout.write("subgraph vertices_set: " + str(sorted(subgraph.vertex_info)))
        average_cov, ave_std = self.estimate_copy_and_depth_by_cov(
            limited_vertices, mode=mode, re_initialize=True, verbose=verbose)
        vertices_peeling_ratios = {}
        checked = set()
        for peel_name in sorted(limited_vertices):
            for peel_end, peel_connection_set in self.vertex_info[peel_name].connections.items():
                if (peel_name, not peel_end) in checked:
                    continue
                else:
                    checked.add((peel_name, not peel_end))
                for (external_v_n, external_v_e) in sorted(peel_connection_set):
                    if external_v_n in subgraph_vertices:
                        continue
                    if self.vertex_to_float_copy[peel_name] > self.vertex_to_copy[peel_name]:
                        # only peel the average part
                        vertices_peeling_ratios[peel_name] = \
                            1 - self.vertex_to_copy[peel_name] / self.vertex_to_float_copy[peel_name]
                        forward_peeling = [(next_n, not next_e)
                                           for next_n, next_e in self.vertex_info[peel_name].connections[not peel_end]
                                           if next_n in limited_vertices and (next_n, not next_e) not in checked]
                        while forward_peeling:
                            next_name, next_end = forward_peeling.pop(0)
                            if self.vertex_to_float_copy[next_name] > self.vertex_to_copy[next_name]:
                                vertices_peeling_ratios[next_name] = \
                                    1 - self.vertex_to_copy[next_name] / self.vertex_to_float_copy[next_name]
                                checked.add((next_name, next_end))
                                forward_peeling.extend(
                                    [(nx_nx_n, not nx_nx_e)
                                     for nx_nx_n, nx_nx_e in self.vertex_info[next_name].connections[next_end]
                                     if nx_nx_n in limited_vertices and (nx_nx_n, not nx_nx_e) not in checked])
        remove_vertices = {del_v for del_v in limited_vertices if del_v not in vertices_peeling_ratios}
        self.remove_vertex(remove_vertices)
        for peel_this_n in sorted(vertices_peeling_ratios):
            self.vertex_info[peel_this_n].cov *= vertices_peeling_ratios[peel_this_n]
            if "weight" in self.vertex_info[peel_this_n].other_attr and \
                    mode in self.vertex_info[peel_this_n].other_attr["weight"]:
                self.vertex_info[peel_this_n].other_attr["weight"][mode] *= vertices_peeling_ratios[peel_this_n]

    def add_gap_nodes_with_spades_res(self, scaffold_fasta, scaffold_paths, min_cov=0., max_cov=inf, log_handler=None,
                                      update_cluster=True, min_identifier_ws=12):
        if not self.__uni_overlap:
            return False
        else:
            ctg_olp = self.__uni_overlap
        spades_scaffolds = SpadesScaffolds(scaffold_fasta, scaffold_paths, assembly_obj=self, log_handler=log_handler,
                                           min_cov=min_cov, max_cov=max_cov, matching_ws=min_identifier_ws)
        gap_added = False
        rename_dict = {}
        for (l_name, l_end, l_trim), (r_name, r_end, r_trim), (gap_name, gap_seq, gap_len) \
                in spades_scaffolds.vertex_jumping_over_gap:
            if l_name in self.vertex_info and r_name in self.vertex_info:
                if (l_name, l_end) in self.vertex_info[r_name].connections[r_end]:
                    if log_handler:
                        log_handler.warning("Connection between " + l_name + ECHO_DIRECTION[l_end] + " and " +
                                            r_name + ECHO_DIRECTION[r_end] + " already existed!")
                    else:
                        sys.stdout.write("Warning: Connection between " + l_name + ECHO_DIRECTION[l_end] + " and " +
                                         r_name + ECHO_DIRECTION[r_end] + " already existed!\n")
                else:
                    if l_trim > 0:
                        self.vertex_info[l_name].seq[l_end] = self.vertex_info[l_name].seq[l_end][:-l_trim]
                        self.vertex_info[l_name].seq[not l_end] = self.vertex_info[l_name].seq[not l_end][l_trim:]
                        assert self.vertex_info[l_name].seq[l_end] == complementary_seq(
                            self.vertex_info[l_name].seq[not l_end])
                        self.vertex_info[l_name].len -= l_trim
                        # rename vertex: LT ~ left trimming, RT ~ right trimming
                        if l_name in rename_dict:
                            temp_n = rename_dict[l_name]
                            new_name = temp_n + "RT" + str(l_trim) if l_end else temp_n + "LT" + str(l_trim)
                        else:
                            new_name = l_name + "RT" + str(l_trim) if l_end else l_name + "LT" + str(l_trim)
                        rename_dict[l_name] = new_name
                    if r_trim > 0:
                        self.vertex_info[r_name].seq[r_end] = self.vertex_info[r_name].seq[r_end][:-r_trim]
                        self.vertex_info[r_name].seq[not r_end] = self.vertex_info[r_name].seq[not r_end][r_trim:]
                        self.vertex_info[r_name].len -= r_trim
                        # rename vertex: LT ~ left trimming, RT ~ right trimming
                        if r_name in rename_dict:
                            temp_n = rename_dict[r_name]
                            new_name = temp_n + "RT" + str(r_trim) if r_end else temp_n + "LT" + str(r_trim)
                        else:
                            new_name = r_name + "RT" + str(r_trim) if r_end else r_name + "LT" + str(r_trim)
                        assert self.vertex_info[r_name].seq[r_end] == complementary_seq(
                            self.vertex_info[r_name].seq[not r_end])
                        rename_dict[r_name] = new_name
                    if gap_len >= 0:
                        new_seq = self.vertex_info[l_name].seq[l_end][self.vertex_info[l_name].len - ctg_olp:] + \
                                  gap_seq + \
                                  self.vertex_info[r_name].seq[not r_end][:ctg_olp]
                    else:
                        new_seq = self.vertex_info[l_name].seq[l_end][self.vertex_info[l_name].len - ctg_olp:] + \
                                  self.vertex_info[r_name].seq[not r_end][-gap_len:ctg_olp]
                    new_average_cov = (self.vertex_info[l_name].cov * self.vertex_info[l_name].len +
                                       self.vertex_info[r_name].cov * self.vertex_info[r_name].len) / \
                                      (self.vertex_info[l_name].len + self.vertex_info[r_name].len)
                    self.vertex_info[gap_name] = Vertex(gap_name,
                                                        length=len(new_seq),
                                                        coverage=new_average_cov,
                                                        forward_seq=new_seq,
                                                        head_connections=OrderedDict([((l_name, l_end), ctg_olp)]),
                                                        tail_connections=OrderedDict([((r_name, r_end), ctg_olp)]))
                    self.vertex_info[l_name].connections[l_end][(gap_name, False)] = ctg_olp
                    self.vertex_info[r_name].connections[r_end][(gap_name, True)] = ctg_olp
                    gap_added = True
        for old_v_name, new_v_name in sorted(rename_dict.items()):
            self.rename_vertex(old_v_name, new_v_name, update_cluster=False)
        if gap_added and update_cluster:
            self.update_vertex_clusters()
        return gap_added

    def get_all_circular_paths(self, mode="embplant_pt",
                               library_info=None,
                               log_handler=None,
                               reverse_start_direction_for_pt=False,
                               max_paths_num=inf):
        # import time
        # count_time = [0.]
        # count_search = [0]
        def circular_directed_graph_solver(ongoing_path, next_connections, vertices_left, check_all_kinds,
                                           palindromic_repeat_vertices):
            # flush_str = "valid/searching: " + str(len(paths)) + "/" + str(count_search[0])
            # sys.stdout.write(flush_str + "\b" * len(flush_str))
            # sys.stdout.flush()
            # count_search[0] += 1
            if not vertices_left:
                new_path = deepcopy(ongoing_path)
                if palindromic_repeat_vertices:
                    new_path = [(this_v, True) if this_v in palindromic_repeat_vertices else (this_v, this_e)
                                for this_v, this_e in new_path]
                if check_all_kinds:
                    if palindromic_repeat_vertices:
                        rev_path = [(this_v, True) if this_v in palindromic_repeat_vertices else (this_v, not this_e)
                                    for this_v, this_e in new_path[::-1]]
                    else:
                        rev_path = [(this_v, not this_e) for this_v, this_e in new_path[::-1]]
                    this_path_derived = [new_path, rev_path]
                    for change_start in range(1, len(new_path)):
                        this_path_derived.append(new_path[change_start:] + new_path[:change_start])
                        this_path_derived.append(rev_path[change_start:] + rev_path[:change_start])
                    standardized_path = tuple(sorted(this_path_derived)[0])
                    if standardized_path not in paths_set:
                        paths_set.add(standardized_path)
                        paths.append(standardized_path)
                else:
                    new_path = tuple(new_path)
                    if new_path not in paths_set:
                        paths_set.add(new_path)
                        paths.append(new_path)
                return

            for next_vertex, next_end in next_connections:
                if len(paths) >= max_paths_num:
                    return
                # print("ongoing_path", ongoing_path)
                # print("next_vertex", next_vertex, next_connections)
                # print("vertices_left", vertices_left)
                # input()
                if next_vertex in vertices_left:
                    # to speed up
                    if vertices_left[next_vertex] == 1 and len(next_connections) >= 2:
                        # len(next_connections) >= 2 actually makes no big difference here according to a single test
                        # maybe add some threshold to do following calculation, e.g. left copies numbers ...
                        # time0 = time.time()
                        # costs very limited time
                        if not self.check_connected(set(vertices_left) - {next_vertex}):
                            # count_time[0] += time.time() - time0
                            continue
                    new_path = deepcopy(ongoing_path)
                    new_left = deepcopy(vertices_left)
                    new_path.append((next_vertex, not next_end))
                    new_left[next_vertex] -= 1
                    if not new_left[next_vertex]:
                        del new_left[next_vertex]
                    new_connections = self.vertex_info[next_vertex].connections[not next_end]
                    if not new_left:
                        if (start_vertex, not start_direction) in new_connections:
                            if palindromic_repeat_vertices:
                                new_path = [
                                    (this_v, True) if this_v in palindromic_repeat_vertices else (this_v, this_e)
                                    for this_v, this_e in new_path]
                            if check_all_kinds:
                                if palindromic_repeat_vertices:
                                    rev_path = [(this_v, True) if this_v in palindromic_repeat_vertices else
                                                (this_v, not this_e)
                                                for this_v, this_e in new_path[::-1]]
                                else:
                                    rev_path = [(this_v, not this_e) for this_v, this_e in new_path[::-1]]
                                this_path_derived = [new_path, rev_path]
                                for change_start in range(1, len(new_path)):
                                    this_path_derived.append(new_path[change_start:] + new_path[:change_start])
                                    this_path_derived.append(rev_path[change_start:] + rev_path[:change_start])
                                standardized_path = tuple(sorted(this_path_derived)[0])
                                if standardized_path not in paths_set:
                                    paths_set.add(standardized_path)
                                    paths.append(standardized_path)
                            else:
                                new_path = tuple(new_path)
                                if new_path not in paths_set:
                                    paths_set.add(new_path)
                                    paths.append(new_path)
                            return
                        else:
                            return
                    else:
                        new_connect_list = sorted(new_connections)
                        # if next_connections is SSC, reorder
                        # if mode == "embplant_pt" and len(new_connect_list) == 2 and new_connect_list[0][0] == \
                        #         new_connect_list[1][0]:
                        #     new_connect_list.sort(
                        #         key=lambda x: -self.vertex_info[x[0]].other_attr["orf"][x[1]]["sum_len"])
                        circular_directed_graph_solver(new_path, new_connect_list, new_left, check_all_kinds,
                                                       palindromic_repeat_vertices)
            return

        # for palindromic repeats
        palindromic_repeats = set()
        log_palindrome = False
        for vertex_n in self.vertex_info:
            if self.vertex_info[vertex_n].seq[True] == self.vertex_info[vertex_n].seq[False]:
                forward_c = deepcopy(self.vertex_info[vertex_n].connections[True])
                reverse_c = deepcopy(self.vertex_info[vertex_n].connections[False])
                # This is heuristic
                # In the rarely-used expression way, a contig connect itself in one end:
                # (vertex_n, True) in forward_c or (vertex_n, False) in reverse_c
                if forward_c and \
                        ((forward_c == reverse_c) or
                         ((vertex_n, True) in forward_c) or
                         ((vertex_n, False) in reverse_c)):
                    log_palindrome = True
                    if len(forward_c) == len(reverse_c) == 2:  # simple palindromic repeats, prune repeated connections
                        for go_d, (nb_vertex, nb_direction) in enumerate(tuple(forward_c)):
                            del self.vertex_info[nb_vertex].connections[nb_direction][(vertex_n, bool(go_d))]
                            del self.vertex_info[vertex_n].connections[bool(go_d)][(nb_vertex, nb_direction)]
                    elif len(forward_c) == len(reverse_c) == 1:  # connect to the same inverted repeat
                        pass
                    else:  # complicated, recorded
                        palindromic_repeats.add(vertex_n)
        if log_palindrome:
            log_handler.info("Palindromic repeats detected. "
                             "Different paths generating identical sequence will be merged.")

        #
        self.update_orf_total_len()
        paths = []
        paths_set = set()
        if 1 not in self.copy_to_vertex:
            do_check_all_start_kinds = True
            start_vertex = sorted(self.vertex_info,
                                  key=lambda x: (-self.vertex_info[x].len,
                                                 -max(self.vertex_info[x].other_attr["orf"][True]["sum_len"],
                                                      self.vertex_info[x].other_attr["orf"][False]["sum_len"]),
                                                 x))[0]
            start_direction = True
        else:
            # start from a single copy vertex, no need to check all kinds of start vertex
            do_check_all_start_kinds = False
            start_vertex = sorted(self.copy_to_vertex[1])[0]
            start_direction = True

        # each contig stored format:
        first_path = [(start_vertex, start_direction)]
        first_connections = sorted(self.vertex_info[start_vertex].connections[start_direction])
        vertex_to_copy = deepcopy(self.vertex_to_copy)
        vertex_to_copy[start_vertex] -= 1
        if vertex_to_copy[start_vertex] <= 0:
            del vertex_to_copy[start_vertex]
        circular_directed_graph_solver(first_path, first_connections, vertex_to_copy, do_check_all_start_kinds,
                                       palindromic_repeats)
        # log_handler.info("check_connected costs: " + str(count_time[0]))

        if not paths:
            raise ProcessingGraphFailed("Detecting path(s) from remaining graph failed!")
        else:
            # modify start_vertex based on the whole path, if starting from a single copy vertex
            def reseed_a_path(input_path, input_unique_vertex):
                if input_unique_vertex not in input_path:
                    new_path = [(element_v, not element_e) for (element_v, element_e) in input_path[::-1]]
                else:
                    new_path = input_path
                reseed_from = new_path.index(input_unique_vertex)
                return new_path[reseed_from:] + new_path[:reseed_from]

            if 1 in self.copy_to_vertex:
                branching_sc_pair_set = set()
                branching_sc_pair_list = list()
                if mode == "embplant_pt" and 2 in self.copy_to_vertex:
                    # find branching points
                    for candidate_name in self.copy_to_vertex[2]:
                        if not bool(self.is_sequential_repeat(candidate_name)):
                            for neighboring_vertices in self.vertex_info[candidate_name].connections.values():
                                if len(neighboring_vertices) == 2:
                                    (left_v, left_e), (right_v, right_e) = sorted(neighboring_vertices)
                                    if left_v in self.copy_to_vertex[1] and right_v in self.copy_to_vertex[1]:
                                        if ((left_v, not left_e), (right_v, right_e)) not in branching_sc_pair_set:
                                            # from v1_e1 to v2_e2;
                                            # v1_e1 and v2_e2 must have the same direction in the path
                                            branching_sc_pair_set.add(((left_v, not left_e), (right_v, right_e)))
                                            branching_sc_pair_set.add(((right_v, not right_e), (left_v, left_e)))
                                            branching_sc_pair_list.append(((left_v, not left_e), (right_v, right_e)))
                                            branching_sc_pair_list.append(((right_v, not right_e), (left_v, left_e)))
                if branching_sc_pair_set:
                    # more orfs found in the reverse direction of LSC of a typical plastome
                    # different paths may have different LSC
                    # picking the sub-path with the longest length with strand of least orfs as the new start point
                    for go_p, each_path in enumerate(paths):
                        reverse_path = [(element_v, not element_e) for (element_v, element_e) in each_path[::-1]]
                        sub_paths_for_checking = []
                        for (left_v, left_e), (right_v, right_e) in branching_sc_pair_list:
                            if (left_v, left_e) in each_path:
                                if (right_v, right_e) in each_path:
                                    left_id = each_path.index((left_v, left_e))
                                    right_id = each_path.index((right_v, right_e))
                                    if left_id <= right_id:
                                        sub_paths_for_checking.append(each_path[left_id: right_id + 1])
                                    else:
                                        sub_paths_for_checking.append(each_path[left_id:] + each_path[:right_id + 1])
                                else:
                                    sub_paths_for_checking.append([])
                            else:
                                if (right_v, right_e) in reverse_path:
                                    left_id = reverse_path.index((left_v, left_e))
                                    right_id = reverse_path.index((right_v, right_e))
                                    if left_id <= right_id:
                                        sub_paths_for_checking.append(reverse_path[left_id: right_id + 1])
                                    else:
                                        sub_paths_for_checking.append(
                                            reverse_path[left_id:] + reverse_path[:right_id + 1])
                                else:
                                    sub_paths_for_checking.append([])
                        # picking the vertex with the longest length with strand of least orfs
                        lsc_pair_id = sorted(range(len(sub_paths_for_checking)),
                                             key=lambda x:
                                             (
                                                 # sum of the contig lengths
                                                 -sum([self.vertex_info[sub_v].len
                                                       for sub_v, sub_e in sub_paths_for_checking[x]]) +
                                                 # sum of the overlap
                                                 # self.__overlap * (len(sub_paths_for_checking[x]) - 1),
                                                 sum([self.vertex_info[l_v].connections[l_e]
                                                      [(sub_paths_for_checking[x][go_s + 1][0],
                                                        not sub_paths_for_checking[x][go_s + 1][1])]
                                                      for go_s, (l_v, l_e) in
                                                      enumerate(sub_paths_for_checking[x][:-1])]),
                                                 # sum of the orf
                                                 sum([self.vertex_info[sub_v].other_attr["orf"][sub_e]["sum_len"]
                                                      for sub_v, sub_e in sub_paths_for_checking[x]]),
                                                 x))[0]
                        if reverse_start_direction_for_pt:
                            if lsc_pair_id % 2 == 0:
                                paths[go_p] = reseed_a_path(each_path, branching_sc_pair_list[lsc_pair_id + 1][0])
                            else:
                                paths[go_p] = reseed_a_path(each_path, branching_sc_pair_list[lsc_pair_id - 1][0])
                        else:
                            paths[go_p] = reseed_a_path(each_path, branching_sc_pair_list[lsc_pair_id][0])
                else:
                    candidate_single_copy_vertices = set()
                    for single_v in self.copy_to_vertex[1]:
                        candidate_single_copy_vertices.add((single_v, True))
                        candidate_single_copy_vertices.add((single_v, False))
                    if mode == "embplant_pt":
                        # more orfs found in the reverse direction of LSC of a typical plastome
                        # picking the vertex with the longest length with strand of least orfs
                        start_vertex, start_direction = sorted(candidate_single_copy_vertices,
                                                               key=lambda x: (-self.vertex_info[x[0]].len,
                                                                              self.vertex_info[x[0]].other_attr["orf"][
                                                                                  x[1]]["sum_len"],
                                                                              x))[0]
                        if reverse_start_direction_for_pt:
                            start_direction = not start_direction
                    else:
                        # picking the vertex with the longest length with strand of most orfs
                        start_vertex, start_direction = sorted(candidate_single_copy_vertices,
                                                               key=lambda x: (-self.vertex_info[x[0]].len,
                                                                              -self.vertex_info[x[0]].other_attr["orf"][
                                                                                  x[1]]["sum_len"],
                                                                              x))[0]
                    for go_p, each_path in enumerate(paths):
                        paths[go_p] = reseed_a_path(each_path, (start_vertex, start_direction))

            # sorting path by average distance among multi-copy loci
            # the highest would be more symmetrical IR, which turns out to be more reasonable
            sorted_paths = []
            total_len = len(list(paths)[0])
            record_pattern = False
            for original_id, this_path in enumerate(paths):
                acc_dist = 0
                for copy_num in self.copy_to_vertex:
                    if copy_num > 2:
                        record_pattern = True
                        for vertex_name in self.copy_to_vertex[copy_num]:
                            loc_ids = [go_to_id for go_to_id, (v, e) in enumerate(this_path) if v == vertex_name]
                            for id_a, id_b in combinations(loc_ids, 2):
                                acc_dist += min((id_a - id_b) % total_len, (id_b - id_a) % total_len)
                sorted_paths.append((this_path, acc_dist, original_id))
            if record_pattern:
                sorted_paths.sort(key=lambda x: (-x[1], x[2]))
                pattern_dict = {acc_distance: ad_id + 1
                                for ad_id, acc_distance in enumerate(sorted(set([x[1] for x in sorted_paths]),
                                                                            reverse=True))}
                if len(pattern_dict) > 1:
                    if mode == "embplant_pt":
                        if log_handler:
                            log_handler.warning("Multiple repeat patterns appeared in your data, "
                                                "a more balanced pattern (always the repeat_pattern1) "
                                                "would be suggested for plastomes with the canonical IR!")
                        else:
                            sys.stdout.write("Warning: Multiple repeat patterns appeared in your data, "
                                             "a more balanced pattern (always the repeat_pattern1) would be suggested "
                                             "for plastomes with the canonical IR!\n")
                    sorted_paths = [(this_path, ".repeat_pattern" + str(pattern_dict[acc_distance]))
                                    for this_path, acc_distance, foo_id in sorted_paths]
                else:
                    sorted_paths = [(this_path, "") for this_path in paths]
            else:
                sorted_paths = [(this_path, "") for this_path in paths]

            if mode == "embplant_pt":
                if len(sorted_paths) > 2 and not (100000 < len(self.export_path(sorted_paths[0][0]).seq) < 200000):
                    if log_handler:
                        log_handler.warning("Multiple circular genome structures with abnormal length produced!")
                        log_handler.warning("Please check the assembly graph and selected graph to confirm.")
                    else:
                        sys.stdout.write(
                            "Warning: Multiple circular genome structures with abnormal length produced!\n")
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
                                            "simply different in SSC direction (two flip-flop configurations)!")
                    else:
                        sys.stdout.write("More than one circular genome structure produced ...\n")
                        sys.stdout.write("Please check the final result to confirm whether they are "
                                         "simply different in SSC direction (two flip-flop configurations)!\n")
            return sorted_paths

    def standardize_paths(self, raw_paths, undirected_vertices={}, only_res=True):
        if undirected_vertices:
            corrected_paths = [[(this_v, True) if this_v in undirected_vertices else (this_v, this_e)
                                for this_v, this_e in path_part]
                               for path_part in raw_paths]
        else:
            corrected_paths = deepcopy(raw_paths)
        here_standardized_path = []
        for part_path in corrected_paths:
            if undirected_vertices:
                rev_part = [(this_v, True) if this_v in undirected_vertices else (this_v, not this_e)
                            for this_v, this_e in part_path[::-1]]
            else:
                rev_part = [(this_v, not this_e) for this_v, this_e in part_path[::-1]]
            if (part_path[0][0], not part_path[0][1]) \
                    in self.vertex_info[part_path[-1][0]].connections[part_path[-1][1]]:
                # circular
                this_part_derived = [part_path, rev_part]
                for change_start in range(1, len(part_path)):
                    this_part_derived.append(part_path[change_start:] + part_path[:change_start])
                    this_part_derived.append(rev_part[change_start:] + rev_part[:change_start])
                # try:
                standard_part = tuple(sorted(this_part_derived, key=lambda x: smart_trans_for_sort(x))[0])
                # except TypeError:
                #     for j in this_part_derived:
                #         print(j)
                #     exit()
            else:
                standard_part = tuple(sorted([part_path, rev_part], key=lambda x: smart_trans_for_sort(x))[0])
            here_standardized_path.append(standard_part)
        if only_res:
            return tuple(sorted(here_standardized_path, key=lambda x: smart_trans_for_sort(x)))
        else:
            return corrected_paths, tuple(sorted(here_standardized_path, key=lambda x: smart_trans_for_sort(x)))

    def get_all_paths(self, mode="embplant_pt", max_paths_num=inf, log_handler=None):

        # def standardize_paths(raw_paths, undirected_vertices):
        #     if undirected_vertices:
        #         corrected_paths = [[(this_v, True) if this_v in undirected_vertices else (this_v, this_e)
        #                             for this_v, this_e in path_part]
        #                            for path_part in raw_paths]
        #     else:
        #         corrected_paths = deepcopy(raw_paths)
        #     here_standardized_path = []
        #     for part_path in corrected_paths:
        #         if undirected_vertices:
        #             rev_part = [(this_v, True) if this_v in undirected_vertices else (this_v, not this_e)
        #                         for this_v, this_e in part_path[::-1]]
        #         else:
        #             rev_part = [(this_v, not this_e) for this_v, this_e in part_path[::-1]]
        #         if (part_path[0][0], not part_path[0][1]) \
        #                 in self.vertex_info[part_path[-1][0]].connections[part_path[-1][1]]:
        #             # circular
        #             this_part_derived = [part_path, rev_part]
        #             for change_start in range(1, len(part_path)):
        #                 this_part_derived.append(part_path[change_start:] + part_path[:change_start])
        #                 this_part_derived.append(rev_part[change_start:] + rev_part[:change_start])
        #             # try:
        #             standard_part = tuple(sorted(this_part_derived, key=lambda x: smart_trans_for_sort(x))[0])
        #             # except TypeError:
        #             #     for j in this_part_derived:
        #             #         print(j)
        #             #     exit()
        #         else:
        #             standard_part = tuple(sorted([part_path, rev_part], key=lambda x: smart_trans_for_sort(x))[0])
        #         here_standardized_path.append(standard_part)
        #     return corrected_paths, tuple(sorted(here_standardized_path, key=lambda x: smart_trans_for_sort(x)))

        def directed_graph_solver(ongoing_paths, next_connections, vertices_left, in_all_start_ve, undirected_vertices):
            # print("-----------------------------")
            # print("ongoing_path", ongoing_path)
            # print("next_connect", next_connections)
            # print("vertices_lef", vertices_left)
            # print("vertices_lef", len(vertices_left))
            if not vertices_left:
                new_paths, new_standardized = self.standardize_paths(ongoing_paths, undirected_vertices, False)
                if new_standardized not in paths_set:
                    paths.append(new_paths)
                    paths_set.add(new_standardized)
                return

            # to speed up under max_paths_num, rank choices by complete, then incomplete
            if len(next_connections) >= 2:
                go_n = 0
                incomplete_choices = []
                while go_n < len(next_connections):
                    next_vertex, next_end = next_connections[go_n]
                    if next_vertex in vertices_left:
                        # maybe add some threshold to do following calculation to speed up, e.g. left copies numbers ...
                        if vertices_left[next_vertex] > 1 or \
                                not self.check_connected(set(vertices_left) - {next_vertex}):
                            incomplete_choices.append(next_connections.pop(go_n))
                        else:
                            go_n += 1
                    else:
                        del next_connections[go_n]
                next_connections.extend(incomplete_choices)

            find_next = False
            for next_vertex, next_end in next_connections:
                # print("next_vertex, next_end: {}, {}".format(next_vertex, next_end))
                # print("vertices_left: {}".format(vertices_left))
                # input()
                if len(paths) >= max_paths_num:
                    return
                # print("next_vertex", next_vertex, next_end)
                if next_vertex in vertices_left:
                    find_next = True
                    new_paths = deepcopy(ongoing_paths)
                    new_left = deepcopy(vertices_left)
                    new_paths[-1].append((next_vertex, not next_end))
                    new_left[next_vertex] -= 1
                    if not new_left[next_vertex]:
                        del new_left[next_vertex]
                    new_connect_list = sorted(self.vertex_info[next_vertex].connections[not next_end])
                    if not new_left:
                        new_paths, new_standardized = self.standardize_paths(new_paths, undirected_vertices, False)
                        if new_standardized not in paths_set:
                            paths.append(new_paths)
                            paths_set.add(new_standardized)
                        return
                    else:
                        # if mode == "embplant_pt" and len(new_connect_list) == 2 and new_connect_list[0][0] == \
                        #         new_connect_list[1][0]:
                        #     new_connect_list.sort(
                        #         key=lambda x: self.vertex_info[x[0]].other_attr["orf"][x[1]]["sum_len"])
                        # to_print = "len(paths)={}".format(len(paths))
                        # sys.stdout.write(to_print + "\b" * len(to_print))
                        # sys.stdout.flush()
                        directed_graph_solver(new_paths, new_connect_list, new_left, in_all_start_ve,
                                              undirected_vertices)
            if not find_next:
                new_all_start_ve = deepcopy(in_all_start_ve)
                while new_all_start_ve:
                    if len(paths) >= max_paths_num:
                        return
                    new_start_vertex, new_start_end = new_all_start_ve.pop(0)
                    if new_start_vertex in vertices_left:
                        new_paths = deepcopy(ongoing_paths)
                        new_left = deepcopy(vertices_left)
                        new_paths.append([(new_start_vertex, new_start_end)])
                        new_left[new_start_vertex] -= 1
                        if not new_left[new_start_vertex]:
                            del new_left[new_start_vertex]
                        new_connect_list = sorted(self.vertex_info[new_start_vertex].connections[new_start_end])
                        if not new_left:
                            new_paths, new_standardized = self.standardize_paths(new_paths, undirected_vertices, False)
                            if new_standardized not in paths_set:
                                paths.append(new_paths)
                                paths_set.add(new_standardized)
                        else:
                            # if mode == "embplant_pt" and len(new_connect_list) == 2 and new_connect_list[0][0] == \
                            #         new_connect_list[1][0]:
                            #     new_connect_list.sort(
                            #         key=lambda x: self.vertex_info[x[0]].other_attr["orf"][x[1]]["sum_len"])
                            # to_print = "len(paths)={}".format(len(paths))
                            # sys.stdout.write(to_print + "\b" * len(to_print))
                            # sys.stdout.flush()
                            # TODO:prioritize those makes the component more connected

                            directed_graph_solver(new_paths, new_connect_list, new_left, new_all_start_ve,
                                                  undirected_vertices)
                            break
                if not new_all_start_ve:
                    return

        paths = list()
        paths_set = set()
        # start from a terminal vertex in an open graph/subgraph
        #         or a single copy vertex in a closed graph/subgraph
        self.update_orf_total_len()

        # 2019-12-28 palindromic repeats
        palindromic_repeats = set()
        log_palindrome = False
        for vertex_n in self.vertex_info:
            if self.vertex_info[vertex_n].seq[True] == self.vertex_info[vertex_n].seq[False]:
                temp_f = self.vertex_info[vertex_n].connections[True]
                temp_r = self.vertex_info[vertex_n].connections[False]
                if temp_f and temp_f == temp_r:
                    log_palindrome = True
                    if len(temp_f) == len(temp_r) == 2:  # simple palindromic repeats, prune repeated connections
                        for go_d, (nb_vertex, nb_direction) in enumerate(tuple(temp_f)):
                            del self.vertex_info[nb_vertex].connections[nb_direction][(vertex_n, bool(go_d))]
                            del self.vertex_info[vertex_n].connections[bool(go_d)][(nb_vertex, nb_direction)]
                    elif len(temp_f) == len(temp_r) == 1:  # connect to the same inverted repeat
                        pass
                    else:  # complicated, recorded
                        palindromic_repeats.add(vertex_n)
        if log_palindrome:
            log_handler.info("Palindromic repeats detected. "
                             "Different paths generating identical sequence will be merged.")

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
        first_connections = sorted(self.vertex_info[start_v_e[0]].connections[start_v_e[1]])
        vertex_to_copy = deepcopy(self.vertex_to_copy)
        vertex_to_copy[start_v_e[0]] -= 1
        if not vertex_to_copy[start_v_e[0]]:
            del vertex_to_copy[start_v_e[0]]
        directed_graph_solver(first_path, first_connections, vertex_to_copy, all_start_v_e,
                              undirected_vertices=palindromic_repeats)

        # standardized_path_unique_set = set([this_path_pair[1] for this_path_pair in path_paris])
        # paths = []
        # for raw_path, standardized_path in path_paris:
        #     if standardized_path in standardized_path_unique_set:
        #         paths.append(raw_path)
        #         standardized_path_unique_set.remove(standardized_path)

        if not paths:
            raise ProcessingGraphFailed("Detecting path(s) from remaining graph failed!")
        else:
            sorted_paths = []
            # total_len = len(list(set(paths))[0])
            record_pattern = False
            for original_id, this_path in enumerate(paths):
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
                sorted_paths.append((this_path, acc_dist, original_id))
            if record_pattern:
                sorted_paths.sort(key=lambda x: (-x[1], x[2]))
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
                                    for this_path, acc_distance, foo_id in sorted_paths]
                else:
                    sorted_paths = [(this_path, "") for this_path in sorted(paths)]
            else:
                sorted_paths = [(this_path, "") for this_path in sorted(paths)]

            if mode == "embplant_pt":
                if len(sorted_paths) > 2 and \
                        not (100000 < sum(
                            [len(self.export_path(part_p).seq) for part_p in sorted_paths[0][0]]) < 200000):
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
                                            "simply different in SSC direction (two flip-flop configurations)!")
                    else:
                        sys.stdout.write("More than one structure (gene order) produced ...\n")
                        sys.stdout.write("Please check the final result to confirm whether they are "
                                         "simply different in SSC direction (two flip-flop configurations)!\n")
            return sorted_paths

    def export_path(self, in_path):
        # overlap = self.__overlap if self.__overlap else 0
        seq_names = []
        seq_segments = []
        len_path = len(in_path)
        for go_p, (this_vertex, this_end) in enumerate(in_path):
            last_v, last_e = in_path[(go_p - 1) % len_path]
            this_overlap = self.vertex_info[last_v].connections[last_e].get((this_vertex, not this_end), 0)
            seq_segments.append(self.vertex_info[this_vertex].seq[this_end][this_overlap:])
            seq_names.append(this_vertex + ("-", "+")[this_end])
        if (in_path[0][0], not in_path[0][1]) in self.vertex_info[in_path[-1][0]].connections[in_path[-1][1]]:
            seq_names[-1] += "(circular)"
        # else:
        #     seq_segments[0] = self.vertex_info[in_path[0][0]].seq[in_path[0][1]][:overlap] + seq_segments[0]
        return Sequence(",".join(seq_names), "".join(seq_segments))


class NaiveKmerNodeGraph(Assembly):
    def __init__(self, fasta_file, kmer_len=55, circular="auto", circular_head_ends="(circular)", single_chain=False):
        """
        :param fasta_file:
        :param kmer_len:
        :param circular: "auto" (default), "yes", "no"
        :param circular_head_ends:
        :return:
        """
        super(NaiveKmerNodeGraph, self).__init__(uni_overlap=kmer_len - 1)
        assert circular in ("auto", "yes", "no")
        assert kmer_len >= 3 and kmer_len % 2 == 1
        self.__kmer = kmer_len  # overlap is actually kmer_len - 1
        recorded_kmers = {}
        count_vertices = 0
        seq_matrix = SequenceList(fasta_file)
        for seq_record in seq_matrix.sequences:
            if len(seq_record.seq) >= kmer_len:
                go_circle = circular == "yes" or (circular == "auto" and seq_record.label.endswith(circular_head_ends))
                kmer_list = SeqKmerIndexer(seq_record.seq, kmer_len, is_circular=go_circle)
                # 1. initial kmer seq
                this_kmer_seq = kmer_list[0]
                if this_kmer_seq in recorded_kmers:
                    this_vertex, this_end = recorded_kmers[this_kmer_seq]
                    self.vertex_info[this_vertex].cov += 1.
                else:
                    # add kmer as a new vertex
                    count_vertices += 1
                    recorded_kmers[this_kmer_seq] = this_vertex, this_end = str(count_vertices), True
                    self.vertex_info[this_vertex] = this_v_info = Vertex(this_vertex, kmer_len, 1., this_kmer_seq)
                    # record the connection as dict() rather than set() for counting
                    self.vertex_info[this_vertex].connections = {True: OrderedDict(), False: OrderedDict()}
                    if not single_chain:
                        recorded_kmers[this_v_info.seq[False]] = this_vertex, not this_end
                if go_circle:
                    # add connection between the first kmer and the last kmer if the seq is circular
                    prev_kmer_seq = kmer_list[- 1]
                    if prev_kmer_seq in recorded_kmers:
                        prev_vertex, prev_end = recorded_kmers[prev_kmer_seq]
                    else:
                        # if the last kmer is not recorded, recorded it with coverage of zero
                        count_vertices += 1
                        recorded_kmers[prev_kmer_seq] = prev_vertex, prev_end = str(count_vertices), True
                        self.vertex_info[prev_vertex] = prev_v_info = Vertex(prev_vertex, kmer_len, 0., prev_kmer_seq)
                        self.vertex_info[prev_vertex].connections = {True: OrderedDict(), False: OrderedDict()}
                        if not single_chain:
                            recorded_kmers[prev_v_info.seq[False]] = prev_vertex, not prev_end
                    # the value here is not the overlap, but the connection strength
                    if (this_vertex, not this_end) not in self.vertex_info[prev_vertex].connections[prev_end]:
                        self.vertex_info[prev_vertex].connections[prev_end][(this_vertex, not this_end)] = 0
                    self.vertex_info[prev_vertex].connections[prev_end][(this_vertex, not this_end)] += 1
                    if (prev_vertex, prev_end) not in self.vertex_info[this_vertex].connections[not this_end]:
                        self.vertex_info[this_vertex].connections[not this_end][(prev_vertex, prev_end)] = 0
                    self.vertex_info[this_vertex].connections[not this_end][(prev_vertex, prev_end)] += 1
                # 2. the remaining kmers
                for go_to in range(1, len(kmer_list)):
                    this_kmer_seq = kmer_list[go_to]
                    if this_kmer_seq in recorded_kmers:
                        this_vertex, this_end = recorded_kmers[this_kmer_seq]
                        self.vertex_info[this_vertex].cov += 1.
                    else:
                        # add kmer as a new vertex
                        count_vertices += 1
                        recorded_kmers[this_kmer_seq] = this_vertex, this_end = str(count_vertices), True
                        self.vertex_info[this_vertex] = this_v_info = Vertex(this_vertex, kmer_len, 1., this_kmer_seq)
                        self.vertex_info[this_vertex].connections = {True: OrderedDict(), False: OrderedDict()}
                        if not single_chain:
                            recorded_kmers[this_v_info.seq[False]] = this_vertex, not this_end
                    # add the connection between this_kmer_seq and prev_kmer_seq
                    prev_kmer_seq = kmer_list[go_to - 1]
                    prev_vertex, prev_end = recorded_kmers[prev_kmer_seq]
                    if (this_vertex, not this_end) not in self.vertex_info[prev_vertex].connections[prev_end]:
                        self.vertex_info[prev_vertex].connections[prev_end][(this_vertex, not this_end)] = 0
                    self.vertex_info[prev_vertex].connections[prev_end][(this_vertex, not this_end)] += 1
                    if (prev_vertex, prev_end) not in self.vertex_info[this_vertex].connections[not this_end]:
                        self.vertex_info[this_vertex].connections[not this_end][(prev_vertex, prev_end)] = 0
                    self.vertex_info[this_vertex].connections[not this_end][(prev_vertex, prev_end)] += 1
            else:
                sys.stderr.write("Warning: " + fasta_file + ":" + seq_record.label +
                                 ":length:" + str(len(seq_record.seq)) + " < kmer:" + str(kmer_len) + " .. skipped!\n")

    def generate_assembly_graph(self):
        def generate_contig(
                this_k_n, this_k_e, next_k_n, next_k_e, terminal_k_quota, in_graph,
                in_kmer_e_to_ctg_e, inter_ks, next_kmers):
            kmer_ids_used_this_round = {this_k_n}
            if this_k_n in terminal_k_quota:
                terminal_k_quota[this_k_n][this_k_e] -= 1
            count_k_len = 1
            count_k_sum = self.vertex_info[this_k_n].cov / len(self.vertex_info[this_k_n].connections[this_k_e])
            forward_seq = self.vertex_info[this_k_n].seq[this_k_e]
            while True:
                if next_k_n in terminal_k_quota:
                    terminal_k_quota[next_k_n][next_k_e] -= 1
                count_k_len += 1
                # avoid repeatedly counting coverage for used kmers
                if next_k_n not in kmer_ids_used_this_round:
                    count_k_sum += self.vertex_info[next_k_n].cov
                forward_seq += self.vertex_info[next_k_n].seq[not next_k_e][-1]
                kmer_ids_used_this_round.add(next_k_n)
                if next_k_n in terminal_k_quota:  # terminal
                    next_kmers.append(next_k_n)
                    break
                elif next_k_n == this_k_n and (not next_k_e) == this_k_e:  # circular
                    break
                else:
                    # update next_k_name, next_k_end
                    next_k_n, next_k_e = \
                        sorted(self.vertex_info[next_k_n].connections[not next_k_e])[0]
            c_name = str(len(in_graph.vertex_info) + 1)
            c_end = False  # start of a contig
            # record connections to be added later
            if this_k_n not in in_kmer_e_to_ctg_e:
                in_kmer_e_to_ctg_e[this_k_n] = {True: set(), False: set()}
            if next_k_n not in in_kmer_e_to_ctg_e:
                in_kmer_e_to_ctg_e[next_k_n] = {True: set(), False: set()}
            in_kmer_e_to_ctg_e[this_k_n][this_k_e].add((c_name, c_end))
            in_kmer_e_to_ctg_e[next_k_n][next_k_e].add((c_name, not c_end))
            # record contig
            c_length = count_k_len + self.__kmer - 1
            c_coverage = count_k_sum / float(count_k_len)
            in_graph.vertex_info[c_name] = this_contig = Vertex(
                c_name, length=c_length, coverage=c_coverage, forward_seq=forward_seq,
                fastg_form_long_name="EDGE_" + c_name + "_length_" + str(c_length) + "_cov_" + str(c_coverage))
            for added_kmer_name in kmer_ids_used_this_round:
                inter_ks.discard(added_kmer_name)

        # 1. profile the graph shape
        # find terminal kmers and how many times a terminal kmer will be used
        # find intermediate kmers
        # find singletons
        point_kmer_link_quota = {}
        intermediate_kmers = set()
        singletons = []
        for this_k_name in sorted(self.vertex_info):
            connect_n_link_tail = len(self.vertex_info[this_k_name].connections[True])
            connect_n_link_head = len(self.vertex_info[this_k_name].connections[False])
            is_terminal = False
            if connect_n_link_tail == 0:
                if connect_n_link_head == 0:
                    singletons.append(this_k_name)
                else:
                    is_terminal = True
            elif connect_n_link_tail == 1:
                if connect_n_link_head == 0:
                    is_terminal = True
                elif connect_n_link_head == 1:
                    intermediate_kmers.add(this_k_name)
                else:
                    is_terminal = True
            else:
                is_terminal = True
            if is_terminal:
                point_kmer_link_quota[this_k_name] = {False: connect_n_link_head, True: connect_n_link_tail}

        # 2. generating
        # 2.1 initialization
        assembly_graph = Assembly(uni_overlap=self.__kmer)  # overlap=self.__kmer
        joint_kmer_end_to_contig_end = {}  # record contig connections using joint_kmer_end_to_contig_end
        # 2.2 add singleton
        for this_k_name in singletons:
            new_contig_name = str(len(assembly_graph.vertex_info) + 1)
            assembly_graph.vertex_info[new_contig_name] = Vertex(
                new_contig_name, length=self.__kmer, coverage=self.vertex_info[this_k_name].cov,
                forward_seq=self.vertex_info[this_k_name].seq[True],
                fastg_form_long_name="EDGE_" + new_contig_name + "_length_" + str(self.__kmer) + "_cov_" + str(
                    self.vertex_info[this_k_name].cov))

        # 2.3 add contigs starting from terminal kmers
        # use waiting_list to record downstream kmers, so that we can label contigs in connecting order
        waiting_terminal_kmers = []
        for this_k_name in sorted(point_kmer_link_quota):
            while True:  # this_k_name in point_kmer_link_quota:
                # if this_k_name in point_kmer_link_quota:
                for this_k_end in (True, False):
                    for next_k_name, next_k_end in sorted(self.vertex_info[this_k_name].connections[this_k_end]):
                        # for every start
                        # check point_kmer_link_quota
                        # check next_k: 1) valid intermediate kmers or 2) valid point kmers
                        if point_kmer_link_quota[this_k_name][this_k_end] > 0 and (
                                next_k_name in intermediate_kmers or (
                                next_k_name in point_kmer_link_quota and
                                point_kmer_link_quota[next_k_name][next_k_end] > 0)):
                            generate_contig(
                                this_k_n=this_k_name, this_k_e=this_k_end,
                                next_k_n=next_k_name, next_k_e=next_k_end, terminal_k_quota=point_kmer_link_quota,
                                in_graph=assembly_graph, in_kmer_e_to_ctg_e=joint_kmer_end_to_contig_end,
                                inter_ks=intermediate_kmers, next_kmers=waiting_terminal_kmers)
                if waiting_terminal_kmers:
                    this_k_name = waiting_terminal_kmers.pop(0)
                else:
                    break
        # 2.4 add contigs starting from remaining kmers (circular sub-graphs)
        # for this_k_name in sorted(point_kmer_link_quota):
        for this_k_name in sorted(intermediate_kmers):
            if this_k_name in intermediate_kmers:
                for this_k_end in (True, False):
                    # clean point_kmer_link_quota
                    # if point_kmer_link_quota[this_k_name][True] <= 0 and \
                    #         point_kmer_link_quota[this_k_name][False] <= 0:
                    #     del point_kmer_link_quota[this_k_name]
                    #     break
                    # else:
                    for next_k_name, next_k_end in sorted(self.vertex_info[this_k_name].connections[this_k_end]):
                        if this_k_name in intermediate_kmers and next_k_name in intermediate_kmers:
                            generate_contig(this_k_n=this_k_name, this_k_e=this_k_end,
                                            next_k_n=next_k_name, next_k_e=next_k_end,
                                            terminal_k_quota=point_kmer_link_quota,
                                            in_graph=assembly_graph, inter_ks=intermediate_kmers,
                                            in_kmer_e_to_ctg_e=joint_kmer_end_to_contig_end, next_kmers=[])
        # 2.5 transfer the connections among kmers to contigs
        for joint_kmer_n in joint_kmer_end_to_contig_end:
            for contig_n_1, contig_e_1 in joint_kmer_end_to_contig_end[joint_kmer_n][True]:
                for contig_n_2, contig_e_2 in joint_kmer_end_to_contig_end[joint_kmer_n][False]:
                    assembly_graph.vertex_info[contig_n_1].connections[contig_e_1][(contig_n_2, contig_e_2)] = \
                        assembly_graph.vertex_info[contig_n_2].connections[contig_e_2][(contig_n_1, contig_e_1)] = \
                        self.__kmer
        return assembly_graph


class SpadesNodes(Vertex):
    def __init__(self, v_name, length=None, coverage=None, forward_seq=None, reverse_seq=None,
                 tail_connections=None, head_connections=None, path_str=None, assembly_obj=None):
        vertices_path = []
        if path_str:
            for this_v_str in path_str.strip().split(","):
                vertices_path.append((this_v_str[:-1], VERTEX_DIRECTION_STR_TO_BOOL[this_v_str[-1]]))
            if assembly_obj:
                path_seq = assembly_obj.export_path(vertices_path).seq
                if forward_seq:
                    if path_seq != forward_seq:
                        raise ProcessingGraphFailed("incompatible assembly_graph and scaffolds.fasta/scaffolds.paths!")
                else:
                    forward_seq = path_seq
                    length = len(forward_seq)
        super(SpadesNodes, self).__init__(v_name=v_name, length=length, coverage=coverage, forward_seq=forward_seq,
                                          reverse_seq=reverse_seq, tail_connections=tail_connections,
                                          head_connections=head_connections)
        self.vertices_path = vertices_path

    def __eq__(self, other):
        return self.name == other.name and self.len == other.len and self.cov == other.cov and self.seq == other.seq


class SpadesScaffolds(object):
    def __init__(self, scaffold_fasta, scaffold_paths, assembly_obj, min_cov=0., max_cov=inf,
                 min_contig_len=11, matching_ws=15, log_handler=None):
        if not os.path.exists(scaffold_fasta):
            raise FileNotFoundError(scaffold_fasta + " not found!")
        if not os.path.exists(scaffold_paths):
            raise FileNotFoundError(scaffold_paths + " not found!")
        self.nodes = {}
        self.vertex_jumping_over_gap = []
        # each_item = [(left_name, left_end, left_trim),
        #              (right_name, right_end, right_trim),
        #              (gap_name, gap_seq, gap_len)]
        sequence_matrix = SequenceList(scaffold_fasta, indexed=True)
        graph_overlap = assembly_obj.uni_overlap()
        with open(scaffold_paths) as path_handler:
            go_patch = 1
            line = path_handler.readline().strip()
            while line:
                # supposing it's recorded in balance by SPAdes, skip the reverse direction record
                if line.startswith("NODE") and not line.endswith("'"):
                    n_tag, node_name, l_tag, node_len, c_tag, node_cov = line.split("_")
                    long_name = line
                    # skip nodes with cov out of bounds
                    node_cov = float(node_cov)
                    if not (min_cov <= node_cov <= max_cov):
                        continue
                    # read paths
                    path_strings = []
                    line = path_handler.readline().strip()
                    path_strings.append(line.strip(";"))
                    while line and line.endswith(";"):
                        line = path_handler.readline().strip()
                        path_strings.append(line.strip(";"))
                    # read sequences
                    if len(path_strings) != 1:
                        full_path_seq = sequence_matrix[long_name].seq
                        last_v_name, last_v_end = "", False
                        # arbitrary_break = False  # cannot find the position of this sub-seq in scaffolds.fasta
                        gap_len = full_path_seq.count("N")
                        last_matching_end_at_scaffold = None
                        last_end_at_scaffold = None
                        for go_sub_path, sub_path in enumerate(path_strings):
                            sub_name = node_name + "S" + str(go_sub_path)
                            try:
                                self.nodes[sub_name] = SpadesNodes(sub_name, coverage=node_cov,
                                                                   path_str=sub_path, assembly_obj=assembly_obj)
                            except KeyError:  # skip paths with pruned nodes
                                break
                            this_sub_seq = str(self.nodes[sub_name].seq[True])
                            len_this_sub = len(this_sub_seq)
                            if len_this_sub <= min_contig_len:
                                break
                            check_ws = min(matching_ws, len_this_sub)
                            # if go_sub_path == 0:
                            #     searching_min = 0
                            #     searching_max = 0
                            # else:
                            #     if arbitrary_break:
                            #         if last_matching_end_at_scaffold is None:
                            #             searching_min = None
                            #         else:
                            #             searching_min = last_matching_end_at_scaffold - check_ws
                            #         searching_max = None
                            #     else:
                            #         # for speeding up
                            #         searching_min = last_matching_end_at_scaffold - check_ws
                            #         searching_max = last_matching_end_at_scaffold - check_ws + gap_len + len_this_sub
                            try:
                                this_start_at_scaffold, matching_start_at_alm, matching_end_at_alm = \
                                    map_contigs_to_scaffolds(scaffold=full_path_seq, contig=this_sub_seq,
                                                             scaffold_start_searching_min=None,  # searching_min,
                                                             scaffold_start_searching_max=None,  # searching_max,
                                                             word_size=check_ws)
                            except AssertionError:  # skip pass with abnormal path_seq_len
                                break
                            # check matching location
                            # print((last_v_name, last_v_end), self.nodes[sub_name].vertices_path[0],
                            #       this_start_at_scaffold, matching_start_at_alm, matching_end_at_alm)
                            # print(this_sub_seq[:20])
                            if this_start_at_scaffold is None or \
                                    (matching_end_at_alm - matching_start_at_alm + 1 < graph_overlap):
                                # arbitrary_break = True
                                continue
                            else:
                                # arbitrary_break = False
                                if last_end_at_scaffold is None or last_matching_end_at_scaffold is None:
                                    # if go_sub_path == 0:
                                    last_v_name, last_v_end = self.nodes[sub_name].vertices_path[-1]
                                    last_matching_end_at_scaffold = this_start_at_scaffold + matching_end_at_alm
                                    last_end_at_scaffold = this_start_at_scaffold + len_this_sub - 1
                                else:
                                    this_matching_start_at_scaffold = this_start_at_scaffold + matching_start_at_alm
                                    this_matching_end_at_scaffold = this_start_at_scaffold + matching_end_at_alm
                                    this_end_at_scaffold = this_start_at_scaffold + len_this_sub - 1
                                    trim_last = last_end_at_scaffold - last_matching_end_at_scaffold
                                    trim_this = this_matching_start_at_scaffold - this_start_at_scaffold
                                    # code checking
                                    assert trim_this >= 0
                                    assert trim_last >= 0
                                    # code checking
                                    seemingly_gap_len = this_start_at_scaffold - last_end_at_scaffold - 1
                                    if seemingly_gap_len >= 0:
                                        # if there's extra bases beyond the original two sub seqs
                                        gap_name = str(go_patch) + "GAP" + str(seemingly_gap_len)
                                    else:
                                        gap_name = str(go_patch) + "OVL" + str(-seemingly_gap_len)
                                    real_gap_len = this_matching_start_at_scaffold - last_matching_end_at_scaffold - 1
                                    this_v_name, this_v_end = self.nodes[sub_name].vertices_path[0]
                                    this_v_name, this_v_end = this_v_name, not this_v_end  # starting vertex
                                    last_fix = bool(assembly_obj.vertex_info[last_v_name].connections[last_v_end])
                                    this_fix = bool(assembly_obj.vertex_info[this_v_name].connections[this_v_end])
                                    go_to_standard_jumping = False
                                    if (last_fix and trim_last) or (this_fix and trim_this):
                                        if log_handler:
                                            log_handler.warning(
                                                str(-gap_len) + "-bp gap/overlap between " + last_v_name + " and " +
                                                this_v_name + " indicated while conflicting connections existed!")
                                        else:
                                            sys.stdout.write("Warning: " + str(-gap_len) + "-bp gap/overlap "
                                                             "between " + last_v_name + " and " + this_v_name +
                                                             " indicated while conflicting connections existed!\n")
                                    elif last_fix and this_fix:
                                        # trim_last == trim_this == 0
                                        # seemingly_gap_len == real_gap_len
                                        if -real_gap_len > graph_overlap:
                                            # if there's >kmer overlap but both of those vertices_set should fix
                                            # leading gap sequence to be illegal
                                            if log_handler:
                                                log_handler.warning(
                                                    str(-gap_len) + "-bp overlap between " + last_v_name + " and " +
                                                    this_v_name + " indicated while conflicting connections existed!")
                                            else:
                                                sys.stdout.write(
                                                    "Warning: " + str(-gap_len) + "-bp overlap between " + last_v_name
                                                    + " and " + this_v_name + " indicated while conflicting connections existed!\n")
                                        else:
                                            go_to_standard_jumping = True
                                    elif this_fix:
                                        # trim_last==0
                                        if -real_gap_len > graph_overlap:
                                            self.vertex_jumping_over_gap.append(
                                                [(last_v_name, last_v_end, trim_last - real_gap_len - graph_overlap),
                                                 (this_v_name, this_v_end, trim_this),
                                                 (gap_name, "", -graph_overlap)])
                                            go_patch += 1
                                        else:
                                            go_to_standard_jumping = True
                                    else:
                                        # last_fix or not (this_fix or last_fix)
                                        # trim_last==0
                                        if -real_gap_len > graph_overlap:
                                            self.vertex_jumping_over_gap.append(
                                                [(last_v_name, last_v_end, trim_last),
                                                 (this_v_name, this_v_end, trim_this - real_gap_len - graph_overlap),
                                                 (gap_name, "", -graph_overlap)])
                                            go_patch += 1
                                        else:
                                            go_to_standard_jumping = True
                                    if go_to_standard_jumping:
                                        this_gap_seq = full_path_seq[last_matching_end_at_scaffold + 1:
                                                                     this_matching_start_at_scaffold]
                                        self.vertex_jumping_over_gap.append(
                                            [(last_v_name, last_v_end, trim_last),
                                             (this_v_name, this_v_end, trim_this),
                                             (gap_name, this_gap_seq, real_gap_len)])
                                        go_patch += 1
                                    last_v_name, last_v_end = self.nodes[sub_name].vertices_path[-1]
                                    last_matching_end_at_scaffold = this_matching_end_at_scaffold
                                    last_end_at_scaffold = this_end_at_scaffold
                line = path_handler.readline().strip()


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
                    all_e[go_e] = hash(all_e[go_e])
        return all_e


def average_np_free(vals):
    return sum(vals) / float(len(vals))


def average_weighted_np_free(vals, weights):
    return sum([val * weights[go_v] for go_v, val in enumerate(vals)]) / float(sum(weights))


def weighted_mean_and_std_np_free(values, weights):
    mean = average_weighted_np_free(values, weights=weights)
    std = average_weighted_np_free([(val - mean) ** 2 for val in values], weights=weights) ** 0.5
    return mean, std


def get_graph_coverages_range_simple(fasta_matrix, drop_low_percent=0.10, drop_high_percent=0.40, drop_ssr=True):
    coverages = []
    lengths = []
    for go_seq, fastg_name in enumerate(fasta_matrix[0]):
        # remove sequences like "ATATATATATAT", "AGAGAGAGAGAGAG"
        if not (drop_ssr and len(set(fasta_matrix[1][go_seq])) < 3):
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
    try:
        cov_mean, cov_std = weighted_mean_and_std_np_free(coverages, lengths)
    except ZeroDivisionError:
        cov_mean, cov_std = 0., 0.
        coverages = [0.]
    return max(cov_mean - cov_std, min(coverages)), cov_mean, min(cov_mean + cov_std, max(coverages))
