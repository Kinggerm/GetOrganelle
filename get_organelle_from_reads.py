#!/usr/bin/env python

import datetime
import shutil
from copy import deepcopy
from math import log, exp
try:
    from math import inf
except ImportError:
    inf = float("inf")
from argparse import ArgumentParser
import GetOrganelleLib
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.pipe_control_func import *
import time
import random
import subprocess
import sys
import os
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
import platform
SYSTEM_NAME = ""
if platform.system() == "Linux":
    SYSTEM_NAME = "linux"
elif platform.system() == "Darwin":
    SYSTEM_NAME = "macOS"
else:
    sys.stdout.write("Error: currently GetOrganelle is not supported for " + platform.system() + "! ")
    exit()
GO_LIB_PATH = os.path.split(GetOrganelleLib.__file__)[0]
GO_DEP_PATH = os.path.realpath(os.path.join(GO_LIB_PATH, "..", "GetOrganelleDep", SYSTEM_NAME))
UTILITY_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "Utilities")

_GO_PATH = GO_PATH
_LBL_DB_PATH = LBL_DB_PATH
_SEQ_DB_PATH = SEQ_DB_PATH

MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    PYTHON_VERSION = "2.7+"
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    PYTHON_VERSION = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)


MAX_RATIO_RL_WS = 0.75
AUTO_MIN_WS = 49
AUTO_MIN_WS_ANIMAL_MT = 41
AUTO_MIN_WS_PLANT_MT = 55
GLOBAL_MIN_WS = 29

BASE_COV_SAMPLING_PERCENT = 0.06
GUESSING_FQ_GZIP_COMPRESSING_RATIO = 3.58
GUESSING_FQ_SEQ_INFLATE_TO_FILE = 3.22

SUPPORTED_ORGANELLE_TYPES = ["embplant_pt", "embplant_mt", "embplant_nr", "other_pt", "animal_mt", "fungus_mt", "fungus_nr"]
ORGANELLE_EXPECTED_GRAPH_SIZES = {"embplant_pt": 130000,
                                  "embplant_mt": 390000,
                                  "embplant_nr": 13000,
                                  "fungus_nr": 13000,
                                  "other_pt": 39000,
                                  "animal_mt": 13000,
                                  "fungus_mt": 65000}

READ_LINE_TO_INF = int(HEAD_MAXIMUM_LINES/4)


def get_options(description, version):
    version = version
    usage = "\n###  Embryophyta plant plastome, 2*(1G raw data, 150 bp) reads\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s cp_seed.fasta -o plastome_output " \
            " -R 15 -k 21,45,65,85,105 -F embplant_pt\n" \
            "###  Embryophyta plant mitogenome\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s mt_seed.fasta -o mitogenome_output " \
            " -R 30 -k 21,45,65,85,105 -F embplant_mt"
    parser = ArgumentParser(usage=usage, description=description, add_help=False)
    # simple help mode
    if "-h" in sys.argv:
        parser.add_argument("-1", dest="fq_file_1", help="Input file with forward paired-end reads (*.fq/.gz/.tar.gz).")
        parser.add_argument("-2", dest="fq_file_2", help="Input file with reverse paired-end reads (*.fq/.gz/.tar.gz).")
        parser.add_argument("-u", dest="unpaired_fq_files", help="Input file(s) with unpaired (single-end) reads. ")
        parser.add_argument("-o", dest="output_base", help="Output directory.")
        parser.add_argument("-s", dest="seed_file", help="Input fasta format file as initial seed. "
                                                         "Default: " + os.path.join(SEQ_DB_PATH, "*.fasta"))
        parser.add_argument("-w", dest="word_size", help="Word size (W) for extension. Default: auto-estimated")
        parser.add_argument("-R", dest="max_rounds", help="Maximum extension rounds (suggested: >=2). "
                                                          "Default: 15 (embplant_pt)")
        parser.add_argument("-F", dest="organelle_type",
                            help="Target organelle genome type(s): "
                                 "embplant_pt/other_pt/embplant_mt/embplant_nr/animal_mt/fungus_mt/fungus_nr/anonym/"
                                 "embplant_pt,embplant_mt/other_pt,embplant_mt,fungus_mt")
        parser.add_argument("--max-reads", type=float,
                            help="Maximum number of reads to be used per file. "
                                 "Default: 1.5E7 (-F embplant_pt/embplant_nr/fungus_mt/fungus_nr); "
                                 "7.5E7 (-F embplant_mt/other_pt/anonym); 3E8 (-F animal_mt)")
        parser.add_argument("--fast", dest="fast_strategy",
                            help="=\"-R 10 -t 4 -J 5 -M 7 --max-n-words 3E7 --larger-auto-ws "
                                 "--disentangle-time-limit 360\"")
        parser.add_argument("-k", dest="spades_kmer", default="21,55,85,115",
                            help="SPAdes kmer settings. Default: %(default)s")
        parser.add_argument("-t", dest="threads", type=int, default=1,
                            help="Maximum threads to use. Default: %(default)s")
        parser.add_argument("-P", dest="pre_grouped", default=int(2E5), help="Pre-grouping value. Default: %(default)s")
        parser.add_argument("-v", "--version", action="version",
                            version="GetOrganelle v{version}".format(version=version))
        parser.add_argument("-h", dest="simple_help", default=False, action="store_true",
                            help="print brief introduction for frequently-used options.")
        parser.add_argument("--help", dest="verbose_help", default=False, action="store_true",
                            help="print verbose introduction for all options.")
        parser.print_help()
        sys.stdout.write("\n")
        exit()
    else:
        # verbose help mode
        # group 1
        group_inout = parser.add_argument_group("IN-OUT OPTIONS", "Options on inputs and outputs")
        # group_inout = OptionGroup(parser, "IN-OUT OPTIONS", "Options on inputs and outputs")
        group_inout.add_argument("-1", dest="fq_file_1",
                                 help="Input file with forward paired-end reads (format: fastq/fastq.gz/fastq.tar.gz).")
        group_inout.add_argument("-2", dest="fq_file_2",
                                 help="Input file with reverse paired-end reads (format: fastq/fastq.gz/fastq.tar.gz).")
        group_inout.add_argument("-u", dest="unpaired_fq_files",
                                 help="Input file(s) with unpaired (single-end) reads (format: fastq/fastq.gz/fastq.tar.gz). "
                                      "files could be comma-separated lists such as 'seq1.fq,seq2.fq'.")
        group_inout.add_argument("-o", dest="output_base",
                                 help="Output directory. Overwriting files if directory exists.")
        group_inout.add_argument("-s", dest="seed_file", default=None,
                                 help="Seed sequence(s). Input fasta format file as initial seed. "
                                      "A seed sequence in GetOrganelle is only used for identifying initial "
                                      "organelle reads. The assembly process is purely de novo. "
                                      "Should be a list of files split by comma(s) on a multi-organelle mode, "
                                      "with the same list length to organelle_type (followed by '-F'). "
                                      "Default: '" + os.path.join(SEQ_DB_PATH, "*.fasta") + "' "
                                      "(* depends on the value followed with flag '-F')")
        group_inout.add_argument("-a", dest="anti_seed",
                                 help="Anti-seed(s). Not suggested unless what you really know what you are doing. "
                                      "Input fasta format file as anti-seed, where the extension process "
                                      "stop. Typically serves as excluding plastid reads when extending mitochondrial "
                                      "reads, or the other way around. You should be cautious about using this option, "
                                      "because if the anti-seed includes some word in the target but not in the seed, "
                                      "the result would have gaps. For example, use the embplant_mt and embplant_pt "
                                      "from the same plant-species as seed and anti-seed.")
        group_inout.add_argument("--max-reads", dest="maximum_n_reads", type=float, default=1.5E7,
                                 help="Hard bound for maximum number of reads to be used per file. "
                                      "A input larger than " + str(
                                     READ_LINE_TO_INF) + " will be treated as infinity (INF). "
                                                         "Default: 1.5E7 (-F embplant_pt/embplant_nr/fungus_mt/fungus_nr); "
                                                         "7.5E7 (-F embplant_mt/other_pt/anonym); 3E8 (-F animal_mt)")
        group_inout.add_argument("--reduce-reads-for-coverage", dest="reduce_reads_for_cov", type=float, default=500,
                                 help="Soft bound for maximum number of reads to be used according to "
                                      "target-hitting base coverage. "
                                      "If the estimated target-hitting base coverage is too high and "
                                      "over this VALUE, GetOrganelle automatically reduce the number of reads to "
                                      "generate a final assembly with base coverage close to this VALUE. "
                                      "This design could greatly save computational resources in many situations. "
                                      "A mean base coverage over 500 is extremely sufficient for most cases. "
                                      "This VALUE must be larger than 10. Set this VALUE to inf to disable reducing. "
                                      "Default: %(default)s.")
        group_inout.add_argument("--max-ignore-percent", dest="maximum_ignore_percent", type=float, default=0.01,
                                 help="The maximum percent of bases to be ignore in extension, due to low quality. "
                                      "Default: %(default)s")
        group_inout.add_argument("--phred-offset", dest="phred_offset", default=-1, type=int,
                                 help="Phred offset for spades-hammer. Default: GetOrganelle-autodetect")
        group_inout.add_argument("--min-quality-score", dest="min_quality_score", type=int, default=1,
                                 help="Minimum quality score in extension. This value would be automatically decreased "
                                      "to prevent ignoring too much raw data (see --max-ignore-percent)."
                                      "Default: %(default)s ('\"' in Phred+33; 'A' in Phred+64/Solexa+64)")
        group_inout.add_argument("--prefix", dest="prefix", default="",
                                 help="Add extra prefix to resulting files under the output directory.")
        group_inout.add_argument("--out-per-round", dest="fg_out_per_round", action="store_true", default=False,
                                 help="Enable output per round. Choose to save memory but cost more time per round.")
        group_inout.add_argument("--zip-files", dest="zip_files", action="store_true", default=False,
                                 help="Choose to compress fq/sam files using gzip.")
        group_inout.add_argument("--keep-temp", dest="keep_temp_files", action="store_true", default=False,
                                 help="Choose to keep the running temp/index files.")
        group_inout.add_argument("--config-dir", dest="get_organelle_path", default=None,
                                 help="The directory where the configuration file and default databases were placed. "
                                      "The default value also can be changed by adding 'export GETORG_PATH=your_favor' "
                                      "to the shell script (e.g. ~/.bash_profile or ~/.bashrc) "
                                      "Default: " + GO_PATH)
        # group 2
        group_scheme = parser.add_argument_group("SCHEME OPTIONS", "Options on running schemes.")
        group_scheme.add_argument("-F", dest="organelle_type",
                                  help="This flag should be followed with embplant_pt (embryophyta plant plastome), "
                                       "other_pt (non-embryophyta plant plastome), embplant_mt "
                                       "(plant mitogenome), embplant_nr (plant nuclear ribosomal RNA), animal_mt "
                                       "(animal mitogenome), fungus_mt (fungus mitogenome), "
                                       "fungus_nr (fungus nuclear ribosomal RNA)"
                                       "or embplant_mt,other_pt,fungus_mt "
                                       "(the combination of any of above organelle genomes split by comma(s), "
                                       "which might be computationally more intensive than separate runs), "
                                       "or anonym (uncertain organelle genome type). "
                                       "The anonym should be used with customized seed and label databases "
                                       "('-s' and '--genes'). "
                                       "For easy usage and compatibility of old versions, following redirection "
                                       "would be automatically fulfilled without warning:\t"
                                       "\nplant_cp->embplant_pt; plant_pt->embplant_pt; "
                                       "\nplant_mt->embplant_mt; plant_nr->embplant_nr")
        group_scheme.add_argument("--fast", dest="fast_strategy", default=False, action="store_true",
                                  help="=\"-R 10 -t 4 -J 5 -M 7 --max-n-words 3E7 --larger-auto-ws "
                                       "--disentangle-time-limit 360\" "
                                       "This option is suggested for homogeneously and highly covered data (very fine data). "
                                       "You can overwrite the value of a specific option listed above by adding "
                                       "that option along with the \"--fast\" flag. "
                                       "You could try GetOrganelle with this option for a list of samples and run a second "
                                       "time without this option for the rest with incomplete results. ")
        group_scheme.add_argument("--memory-save", dest="memory_save", default=False, action="store_true",
                                  help="=\"--out-per-round -P 0 --remove-duplicates 0\" "
                                       "You can overwrite the value of a specific option listed above by adding "
                                       "that option along with the \"--memory-save\" flag. A larger '-R' value is suggested "
                                       "when \"--memory-save\" is chosen.")
        group_scheme.add_argument("--memory-unlimited", dest="memory_unlimited", default=False, action="store_true",
                                  help="=\"-P 1E7 --index-in-memory --remove-duplicates 2E8 "
                                       "--min-quality-score -5 --max-ignore-percent 0\" "
                                       "You can overwrite the value of a specific option listed above by adding "
                                       "that option along with the \"--memory-unlimited\" flag. ")
        # group 3
        group_extending = parser.add_argument_group("EXTENDING OPTIONS",
                                                    "Options on the performance of extending process")
        group_extending.add_argument("-w", dest="word_size", type=float,
                                     help="Word size (W) for pre-grouping (if not assigned by '--pre-w') and extending "
                                          "process. This script would try to guess (auto-estimate) a proper W "
                                          "using an empirical function based on average read length, reads quality, "
                                          "target genome coverage, and other variables that might influence the extending "
                                          "process. You could assign the ratio (1>input>0) of W to "
                                          "read_length, based on which this script would estimate the W for you; "
                                          "or assign an absolute W value (read length>input>=35). Default: auto-estimated.")
        group_extending.add_argument("--pre-w", dest="pregroup_word_size", type=float,
                                     help="Word size (W) for pre-grouping. Used to reproduce result when word size is "
                                          "a certain value during pregrouping process and later changed during reads "
                                          "extending process. Similar to word size. Default: the same to word size.")
        group_extending.add_argument("-R", "--max-rounds", dest="max_rounds", type=int,  # default=inf,
                                     help="Maximum number of extending rounds (suggested: >=2). "
                                          "Default: 15 (-F embplant_pt), 30 (-F embplant_mt/other_pt), "
                                          "10 (-F embplant_nr/animal_mt/fungus_mt/fungus_nr), inf (-P 0).")
        group_extending.add_argument("--max-n-words", dest="maximum_n_words", type=float, default=4E8,
                                     help="Maximum number of words to be used in total."
                                          "Default: 4E8 (-F embplant_pt), 2E8 (-F embplant_nr/fungus_mt/fungus_nr/animal_mt), "
                                          "2E9 (-F embplant_mt/other_pt)")
        group_extending.add_argument("-J", dest="jump_step", type=int, default=3,
                                     help="The length of step for checking words in reads during extending process "
                                          "(integer >= 1). When you have reads of high quality, the larger the number is, "
                                          "the faster the extension will be, "
                                          "the more risk of missing reads in low coverage area. "
                                          "Choose 1 to choose the slowest but safest extension strategy. Default: %(default)s")
        group_extending.add_argument("-M", dest="mesh_size", type=int, default=2,
                                     help="(Beta parameter) "
                                          "The length of step for building words from seeds during extending process "
                                          "(integer >= 1). When you have reads of high quality, the larger the number is, "
                                          "the faster the extension will be, "
                                          "the more risk of missing reads in low coverage area. "
                                          "Another usage of this mesh size is to choose a larger mesh size coupled with a "
                                          "smaller word size, which makes smaller word size feasible when memory is limited."
                                          "Choose 1 to choose the slowest but safest extension strategy. Default: %(default)s")
        group_extending.add_argument("--bowtie2-options", dest="bowtie2_options", default="--very-fast -t",
                                     help="Bowtie2 options, such as '--ma 3 --mp 5,2 --very-fast -t'. Default: %(default)s.")
        group_extending.add_argument("--larger-auto-ws", dest="larger_auto_ws", default=False, action="store_true",
                                     help="By using this flag, the empirical function for estimating W would tend to "
                                          "produce a relative larger W, which would speed up the matching in extending, "
                                          "reduce the memory cost in extending, but increase the risk of broken final "
                                          "graph. Suggested when the data is good with high and homogenous coverage.")
        mixed_organelles = ("other_pt", "embplant_mt", "fungus_mt")
        group_extending.add_argument("--target-genome-size", dest="target_genome_size", default='130000', type=str,
                                     help="Hypothetical value(s) of target genome size. This is only used for estimating "
                                          "word size when no '-w word_size' is given. "
                                          "Should be a list of INTEGER numbers split by comma(s) on a multi-organelle mode, "
                                          "with the same list length to organelle_type (followed by '-F'). "
                                          "Default: " +
                                          " or ".join(
                                              [str(
                                                  ORGANELLE_EXPECTED_GRAPH_SIZES[this_type]) + " (-F " + this_type + ")"
                                               for this_type in SUPPORTED_ORGANELLE_TYPES]) + " or " +
                                          ",".join([str(ORGANELLE_EXPECTED_GRAPH_SIZES[this_type])
                                                    for this_type in mixed_organelles]) +
                                          " (-F " + ",".join(mixed_organelles) + ")")
        group_extending.add_argument("--max-extending-len", dest="max_extending_len", type=str,
                                     help="Maximum extending length(s) derived from the seed(s). "
                                          "A single value could be a non-negative number, or inf (infinite) "
                                          "or auto (automatic estimation). "
                                          "This is designed for properly stopping the extending from getting too long and "
                                          "saving computational resources. However, empirically, a maximum extending length "
                                          "value larger than 6000 would not be helpful for saving computational resources. "
                                          "This value would not be precise in controlling output size, especially "
                                          "when pre-group (followed by '-P') is turn on."
                                          "In the auto mode, the maximum extending length is estimated based on the sizes of "
                                          "the gap regions that not covered in the seed sequences. A sequence of a closely "
                                          "related species would be preferred for estimating a better maximum extending "
                                          "length value. If you are using limited loci, e.g. rbcL gene as the seed for "
                                          "assembling the whole plastome (with extending length ca. 75000 >> 6000), "
                                          "you should set maximum extending length to inf. "
                                          "Should be a list of numbers/auto/inf split by comma(s) on a multi-organelle mode, "
                                          "with the same list length to organelle_type (followed by '-F'). "
                                          "Default: inf. ")
        # group 4
        group_assembly = parser.add_argument_group("ASSEMBLY OPTIONS", "These options are about the assembly and "
                                                                       "graph disentangling")
        group_assembly.add_argument("-k", dest="spades_kmer", default="21,55,85,115",
                                    help="SPAdes kmer settings. Use the same format as in SPAdes. illegal kmer values "
                                         "would be automatically discarded by GetOrganelle. "
                                         "Default: %(default)s")
        group_assembly.add_argument("--spades-options", dest="other_spades_options", default="",
                                    help="Other SPAdes options. Use double quotation marks to include all "
                                         "the arguments and parameters.")
        group_assembly.add_argument("--no-spades", dest="run_spades", action="store_false", default=True,
                                    help="Disable SPAdes.")
        group_assembly.add_argument("--ignore-k", dest="ignore_kmer_res", default=40, type=int,
                                    help="A kmer threshold below which, no slimming/disentangling would be executed"
                                         " on the result. Default: %(default)s")
        group_assembly.add_argument("--genes", dest="genes_fasta",
                                    help="Followed with a customized database (a fasta file or the base name of a "
                                         "blast database) containing or made of ONE set of protein coding genes "
                                         "and ribosomal RNAs extracted from ONE reference genome that you want to assemble. "
                                         "Should be a list of databases split by comma(s) on a multi-organelle mode, "
                                         "with the same list length to organelle_type (followed by '-F'). "
                                         "This is optional for any organelle mentioned in '-F' but required for 'anonym'. "
                                         "By default, certain database(s) in " + str(LBL_DB_PATH) + " would be used "
                                         "contingent on the organelle types chosen (-F). "
                                         "The default value become invalid when '--genes' or '--ex-genes' is used.")
        group_assembly.add_argument("--ex-genes", dest="exclude_genes",
                                    help="This is optional and Not suggested, since non-target contigs could contribute "
                                         "information for better downstream coverage-based clustering. "
                                         "Followed with a customized database (a fasta file or the base name of a "
                                         "blast database) containing or made of protein coding genes "
                                         "and ribosomal RNAs extracted from reference genome(s) that you want to exclude. "
                                         "Could be a list of databases split by comma(s) but "
                                         "NOT required to have the same list length to organelle_type (followed by '-F'). "
                                         "The default value will become invalid when '--genes' or '--ex-genes' is used.")
        group_assembly.add_argument("--disentangle-df", dest="disentangle_depth_factor", default=5.0, type=float,
                                    help="Depth factor for differentiate genome type of contigs. "
                                         "The genome type of contigs are determined by blast. "
                                         "Default: %(default)s")
        group_assembly.add_argument("--disentangle-tf", dest="disentangle_type_factor", type=float, default=3.,
                                    help="Type factor for identifying contig type tag when multiple tags exist in one contig. "
                                         "Default:%(default)s")
        group_assembly.add_argument("--contamination-depth", dest="contamination_depth", default=3., type=float,
                                    help="Depth factor for confirming contamination in parallel contigs. Default: %(default)s")
        group_assembly.add_argument("--contamination-similarity", dest="contamination_similarity", default=0.9,
                                    type=float,
                                    help="Similarity threshold for confirming contaminating contigs. Default: %(default)s")
        group_assembly.add_argument("--no-degenerate", dest="degenerate", default=True, action="store_false",
                                    help="Disable making consensus from parallel contig based on nucleotide degenerate table.")
        group_assembly.add_argument("--degenerate-depth", dest="degenerate_depth", default=1.5, type=float,
                                    help="Depth factor for confirming parallel contigs. Default: %(default)s")
        group_assembly.add_argument("--degenerate-similarity", dest="degenerate_similarity", default=0.98, type=float,
                                    help="Similarity threshold for confirming parallel contigs. Default: %(default)s")
        group_assembly.add_argument("--disentangle-time-limit", dest="disentangle_time_limit", default=1800, type=int,
                                    help="Time limit (second) for each try of disentangling a graph file as a circular "
                                         "genome. Disentangling a graph as contigs is not limited. Default: %(default)s")
        group_assembly.add_argument("--expected-max-size", dest="expected_max_size", default='250000', type=str,
                                    help="Expected maximum target genome size(s) for disentangling. "
                                         "Should be a list of INTEGER numbers split by comma(s) on a multi-organelle mode, "
                                         "with the same list length to organelle_type (followed by '-F'). "
                                         "Default: 250000 (-F embplant_pt/fungus_mt), "
                                         "25000 (-F embplant_nr/animal_mt/fungus_nr), 1000000 (-F embplant_mt/other_pt),"
                                         "1000000,1000000,250000 (-F other_pt,embplant_mt,fungus_mt)")
        group_assembly.add_argument("--expected-min-size", dest="expected_min_size", default=10000, type=str,
                                    help="Expected minimum target genome size(s) for disentangling. "
                                         "Should be a list of INTEGER numbers split by comma(s) on a multi-organelle mode, "
                                         "with the same list length to organelle_type (followed by '-F'). "
                                         "Default: %(default)s for all.")
        group_assembly.add_argument("--reverse-lsc", dest="reverse_lsc", default=False, action="store_true",
                                    help="For '-F embplant_pt' with complete circular result, "
                                         "by default, the direction of the starting contig (usually "
                                         "the LSC region) is determined as the direction with less ORFs. Choose this option "
                                         "to reverse the direction of the starting contig when result is circular. "
                                         "Actually, both directions are biologically equivalent to each other. The "
                                         "reordering of the direction is only for easier downstream analysis.")
        group_assembly.add_argument("--max-paths-num", dest="max_paths_num", default=1000, type=int,
                                    help="Repeats would dramatically increase the number of potential isomers (paths). "
                                         "This option was used to export a certain amount of paths out of all possible paths "
                                         "per assembly graph. Default: %(default)s")
        # group 5
        group_computational = parser.add_argument_group("ADDITIONAL OPTIONS", "")
        group_computational.add_argument("-t", dest="threads", type=int, default=1,
                                         help="Maximum threads to use.")
        group_computational.add_argument("-P", dest="pre_grouped", type=float, default=2E5,
                                         help="The maximum number (integer) of high-covered reads to be pre-grouped "
                                              "before extending process. pre_grouping is suggested when the whole genome "
                                              "coverage is shallow but the organ genome coverage is deep. "
                                              "The default value is 2E5. "
                                              "For personal computer with 8G memory, we suggest no more than 3E5. "
                                              "A larger number (ex. 6E5) would run faster but exhaust memory "
                                              "in the first few minutes. Choose 0 to disable this process.")
        group_computational.add_argument("--which-blast", dest="which_blast", default="",
                                         help="Assign the path to BLAST binary files if not added to the path. "
                                              "Default: try \"" + os.path.realpath(GO_DEP_PATH) +
                                              "/ncbi-blast\" first, then $PATH")
        group_computational.add_argument("--which-bowtie2", dest="which_bowtie2", default="",
                                         help="Assign the path to Bowtie2 binary files if not added to the path. "
                                              "Default: try \"" + os.path.realpath(GO_DEP_PATH) +
                                              "/bowtie2\" first, then $PATH")
        group_computational.add_argument("--which-spades", dest="which_spades", default="",
                                         help="Assign the path to SPAdes binary files if not added to the path. "
                                              "Default: try \"" + os.path.realpath(GO_DEP_PATH) +
                                              "/SPAdes\" first, then $PATH")
        group_computational.add_argument("--which-bandage", dest="which_bandage", default="",
                                         help="Assign the path to bandage binary file if not added to the path. "
                                              "Default: try $PATH")
        group_computational.add_argument("--continue", dest="script_resume", default=False, action="store_true",
                                         help="Several check points based on files produced, rather than on the log file, "
                                              "so keep in mind that this script will NOT detect the difference "
                                              "between this input parameters and the previous ones.")
        group_computational.add_argument("--overwrite", dest="script_overwrite", default=False, action="store_true",
                                         help="Overwrite previous file if existed. ")
        group_computational.add_argument("--index-in-memory", dest="index_in_memory", action="store_true",
                                         default=False,
                                         help="Keep index in memory. Choose save index in memory than in disk.")
        group_computational.add_argument("--remove-duplicates", dest="rm_duplicates", default=1E7, type=float,
                                         help="By default this script use unique reads to extend. Choose the number of "
                                              "duplicates (integer) to be saved in memory. A larger number (ex. 2E7) would "
                                              "run faster but exhaust memory in the first few minutes. "
                                              "Choose 0 to disable this process. "
                                              "Note that whether choose or not will not disable "
                                              "the calling of replicate reads. Default: %(default)s.")
        group_computational.add_argument("--flush-step", dest="echo_step", default=54321,
                                         help="Flush step (INTEGER OR INF) for presenting progress. "
                                              "For running in the background, you could set this to inf, "
                                              "which would disable this. Default: %(default)s")
        group_computational.add_argument("--random-seed", dest="random_seed", default=12345, type=int,
                                         help="Default: %(default)s")
        group_computational.add_argument("--verbose", dest="verbose_log", action="store_true", default=False,
                                         help="Verbose output. Choose to enable verbose running log_handler.")

        parser.add_argument("-v", "--version", action="version",
                            version="GetOrganelle v{version}".format(version=version))
        parser.add_argument("-h", dest="simple_help", default=False, action="store_true",
                            help="print brief introduction for frequently-used options.")
        parser.add_argument("--help", dest="verbose_help", default=False, action="store_true",
                            help="print verbose introduction for all options.")
        if "--help" in sys.argv:
            parser.print_help()
            exit()

    # if "--help" in sys.argv:
    #     parser.add_option_group(group_inout)
    #     parser.add_option_group(group_scheme)
    #     parser.add_option_group(group_extending)
    #     parser.add_option_group(group_assembly)
    #     parser.add_option_group(group_computational)
    #
    # elif "-h" in sys.argv:
    #     for not_often_used in ("-a", "--max-ignore-percent", "--reduce-reads-for-coverage", "--phred-offset",
    #                            "--min-quality-score", "--prefix", "--out-per-round", "--zip-files", "--keep-temp",
    #                            "--config-dir",
    #                            "--memory-save", "--memory-unlimited", "--pre-w", "--max-n-words",
    #                            "-J", "-M", "--bowtie2-options",
    #                            "--larger-auto-ws", "--target-genome-size", "--spades-options", "--no-spades",
    #                            "--ignore-k", "--genes", "--ex-genes", "--disentangle-df",
    #                            "--contamination-depth", "--contamination-similarity", "--no-degenerate",
    #                            "--degenerate-depth", "--degenerate-similarity", "--disentangle-time-limit",
    #                            "--expected-max-size", "--expected-min-size", "--reverse-lsc", "--max-paths-num",
    #                            "--which-blast", "--which-bowtie2", "--which-spades", "--which-bandage",
    #                            "--continue", "--overwrite", "--index-in-memory",
    #                            "--remove-duplicates", "--flush-step", "--verbose"):
    #         parser.remove_option(not_often_used)
    #
    # else:
    #     parser.add_option_group(group_inout)
    #     parser.add_option_group(group_scheme)
    #     parser.add_option_group(group_extending)
    #     parser.add_option_group(group_assembly)
    #     parser.add_option_group(group_computational)
    # redirect organelle types before parsing arguments
    redirect_organelle_types = {"plant_cp": "embplant_pt",
                                "plant_pt": "embplant_pt",
                                "plant_mt": "embplant_mt",
                                "plant_nr": "embplant_nr"}
    for go_arg, candidate_arg in enumerate(sys.argv):
        if go_arg > 1 and sys.argv[go_arg - 1] in {"-F", "-E"}:
            if candidate_arg in redirect_organelle_types:
                sys.argv[go_arg] = redirect_organelle_types[candidate_arg]
            elif "," in candidate_arg:
                new_arg = []
                for sub_arg in candidate_arg.split(","):
                    if sub_arg in redirect_organelle_types:
                        new_arg.append(redirect_organelle_types[sub_arg])
                    else:
                        new_arg.append(sub_arg)
                sys.argv[go_arg] = ",".join(new_arg)
    #
    try:
        options = parser.parse_args()
    except Exception as e:
        sys.stderr.write("\n############################################################################\n" + str(e))
        sys.stderr.write("\n\"-h\" for more usage\n")
        exit()
    else:
        # if pos_args:
        #     sys.stderr.write("\n############################################################################"
        #                      "\nUnrecognized options: " + "\", \"".join(pos_args) + "\n")
        #     exit()
        if not ((options.fq_file_1 and options.fq_file_2) or options.unpaired_fq_files):
            sys.stderr.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\n")
            sys.stderr.write("Missing/Illegal input reads file(s) (followed after '-1&-2' and/or '-u')!\n")
            exit()
        if not options.output_base:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\n")
            sys.stderr.write("Missing option: output directory (followed after '-o')!\n")
            exit()
        if not options.organelle_type:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\n")
            sys.stderr.write("Missing option: organelle type (followed after '-F')!\n")
            exit()
        else:
            options.organelle_type = options.organelle_type.split(",")
        if int(bool(options.fq_file_1)) + int(bool(options.fq_file_2)) == 1:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: unbalanced paired reads!\n\n")
            exit()

        global _GO_PATH, _LBL_DB_PATH, _SEQ_DB_PATH
        if options.get_organelle_path:
            _GO_PATH = os.path.expanduser(options.get_organelle_path)
            if os.path.isdir(_GO_PATH):
                _LBL_DB_PATH = os.path.join(_GO_PATH, LBL_NAME)
                _SEQ_DB_PATH = os.path.join(_GO_PATH, SEQ_NAME)
            else:
                sys.stderr.write("\n############################################################################"
                                 "\nERROR: path " + _GO_PATH + " invalid!\n")
                exit()

        def _check_default_db(this_sub_organelle, extra_type=""):
            if not ((os.path.isfile(os.path.join(_LBL_DB_PATH, this_sub_organelle + ".fasta")) or options.genes_fasta)
                    and
                    (os.path.isfile(os.path.join(_SEQ_DB_PATH, this_sub_organelle + ".fasta")) or options.seed_file)):
                sys.stderr.write("\n############################################################################"
                                 "\nERROR: default " + this_sub_organelle + "," * int(bool(extra_type)) + extra_type +
                                 " database not added yet!\n"
                                 "\nInstall it by: get_organelle_config.py -a " + this_sub_organelle +
                                 "," * int(bool(extra_type)) + extra_type +
                                 "\nor\nInstall all types by: get_organelle_config.py -a all\n")
                exit()
        for sub_organelle_t in options.organelle_type:
            if sub_organelle_t not in {"embplant_pt", "other_pt", "embplant_mt", "embplant_nr", "animal_mt",
                                       "fungus_mt", "fungus_nr", "anonym"}:
                sys.stderr.write("\n############################################################################"
                                 "\nERROR: \"-F\" MUST be one of 'embplant_pt', 'other_pt', 'embplant_mt', "
                                 "'embplant_nr', 'animal_mt', 'fungus_mt', 'fungus_nr', 'anonym', "
                                 "or a combination of above split by comma(s)!\n\n")
                exit()
            elif sub_organelle_t == "anonym":
                if not options.seed_file or not options.genes_fasta:
                    sys.stderr.write("\n############################################################################"
                                     "\nERROR: \"-s\" and \"--genes\" must be specified when \"-F anonym\"!\n\n")
                    exit()
            else:
                if sub_organelle_t in ("embplant_pt", "embplant_mt"):
                    for go_t, check_sub in enumerate(["embplant_pt", "embplant_mt"]):
                        _check_default_db(check_sub, ["embplant_pt", "embplant_mt"][not go_t])
                else:
                    _check_default_db(sub_organelle_t)

        organelle_type_len = len(options.organelle_type)

        if not options.seed_file:
            use_default_seed = True
            options.seed_file = [os.path.join(_SEQ_DB_PATH, sub_o + ".fasta") for sub_o in options.organelle_type]
        else:
            use_default_seed = False
            options.seed_file = str(options.seed_file).split(",")
            if len(options.seed_file) != organelle_type_len:
                sys.stderr.write("\n############################################################################"
                                 "\nERROR: -F is followed with " + str(organelle_type_len) + " organelle types, " +
                                 "while -s is followed with " + str(len(options.seed_file)) + " file(s)!\n")
                exit()
        for check_file in [options.fq_file_1, options.fq_file_2, options.anti_seed] + options.seed_file:
            if check_file:
                if not os.path.exists(check_file):
                    sys.stderr.write("\n############################################################################"
                                     "\nERROR: " + check_file + " not found!\n\n")
                    exit()
                if os.path.getsize(check_file) == 0:
                    sys.stderr.write("\n############################################################################"
                                     "\nERROR: " + check_file + " is empty!\n\n")
                    exit()
        if options.unpaired_fq_files:
            options.unpaired_fq_files = options.unpaired_fq_files.split(",")
            for fastq_file in options.unpaired_fq_files:
                if not os.path.exists(fastq_file):
                    sys.stderr.write("\n############################################################################"
                                     "\nERROR: " + fastq_file + " not found!\n\n")
                    exit()
        else:
            options.unpaired_fq_files = []
        if options.jump_step < 1:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: Jump step MUST be an integer that >= 1\n")
            exit()
        if options.mesh_size < 1:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: Mesh size MUST be an integer that >= 1\n")
            exit()
        if options.fq_file_1 == options.fq_file_2 and options.fq_file_1:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: 1st fastq file is the same with 2nd fastq file!\n")
            exit()
        if options.memory_save and options.memory_unlimited:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: \"--memory-save\" and \"--memory-unlimited\" are not compatible!\n")
        assert options.threads > 0
        if options.reduce_reads_for_cov < 10:
            sys.stderr.write("\n############################################################################"
                             "\nERROR: value after \"--reduce-reads-for-coverage\" must be larger than 10!\n")
            exit()
        if options.echo_step == "inf":
            options.echo_step = inf
        elif type(options.echo_step) == int:
            pass
        elif type(options.echo_step) == str:
            try:
                options.echo_step = int(float(options.echo_step))
            except ValueError:
                sys.stderr.write("\n############################################################################"
                                 "\n--flush-step should be followed by positive integer or inf!\n")
                exit()
        assert options.echo_step > 0
        assert options.max_paths_num > 0
        assert options.phred_offset in (-1, 64, 33)
        assert options.script_resume + options.script_overwrite < 2, "'--overwrite' conflicts with '--continue'"
        options.prefix = os.path.basename(options.prefix)
        if os.path.isdir(options.output_base):
            if options.script_resume:
                previous_attributes = LogInfo(options.output_base, options.prefix).__dict__
            else:
                if options.script_overwrite:
                    try:
                        shutil.rmtree(options.output_base)
                    except OSError as e:
                        sys.stderr.write(
                            "\n############################################################################"
                            "\nRemoving existed " + options.output_base + " failed! "
                            "\nPlease manually remove it or use a new output directory!\n")
                    os.mkdir(options.output_base)
                else:
                    sys.stderr.write("\n############################################################################"
                                     "\n" + options.output_base + " existed! "
                                     "\nPlease use a new output directory, or use '--continue'/'--overwrite'\n")
                    exit()
                previous_attributes = {}
        else:
            options.script_resume = False
            os.mkdir(options.output_base)
            previous_attributes = {}
        # if options.script_resume and os.path.isdir(options.output_base):
        #     previous_attributes = LogInfo(options.output_base, options.prefix).__dict__
        # else:
        #     previous_attributes = {}
        #     options.script_resume = False
        # if not os.path.isdir(options.output_base):
        #     os.mkdir(options.output_base)
        #     options.script_resume = False
        log_handler = simple_log(logging.getLogger(), options.output_base, options.prefix + "get_org.")
        log_handler.info("")
        log_handler.info(description)
        log_handler.info("Python " + str(sys.version).replace("\n", " "))
        log_handler.info("PLATFORM: " + " ".join(platform.uname()))
        # log versions of dependencies
        lib_versions_info = []
        lib_not_available = []
        lib_versions_info.append("GetOrganelleLib " + GetOrganelleLib.__version__)
        try:
            import numpy as np
        except ImportError:
            lib_not_available.append("numpy")
        else:
            lib_versions_info.append("numpy " + np.__version__)
        # try:
        #     import sympy
        # except ImportError:
        #     lib_not_available.append("sympy")
        # else:
        #     lib_versions_info.append("sympy " + sympy.__version__)
        try:
            import gekko
        except ImportError:
            lib_not_available.append("gekko")
        else:
            lib_versions_info.append("gekko " + gekko.__version__)
        try:
            import psutil
        except ImportError:
            pass
        else:
            lib_versions_info.append("psutil " + psutil.__version__)
        log_handler.info("PYTHON LIBS: " + "; ".join(lib_versions_info))
        options.which_bowtie2 = detect_bowtie2_path(options.which_bowtie2, GO_DEP_PATH)
        if options.run_spades:
            options.which_spades = detect_spades_path(options.which_spades, GO_DEP_PATH)
        options.which_blast = detect_blast_path(options.which_blast, GO_DEP_PATH)
        dep_versions_info = []
        dep_versions_info.append(detect_bowtie2_version(options.which_bowtie2))
        if options.run_spades:
            dep_versions_info.append(detect_spades_version(options.which_spades))
        dep_versions_info.append(detect_blast_version(options.which_blast))
        if executable(os.path.join(options.which_bandage, "Bandage -v")):
            dep_versions_info.append(detect_bandage_version(options.which_bandage))
        log_handler.info("DEPENDENCIES: " + "; ".join(dep_versions_info))
        # log database
        log_handler.info("GETORG_PATH=" + _GO_PATH)
        existing_seed_db, existing_label_db = get_current_db_versions(db_type="both", seq_db_path=_SEQ_DB_PATH,
                                                                      lbl_db_path=_LBL_DB_PATH, silent=True)
        if use_default_seed:
            log_seed_types = deepcopy(options.organelle_type)
            if "embplant_pt" in log_seed_types and "embplant_mt" not in log_seed_types:
                log_seed_types.append("embplant_mt")
            if "embplant_mt" in log_seed_types and "embplant_pt" not in log_seed_types:
                log_seed_types.append("embplant_pt")
            log_handler.info("SEED  DB: " + single_line_db_versions(existing_seed_db, log_seed_types))
        if not options.genes_fasta:
            log_label_types = deepcopy(options.organelle_type)
            if "embplant_pt" in log_label_types and "embplant_mt" not in log_label_types:
                log_label_types.append("embplant_mt")
            if "embplant_mt" in log_label_types and "embplant_pt" not in log_label_types:
                log_label_types.append("embplant_pt")
            log_handler.info("LABEL DB: " + single_line_db_versions(existing_label_db, log_label_types))
        # working directory
        log_handler.info("WORKING DIR: " + os.getcwd())
        log_handler.info(" ".join(["\"" + arg + "\"" if " " in arg else arg for arg in sys.argv]) + "\n")

        # if options.run_spades:
        # space is forbidden for both spades and blast
        for fq_file in [options.fq_file_1, options.fq_file_2] * int(bool(options.fq_file_1 and options.fq_file_2))\
                       + options.unpaired_fq_files:
            assert is_valid_path(os.path.basename(fq_file)), \
                "Invalid characters (e.g. space, non-ascii) for SPAdes in file name: " + os.path.basename(fq_file)
        for fq_file in [options.output_base, options.prefix]:
            assert is_valid_path(os.path.realpath(fq_file)), \
                "Invalid characters (e.g. space, non-ascii) for SPAdes in path: " + os.path.realpath(fq_file)

        log_handler = timed_log(log_handler, options.output_base, options.prefix + "get_org.")
        if options.word_size is None:
            pass
        elif 0 < options.word_size < 1:
            pass
        elif options.word_size >= GLOBAL_MIN_WS:
            options.word_size = int(options.word_size)
        else:
            log_handler.error("Illegal word size (\"-w\") value!")
            exit()

        if options.pregroup_word_size:
            if 0 < options.pregroup_word_size < 1:
                pass
            elif options.pregroup_word_size >= GLOBAL_MIN_WS:
                options.pregroup_word_size = int(options.pregroup_word_size)
            else:
                log_handler.error("Illegal word size (\"--pre-w\") value!")
                exit()

        if options.fast_strategy:
            if "-R" not in sys.argv and "--max-rounds" not in sys.argv:
                options.max_rounds = 10
            if "-t" not in sys.argv:
                options.threads = 4
            if "-J" not in sys.argv:
                options.jump_step = 5
            if "-M" not in sys.argv:
                options.mesh_size = 7
            if "--max-n-words" not in sys.argv:
                options.maximum_n_words = 3E7
            options.larger_auto_ws = True
            if "--disentangle-time-limit" not in sys.argv:
                options.disentangle_time_limit = 360

        if options.memory_save:
            if "-P" not in sys.argv:
                options.pre_grouped = 0
            if "--remove-duplicates" not in sys.argv:
                options.rm_duplicates = 0

        if options.memory_unlimited:
            if "-P" not in sys.argv:
                options.pre_grouped = 1E7
            if "--remove-duplicates" not in sys.argv:
                options.rm_duplicates = 2E8
            if "--min-quality-score" not in sys.argv:
                options.min_quality_score = -5
            if "--max-ignore-percent" not in sys.argv:
                options.maximum_ignore_percent = 0

        # using the default
        if "--max-reads" not in sys.argv:
            if "embplant_mt" in options.organelle_type or "anonym" in options.organelle_type:
                options.maximum_n_reads *= 5
            elif "animal_mt" in options.organelle_type:
                options.maximum_n_reads *= 20
        if options.maximum_n_reads > READ_LINE_TO_INF:
            options.maximum_n_reads = inf
        else:
            options.maximum_n_reads = int(options.maximum_n_reads)

        if "--max-n-words" not in sys.argv:
            if "embplant_mt" in options.organelle_type or "anonym" in options.organelle_type:
                options.maximum_n_words *= 5
            elif "embplant_nr" in options.organelle_type or "fungus_mt" in options.organelle_type or\
                    "fungus_nr" in options.organelle_type:
                options.maximum_n_words /= 2
            elif "animal_mt" in options.organelle_type:
                options.maximum_n_words /= 2
        if "--genes" not in sys.argv:
            options.genes_fasta = []  #  None] * organelle_type_len
        else:
            temp_val_len = len(str(options.genes_fasta).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --genes is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            temp_vals = []
            for sub_genes in str(options.genes_fasta).split(","):
                # if sub_genes == "":
                #     temp_vals.append(sub_genes)
                if not os.path.exists(sub_genes):
                    log_handler.error(sub_genes + " not found!")
                    exit()
                else:
                    temp_vals.append(sub_genes)
            options.genes_fasta = temp_vals
        if "--ex-genes" not in sys.argv:
            options.exclude_genes = []
        else:
            temp_vals = []
            for sub_genes in str(options.exclude_genes).split(","):
                if not (os.path.exists(sub_genes) or os.path.exists(remove_db_postfix(sub_genes) + ".nhr")):
                    log_handler.error(sub_genes + " not found!")
                    exit()
                else:
                    temp_vals.append(sub_genes)
            options.exclude_genes = temp_vals

        if "--target-genome-size" not in sys.argv:
            raw_default_value = int(str(options.target_genome_size))
            options.target_genome_size = []
            for go_t, sub_organelle_t in enumerate(options.organelle_type):
                if sub_organelle_t == "embplant_mt":
                    options.target_genome_size.append(int(raw_default_value * 3))
                elif sub_organelle_t == "fungus_mt":
                    options.target_genome_size.append(int(raw_default_value / 2))
                elif sub_organelle_t in ("embplant_nr", "animal_mt", "fungus_nr"):
                    options.target_genome_size.append(int(raw_default_value / 10))
                elif sub_organelle_t == "anonym":
                    ref_seqs = read_fasta(options.genes_fasta[go_t])[1]
                    options.target_genome_size.append(2 * sum([len(this_seq) for this_seq in ref_seqs]))
                    log_handler.info(
                        "Setting '--target-genome-size " + ",".join([str(t_s) for t_s in options.target_genome_size]) +
                        "' for estimating the word size value for anonym type.")
                else:
                    options.target_genome_size.append(raw_default_value)
        else:
            temp_val_len = len(str(options.target_genome_size).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --target-genome-size is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            try:
                options.target_genome_size = [int(sub_size) for sub_size in str(options.target_genome_size).split(",")]
            except ValueError:
                log_handler.error("Invalid --target-genome-size value(s): " + str(options.target_genome_size))
                exit()

        if "--expected-max-size" not in sys.argv:
            raw_default_value = int(str(options.expected_max_size))
            options.expected_max_size = []
            for got_t, sub_organelle_t in enumerate(options.organelle_type):
                if sub_organelle_t == "embplant_pt":
                    options.expected_max_size.append(raw_default_value)
                elif sub_organelle_t in ("embplant_mt", "other_pt"):
                    options.expected_max_size.append(int(raw_default_value * 4))
                elif sub_organelle_t == "fungus_mt":
                    options.expected_max_size.append(raw_default_value)
                elif sub_organelle_t in ("embplant_nr", "fungus_nr", "animal_mt"):
                    options.expected_max_size.append(int(raw_default_value / 10))
                elif sub_organelle_t == "anonym":
                    ref_seqs = read_fasta(options.genes_fasta[got_t])[1]
                    options.expected_max_size.append(10 * sum([len(this_seq) for this_seq in ref_seqs]))
                    log_handler.info(
                        "Setting '--expected-max-size " + ",".join([str(t_s) for t_s in options.expected_max_size]) +
                        "' for estimating the word size value for anonym type.")
        else:
            temp_val_len = len(str(options.expected_max_size).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --expected-max-size is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            try:
                options.expected_max_size = [int(sub_size) for sub_size in str(options.expected_max_size).split(",")]
            except ValueError:
                log_handler.error("Invalid --expected-max-size value(s): " + str(options.expected_max_size))
                exit()

        if "--expected-min-size" not in sys.argv:
            raw_default_value = int(str(options.expected_min_size))
            options.expected_min_size = []
            for sub_organelle_t in options.organelle_type:
                options.expected_min_size.append(raw_default_value)
        else:
            temp_val_len = len(str(options.expected_min_size).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --expected-min-size is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            try:
                options.expected_min_size = [int(sub_size) for sub_size in str(options.expected_min_size).split(",")]
            except ValueError:
                log_handler.error("Invalid --expected-min-size value(s): " + str(options.expected_min_size))
                exit()

        if "--max-extending-len" not in sys.argv:
            options.max_extending_len = []  # -1 means auto
            for go_t, seed_f in enumerate(options.seed_file):
                # using auto as the default when using default seed files
                # if os.path.realpath(seed_f) == os.path.join(_SEQ_DB_PATH, options.organelle_type[go_t] + ".fasta"):
                #     options.max_extending_len.append(-1)
                # else:
                #     options.max_extending_len.append(inf)
                options.max_extending_len.append(inf)
        else:
            temp_val_len = len(str(options.max_extending_len).split(","))
            if temp_val_len != organelle_type_len:
                log_handler.error("-F is followed with " + str(organelle_type_len) + " organelle types, " +
                                  "while --max-extending-len is followed with " + str(temp_val_len) + " value(s)!\n")
                exit()
            try:
                options.max_extending_len = [-1 if sub_size == "auto" else float(sub_size)
                                             for sub_size in str(options.max_extending_len).split(",")]
            except ValueError:
                log_handler.error("Invalid --max-extending-len value(s): " + str(options.max_extending_len))
                exit()

        for sub_organelle_t in options.organelle_type:
            if sub_organelle_t in ("fungus_mt", "animal_mt", "anonym"):
                global MAX_RATIO_RL_WS
                MAX_RATIO_RL_WS = 0.8
                break

        if not executable(os.path.join(options.which_bowtie2, "bowtie2")):
            log_handler.error(os.path.join(options.which_bowtie2, "bowtie2") + " not accessible!")
            exit()
        if not executable(os.path.join(options.which_bowtie2, "bowtie2-build") + " --large-index"):
            log_handler.error(os.path.join(options.which_bowtie2, "bowtie2-build") + " not accessible!")
            exit()
        # if not executable(os.path.join(options.which_bowtie2, "bowtie2-build-l")):
        #     log_handler.error(os.path.join(options.which_bowtie2, "bowtie2-build-l") + " not accessible!")
        #     exit()
        run_slim = False
        run_disentangle = False
        if options.run_spades:
            if options.which_spades:
                if not executable(os.path.join(options.which_spades, "spades.py -h")):
                    raise Exception("spades.py not found/executable in " + options.which_spades + "!")
                else:
                    run_slim = True
                    run_disentangle = True
            else:
                options.which_spades = ""
                if not executable("spades.py -h"):
                    log_handler.error("spades.py not found in the PATH. "
                                        "Adding SPAdes binary dir to the PATH or using \"--which-spades\" to fix this. "
                                        "Now only get the reads and skip assembly.")
                    options.run_spades = False
                else:
                    run_slim = True
                    run_disentangle = True
            if not executable(os.path.join(options.which_blast, "blastn")):
                log_handler.error(os.path.join(options.which_blast, "blastn") +
                                    " not accessible! Slimming/Disentangling disabled!!\n")
                run_slim = False
                run_disentangle = False
            if options.genes_fasta and not executable(os.path.join(options.which_blast, "makeblastdb")):
                log_handler.error(os.path.join(options.which_blast, "makeblastdb") +
                                    " not accessible! Slimming/Disentangling disabled!!\n")
                run_slim = False
                run_disentangle = False
            if lib_not_available:
                log_handler.error("/".join(lib_not_available) + " not available! Disentangling disabled!!\n")
                run_disentangle = False
        options.rm_duplicates = int(options.rm_duplicates)
        options.pre_grouped = int(options.pre_grouped)
        if not options.rm_duplicates and options.pre_grouped:
            log_handler.warning("removing duplicates was inactive, so that the pre-grouping was disabled.")
            options.pre_grouped = False
        if options.max_rounds and options.max_rounds < 1:
            log_handler.warning("illegal maximum rounds! Set to infinite")
            options.max_rounds = inf
        if not options.max_rounds:
            if not options.pre_grouped:
                options.max_rounds = inf
            else:
                options.max_rounds = 1
            for sub_organelle_t in options.organelle_type:
                if sub_organelle_t in {"embplant_mt", "other_pt"}:
                    options.max_rounds = max(options.max_rounds, 30)
                elif sub_organelle_t in {"embplant_nr", "animal_mt", "fungus_mt", "fungus_nr"}:
                    options.max_rounds = max(options.max_rounds, 10)
                elif sub_organelle_t == "embplant_pt":
                    options.max_rounds = max(options.max_rounds, 15)
        random.seed(options.random_seed)
        try:
            import numpy as np
        except ImportError:
            pass
        else:
            np.random.seed(options.random_seed)
        return options, log_handler, previous_attributes, run_slim, run_disentangle


def estimate_maximum_n_reads_using_mapping(
        twice_max_coverage, check_dir, original_fq_list, reads_paired,
        maximum_n_reads_hard_bound, seed_files, organelle_types, in_customs, ex_customs, target_genome_sizes,
        keep_temp, resume, other_spades_opts,
        which_blast, which_spades, which_bowtie2, threads, random_seed, verbose_log, log_handler):
    from GetOrganelleLib.sam_parser import MapRecords, get_cover_range
    if executable(os.path.join(UTILITY_PATH, "slim_graph.py -h")):
        which_slim = UTILITY_PATH
    elif executable(os.path.join(PATH_OF_THIS_SCRIPT, "slim_graph.py -h")):
        which_slim = PATH_OF_THIS_SCRIPT
    elif executable("slim_graph.py -h"):
        which_slim = ""
    else:
        which_slim = None
    result_n_reads = [maximum_n_reads_hard_bound] * len(original_fq_list)
    data_maximum_n_reads = inf
    if not os.path.exists(check_dir):
        os.mkdir(check_dir)
    check_num_line = 100000
    increase_checking_reads_by = 5
    min_valid_cov_to_estimate = 5.0
    maximum_percent_worth_estimating = 0.1
    previous_file_sizes = [0] * len(original_fq_list)
    no_more_new_reads = [False] * len(original_fq_list)
    estimated_maximum_n_reads_list = [inf] * len(original_fq_list)
    original_fq_sizes = [os.path.getsize(raw_fq) * GUESSING_FQ_GZIP_COMPRESSING_RATIO
                         if raw_fq.endswith(".gz") else os.path.getsize(raw_fq)
                         for raw_fq in original_fq_list]
    # make paired equal size estimation if compressed
    if reads_paired and original_fq_list[0].endswith(".gz") and original_fq_list[1].endswith(".gz") and \
            abs(log(float(original_fq_sizes[0])/original_fq_sizes[1])) < log(1.3):
        original_fq_sizes[0] = original_fq_sizes[1] = (original_fq_sizes[0] + original_fq_sizes[1]) /2.
    # if the original data sizes is too small, no need to reduce
    max_organelle_base_percent = 0.2
    for go_t, organelle_type in enumerate(organelle_types):
        # temporary treat: compatible with previous
        if organelle_type in ORGANELLE_EXPECTED_GRAPH_SIZES:
            min_file_size = ORGANELLE_EXPECTED_GRAPH_SIZES[organelle_type] * twice_max_coverage \
                            / max_organelle_base_percent * GUESSING_FQ_SEQ_INFLATE_TO_FILE
        else:
            min_file_size = target_genome_sizes[go_t] * twice_max_coverage \
                            / max_organelle_base_percent * GUESSING_FQ_SEQ_INFLATE_TO_FILE
        if sum(original_fq_sizes) < min_file_size:
            if not keep_temp:
                try:
                    shutil.rmtree(check_dir)
                except OSError:
                    log_handler.warning("Removing temporary directory " + check_dir + " failed.")
            return result_n_reads
    #
    count_round = 1
    while count_round == 1 or check_num_line < min(maximum_n_reads_hard_bound, data_maximum_n_reads):
        if check_num_line > READ_LINE_TO_INF:
            return [inf] * len(original_fq_list)
        log_handler.info("Tasting " + "+".join([str(check_num_line)] * len(original_fq_list)) + " reads ...")
        this_check_dir = os.path.join(check_dir, str(count_round))
        if not os.path.exists(this_check_dir):
            os.mkdir(this_check_dir)
        check_fq_files = []
        check_percents = []
        for f_id, r_file in enumerate(original_fq_list):
            check_fq = os.path.join(this_check_dir, "check_" + str(f_id + 1))
            if not (os.path.exists(check_fq) and resume):
                if r_file.endswith(".gz"):
                    unzip(r_file, check_fq, 4 * check_num_line, verbose_log, log_handler if verbose_log else None)
                else:
                    os.system("head -n " + str(int(4 * check_num_line)) + " " + r_file + " > " + check_fq + ".temp")
                    os.rename(check_fq + ".temp", check_fq)
            check_f_size = os.path.getsize(check_fq)
            if check_f_size == 0:
                raise ValueError("Empty file" + check_fq + "\n"
                                 "Please check the legality and integrity of your input reads!\n")
            if check_f_size == previous_file_sizes[f_id]:
                no_more_new_reads[f_id] = True
                check_percents.append(1)
                tmp_line = 0
                with open(check_fq) as counter:
                    for foo in counter:
                        tmp_line += 1
                estimated_maximum_n_reads_list[f_id] = int(tmp_line / 4)
            else:
                check_percents.append(min(float(check_f_size) / original_fq_sizes[f_id], 1))
                estimated_maximum_n_reads_list[f_id] = int(check_num_line / check_percents[-1])
            check_fq_files.append(check_fq)
        count_round += 1
        data_maximum_n_reads = max(estimated_maximum_n_reads_list)

        go_next_run = False
        base_cov_of_all_organelles = []
        for go_t, seed_f in enumerate(seed_files):
            organelle_type = organelle_types[go_t]
            if sum([os.path.exists(remove_db_postfix(seed_f) + ".index" + postfix)
                    for postfix in
                    (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")]) != 6:
                new_seed_file = os.path.join(this_check_dir, os.path.basename(seed_f))
                check_fasta_seq_names(seed_f, new_seed_file)
                seed_f = new_seed_file
            bowtie_out_base = os.path.join(this_check_dir, organelle_type + ".check")
            mapped_fq = bowtie_out_base + ".fq"
            mapped_sam = bowtie_out_base + ".sam"
            map_with_bowtie2(
                seed_file=seed_f, original_fq_files=check_fq_files, bowtie_out=bowtie_out_base, resume=resume,
                threads=threads, random_seed=random_seed, generate_fq=True, target_echo_name="seed",
                log_handler=log_handler if verbose_log else None, verbose_log=verbose_log, silent=not verbose_log,
                which_bowtie2=which_bowtie2)
            seed_fq_size = os.path.getsize(mapped_fq)
            if not seed_fq_size:
                if sum(no_more_new_reads) == len(no_more_new_reads):
                    if log_handler:
                        log_handler.error("No " + str(organelle_type) + " seed reads found!")
                        log_handler.error("Please check your raw data or change your " + str(organelle_type) + " seed!")
                else:
                    data_size_checked = [check_percents[go_f] * fq_size
                                         for go_f, fq_size in enumerate(original_fq_sizes)]
                    data_checked_percent = sum(data_size_checked) / float(sum(original_fq_sizes))
                    if data_checked_percent > maximum_percent_worth_estimating:
                        base_cov_of_all_organelles.append(0.)
                        break
                    else:
                        # another run with more reads
                        go_next_run = True
                        break
            mapping_records = MapRecords(mapped_sam)
            mapping_records.update_coverages()
            coverage_info = mapping_records.coverages
            coverages_2 = [pos for ref in coverage_info for pos in coverage_info[ref] if pos > 0]
            base_cov_values = get_cover_range(coverages_2, guessing_percent=BASE_COV_SAMPLING_PERCENT)
            mean_read_len, max_read_len, all_read_nums = \
                get_read_len_mean_max_count(mapped_fq, maximum_n_reads_hard_bound, n_process=1)
                # get_read_len_mean_max_count(mapped_fq, maximum_n_reads_hard_bound, n_process=threads)
            if executable(os.path.join(which_spades, "spades.py -h")) and \
                    executable(os.path.join(which_bowtie2, "bowtie2")):
                try:
                    this_in = "" if not in_customs else in_customs[go_t]
                    this_ex = "" if not ex_customs else ex_customs[go_t]
                    base_cov_values = pre_assembly_mapped_reads_for_base_cov(
                        original_fq_files=check_fq_files, mapped_fq_file=mapped_fq, seed_fs_file=seed_f,
                        mean_read_len=mean_read_len, organelle_type=organelle_type,
                        in_custom=this_in, ex_custom=this_ex, threads=threads, resume=resume,
                        other_spades_opts=other_spades_opts,
                        which_spades=which_spades, which_slim=which_slim, which_blast=which_blast,
                        log_handler=log_handler if verbose_log else None, verbose_log=verbose_log)
                except NotImplementedError:
                    pass
            if base_cov_values[1] < min_valid_cov_to_estimate:
                data_size_checked = [check_percents[go_f] * fq_size
                                     for go_f, fq_size in enumerate(original_fq_sizes)]
                data_checked_percent = sum(data_size_checked) / float(sum(original_fq_sizes))
                if data_checked_percent > maximum_percent_worth_estimating:
                    base_cov_of_all_organelles.append(0.)
                    break
                else:
                    # another run with more reads
                    go_next_run = True
                    break
            else:
                base_cov_of_all_organelles.append(base_cov_values[1])
        if go_next_run:
            check_num_line *= increase_checking_reads_by
            continue
        data_all_size = sum(original_fq_sizes)
        data_size_checked = [check_percents[go_f] * fq_size for go_f, fq_size in enumerate(original_fq_sizes)]
        data_checked_percent = sum(data_size_checked) / float(data_all_size)
        the_check_base_cov = min(base_cov_of_all_organelles)
        the_real_base_cov = the_check_base_cov / data_checked_percent
        if the_real_base_cov > twice_max_coverage:
            reduce_ratio = twice_max_coverage / the_real_base_cov
            result_n_reads = [min(maximum_n_reads_hard_bound, math.ceil(real_num * reduce_ratio))
                              if real_num * reduce_ratio <= READ_LINE_TO_INF else inf
                              for real_num in estimated_maximum_n_reads_list]
        else:
            result_n_reads = [maximum_n_reads_hard_bound] * len(original_fq_list)
        break
    if not keep_temp:
        try:
            shutil.rmtree(check_dir)
        except OSError:
            log_handler.warning("Removing temporary directory " + check_dir + " failed.")
    return result_n_reads


def combination_res_log(all_choices_num, chosen_num):
    res = 0.
    for ch_n in range(chosen_num, 0, -1):
        res += log(all_choices_num - ch_n + 1) - log(ch_n)
    return res


def trans_word_cov(word_cov, base_cov, mean_base_error_rate, read_length):
    if mean_base_error_rate == 0.:
        return word_cov
    wrong_words_percent = 0
    for error_site_num in range(1, int(min(read_length * mean_base_error_rate * 10, read_length))):
        prob_of_err_site_num = combination_res_log(read_length, error_site_num) \
                               + error_site_num * log(mean_base_error_rate) \
                               + (read_length - error_site_num) * log(1 - mean_base_error_rate)
        wrong_words_percent += (1 - 2 ** (-error_site_num)) * exp(prob_of_err_site_num)
    # if word size < read_len/2, wrong words percent decreases
    increase_word_cov = word_cov / (1 - wrong_words_percent) - word_cov
    if word_cov > 0.5 * base_cov:
        word_cov += increase_word_cov ** 0.34
    elif word_cov + increase_word_cov > 0.5 * base_cov:
        word_cov = 0.5 * base_cov + (word_cov + increase_word_cov - 0.5 * base_cov) ** 0.34
    else:
        word_cov += increase_word_cov
    return word_cov


def estimate_word_size(base_cov, base_cov_deviation, read_length, target_size, mean_error_rate=0.015, log_handler=None,
                       max_discontinuous_prob=0.01, min_word_size=AUTO_MIN_WS, max_effective_word_cov=60,
                       wc_bc_ratio_constant=0.35, organelle_type=""):
    # base_cov_deviation cannot be well estimated and thus excluded from the estimation
    echo_problem = False
    # G: genome size, N: Number of reads from data, L: read length,
    #    ## Poisson distribution
    #    mean read cov = N/(G-L+1)
    #    expected # reads starting within any specific interval of C consecutive nucleotides = (N/(G-L+1))*C
    #    P(no read starts in the interval) = e^(-C*N/(G-L+1))
    #    P(>=1 reads start in the interval) = 1-e^(-C*N/(G-L+1))
    #    P(the interval is not continuous) = 1-(1-e^(-N/(G-L+1)))^C
    #
    # 1. The higher the base coverage is, the larger the word size should be. # to exclude unnecessary contigs.
    # 2. The longer the read length is, the larger the word size should be
    # 3. The higher the error rate is, the smaller the word size should be

    # empirical functions:
    word_cov = min(max_effective_word_cov, base_cov * wc_bc_ratio_constant)
    # min_word_cov = log(-1/(max_discontinuous_prob**(1/target_size) - 1))
    min_word_cov = 5
    while 1 - (1 - math.e ** (-min_word_cov)) ** target_size > max_discontinuous_prob:
        min_word_cov += 0.05
    # print(min_word_cov)
    #
    wc_bc_ratio_max = 1 - (min_word_size - 1) / read_length
    if base_cov * wc_bc_ratio_max < min_word_cov:
        min_word_cov = base_cov * wc_bc_ratio_max
        echo_problem = True
    word_cov = max(min_word_cov, word_cov)
    word_cov = trans_word_cov(word_cov, base_cov, mean_error_rate / 2., read_length)
    # 1. relationship between kmer coverage and base coverage, k_cov = base_cov * (read_len - k_len + 1) / read_len
    estimated_word_size = int(read_length * (1 - word_cov / base_cov)) + 1
    # print(estimated_word_size)
    estimated_word_size = min(int(read_length * MAX_RATIO_RL_WS), max(min_word_size, estimated_word_size))
    if echo_problem:
        if log_handler:
            log_handler.warning("Guessing that you are using too few data for assembling " + organelle_type + "!")
            log_handler.warning("GetOrganelle is still trying ...")
        else:
            sys.stdout.write("Guessing that you are using too few data for assembling " + organelle_type + "!\n")
            sys.stdout.write("GetOrganelle is still trying ...\n")
    return int(round(estimated_word_size, 0))


def calculate_word_size_according_to_ratio(word_size_ratio, mean_read_len, log_handler):
    if word_size_ratio < 1:
        new_word_size = int(round(word_size_ratio * mean_read_len, 0))
        if new_word_size < GLOBAL_MIN_WS:
            new_word_size = GLOBAL_MIN_WS
            log_handler.warning("Too small ratio " + str(new_word_size) + ", setting '-w " + str(GLOBAL_MIN_WS) + "'")
        else:
            log_handler.info("Setting '-w " + str(new_word_size) + "'")
        return new_word_size
    else:
        max_ws = int(round(mean_read_len * 0.9))
        if word_size_ratio > max_ws:
            word_size_ratio = max_ws
            log_handler.warning("Too large word size for mean read length " + str(mean_read_len) +
                                ", setting '-w " + str(word_size_ratio) + "'")
        return word_size_ratio


def extend_with_constant_words(baits_pool, raw_fq_files, word_size, output, jump_step=3):
    output_handler = open(output + ".Temp", "w")
    for fq_f in raw_fq_files:
        with open(fq_f) as fq_f_input_handler:
            head_line = fq_f_input_handler.readline()
            while head_line:
                seq_line = fq_f_input_handler.readline()
                seq_len = len(seq_line.strip())
                accepted = False
                for i in range(0, seq_len, jump_step):
                    if seq_line[i:i + word_size] in baits_pool:
                        accepted = True
                        break
                if accepted:
                    output_handler.write(head_line)
                    output_handler.write(seq_line)
                    output_handler.write(fq_f_input_handler.readline())
                    output_handler.write(fq_f_input_handler.readline())
                else:
                    fq_f_input_handler.readline()
                    fq_f_input_handler.readline()
                head_line = fq_f_input_handler.readline()
    output_handler.close()
    os.rename(output + ".Temp", output)


def pre_assembly_mapped_reads_for_base_cov(
        original_fq_files, mapped_fq_file, seed_fs_file, mean_read_len, organelle_type, in_custom, ex_custom,
        threads, resume, other_spades_opts, which_spades, which_slim, which_blast,
        log_handler=None, verbose_log=False, keep_temp=False):
    from GetOrganelleLib.assembly_parser import get_graph_coverages_range_simple
    draft_kmer = min(45, int(mean_read_len / 2) * 2 - 3)
    this_modified_dir = os.path.realpath(mapped_fq_file) + ".spades"
    this_original_graph = os.path.join(this_modified_dir, "assembly_graph.fastg")
    this_modified_base = "assembly_graph.fastg.modified"
    this_modified_graph = this_original_graph + ".modified.fastg"
    more_fq_file = os.path.realpath(mapped_fq_file) + ".more.fq"
    more_modified_dir = more_fq_file + ".spades"
    more_original_graph = os.path.join(more_modified_dir, "assembly_graph.fastg")
    more_modified_base = "assembly_graph.fastg.modified"
    more_modified_graph = more_original_graph + ".modified.fastg"

    if in_custom or ex_custom:
        include_priority_db = in_custom
        exclude_db = ex_custom
    else:
        include_priority_db = os.path.join(_LBL_DB_PATH, organelle_type + ".fasta")
        exclude_db = ""
    db_command = ""
    if include_priority_db:
        db_command += " --include-priority " + include_priority_db
    if exclude_db:
        db_command += " --exclude " + exclude_db

    if resume and (os.path.exists(this_modified_graph) or os.path.exists(more_modified_graph)):
        if os.path.exists(more_modified_graph) and os.path.getsize(more_modified_graph) > 0:
            kmer_cov_values = get_graph_coverages_range_simple(read_fasta(more_modified_graph))
            base_cov_values = [this_word_cov * mean_read_len / (mean_read_len - draft_kmer + 1)
                               for this_word_cov in kmer_cov_values]
        elif os.path.exists(this_modified_graph) and os.path.getsize(this_modified_graph) > 0:
            kmer_cov_values = get_graph_coverages_range_simple(read_fasta(this_modified_graph))
            base_cov_values = [this_word_cov * mean_read_len / (mean_read_len - draft_kmer + 1)
                               for this_word_cov in kmer_cov_values]
        else:
            base_cov_values = [0.0, 0.0, 0.0]
    else:
        try:
            # log_handler.info(" ...")
            this_command = os.path.join(which_spades, "spades.py") + " -t " + str(threads) + \
                           " -s " + mapped_fq_file + " " + other_spades_opts + \
                           " -k " + str(draft_kmer) + " --only-assembler -o " + this_modified_dir
            pre_assembly = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if verbose_log and log_handler:
                log_handler.info(this_command)
            output = monitor_spades_log(pre_assembly, log_handler, sensitively_stop=True, silent=True)
            if not os.path.exists(this_original_graph) or os.path.getsize(this_original_graph) == 0:
                raise OSError("original graph")
            if "== Error ==" in output:
                if verbose_log and log_handler:
                    log_handler.error('\n' + output.strip())
                raise NotImplementedError
            if which_slim is None:
                shutil.copy(this_original_graph, this_modified_graph)
            else:
                which_bl_str = " --which-blast " + which_blast if which_blast else ""
                slim_command = os.path.join(which_slim, "slim_graph.py") + \
                               " --verbose " * int(bool(verbose_log)) + which_bl_str + \
                               " --log -t " + str(threads) + " --wrapper " + this_original_graph + \
                               " -o " + this_modified_dir + " --out-base " + this_modified_base + \
                               " " + db_command + " --keep-temp " * int(bool(keep_temp))
                do_slim = subprocess.Popen(slim_command,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                if verbose_log and log_handler:
                    log_handler.info(slim_command)
                output, err = do_slim.communicate()
                if not os.path.exists(this_modified_graph):
                    if log_handler:
                        log_handler.error("slimming the pre-assembled graph failed.")
                    if verbose_log and log_handler:
                        log_handler.error("\n" + output.decode("utf8").strip())
                    shutil.copy(this_original_graph, this_modified_graph)
                elif os.path.getsize(this_modified_graph) == 0:
                    raise OSError("modified graph")
            kmer_cov_values = get_graph_coverages_range_simple(read_fasta(this_modified_graph))
            base_cov_values = [this_word_cov * mean_read_len / (mean_read_len - draft_kmer + 1)
                               for this_word_cov in kmer_cov_values]
        except OSError:
            # if os.path.exists(mapped_fq_file + ".spades"):
            #     shutil.rmtree(mapped_fq_file + ".spades")
            # using words to recruit more reads for word size estimation
            # gathering_word_size = min(auto_min_word_size, 2 * int(mean_read_len * auto_min_word_size/100.) - 1)
            if log_handler:
                log_handler.info("Retrying with more reads ..")
            gathering_word_size = 25
            if resume and os.path.exists(more_fq_file):
                pass
            else:
                theses_words = chop_seqs(
                    fq_simple_generator(mapped_fq_file), word_size=gathering_word_size)
                theses_words |= chop_seqs(
                    read_fasta(seed_fs_file)[1], word_size=gathering_word_size)
                extend_with_constant_words(
                    theses_words, original_fq_files, word_size=gathering_word_size, output=more_fq_file)
            more_command = os.path.join(which_spades, "spades.py") + " -t " + str(threads) + " -s " + \
                           more_fq_file + " " + other_spades_opts + " -k " + str(draft_kmer) + \
                           " --only-assembler -o " + this_modified_dir
            pre_assembly = subprocess.Popen(
                more_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            if verbose_log and log_handler:
                log_handler.info(more_command)
            output = monitor_spades_log(pre_assembly, log_handler, sensitively_stop=True)
            if not os.path.exists(more_original_graph) or os.path.getsize(more_original_graph) == 0:
                if verbose_log and log_handler:
                    log_handler.error(more_original_graph + " not found/valid!")
                raise NotImplementedError
            elif "== Error ==" in output:
                if verbose_log and log_handler:
                    log_handler.error('\n' + output.strip())
                raise NotImplementedError
            else:
                if which_slim is None or not executable(os.path.join(which_blast, "blastn")):
                    shutil.copy(more_original_graph, more_modified_graph)
                else:
                    which_bl_str = " --which-blast " + which_blast if which_blast else ""
                    slim_command = os.path.join(which_slim, "slim_graph.py") + \
                                   " --verbose " * int(bool(verbose_log)) + which_bl_str + \
                                   " --log -t " + str(threads) + " --wrapper " + more_original_graph + \
                                   " -o " + more_modified_dir + " --out-base " + more_modified_base + \
                                   " " + db_command + " --keep-temp " * int(bool(keep_temp))
                    do_slim = subprocess.Popen(slim_command,
                                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                    if verbose_log and log_handler:
                        log_handler.info(slim_command)
                    output, err = do_slim.communicate()
                    if not os.path.exists(more_modified_graph):
                        if log_handler:
                            log_handler.error("slimming the pre-assembled graph failed.")
                        if verbose_log and log_handler:
                            log_handler.error("\n" + output.decode("utf8").strip())
                        shutil.copy(more_original_graph, more_modified_graph)
                    elif os.path.getsize(more_modified_graph) == 0:
                        if log_handler:
                            log_handler.warning("No target found in the pre-assembled graph/seed. "
                                                "GetOrganelle is still trying ..")
                        shutil.copy(more_original_graph, more_modified_graph)
                kmer_cov_values = get_graph_coverages_range_simple(read_fasta(more_modified_graph))
                base_cov_values = [this_word_cov * mean_read_len / (mean_read_len - draft_kmer + 1)
                                   for this_word_cov in kmer_cov_values]
    if not keep_temp and os.path.exists(this_modified_dir):
        shutil.rmtree(this_modified_dir)
    return base_cov_values


def check_parameters(word_size, original_fq_files, seed_fs_files, seed_fq_files, seed_sam_files,
                     organelle_types, in_custom_list, ex_custom_list, mean_error_rate, target_genome_sizes,
                     max_extending_len, mean_read_len, max_read_len, low_quality_pattern,
                     all_read_nums, reduce_reads_for_cov,
                     log_handler, other_spades_opts, which_spades, which_blast, which_bowtie2,
                     wc_bc_ratio_constant=0.35, larger_auto_ws=False,
                     threads=1, resume=False, random_seed=12345, verbose_log=False, zip_files=False):
    from GetOrganelleLib.sam_parser import MapRecords, get_cover_range, mapping_gap_info_from_coverage_dict
    from itertools import combinations

    if word_size is None or -1 in max_extending_len:
        log_handler.info("The automatically-estimated parameter(s) do not ensure the best choice(s).")
        log_handler.info("If the result graph is not a circular organelle genome, ")
        log_handler.info("  you could adjust the value(s) of "
                         "'-w'" +
                         "/'--max-extending-len'" * int(bool(-1 in max_extending_len)) +
                         "/'-R' for another new run.")
    auto_max_extending_len = [m_e_l == -1 for m_e_l in max_extending_len]
    if "animal_mt" in organelle_types:
        auto_min_word_size = AUTO_MIN_WS_ANIMAL_MT
    elif "embplant_mt" not in organelle_types:
        auto_min_word_size = AUTO_MIN_WS
    else:
        auto_min_word_size = AUTO_MIN_WS_PLANT_MT

    if executable(os.path.join(UTILITY_PATH, "slim_graph.py -h")):
        which_slim = UTILITY_PATH
    elif executable(os.path.join(PATH_OF_THIS_SCRIPT, "slim_graph.py -h")):
        which_slim = PATH_OF_THIS_SCRIPT
    elif executable("slim_graph.py -h"):
        which_slim = ""
    else:
        which_slim = None

    base_coverages_by_organelles = []
    for go_t, this_sam_f in enumerate(seed_sam_files):
        gathering_word_size = None
        mapping_records = MapRecords(this_sam_f)
        mapping_records.update_coverages()
        coverage_info = mapping_records.coverages
        # multiple ref ?
        coverages_2 = [pos for ref in coverage_info for pos in coverage_info[ref] if pos > 2]
        if not coverages_2:
            coverages_2 = [pos for ref in coverage_info for pos in coverage_info[ref] if pos > 0]
            if not coverages_2:
                if log_handler:
                    log_handler.error("No " + organelle_types[go_t] + " seed reads found!")
                    log_handler.error("Please check your raw data or change your " + organelle_types[go_t] + " seed!")
                exit()
        # top BASE_COV_SAMPLING_PERCENT from mapped reads
        base_cov_values = get_cover_range(coverages_2, guessing_percent=BASE_COV_SAMPLING_PERCENT)
        # log_handler.info(
        #     "Estimated " + organelle_types[go_t] + "-hitting base-coverage = " + "%.2f" % base_cov_values[1])
        #     # + "~".join(["%.2f" % base_c for base_c in base_cov_values]))

        this_modified_dir = seed_fq_files[go_t] + ".spades"
        this_modified_graph = os.path.join(this_modified_dir, "assembly_graph.fastg.modified.fastg")
        # if base_cov_values[0] < 100 and set(organelle_types) != {"embplant_pt"}:
        # if word_size is None:
        if word_size is None or max_extending_len[go_t] == -1:
            if executable(os.path.join(which_spades, "spades.py -h")) and \
                    executable(os.path.join(which_bowtie2, "bowtie2")):
                log_handler.info("Pre-assembling mapped reads ...")
                try:
                    this_in = "" if not in_custom_list else in_custom_list[go_t]
                    this_ex = "" if not ex_custom_list else ex_custom_list[go_t]
                    base_cov_values = pre_assembly_mapped_reads_for_base_cov(
                        original_fq_files=original_fq_files, mapped_fq_file=seed_fq_files[go_t],
                        seed_fs_file=seed_fs_files[go_t], mean_read_len=mean_read_len,
                        # TODO check in_customs lengths
                        organelle_type=organelle_types[go_t], in_custom=this_in, ex_custom=this_ex,
                        threads=threads, resume=resume, log_handler=log_handler, verbose_log=verbose_log,
                        other_spades_opts=other_spades_opts,
                        which_spades=which_spades, which_slim=which_slim, which_blast=which_blast)
                except NotImplementedError:
                    if max_extending_len[go_t] == -1:
                        log_handler.warning(
                            "Pre-assembling failed. The estimations for " + organelle_types[go_t] + "-hitting "
                            "base-coverage, -w, --max-extending-len may be misleading.")
                    else:
                        log_handler.warning(
                            "Pre-assembling failed. "
                            "The estimations for " + organelle_types[go_t] + "-hitting base-coverage "
                            "and word size may be misleading.")
                    pass
                else:
                    log_handler.info("Pre-assembling mapped reads finished.")
            else:
                log_handler.warning(
                    "No pre-assembling due to insufficient dependencies! "
                    "The estimations for " + organelle_types[go_t] +
                    "-hitting base-coverage and word size may be misleading.")
        base_coverages_by_organelles.append((base_cov_values[1], (base_cov_values[2] - base_cov_values[0]) / 2))
        log_handler.info(
            "Estimated " + organelle_types[go_t] + "-hitting base-coverage = " + "%.2f" % base_cov_values[1])
        #
        if executable(os.path.join(which_spades, "spades.py -h")) and \
                executable(os.path.join(which_bowtie2, "bowtie2")):
            if max_extending_len[go_t] == -1:  # auto
                best_seed, gap_percent, largest_gap_lens = mapping_gap_info_from_coverage_dict(coverage_info)
                log_handler.info("Closest " + organelle_types[go_t] + " seed sequence: " + str(best_seed))
                # redo quick-mapping with the closest seed
                if os.path.exists(this_modified_graph):
                    simulated_fq_f = os.path.join(seed_fq_files[go_t] + ".spades",
                                                  "get_org.assembly_graph.simulated.fq")
                    simulate_fq_simple(from_fasta_file=this_modified_graph,
                                       out_dir=seed_fq_files[go_t] + ".spades",
                                       out_name="get_org.assembly_graph.simulated.fq",
                                       sim_read_jump_size=7, resume=resume, random_obj=random)
                    closest_seed_f = os.path.join(seed_fq_files[go_t] + ".spades", "get_org.closest_seed.fasta")
                    seed_seq_list = SequenceList(seed_fs_files[go_t])
                    for seq_record in seed_seq_list:
                        if seq_record.label.startswith(best_seed):
                            with open(closest_seed_f + ".Temp", "w") as out_closest:
                                out_closest.write(seq_record.fasta_str() + "\n")
                            os.rename(closest_seed_f + ".Temp", closest_seed_f)
                            break
                    bowtie_out_base = os.path.join(seed_fq_files[go_t] + ".spades", "get_org.map_to_closest")
                    map_with_bowtie2(seed_file=closest_seed_f, original_fq_files=[simulated_fq_f],
                                     bowtie_out=bowtie_out_base, resume=resume, threads=threads,
                                     random_seed=random_seed, target_echo_name=organelle_types[go_t],
                                     log_handler=log_handler, generate_fq=False, silent=verbose_log,
                                     which_bowtie2=which_bowtie2, bowtie2_mode="--very-fast-local")
                    mapping_records = MapRecords(bowtie_out_base + ".sam")
                    mapping_records.update_coverages()
                    coverage_info = mapping_records.coverages
                    best_seed, gap_percent, largest_gap_lens = mapping_gap_info_from_coverage_dict(coverage_info)
                    # if not keep_temp:
                    #     os.remove(simulated_fq_f)
                    if zip_files:
                        zip_file(source=bowtie_out_base + ".sam", target=bowtie_out_base + ".sam.tar.gz",
                                 verbose_log=verbose_log, log_handler=log_handler, remove_source=True)
                        zip_file(source=simulated_fq_f, target=simulated_fq_f + ".tar.gz",
                                 verbose_log=verbose_log, log_handler=log_handler, remove_source=True)
                log_handler.info("Unmapped percentage " + "%1.4f" % gap_percent + " and unmapped lengths " +
                                 " ".join([str(g_l) for g_l in largest_gap_lens[:5]]) + " ..")
                cov_dev_percent = (base_cov_values[2] - base_cov_values[0]) / 2 / base_cov_values[1]
                # if organelle_types[go_t] == "animal_mt":
                #     # empirical function
                #     max_extending_len[go_t] = largest_gap_lens[0] / 2. * (1 + gap_percent ** 2) * (1 + cov_dev_percent ** 2)
                #     max_extending_len[go_t] = min(int(math.ceil(max_extending_len[go_t] / 100)) * 100, 15000)
                # else:
                if len(coverage_info[best_seed]) < target_genome_sizes[go_t] / 10. or gap_percent > 0.4:
                    max_extending_len[go_t] = inf
                else:
                    if largest_gap_lens:
                        # empirical function
                        max_extending_len[go_t] = largest_gap_lens[0] / 2. * (1 + gap_percent ** 0.5) \
                                                  * (1 + cov_dev_percent ** 0.5)
                    else:
                        max_extending_len[go_t] = 1
                    # if more.fq was used,
                    # previous empirical formula is not estimating the gap based on the initial mapped fq
                    # so the gap could be actually larger by 2 * (max_read_len - gathering_word_size)
                    if gathering_word_size is not None:
                        max_extending_len[go_t] += 2 * (max_read_len - gathering_word_size)
                    max_extending_len[go_t] = int(math.ceil(max_extending_len[go_t]/100)) * 100
        else:
            max_extending_len[go_t] = inf

    # check the divergence of coverages of different organelle genomes
    for estimated_a, estimated_b in combinations(base_coverages_by_organelles, 2):
        if abs(log(estimated_a[0]) - log(estimated_b[0])) > log(10):
            log_handler.warning("Multi-organelle mode (with the same data size and word size) is not suggested "
                                "for organelles with divergent base-coverages.")
            log_handler.warning("Please try to get different organelles in separate runs, "
                                "or to use other seeds to get a better estimation of coverage values.")
            break

    # check the base coverage to ensure not using too much data
    this_minimum_base_cov = min([value_set[0] for value_set in base_coverages_by_organelles])
    if this_minimum_base_cov > reduce_reads_for_cov:
        reduce_ratio = reduce_reads_for_cov / this_minimum_base_cov
        for go_r_n, read_num in enumerate(all_read_nums):
            all_read_nums[go_r_n] = int(read_num * reduce_ratio)
        log_handler.info("Reads reduced to = " + "+".join([str(sub_num) for sub_num in all_read_nums]))
        for go_t, (t_base_cov, t_base_sd) in enumerate(base_coverages_by_organelles):
            base_coverages_by_organelles[go_t] = t_base_cov * reduce_ratio, t_base_sd * reduce_ratio
            log_handler.info("Adjusting expected " + organelle_types[go_t] + " base coverage to " +
                             "%.2f" % (t_base_cov * reduce_ratio))

    if word_size is None:
        all_ws_values = []
        for go_type, (this_base_cov, cov_dev) in enumerate(base_coverages_by_organelles):
            if larger_auto_ws:
                word_size = estimate_word_size(
                    base_cov=this_base_cov, base_cov_deviation=cov_dev,
                    read_length=mean_read_len, target_size=target_genome_sizes[go_type],
                    max_discontinuous_prob=0.05, min_word_size=69, mean_error_rate=mean_error_rate,
                    log_handler=log_handler, wc_bc_ratio_constant=wc_bc_ratio_constant - 0.03,
                    organelle_type=organelle_types[go_type])
            else:
                word_size = estimate_word_size(
                    base_cov=this_base_cov, base_cov_deviation=cov_dev,
                    read_length=mean_read_len, target_size=target_genome_sizes[go_type],
                    max_discontinuous_prob=0.01, min_word_size=auto_min_word_size,
                    mean_error_rate=mean_error_rate, log_handler=log_handler,
                    wc_bc_ratio_constant=wc_bc_ratio_constant, organelle_type=organelle_types[go_type])
            all_ws_values.append(word_size)
        word_size = min(all_ws_values)
        log_handler.info("Estimated word size(s): " + ",".join([str(here_w) for here_w in all_ws_values]))
        log_handler.info("Setting '-w " + str(word_size) + "'")
    elif float(str(word_size)) < 1:
        new_word_size = int(round(word_size * mean_read_len, 0))
        if new_word_size < GLOBAL_MIN_WS:
            word_size = GLOBAL_MIN_WS
            log_handler.warning("Too small ratio " + str(word_size) + ", setting '-w " + str(GLOBAL_MIN_WS) + "'")
        else:
            word_size = new_word_size
            log_handler.info("Setting '-w " + str(word_size) + "'")

    all_infinite = True
    for go_t, max_ex_len in enumerate(max_extending_len):
        if not auto_max_extending_len[go_t]:  # not user defined
            all_infinite = False
            break
        # if organelle_types[go_t] == "animal_mt":
        #     if max_extending_len[go_t] < 15000:
        #         all_infinite = False
        #         break
        # else:
        if max_extending_len[go_t] < 6000:  # empirically not efficient for max_extending_len > 6000
            all_infinite = False
            break
    if all_infinite:
        for go_t in range(len(max_extending_len)):
            max_extending_len[go_t] = inf
    log_handler.info(
        "Setting '--max-extending-len " + ",".join([str(max_ex_l) for max_ex_l in max_extending_len]) + "'")

    if float(word_size) / max_read_len <= 0.5 and len(low_quality_pattern) > 2:
        keep_seq_parts = True
    else:
        keep_seq_parts = False

    return word_size, keep_seq_parts, base_coverages_by_organelles, max_extending_len, all_read_nums


def check_kmers(kmer_str, word_s, max_r_len, log_handler):
    if kmer_str:
        # delete illegal kmer
        try:
            kmer_values = [int(kmer_v) for kmer_v in kmer_str.split(",")]
        except ValueError:
            raise ValueError("Invalid kmer value string: " + kmer_str)
        else:
            for kmer_v_out in kmer_values:
                if kmer_v_out % 2 == 0:
                    raise ValueError("Invalid kmer value: " + str(kmer_v_out) + "! kmer values must be odd numbers!")
            # delete illegal kmer
            kmer_values = [kmer_v for kmer_v in kmer_values if 21 <= kmer_v <= min(max_r_len, 127)]

            spades_kmer = ",".join([str(kmer_v) for kmer_v in sorted(kmer_values)])
            log_handler.info("Setting '-k " + spades_kmer + "'")
            return spades_kmer
    else:
        return None


try:
    import psutil
except ImportError:
    this_process = None
else:
    this_process = psutil.Process(os.getpid())


def write_fq_results(original_fq_files, accepted_contig_id, out_file_name, temp2_clusters_dir, fq_info_in_memory,
                     all_read_limits, echo_step, verbose, index_in_memory, log_handler, extra_accepted_lines=set()):
    if verbose:
        if echo_step != inf:
            sys.stdout.write(' ' * 100 + '\b' * 100)
            sys.stdout.flush()
        log_handler.info("Producing output ...")
        log_handler.info("reading indices ...")
    accepted_lines = []
    if index_in_memory:
        # read cluster indices
        for this_index in accepted_contig_id:
            accepted_lines += fq_info_in_memory[1][this_index]
        # produce the pair-end output
        accepted_lines = set(accepted_lines)
    else:
        # read cluster indices
        temp2_indices_file_in = open(temp2_clusters_dir, 'r')
        this_index = 0
        for line in temp2_indices_file_in:
            if this_index in accepted_contig_id:
                accepted_lines += [int(x) for x in line.strip().split('\t')]
            this_index += 1
        accepted_lines = set(accepted_lines)
    # add initial mapped read ids
    for line_id in extra_accepted_lines:
        accepted_lines.add(line_id)

    # write by line
    if verbose:
        log_handler.info("writing fastq lines ...")
    post_reading = [open(fq_file, 'r') for fq_file in original_fq_files]
    files_out = [open(out_file_name + '_' + str(i + 1) + '.temp', 'w') for i in range(len(original_fq_files))]
    line_count = 0
    for i in range(len(original_fq_files)):
        count_r = 0
        line = post_reading[i].readline()
        while line:
            count_r += 1
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
            if count_r >= all_read_limits[i]:
                break
        files_out[i].close()
        post_reading[i].close()
    del accepted_lines
    for i in range(len(original_fq_files)):
        os.rename(out_file_name + '_' + str(i + 1) + '.temp', out_file_name + '_' + str(i + 1) + '.fq')
    if verbose:
        log_handler.info("writing fastq lines finished.")


def make_read_index(original_fq_files, direction_according_to_user_input, all_read_limits, rm_duplicates, output_base,
                    word_size, anti_lines, pre_grouped, index_in_memory, anti_seed, keep_seq_parts,
                    low_quality, echo_step, resume, log_handler):
    # read original reads
    # line_cluster (list) ~ forward_reverse_reads
    line_clusters = []
    seq_duplicates = {}
    forward_reverse_reads = []
    line_count = 0
    this_index = 0
    do_split_low_quality = len(low_quality) > 2
    #
    name_to_line = {}
    #
    temp1_contig_dir = [os.path.join(output_base, k + 'temp.indices.1') for k in ("_", "")]
    temp2_clusters_dir = [os.path.join(output_base, k + 'temp.indices.2') for k in ("_", "")]
    cancel_seq_parts = True
    if resume and os.path.exists(temp1_contig_dir[1]) and os.path.exists(temp2_clusters_dir[1]):
        if pre_grouped or index_in_memory:
            log_handler.info("Reading existed indices for fastq ...")
            #
            if keep_seq_parts:
                forward_reverse_reads = [x.strip().split("\t") for x in open(temp1_contig_dir[1], 'r')]
                cancel_seq_parts = True if max([len(x) for x in forward_reverse_reads]) == 1 else False
            else:
                forward_reverse_reads = [x.strip() for x in open(temp1_contig_dir[1], 'r')]
            #
            line_clusters = [[int(x) for x in y.split('\t')] for y in open(temp2_clusters_dir[1], 'r')]
            if rm_duplicates:
                line_count = sum([len(x) for x in line_clusters]) * 4
            # log
            len_indices = len(line_clusters)
            if this_process:
                memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
            else:
                memory_usage = ''
            if rm_duplicates:
                log_handler.info(memory_usage + str(len_indices) + " unique reads in all " +
                                 str(line_count // 4) + " reads")
            else:
                log_handler.info(memory_usage + str(len_indices) + " reads")
        else:
            log_handler.info("indices for fastq existed!")
            len_indices = len([x for x in open(temp2_clusters_dir[1], 'r')])
    elif resume and os.path.exists(temp1_contig_dir[1]) and not rm_duplicates:
        if index_in_memory:
            log_handler.info("Reading existed indices for fastq ...")
            #
            if keep_seq_parts:
                forward_reverse_reads = [x.strip().split("\t") for x in open(temp1_contig_dir[1], 'r')]
                cancel_seq_parts = True if max([len(x) for x in forward_reverse_reads]) == 1 else False
            else:
                forward_reverse_reads = [x.strip() for x in open(temp1_contig_dir[1], 'r')]

        # lengths = []
        use_user_direction = False
        for id_file, file_name in enumerate(original_fq_files):
            file_in = open(file_name, "r")
            count_this_read_n = 0
            line = file_in.readline()
            # if anti seed input, name & direction should be recognized
            if anti_seed:
                while line and count_this_read_n < all_read_limits[id_file]:
                    if line.startswith("@"):
                        count_this_read_n += 1
                        # parsing name & direction
                        if use_user_direction:
                            this_name = line[1:].strip()
                            direction = direction_according_to_user_input[id_file]
                        else:
                            try:
                                if ' ' in line:
                                    this_head = line[1:].split(' ')
                                    this_name, direction = this_head[0], int(this_head[1][0])
                                elif '#' in line:
                                    this_head = line[1:].split('#')
                                    this_name, direction = this_head[0], int(this_head[1].strip("/")[0])
                                elif line[-3] == "/" and line[-2].isdigit():  # 2019-04-22 added
                                    this_name, direction = line[1:-3], int(line[-2])
                                elif line[1:].strip().isdigit():
                                    log_handler.info("Using user-defined read directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                                else:
                                    log_handler.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                    log_handler.info("Using user-defined read directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                            except (ValueError, IndexError):
                                log_handler.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                log_handler.info("Using user-defined read directions. ")
                                use_user_direction = True
                                this_name = line[1:].strip()
                                direction = direction_according_to_user_input[id_file]
                        if (this_name, direction) in anti_lines:
                            line_count += 4
                            for i in range(4):
                                line = file_in.readline()
                            continue
                        this_seq = file_in.readline().strip()
                        # drop nonsense reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue
                        file_in.readline()
                        quality_str = file_in.readline()
                        if do_split_low_quality:
                            this_seq = split_seq_by_quality_pattern(this_seq, quality_str, low_quality, word_size)
                            # drop nonsense reads
                            if not this_seq:
                                line_count += 4
                                line = file_in.readline()
                                continue
                        line_clusters.append([line_count])
                    else:
                        log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if echo_step != inf and line_count % echo_step == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            else:
                while line and count_this_read_n < all_read_limits[id_file]:
                    if line.startswith("@"):
                        count_this_read_n += 1
                        this_seq = file_in.readline().strip()

                        # drop nonsense reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue

                        file_in.readline()
                        quality_str = file_in.readline()
                        if do_split_low_quality:
                            this_seq = split_seq_by_quality_pattern(this_seq, quality_str, low_quality, word_size)
                            # drop nonsense reads
                            if not this_seq:
                                line_count += 4
                                line = file_in.readline()
                                continue
                        line_clusters.append([line_count])
                    else:
                        log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if echo_step != inf and line_count % echo_step == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            line = file_in.readline()
            file_in.close()
            if line:
                log_handler.info("For " + file_name + ", only top " + str(int(all_read_limits[id_file])) +
                                 " reads are used in downstream analysis.")
        if this_process:
            memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
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
        len_indices = len(line_clusters)
        log_handler.info(memory_usage + str(len_indices) + " reads")
    elif resume and os.path.exists(temp1_contig_dir[1]):
        # lengths = []
        use_user_direction = False
        for id_file, file_name in enumerate(original_fq_files):
            file_in = open(file_name, "r")
            count_this_read_n = 0
            line = file_in.readline()
            # if anti seed input, name & direction should be recognized
            if anti_seed:
                while line and count_this_read_n < all_read_limits[id_file]:
                    if line.startswith("@"):
                        count_this_read_n += 1
                        # parsing name & direction
                        if use_user_direction:
                            this_name = line[1:].strip()
                            direction = direction_according_to_user_input[id_file]
                        else:
                            try:
                                if ' ' in line:
                                    this_head = line[1:].split(' ')
                                    this_name, direction = this_head[0], int(this_head[1][0])
                                elif '#' in line:
                                    this_head = line[1:].split('#')
                                    this_name, direction = this_head[0], int(this_head[1].strip("/")[0])
                                elif line[-3] == "/" and line[-2].isdigit():  # 2019-04-22 added
                                    this_name, direction = line[1:-3], int(line[-2])
                                elif line[1:].strip().isdigit():
                                    log_handler.info("Using user-defined read directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                                else:
                                    log_handler.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                    log_handler.info("Using user-defined read directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                            except (ValueError, IndexError):
                                log_handler.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                log_handler.info("Using user-defined read directions. ")
                                use_user_direction = True
                                this_name = line[1:].strip()
                                direction = direction_according_to_user_input[id_file]

                        if (this_name, direction) in anti_lines:
                            line_count += 4
                            for i in range(4):
                                line = file_in.readline()
                            continue
                        this_seq = file_in.readline().strip()
                        # drop nonsense reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue

                        file_in.readline()
                        quality_str = file_in.readline()
                        if do_split_low_quality:
                            this_seq = split_seq_by_quality_pattern(this_seq, quality_str, low_quality, word_size)
                            # drop nonsense reads
                            if not this_seq:
                                line_count += 4
                                line = file_in.readline()
                                continue

                            if keep_seq_parts:
                                if cancel_seq_parts and len(this_seq) > 1:
                                    cancel_seq_parts = False
                                this_c_seq = complementary_seqs(this_seq)
                                # lengths.extend([len(seq_part) for seq_part in this_seq])
                            else:
                                this_seq = this_seq[0]
                                this_c_seq = complementary_seq(this_seq)
                                # lengths.append(len(this_seq))
                        else:
                            this_c_seq = complementary_seq(this_seq)
                            # lengths.append(len(this_seq))
                        if this_seq in seq_duplicates:
                            line_clusters[seq_duplicates[this_seq]].append(line_count)
                        elif this_c_seq in seq_duplicates:
                            line_clusters[seq_duplicates[this_c_seq]].append(line_count)
                        else:
                            if index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            seq_duplicates[this_seq] = this_index
                            line_clusters.append([line_count])
                            this_index += 1
                        if len(seq_duplicates) > rm_duplicates:
                            seq_duplicates = {}
                    else:
                        log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if echo_step != inf and line_count % echo_step == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            else:
                while line and count_this_read_n < all_read_limits[id_file]:
                    if line.startswith("@"):
                        count_this_read_n += 1
                        this_seq = file_in.readline().strip()

                        # drop nonsense reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue

                        file_in.readline()
                        quality_str = file_in.readline()
                        if do_split_low_quality:
                            this_seq = split_seq_by_quality_pattern(this_seq, quality_str, low_quality, word_size)
                            # drop nonsense reads
                            if not this_seq:
                                line_count += 4
                                line = file_in.readline()
                                continue
                            if keep_seq_parts:
                                if cancel_seq_parts and len(this_seq) > 1:
                                    cancel_seq_parts = False
                                this_c_seq = complementary_seqs(this_seq)
                                # lengths.extend([len(seq_part) for seq_part in this_seq])
                            else:
                                this_seq = this_seq[0]
                                this_c_seq = complementary_seq(this_seq)
                                # lengths.append(len(this_seq))
                        else:
                            this_c_seq = complementary_seq(this_seq)
                            # lengths.append(len(this_seq))
                        if this_seq in seq_duplicates:
                            line_clusters[seq_duplicates[this_seq]].append(line_count)
                        elif this_c_seq in seq_duplicates:
                            line_clusters[seq_duplicates[this_c_seq]].append(line_count)
                        else:
                            if index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            seq_duplicates[this_seq] = this_index
                            line_clusters.append([line_count])
                            this_index += 1
                        if len(seq_duplicates) > rm_duplicates:
                            seq_duplicates = {}
                    else:
                        log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if echo_step != inf and line_count % echo_step == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            line = file_in.readline()
            file_in.close()
            if line:
                log_handler.info("For " + file_name + ", only top " + str(int(all_read_limits[id_file])) +
                                 " reads are used in downstream analysis.")
        if this_process:
            memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
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
        len_indices = len(line_clusters)
        if len_indices == 0 and line_count // 4 > 0:
            log_handler.error("No qualified reads found!")
            log_handler.error("Word size (" + str(word_size) + ") CANNOT be larger than your "
                              "post-trimmed maximum read length!")
            exit()
        log_handler.info(memory_usage + str(len_indices) + " candidates in all " + str(line_count // 4) + " reads")
    else:
        if not index_in_memory:
            temp1_contig_out = open(temp1_contig_dir[0], 'w')
        # lengths = []
        use_user_direction = False
        for id_file, file_name in enumerate(original_fq_files):
            file_in = open(file_name, "r")
            count_this_read_n = 0
            line = file_in.readline()
            # if anti seed input, name & direction should be recognized
            if anti_seed:
                while line and count_this_read_n < all_read_limits[id_file]:
                    if line.startswith("@"):
                        count_this_read_n += 1
                        # parsing name & direction
                        if use_user_direction:
                            this_name = line[1:].strip()
                            direction = direction_according_to_user_input[id_file]
                        else:
                            try:
                                if ' ' in line:
                                    this_head = line[1:].split(' ')
                                    this_name, direction = this_head[0], int(this_head[1][0])
                                elif '#' in line:
                                    this_head = line[1:].split('#')
                                    this_name, direction = this_head[0], int(this_head[1].strip("/")[0])
                                elif line[-3] == "/" and line[-2].isdigit():  # 2019-04-22 added
                                    this_name, direction = line[1:-3], int(line[-2])
                                elif line[1:].strip().isdigit():
                                    log_handler.info("Using user-defined read directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                                else:
                                    log_handler.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                    log_handler.info("Using user-defined read directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                            except (ValueError, IndexError):
                                log_handler.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                log_handler.info("Using user-defined read directions. ")
                                use_user_direction = True
                                this_name = line[1:].strip()
                                direction = direction_according_to_user_input[id_file]

                        if (this_name, direction) in anti_lines:
                            line_count += 4
                            for i in range(4):
                                line = file_in.readline()
                            continue
                        this_seq = file_in.readline().strip()
                        # drop nonsense reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue

                        file_in.readline()
                        quality_str = file_in.readline()
                        if do_split_low_quality:
                            this_seq = split_seq_by_quality_pattern(this_seq, quality_str, low_quality, word_size)
                            # drop nonsense reads
                            if not this_seq:
                                line_count += 4
                                line = file_in.readline()
                                continue

                            if keep_seq_parts:
                                if cancel_seq_parts and len(this_seq) > 1:
                                    cancel_seq_parts = False
                                this_c_seq = complementary_seqs(this_seq)
                                # lengths.extend([len(seq_part) for seq_part in this_seq])
                            else:
                                this_seq = this_seq[0]
                                this_c_seq = complementary_seq(this_seq)
                                # lengths.append(len(this_seq))
                        else:
                            this_c_seq = complementary_seq(this_seq)
                            # lengths.append(len(this_seq))
                        if rm_duplicates:
                            if this_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_seq]].append(line_count)
                            elif this_c_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_c_seq]].append(line_count)
                            else:
                                if index_in_memory:
                                    forward_reverse_reads.append(this_seq)
                                    forward_reverse_reads.append(this_c_seq)
                                else:
                                    if do_split_low_quality and keep_seq_parts:
                                        temp1_contig_out.write(
                                            "\t".join(this_seq) + '\n' + "\t".join(this_c_seq) + '\n')
                                    else:
                                        temp1_contig_out.write(this_seq + '\n' + this_c_seq + '\n')
                                seq_duplicates[this_seq] = this_index
                                line_clusters.append([line_count])
                                this_index += 1
                            if len(seq_duplicates) > rm_duplicates:
                                seq_duplicates = {}
                        else:
                            line_clusters.append([line_count])
                            if index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            else:
                                if do_split_low_quality and keep_seq_parts:
                                    temp1_contig_out.write("\t".join(this_seq) + '\n' + "\t".join(this_c_seq) + '\n')
                                else:
                                    temp1_contig_out.write(this_seq + '\n' + this_c_seq + '\n')
                    else:
                        log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if echo_step != inf and line_count % echo_step == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            else:
                while line and count_this_read_n < all_read_limits[id_file]:
                    if line.startswith("@"):
                        count_this_read_n += 1
                        this_seq = file_in.readline().strip()

                        # drop nonsense reads
                        if len(this_seq) < word_size:
                            line_count += 4
                            for i in range(3):
                                line = file_in.readline()
                            continue

                        file_in.readline()
                        quality_str = file_in.readline()
                        if do_split_low_quality:
                            this_seq = split_seq_by_quality_pattern(this_seq, quality_str, low_quality, word_size)
                            # drop nonsense reads
                            if not this_seq:
                                line_count += 4
                                line = file_in.readline()
                                continue
                            if keep_seq_parts:
                                if cancel_seq_parts and len(this_seq) > 1:
                                    cancel_seq_parts = False
                                this_c_seq = complementary_seqs(this_seq)
                                # lengths.extend([len(seq_part) for seq_part in this_seq])
                            else:
                                this_seq = this_seq[0]
                                this_c_seq = complementary_seq(this_seq)
                                # lengths.append(len(this_seq))
                        else:
                            this_c_seq = complementary_seq(this_seq)
                            # lengths.append(len(this_seq))
                        if rm_duplicates:
                            if this_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_seq]].append(line_count)
                            elif this_c_seq in seq_duplicates:
                                line_clusters[seq_duplicates[this_c_seq]].append(line_count)
                            else:
                                if index_in_memory:
                                    forward_reverse_reads.append(this_seq)
                                    forward_reverse_reads.append(this_c_seq)
                                else:
                                    if do_split_low_quality and keep_seq_parts:
                                        temp1_contig_out.write(
                                            "\t".join(this_seq) + '\n' + "\t".join(this_c_seq) + '\n')
                                    else:
                                        temp1_contig_out.write(this_seq + '\n' + this_c_seq + '\n')
                                seq_duplicates[this_seq] = this_index
                                line_clusters.append([line_count])
                                this_index += 1
                            if len(seq_duplicates) > rm_duplicates:
                                seq_duplicates = {}
                        else:
                            line_clusters.append([line_count])
                            if index_in_memory:
                                forward_reverse_reads.append(this_seq)
                                forward_reverse_reads.append(this_c_seq)
                            else:
                                if do_split_low_quality and keep_seq_parts:
                                    temp1_contig_out.write("\t".join(this_seq) + '\n' + "\t".join(this_c_seq) + '\n')
                                else:
                                    temp1_contig_out.write(this_seq + '\n' + this_c_seq + '\n')
                    else:
                        log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if echo_step != inf and line_count % echo_step == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            line = file_in.readline()
            file_in.close()
            if line:
                log_handler.info("For " + file_name + ", only top " + str(int(all_read_limits[id_file])) +
                                 " reads are used in downstream analysis.")
        if not index_in_memory:
            temp1_contig_out.close()
            os.rename(temp1_contig_dir[0], temp1_contig_dir[1])

        if this_process:
            memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
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
        len_indices = len(line_clusters)
        if rm_duplicates:
            if len_indices == 0 and line_count // 4 > 0:
                log_handler.error("No qualified reads found!")
                log_handler.error("Word size (" + str(word_size) + ") CANNOT be larger than your "
                                  "post-trimmed maximum read length!")
                exit()
            log_handler.info(memory_usage + str(len_indices) + " candidates in all " + str(line_count // 4) + " reads")
        else:
            # del lengths
            log_handler.info(memory_usage + str(len_indices) + " reads")
    if keep_seq_parts and cancel_seq_parts:
        keep_seq_parts = False
        for go_to, all_seq_parts in enumerate(forward_reverse_reads):
            forward_reverse_reads[go_to] = all_seq_parts[0]
    return forward_reverse_reads, line_clusters, len_indices, keep_seq_parts


def pre_grouping(fastq_indices_in_memory, dupli_threshold, out_base, index_in_memory, preg_word_size, log_handler):
    forward_and_reverse_reads, line_clusters, len_indices, keep_seq_parts = fastq_indices_in_memory
    log_handler.info("Pre-grouping reads ...")
    log_handler.info("Setting '--pre-w " + str(preg_word_size) + "'")
    lines_with_duplicates = {}
    count_dupli = 0
    for j in range(len(line_clusters)):
        if len(line_clusters[j]) >= 2:
            if count_dupli < dupli_threshold:
                lines_with_duplicates[j] = int
            count_dupli += 1
    if this_process:
        memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
    else:
        memory_usage = ''
    log_handler.info(memory_usage + str(len(lines_with_duplicates)) + "/" + str(count_dupli) + " used/duplicated")

    groups_of_duplicate_lines = {}
    count_groups = 0
    these_words = {}

    if index_in_memory:

        def generate_forward_and_reverse(here_unique_id):
            return forward_and_reverse_reads[2 * here_unique_id], forward_and_reverse_reads[2 * here_unique_id + 1]
    else:
        # variable outside the function
        here_go_to = [0]
        temp_seq_file = open(os.path.join(out_base, 'temp.indices.1'))
        if keep_seq_parts:
            def generate_forward_and_reverse(here_unique_id):
                forward_seq_line = temp_seq_file.readline()
                reverse_seq_line = temp_seq_file.readline()
                # skip those reads that are not unique/represented by others
                while here_go_to[0] < 2 * here_unique_id:
                    forward_seq_line = temp_seq_file.readline()
                    reverse_seq_line = temp_seq_file.readline()
                    here_go_to[0] += 2
                here_go_to[0] += 2
                return forward_seq_line.strip().split("\t"), reverse_seq_line.strip().split("\t")
        else:
            def generate_forward_and_reverse(here_unique_id):
                forward_seq_line = temp_seq_file.readline()
                reverse_seq_line = temp_seq_file.readline()
                # skip those reads that are not unique/represented by others
                while here_go_to[0] < 2 * here_unique_id:
                    forward_seq_line = temp_seq_file.readline()
                    reverse_seq_line = temp_seq_file.readline()
                    here_go_to[0] += 2
                here_go_to[0] += 2
                return forward_seq_line.strip(), reverse_seq_line.strip()

    for this_unique_read_id in sorted(lines_with_duplicates):
        this_seq, this_c_seq = generate_forward_and_reverse(this_unique_read_id)
        these_group_id = set()
        this_words = []
        if keep_seq_parts:
            for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                seq_len = len(this_seq_part)
                temp_length = seq_len - preg_word_size
                for i in range(0, temp_length + 1):
                    forward = this_seq_part[i:i + preg_word_size]
                    reverse = this_c_seq_part[temp_length - i:seq_len - i]
                    if forward in these_words:
                        these_group_id.add(these_words[forward])
                    else:
                        this_words.append(forward)
                        this_words.append(reverse)
        else:
            seq_len = len(this_seq)
            temp_length = seq_len - preg_word_size
            for i in range(0, temp_length + 1):
                forward = this_seq[i:i + preg_word_size]
                reverse = this_c_seq[temp_length - i:seq_len - i]
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
    if this_process:
        memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
    else:
        memory_usage = ''
    del these_words
    group_id_to_read_counts = {}
    for cal_copy_group_id in groups_of_duplicate_lines:
        group_id_to_read_counts[cal_copy_group_id] = sum([len(line_clusters[line_id])
                                                          for line_id in groups_of_duplicate_lines[cal_copy_group_id]])
    log_handler.info(memory_usage + str(len(groups_of_duplicate_lines)) + " groups made.")
    return groups_of_duplicate_lines, lines_with_duplicates, group_id_to_read_counts


class RoundLimitException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class WordsLimitException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NoMoreReads(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def extending_no_lim(word_size, seed_file, original_fq_files, len_indices, pre_grouped,
                     groups_of_duplicate_lines, lines_with_duplicates, fq_info_in_memory, output_base,
                     max_rounds, min_rounds, fg_out_per_round, jump_step, mesh_size, verbose, resume,
                     all_read_limits, maximum_n_words, keep_seq_parts, low_qual_pattern, echo_step,
                     log_handler):
    # adding initial word
    log_handler.info("Adding initial words ...")
    if keep_seq_parts:
        accepted_words = chop_seq_list(
            fq_simple_generator(seed_file[0], split_pattern=low_qual_pattern, min_sub_seq=word_size),
            word_size)
        for go_type in range(1, len(seed_file)):
            chop_seq_list(
                fq_simple_generator(seed_file[go_type], split_pattern=low_qual_pattern, min_sub_seq=word_size),
                word_size, previous_words=accepted_words)
    else:
        accepted_words = chop_seqs(fq_simple_generator(seed_file[0]), word_size)
        for go_type in range(1, len(seed_file)):
            chop_seqs(fq_simple_generator(seed_file[go_type]), word_size, previous_words=accepted_words)
    log_handler.info("AW " + str(len(accepted_words)))

    accepted_rd_id = set()
    accepted_rd_id_this_round = set()
    g_duplicate_lines = deepcopy(groups_of_duplicate_lines)
    l_with_duplicates = deepcopy(lines_with_duplicates)
    line_to_accept = set()
    round_count = 1
    initial_aw_count = len(accepted_words)
    prev_aw_count = initial_aw_count
    accumulated_num_words = initial_aw_count
    check_times = 1000
    check_step = max(int(len_indices / check_times), 1)
    if fg_out_per_round:
        round_dir = os.path.join(output_base, "intermediate_reads")
        if not os.path.exists(round_dir):
            os.mkdir(round_dir)
    if not this_process and verbose:
        log_handler.warning("Package psutil is not installed, so that memory usage will not be logged\n"
                            "Don't worry. This will not affect the result.")
    try:
        def summarise_round(acc_words, acc_contig_id_this_round, pre_aw, r_count, acc_num_words, unique_id):
            len_aw = len(acc_words)
            len_al = len(acc_contig_id_this_round)
            # for check words limit; memory control
            acc_num_words += len_aw - pre_aw
            if this_process:
                inside_memory_usage = " Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
            else:
                inside_memory_usage = ''
            if fg_out_per_round:
                write_fq_results(original_fq_files, acc_contig_id_this_round,
                                 os.path.join(round_dir, "Round." + str(r_count)),
                                 os.path.join(output_base, 'temp.indices.2'), fq_info_in_memory, all_read_limits,
                                 echo_step, verbose, bool(fq_info_in_memory), log_handler)
                # clear former accepted words from memory
                del acc_words
                # then add new accepted words into memory
                if keep_seq_parts:
                    acc_words = chop_seq_list(
                        fq_simple_generator(
                            [os.path.join(round_dir, "Round." + str(r_count) + '_' + str(x + 1) + '.fq') for x in
                             range(len(original_fq_files))],
                            split_pattern=low_qual_pattern, min_sub_seq=word_size),
                        word_size, mesh_size)
                else:
                    acc_words = chop_seqs(
                        fq_simple_generator(
                            [os.path.join(round_dir, "Round." + str(r_count) + '_' + str(x + 1) + '.fq') for x in
                             range(len(original_fq_files))]),
                        word_size, mesh_size)
                acc_contig_id_this_round = set()
            log_handler.info("Round " + str(r_count) + ': ' + str(unique_id + 1) + '/' + str(len_indices) + " AI " + str(
                len_al) + " AW " + str(len_aw) + inside_memory_usage)
            #
            if len_aw == pre_aw:
                raise NoMoreReads('')
            pre_aw = len(acc_words)
            #
            if r_count == max_rounds:
                raise RoundLimitException(r_count)
            r_count += 1
            return acc_words, acc_contig_id_this_round, pre_aw, r_count, acc_num_words

        def echo_to_screen():
            inside_this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round " \
                                + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                " AI " + str(len(accepted_rd_id_this_round)) + " AW " + str(len(accepted_words))
            sys.stdout.write(inside_this_print + '\b' * len(inside_this_print))
            sys.stdout.flush()

        def check_words_limit(inside_max_n_words):
            if accumulated_num_words + len(accepted_words) - prev_aw_count > inside_max_n_words:
                if this_process:
                    inside_memory_usage = " Mem " + str(
                        round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
                else:
                    inside_memory_usage = ''
                log_handler.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                 str(len_indices) + " AI " + str(len(accepted_rd_id_this_round)) +
                                 " AW " + str(len(accepted_words)) + inside_memory_usage)
                raise WordsLimitException("")

        # core extending code
        # here efficiency is more important than code conciseness,
        # so there are four similar structure with minor differences
        reads_generator = tuple()
        while True:
            # if verbose:
            #     log_handler.info("Round " + str(round_count) + ": Start ...")

            if fq_info_in_memory:
                reads_generator = (this_read for this_read in fq_info_in_memory[0])
            else:
                if keep_seq_parts:
                    reads_generator = (this_read.strip().split("\t") for this_read in
                                       open(os.path.join(output_base, 'temp.indices.1'), 'r'))
                else:
                    reads_generator = (this_read.strip() for this_read in
                                       open(os.path.join(output_base, 'temp.indices.1'), 'r'))
            unique_read_id = 0
            if keep_seq_parts:
                if pre_grouped and g_duplicate_lines:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            if unique_read_id in line_to_accept:
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                                line_to_accept.remove(unique_read_id)
                                for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                                    seq_len = len(this_seq_part)
                                    temp_length = seq_len - word_size
                                    for i in range(0, temp_length + 1, mesh_size):
                                        # add forward
                                        accepted_words.add(this_seq_part[i:i + word_size])
                                        # add reverse
                                        accepted_words.add(this_c_seq_part[temp_length - i:seq_len - i])
                            else:
                                accepted = False
                                for this_seq_part in this_seq:
                                    seq_len = len(this_seq_part)
                                    temp_length = seq_len - word_size
                                    for i in range(0, (temp_length + 1) // 2, jump_step):
                                        # from first kmer to the middle
                                        if this_seq_part[i:i + word_size] in accepted_words:
                                            accepted = True
                                            break
                                        # from last kmer to the middle
                                        if this_seq_part[temp_length - i:seq_len - i] in accepted_words:
                                            accepted = True
                                            break
                                    if accepted:
                                        break
                                if accepted:
                                    for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                                        seq_len = len(this_seq_part)
                                        temp_length = seq_len - word_size
                                        for i in range(0, temp_length + 1, mesh_size):
                                            # add forward
                                            accepted_words.add(this_seq_part[i:i + word_size])
                                            # add reverse
                                            accepted_words.add(this_c_seq_part[temp_length - i:seq_len - i])
                                    accepted_rd_id.add(unique_read_id)
                                    accepted_rd_id_this_round.add(unique_read_id)
                                    if unique_read_id in l_with_duplicates:
                                        which_group = l_with_duplicates[unique_read_id]
                                        for id_to_accept in g_duplicate_lines[which_group]:
                                            line_to_accept.add(id_to_accept)
                                            del l_with_duplicates[id_to_accept]
                                        line_to_accept.remove(unique_read_id)
                                        del g_duplicate_lines[which_group]
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
                else:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            accepted = False
                            for this_seq_part in this_seq:
                                seq_len = len(this_seq_part)
                                temp_length = seq_len - word_size
                                for i in range(0, (temp_length + 1) // 2, jump_step):
                                    # from first kmer to the middle
                                    if this_seq_part[i:i + word_size] in accepted_words:
                                        accepted = True
                                        break
                                    # from last kmer to the middle
                                    if this_seq_part[temp_length - i:seq_len - i] in accepted_words:
                                        accepted = True
                                        break
                                if accepted:
                                    break
                            if accepted:
                                for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                                    seq_len = len(this_seq_part)
                                    temp_length = seq_len - word_size
                                    for i in range(0, temp_length + 1, mesh_size):
                                        accepted_words.add(this_seq_part[i:i + word_size])
                                        accepted_words.add(this_c_seq_part[temp_length - i:seq_len - i])
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
            else:
                if pre_grouped and g_duplicate_lines:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            seq_len = len(this_seq)
                            temp_length = seq_len - word_size
                            if unique_read_id in line_to_accept:
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                                line_to_accept.remove(unique_read_id)
                                for i in range(0, temp_length + 1, mesh_size):
                                    # add forward
                                    accepted_words.add(this_seq[i:i + word_size])
                                    # add reverse
                                    accepted_words.add(this_c_seq[temp_length - i:seq_len - i])
                            else:
                                accepted = False
                                for i in range(0, (temp_length + 1) // 2, jump_step):
                                    # from first kmer to the middle
                                    if this_seq[i:i + word_size] in accepted_words:
                                        accepted = True
                                        break
                                    # from last kmer to the middle
                                    if this_seq[temp_length - i:seq_len - i] in accepted_words:
                                        accepted = True
                                        break
                                if accepted:
                                    for i in range(0, temp_length + 1, mesh_size):
                                        # add forward
                                        accepted_words.add(this_seq[i:i + word_size])
                                        # add reverse
                                        accepted_words.add(this_c_seq[temp_length - i:seq_len - i])
                                    accepted_rd_id.add(unique_read_id)
                                    accepted_rd_id_this_round.add(unique_read_id)
                                    if unique_read_id in l_with_duplicates:
                                        which_group = l_with_duplicates[unique_read_id]
                                        for id_to_accept in g_duplicate_lines[which_group]:
                                            line_to_accept.add(id_to_accept)
                                            del l_with_duplicates[id_to_accept]
                                        line_to_accept.remove(unique_read_id)
                                        del g_duplicate_lines[which_group]
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
                else:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            accepted = False
                            seq_len = len(this_seq)
                            temp_length = seq_len - word_size
                            for i in range(0, (temp_length + 1) // 2, jump_step):
                                # from first kmer to the middle
                                if this_seq[i:i + word_size] in accepted_words:
                                    accepted = True
                                    break
                                # from last kmer to the middle
                                if this_seq[temp_length - i:seq_len - i] in accepted_words:
                                    accepted = True
                                    break
                            if accepted:
                                for i in range(0, temp_length + 1, mesh_size):
                                    accepted_words.add(this_seq[i:i + word_size])
                                    accepted_words.add(this_c_seq[temp_length - i:seq_len - i])
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
            reads_generator.close()
    except KeyboardInterrupt:
        reads_generator.close()
        if echo_step != inf:
            sys.stdout.write(' ' * 100 + '\b' * 100)
            sys.stdout.flush()
        log_handler.info(
            "Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + " AI " + str(
                len(accepted_rd_id_this_round)) + " AW " + str(len(accepted_words)))
        log_handler.info("KeyboardInterrupt")
    except NoMoreReads:
        reads_generator.close()
        if round_count < min_rounds:
            log_handler.info("No more reads found and terminated ...")
            log_handler.warning("Terminated at an insufficient number of rounds. "
                                "Try decrease '-w' if failed in the end.")
        else:
            log_handler.info("No more reads found and terminated ...")
    except WordsLimitException:
        reads_generator.close()
        if round_count <= min_rounds:
            log_handler.info("Hit the words limit and terminated ...")
            log_handler.warning("Terminated at an insufficient number of rounds, see '--max-n-words'/'--max-extending-len' for more.")
        else:
            log_handler.info("Hit the words limit and terminated ...")
    except RoundLimitException as r_lim:
        reads_generator.close()
        log_handler.info("Hit the round limit " + str(r_lim) + " and terminated ...")
    del reads_generator
    accepted_words = set()
    accepted_rd_id_this_round = set()
    del l_with_duplicates
    return accepted_rd_id


def extending_with_lim(word_size, seed_file, original_fq_files, len_indices, pre_grouped,
                       groups_of_duplicate_lines, lines_with_duplicates, group_id_to_read_counts, fq_info_in_memory,
                       output_base, max_rounds, min_rounds, fg_out_per_round, jump_step, mesh_size, verbose, resume,
                       all_read_limits, extending_dist_limit, maximum_n_words, keep_seq_parts, low_qual_pattern,
                       mean_read_len, mean_base_cov,
                       echo_step, log_handler):
    # adding initial word
    log_handler.info("Adding initial words ...")
    if keep_seq_parts:
        accepted_words = chop_seq_list_as_empty_dict(
            seq_iter=fq_simple_generator(seed_file[0], split_pattern=low_qual_pattern, min_sub_seq=word_size),
            word_size=word_size, val_len=extending_dist_limit[0])
        for go_type in range(1, len(seed_file)):
            chop_seq_list_as_empty_dict(
                seq_iter=fq_simple_generator(seed_file[go_type], split_pattern=low_qual_pattern,
                                             min_sub_seq=word_size),
                word_size=word_size, val_len=extending_dist_limit[go_type], previous_words=accepted_words)
    else:
        accepted_words = chop_seqs_as_empty_dict(
            seq_iter=fq_simple_generator(seed_file[0]),
            word_size=word_size, val_len=extending_dist_limit[0])
        for go_type in range(1, len(seed_file)):
            chop_seqs_as_empty_dict(
                seq_iter=fq_simple_generator(seed_file[go_type]),
                word_size=word_size, val_len=extending_dist_limit[go_type], previous_words=accepted_words)
    log_handler.info("AW " + str(len(accepted_words)))

    accepted_rd_id = set()
    accepted_rd_id_this_round = set()
    g_duplicate_lines = deepcopy(groups_of_duplicate_lines)
    l_with_duplicates = deepcopy(lines_with_duplicates)
    line_to_accept = {}
    round_count = 1
    initial_aw_count = len(accepted_words)
    prev_aw_count = initial_aw_count
    accumulated_num_words = initial_aw_count
    check_times = 1000
    check_step = max(int(len_indices / check_times), 1)
    if fg_out_per_round:
        round_dir = os.path.join(output_base, "intermediate_reads")
        if not os.path.exists(round_dir):
            os.mkdir(round_dir)
    if this_process and verbose:
        log_handler.warning("Package psutil is not installed, so that memory usage will not be logged\n"
                            "Don't worry. This will not affect the result.")
    try:
        def summarise_round(acc_words, acc_contig_id_this_round, pre_aw, r_count, acc_num_words, unique_id):
            len_aw = len(acc_words)
            len_al = len(acc_contig_id_this_round)
            # for check words limit; memory control
            acc_num_words += len_aw - pre_aw
            if this_process:
                inside_memory_usage = " Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
            else:
                inside_memory_usage = ''
            if fg_out_per_round:
                write_fq_results(original_fq_files, acc_contig_id_this_round,
                                 os.path.join(round_dir, "Round." + str(r_count)),
                                 os.path.join(output_base, 'temp.indices.2'), fq_info_in_memory, all_read_limits,
                                 echo_step, verbose, bool(fq_info_in_memory), log_handler)
                acc_contig_id_this_round = set()

            log_handler.info("Round " + str(r_count) + ': ' + str(unique_id + 1) + '/' + str(len_indices) + " AI " + str(
                len_al) + " AW " + str(len_aw) + inside_memory_usage)

            # cost too much time
            # acc_words = {in_k: in_v for in_k, in_v in acc_words.items() if in_v < extending_dist_limit}
            #
            if len_aw == pre_aw:
                raise NoMoreReads('')
            pre_aw = len(acc_words)
            #
            if r_count == max_rounds:
                raise RoundLimitException(r_count)
            r_count += 1
            return acc_words, acc_contig_id_this_round, pre_aw, r_count, acc_num_words

        def echo_to_screen():
            inside_this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round " \
                                + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                " AI " + str(len(accepted_rd_id_this_round)) + " AW " + str(len(accepted_words))
            sys.stdout.write(inside_this_print + '\b' * len(inside_this_print))
            sys.stdout.flush()

        def check_words_limit(inside_max_n_words):
            if accumulated_num_words + len(accepted_words) - prev_aw_count > inside_max_n_words:
                if this_process:
                    inside_memory_usage = " Mem " + str(
                        round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
                else:
                    inside_memory_usage = ''
                log_handler.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                 str(len_indices) + " AI " + str(len(accepted_rd_id_this_round)) +
                                 " AW " + str(len(accepted_words)) + inside_memory_usage)
                raise WordsLimitException("")

        # core extending code
        # here efficiency is more important than code conciseness,
        # so there are four similar structure with minor differences
        while True:
            if verbose:
                log_handler.info("Round " + str(round_count) + ": Start ...")

            if fq_info_in_memory:
                reads_generator = (this_read for this_read in fq_info_in_memory[0])
            else:
                if keep_seq_parts:
                    reads_generator = (this_read.strip().split("\t") for this_read in
                                       open(os.path.join(output_base, 'temp.indices.1'), 'r'))
                else:
                    reads_generator = (this_read.strip() for this_read in
                                       open(os.path.join(output_base, 'temp.indices.1'), 'r'))
            unique_read_id = 0
            if keep_seq_parts:
                if pre_grouped and g_duplicate_lines:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            if unique_read_id in line_to_accept:
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                                group_left = line_to_accept.pop(unique_read_id)
                                for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                                    seq_len = len(this_seq_part)
                                    temp_length = seq_len - word_size
                                    for i in range(0, temp_length + 1, mesh_size):
                                        # add forward & reverse
                                        this_w = this_seq_part[i:i + word_size]
                                        accepted_words[this_w] \
                                            = accepted_words[this_c_seq_part[temp_length - i:seq_len - i]] \
                                            = max(group_left, accepted_words.get(this_w, 0))
                            else:
                                accepted = False
                                accept_go_to_word = 0
                                part_accumulated_go_to = 0
                                accept_dist = 0
                                for this_seq_part in this_seq:
                                    seq_len = len(this_seq_part)
                                    temp_length = seq_len - word_size
                                    for i in range(0, (temp_length + 1) // 2, jump_step):
                                        # from first kmer to the middle
                                        this_w = this_seq_part[i:i + word_size]
                                        if this_w in accepted_words:
                                            accepted = True
                                            accept_go_to_word = part_accumulated_go_to + i
                                            accept_dist = accepted_words[this_w]
                                            break
                                        # from last kmer to the middle
                                        this_w = this_seq_part[temp_length - i:seq_len - i]
                                        if this_w in accepted_words:
                                            accepted = True
                                            accept_go_to_word = part_accumulated_go_to + temp_length - i
                                            accept_dist = accepted_words[this_w]
                                            break
                                    if accepted:
                                        break
                                    part_accumulated_go_to += seq_len
                                if accepted:
                                    if accept_dist - mean_read_len > 0:
                                        part_accumulated_go_to = 0
                                        for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                                            seq_len = len(this_seq_part)
                                            temp_length = seq_len - word_size
                                            for i in range(0, temp_length + 1, mesh_size):
                                                this_dist = accept_dist - \
                                                            abs(accept_go_to_word - (part_accumulated_go_to + i))
                                                this_w = this_seq_part[i:i + word_size]
                                                # add forward & reverse
                                                accepted_words[this_w] \
                                                    = accepted_words[this_c_seq_part[temp_length - i:seq_len - i]] \
                                                    = max(this_dist, accepted_words.get(this_w, 0))
                                            part_accumulated_go_to += seq_len
                                    accepted_rd_id.add(unique_read_id)
                                    accepted_rd_id_this_round.add(unique_read_id)
                                    if unique_read_id in l_with_duplicates:
                                        which_group = l_with_duplicates[unique_read_id]
                                        # N_reads = (contig_len - read_len) * (base_cov / read_len)
                                        expected_contig_len = \
                                            group_id_to_read_counts[which_group] * mean_read_len / mean_base_cov + \
                                            mean_read_len
                                        group_left = accept_dist - (expected_contig_len - word_size + 1)
                                        if group_left < 0:
                                            for id_to_accept in g_duplicate_lines[which_group]:
                                                line_to_accept[id_to_accept] = group_left
                                                del l_with_duplicates[id_to_accept]
                                            del line_to_accept[unique_read_id]
                                            del g_duplicate_lines[which_group]
                                        else:
                                            g_duplicate_lines[which_group].remove(unique_read_id)
                                            del l_with_duplicates[unique_read_id]
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count, accumulated_num_words, unique_read_id)
                else:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            accepted = False
                            accept_go_to_word = 0
                            part_accumulated_go_to = 0
                            accept_dist = 0
                            for this_seq_part in this_seq:
                                seq_len = len(this_seq_part)
                                temp_length = seq_len - word_size
                                for i in range(0, (temp_length + 1) // 2, jump_step):
                                    # from first kmer to the middle
                                    this_w = this_seq_part[i:i + word_size]
                                    if this_w in accepted_words:
                                        accepted = True
                                        accept_go_to_word = part_accumulated_go_to + i
                                        accept_dist = accepted_words[this_w]
                                        break
                                    # from last kmer to the middle
                                    this_w = this_seq_part[temp_length - i:seq_len - i]
                                    if this_w in accepted_words:
                                        accepted = True
                                        accept_go_to_word = part_accumulated_go_to + temp_length - i
                                        accept_dist = accepted_words[this_w]
                                        break
                                if accepted:
                                    break
                                part_accumulated_go_to += seq_len
                            if accepted:
                                if accept_dist - mean_read_len > 0:
                                    part_accumulated_go_to = 0
                                    for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                                        seq_len = len(this_seq_part)
                                        temp_length = seq_len - word_size
                                        for i in range(0, temp_length + 1, mesh_size):
                                            this_dist = accept_dist - \
                                                        abs(accept_go_to_word - (part_accumulated_go_to + i))
                                            # if this_dist < extending_dist_limit:
                                            this_w = this_seq_part[i:i + word_size]
                                            accepted_words[this_w] \
                                                = accepted_words[this_c_seq_part[temp_length - i:seq_len - i]] \
                                                = max(this_dist, accepted_words.get(this_w, 0))
                                        part_accumulated_go_to += seq_len
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
            else:
                if pre_grouped and g_duplicate_lines:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            seq_len = len(this_seq)
                            temp_length = seq_len - word_size
                            if unique_read_id in line_to_accept:
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                                group_left = line_to_accept.pop(unique_read_id)
                                for i in range(0, temp_length + 1, mesh_size):
                                    # add forward & reverse
                                    this_w = this_seq[i:i + word_size]
                                    accepted_words[this_w] \
                                        = accepted_words[this_c_seq[temp_length - i:seq_len - i]] \
                                        = max(group_left, accepted_words.get(this_w, 0))
                            else:
                                accepted = False
                                accept_go_to_word = 0
                                accept_dist = 0
                                for i in range(0, (temp_length + 1) // 2, jump_step):
                                    # from first kmer to the middle
                                    this_w = this_seq[i:i + word_size]
                                    if this_w in accepted_words:
                                        accepted = True
                                        accept_go_to_word = i
                                        accept_dist = accepted_words[this_w]
                                        break
                                    # from last kmer to the middle
                                    this_w = this_seq[temp_length - i:seq_len - i]
                                    if this_w in accepted_words:
                                        accepted = True
                                        accept_go_to_word = temp_length - i
                                        accept_dist = accepted_words[this_w]
                                        break
                                if accepted:
                                    if accept_dist - mean_read_len > 0:
                                        for i in range(0, temp_length + 1, mesh_size):
                                            this_dist = accept_dist - abs(accept_go_to_word - i)
                                            # if this_dist < extending_dist_limit:
                                            # add forward & reverse
                                            this_w = this_seq[i:i + word_size]
                                            accepted_words[this_w] \
                                                = accepted_words[this_c_seq[temp_length - i:seq_len - i]] \
                                                = max(this_dist, accepted_words.get(this_w, 0))
                                    accepted_rd_id.add(unique_read_id)
                                    accepted_rd_id_this_round.add(unique_read_id)
                                    if unique_read_id in l_with_duplicates:
                                        which_group = l_with_duplicates[unique_read_id]
                                        # using unique reads
                                        expected_contig_len = \
                                            group_id_to_read_counts[which_group] * mean_read_len / mean_base_cov + \
                                            mean_read_len
                                        # print(group_id_to_read_counts[which_group], expected_contig_len)
                                        group_left = accept_dist - (expected_contig_len - word_size + 1)
                                        if group_left < 0:
                                            for id_to_accept in g_duplicate_lines[which_group]:
                                                line_to_accept[id_to_accept] = group_left
                                                del l_with_duplicates[id_to_accept]
                                            del line_to_accept[unique_read_id]
                                            del g_duplicate_lines[which_group]
                                        else:
                                            g_duplicate_lines[which_group].remove(unique_read_id)
                                            del l_with_duplicates[unique_read_id]
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
                else:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_rd_id:
                            accepted = False
                            accept_go_to_word = 0
                            accept_dist = 0
                            seq_len = len(this_seq)
                            temp_length = seq_len - word_size
                            for i in range(0, (temp_length + 1) // 2, jump_step):
                                # from first kmer to the middle
                                this_w = this_seq[i:i + word_size]
                                if this_w in accepted_words:
                                    accepted = True
                                    accept_go_to_word = i
                                    accept_dist = accepted_words[this_w]
                                    break
                                # from last kmer to the middle
                                this_w = this_seq[temp_length - i:seq_len - i]
                                if this_w in accepted_words:
                                    accepted = True
                                    accept_dist = accepted_words[this_w]
                                    break
                            if accepted:
                                if accept_dist - mean_read_len > 0:
                                    for i in range(0, temp_length + 1, mesh_size):
                                        this_dist = accept_dist - abs(accept_go_to_word - i)
                                        # if this_dist < extending_dist_limit:
                                        this_w = this_seq[i:i + word_size]
                                        accepted_words[this_w] \
                                            = accepted_words[this_c_seq[temp_length - i:seq_len - i]] \
                                            = max(this_dist, accepted_words.get(this_w, 0))
                                accepted_rd_id.add(unique_read_id)
                                accepted_rd_id_this_round.add(unique_read_id)
                        if echo_step != inf and unique_read_id % echo_step == 0:
                            echo_to_screen()
                        if unique_read_id % check_step == 0:
                            check_words_limit(maximum_n_words)
                    accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                        = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                          round_count,
                                          accumulated_num_words, unique_read_id)
            reads_generator.close()
    except KeyboardInterrupt:
        reads_generator.close()
        if echo_step != inf:
            sys.stdout.write(' ' * 100 + '\b' * 100)
            sys.stdout.flush()
        log_handler.info(
            "Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + " AI " + str(
                len(accepted_rd_id_this_round)) + " AW " + str(len(accepted_words)))
        log_handler.info("KeyboardInterrupt")
    except NoMoreReads:
        reads_generator.close()
        if round_count < min_rounds:
            log_handler.info("No more reads found and terminated ...")
            log_handler.warning("Terminated at an insufficient number of rounds. "
                                "Try decrease '-w' if failed in the end.")
        else:
            log_handler.info("No more reads found and terminated ...")
    except WordsLimitException:
        reads_generator.close()
        if round_count <= min_rounds:
            log_handler.info("Hit the words limit and terminated ...")
            log_handler.warning("Terminated at an insufficient number of rounds, see '--max-n-words'/'--max-extending-len' for more.")
        else:
            log_handler.info("Hit the words limit and terminated ...")
    except RoundLimitException as r_lim:
        reads_generator.close()
        log_handler.info("Hit the round limit " + str(r_lim) + " and terminated ...")
    del reads_generator
    accepted_words = set()
    accepted_rd_id_this_round = set()
    del l_with_duplicates
    return accepted_rd_id


def get_anti_with_fas(word_size, anti_words, anti_input, original_fq_files, log_handler):
    anti_lines = set()
    pre_reading_handler = [open(fq_file, 'r') for fq_file in original_fq_files]
    line_count = 0

    def add_to_anti_lines(here_head):
        try:
            if ' ' in here_head:
                here_head_split = here_head.split(' ')
                this_name, direction = here_head_split[0], int(here_head_split[1][0])
            elif '#' in here_head:
                here_head_split = here_head.split('#')
                this_name, direction = here_head_split[0], int(here_head_split[1].strip("/")[0])
            elif here_head[-2] == "/" and here_head[-1].isdigit():    # 2019-04-22 added
                this_name, direction = here_head[:-2], int(here_head[-1])
            else:
                this_name, direction = here_head, 1
        except (ValueError, IndexError):
            log_handler.error('Unrecognized fq format in ' + str(line_count))
            exit()
        else:
            anti_lines.add((this_name, direction))

    for file_in in pre_reading_handler:
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
                    for i in range(0, temp_length + 1):
                        if this_seq[i:i + word_size] in anti_words:
                            add_to_anti_lines(this_head)
                            break
                        if this_c_seq[i:i + word_size] in anti_words:
                            add_to_anti_lines(this_head)
                            break
                else:
                    log_handler.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                    exit()
                line_count += 1
                for i in range(3):
                    line = file_in.readline()
                    line_count += 1
        file_in.close()
    return anti_lines


def making_seed_reads_using_mapping(seed_file, original_fq_files,
                                    out_base, resume, verbose_log, threads, random_seed, organelle_type, prefix,
                                    keep_temp, bowtie2_other_options, log_handler, which_bowtie2=""):
    seed_dir = os.path.join(out_base, prefix + "seed")
    if not os.path.exists(seed_dir):
        os.mkdir(seed_dir)
    if sum([os.path.exists(remove_db_postfix(seed_file) + ".index" + postfix)
            for postfix in
            (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")]) != 6:
        new_seed_file = os.path.join(seed_dir, os.path.basename(seed_file))
        check_fasta_seq_names(seed_file, new_seed_file, log_handler)
        seed_file = new_seed_file
    bowtie_out_base = os.path.join(seed_dir, prefix + organelle_type + ".initial")
    total_seed_fq = bowtie_out_base + ".fq"
    total_seed_sam = bowtie_out_base + ".sam"
    seed_index_base = seed_file + '.index'
    map_with_bowtie2(seed_file=seed_file, original_fq_files=original_fq_files,
                     bowtie_out=bowtie_out_base, resume=resume, threads=threads, random_seed=random_seed,
                     generate_fq=True, target_echo_name="seed", log_handler=log_handler, verbose_log=verbose_log,
                     which_bowtie2=which_bowtie2, bowtie2_mode="", bowtie2_other_options=bowtie2_other_options)
    if not keep_temp:
        for seed_index_file in [x for x in os.listdir(seed_dir) if x.startswith(os.path.basename(seed_index_base))]:
            os.remove(os.path.join(seed_dir, seed_index_file))
    seed_fq_size = os.path.getsize(total_seed_fq)
    if not seed_fq_size:
        if log_handler:
            log_handler.error("No " + str(organelle_type) + " seed reads found!")
            log_handler.error("Please check your raw data or change your " + str(organelle_type) + " seed!")
        exit()
    log_handler.info("Seed reads made: " + total_seed_fq + " (" + str(int(seed_fq_size)) + " bytes)")
    if seed_fq_size < 10000:
        log_handler.error("Too few seed reads found! "
                          "Please change your seed file (-s) or "
                          "increase your data input (--max-reads/--reduce-reads-for-coverage)!")
        exit()
    return total_seed_fq, total_seed_sam, seed_file


def get_anti_lines_using_mapping(anti_seed, seed_sam_files, original_fq_files,
                                 out_base, resume, verbose_log, threads,
                                 random_seed, prefix, keep_temp, bowtie2_other_options, log_handler, which_bowtie2=""):
    from GetOrganelleLib.sam_parser import get_heads_from_sam_fast
    seed_dir = os.path.join(out_base, prefix + "seed")
    if not os.path.exists(seed_dir):
        os.mkdir(seed_dir)
    if anti_seed:
        new_anti_seed = os.path.join(seed_dir, os.path.basename(anti_seed))
        check_fasta_seq_names(anti_seed, new_anti_seed, log_handler)
        anti_seed = new_anti_seed
    else:
        anti_seed = ""
    anti_index_base = anti_seed + '.index'
    bowtie_out_base = os.path.join(out_base, prefix + "anti_seed_bowtie")
    anti_seed_sam = [os.path.join(out_base, x + prefix + "anti_seed_bowtie.sam") for x in ("temp.", "")]

    if anti_seed:
        map_with_bowtie2(seed_file=anti_seed, original_fq_files=original_fq_files,
                         bowtie_out=bowtie_out_base, resume=resume, threads=threads, random_seed=random_seed,
                         log_handler=log_handler, target_echo_name="anti-seed", generate_fq=False, verbose_log=verbose_log,
                         which_bowtie2=which_bowtie2, bowtie2_mode="", bowtie2_other_options=bowtie2_other_options)

        log_handler.info("Parsing bowtie2 result ...")
        anti_lines = get_heads_from_sam_fast(anti_seed_sam[1]) - get_heads_from_sam_fast(*seed_sam_files)
        log_handler.info("Parsing bowtie2 result finished ...")
    else:
        anti_lines = set()

    if not keep_temp:
        for anti_index_file in [x for x in os.listdir(seed_dir) if x.startswith(os.path.basename(anti_index_base))]:
            os.remove(os.path.join(seed_dir, anti_index_file))

    return anti_lines


def assembly_with_spades(spades_kmer, spades_out_put, parameters, out_base, prefix, original_fq_files, reads_paired,
                         which_spades, verbose_log, resume, threads, log_handler):
    if '-k' in parameters or not spades_kmer:
        kmer = ''
    else:
        kmer = '-k ' + spades_kmer
    if resume and os.path.exists(spades_out_put):
        spades_command = os.path.join(which_spades, "spades.py") + " --continue -o " + spades_out_put
    else:
        spades_out_command = '-o ' + spades_out_put
        if reads_paired['input'] and reads_paired['pair_out']:
            all_unpaired = []
            # spades does not accept empty files
            if os.path.getsize(os.path.join(out_base, prefix + "extended_1_unpaired.fq")):
                all_unpaired.append(os.path.join(out_base, prefix + "extended_1_unpaired.fq"))
            if os.path.getsize(os.path.join(out_base, prefix + "extended_2_unpaired.fq")):
                all_unpaired.append(os.path.join(out_base, prefix + "extended_2_unpaired.fq"))
            for iter_unpaired in range(len(original_fq_files) - 2):
                if os.path.getsize(str(os.path.join(out_base, prefix + "extended_" + str(iter_unpaired + 3) + ".fq"))):
                    all_unpaired.append(
                        str(os.path.join(out_base, prefix + "extended_" + str(iter_unpaired + 3) + ".fq")))
            if os.path.getsize(os.path.join(out_base, prefix + "extended_1_paired.fq")):
                spades_command = ' '.join(
                    [os.path.join(which_spades, "spades.py"), '-t', str(threads), parameters, '-1',
                     os.path.join(out_base, prefix + "extended_1_paired.fq"), '-2',
                     os.path.join(out_base, prefix + "extended_2_paired.fq")] +
                    ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                    [kmer, spades_out_command]).strip()
            else:
                # log_handler.warning("No paired reads found for the target!?")
                spades_command = ' '.join(
                    [os.path.join(which_spades, "spades.py"), '-t', str(threads), parameters] +
                    ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                    [kmer, spades_out_command]).strip()
        else:
            all_unpaired = []
            for iter_unpaired in range(len(original_fq_files)):
                if os.path.getsize(str(os.path.join(out_base, prefix + "extended_" + str(iter_unpaired + 1) + ".fq"))):
                    all_unpaired.append(
                        str(os.path.join(out_base, prefix + "extended_" + str(iter_unpaired + 1) + ".fq")))
            spades_command = ' '.join(
                [os.path.join(which_spades, "spades.py"), '-t', str(threads), parameters] +
                ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                [kmer, spades_out_command]).strip()
    log_handler.info(spades_command)
    spades_running = subprocess.Popen(spades_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    # output, err = spades_running.communicate()
    output = monitor_spades_log(spades_running, log_handler)
    if not os.path.exists(spades_out_put):
        log_handler.error("Assembling failed.")
        return False
    elif "== Error ==" in output or "terminated by segmentation fault" in output:
        # check when other kmer assembly results were produced
        real_kmer_values = sorted([int(kmer_d[1:])
                                   for kmer_d in os.listdir(spades_out_put)
                                   if os.path.isdir(os.path.join(spades_out_put, kmer_d))
                                   and kmer_d.startswith("K")])
        real_kmer_values = [str(k_val) for k_val in real_kmer_values]
        temp_res = False
        failed_at_k = None
        for kmer_val in real_kmer_values:
            this_k_path = os.path.join(spades_out_put, "K" + kmer_val)
            if os.path.exists(os.path.join(this_k_path, "assembly_graph.fastg")):
                temp_res = True
            else:
                failed_at_k = kmer_val
        if temp_res:
            if failed_at_k:
                log_handler.warning("SPAdes failed for '-k " + failed_at_k + "'!")
                log_handler.warning("If you need result based on kmer=" + failed_at_k + " urgently, "
                    "please check " + os.path.join(spades_out_put, "spades.log"))
                del real_kmer_values[real_kmer_values.index(failed_at_k)]
                log_handler.warning("GetOrganelle would continue to process results based on "
                                    "kmer=" + ",".join(real_kmer_values) + ".")
                # os.system("cp " + os.path.join(spades_out_put, "K" + real_kmer_values[-1], "assembly_graph.fastg")
                #           + " " + spades_out_put)
                log_handler.info('Assembling finished with warnings.\n')
                return True
            else:
                log_handler.warning("SPAdes failed with unknown errors!")
                log_handler.warning("If you need to know more details, please check " +
                                    os.path.join(spades_out_put, "spades.log") + " and contact SPAdes developers.")
                log_handler.warning("GetOrganelle would continue to process results based on "
                                    "kmer=" + ",".join(real_kmer_values) + ".")
                # os.system("cp " + os.path.join(spades_out_put, "K" + real_kmer_values[-1], "assembly_graph.fastg")
                #           + " " + spades_out_put)
                log_handler.info("Assembling finished with warnings.\n")
                return True
        else:
            if "mmap(2) failed" in output:
                # https://github.com/ablab/spades/issues/91
                log_handler.error("Guessing your output directory is inside a VirtualBox shared folder!")
                log_handler.error("Assembling failed.")
            else:
                log_handler.error("Assembling failed.")
            return False
    elif not os.path.exists(os.path.join(spades_out_put, "assembly_graph.fastg")):
        if verbose_log:
            log_handler.info(output)
        log_handler.warning("Assembling exited halfway.\n")
        return True
    else:
        spades_log = output.split("\n")
        if verbose_log:
            log_handler.info(output)
        for line in spades_log:
            line = line.strip()
            if line.count(":") > 2 and "Insert size = " in line and \
                    line.split()[0].replace(":", "").replace(".", "").isdigit():
                try:
                    log_handler.info(line.split("   ")[-1].split(", read length =")[0].strip())
                except IndexError:
                    pass
        log_handler.info('Assembling finished.\n')
        return True


def slim_spades_result(organelle_types, in_custom, ex_custom, spades_output, ignore_kmer_res, max_slim_extending_len,
                       verbose_log, log_handler, threads, which_blast="", resume=False, keep_temp=False):
    if executable(os.path.join(UTILITY_PATH, "slim_graph.py -h")):
        which_slim = UTILITY_PATH
    elif executable(os.path.join(PATH_OF_THIS_SCRIPT, "slim_graph.py -h")):
        which_slim = PATH_OF_THIS_SCRIPT
    elif executable("slim_graph.py -h"):
        which_slim = ""
    else:
        raise Exception("slim_graph.py not found!")
    slim_stat_list = []
    if not executable(os.path.join(which_blast, "blastn")):
        if log_handler:
            log_handler.warning(
                os.path.join(which_blast, "blastn") + " not accessible! Skip slimming assembly result ...")
        slim_stat_list.append((1, None))
        return slim_stat_list
    if not executable(os.path.join(which_blast, "makeblastdb")):
        if log_handler:
            log_handler.warning(
                os.path.join(which_blast, "makeblastdb") + " not accessible! Skip slimming assembly result ...")
        slim_stat_list.append((1, None))
        return slim_stat_list
    include_priority_db = []
    exclude_db = []
    if in_custom or ex_custom:
        include_priority_db = in_custom
        exclude_db = ex_custom
    else:
        if organelle_types == ["embplant_pt"]:
            include_priority_db = [os.path.join(_LBL_DB_PATH, "embplant_pt.fasta"),
                                   os.path.join(_LBL_DB_PATH, "embplant_mt.fasta")]
            max_slim_extending_len = \
                max_slim_extending_len if max_slim_extending_len else MAX_SLIM_EXTENDING_LENS[organelle_types[0]]
        elif organelle_types == ["embplant_mt"]:
            include_priority_db = [os.path.join(_LBL_DB_PATH, "embplant_mt.fasta"),
                                   os.path.join(_LBL_DB_PATH, "embplant_pt.fasta")]
            max_slim_extending_len = \
                max_slim_extending_len if max_slim_extending_len else MAX_SLIM_EXTENDING_LENS[organelle_types[0]]
        else:
            include_priority_db = [os.path.join(_LBL_DB_PATH, sub_organelle_t + ".fasta")
                                   for sub_organelle_t in organelle_types]
            if max_slim_extending_len is None:
                max_slim_extending_len = max([MAX_SLIM_EXTENDING_LENS[sub_organelle_t]
                                              for sub_organelle_t in organelle_types])
    kmer_values = sorted([int(kmer_d[1:])
                          for kmer_d in os.listdir(spades_output)
                          if os.path.isdir(os.path.join(spades_output, kmer_d))
                          and kmer_d.startswith("K")
                          and os.path.exists(os.path.join(spades_output, kmer_d, "assembly_graph.fastg"))],
                         reverse=True)
    if not kmer_values:
        return [], ignore_kmer_res  # to avoid "ValueError: max() arg is an empty sequence"
    if max(kmer_values) <= ignore_kmer_res:
        log_handler.info("Small kmer values, resetting \"--ignore-k -1\"")
        ignore_kmer_res = -1
    kmer_dirs = [os.path.join(spades_output, "K" + str(kmer_val))
                 for kmer_val in kmer_values if kmer_val > ignore_kmer_res]
    in_ex_info = generate_in_ex_info_name(include_indices=include_priority_db, exclude_indices=exclude_db)
    for kmer_dir in kmer_dirs:
        graph_file = os.path.join(kmer_dir, "assembly_graph.fastg")
        this_fastg_file_out = os.path.join(kmer_dir, "assembly_graph.fastg" + in_ex_info + ".fastg")
        if resume:
            if os.path.exists(this_fastg_file_out):
                if log_handler:
                    log_handler.info("Slimming " + graph_file + " ... skipped.")
                slim_stat_list.append((0, this_fastg_file_out))
                continue
        run_command = ""
        if include_priority_db:
            run_command += " --include-priority " + ",".join(include_priority_db)
        if exclude_db:
            run_command += " --exclude " + ",".join(exclude_db)
        which_bl_str = " --which-blast " + which_blast if which_blast else ""
        run_command = os.path.join(which_slim, "slim_graph.py") + " --verbose " * int(bool(verbose_log)) + \
                      " --log --wrapper -t " + str(threads) + " --keep-temp " * int(bool(keep_temp)) + \
                      (" --max-slim-extending-len " +
                       str(max_slim_extending_len) + " ") * int(bool(max_slim_extending_len)) + \
                      which_bl_str + " " + graph_file + run_command  # \
        slim_spades = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log and log_handler:
            log_handler.info(run_command)
        output, err = slim_spades.communicate()
        output_file_list = [os.path.join(kmer_dir, x) for x in os.listdir(kmer_dir) if x.count(".fastg") == 2]
        if " failed" in output.decode("utf8") or "- ERROR:" in output.decode("utf8"):
            if log_handler:
                if verbose_log:
                    log_handler.error(output.decode("utf8"))
                log_handler.error("Slimming " + graph_file + " failed. "
                                  "Please check " + os.path.join(kmer_dir, "slim.log.txt") + " for details. ")
            slim_stat_list.append((1, None))
        elif output_file_list and os.path.getsize(output_file_list[0]) == 0:
            if log_handler:
                log_handler.warning("Slimming " + graph_file + " finished with no target organelle contigs found!")
            slim_stat_list.append((2, None))
        elif output_file_list:
            if log_handler:
                if verbose_log:
                    log_handler.info(output.decode("utf8"))
                log_handler.info("Slimming " + graph_file + " finished!")
            slim_stat_list.append((0, this_fastg_file_out))
        else:
            slim_stat_list.append((1, None))
    return slim_stat_list, ignore_kmer_res


def separate_fq_by_pair(out_base, prefix, verbose_log, log_handler):
    log_handler.info("Separating extended fastq file ... ")
    out_paired_1 = os.path.join(out_base, prefix + "extended_1_paired.fq")
    out_paired_2 = os.path.join(out_base, prefix + "extended_2_paired.fq")
    out_unpaired_1 = os.path.join(out_base, prefix + "extended_1_unpaired.fq")
    out_unpaired_2 = os.path.join(out_base, prefix + "extended_2_unpaired.fq")
    get_paired_and_unpaired_reads(input_fq_1=os.path.join(out_base, prefix + "extended_1.fq"),
                                  input_fq_2=os.path.join(out_base, prefix + "extended_2.fq"),
                                  output_p_1=out_paired_1,
                                  output_p_2=out_paired_2,
                                  output_u_1=out_unpaired_1,
                                  output_u_2=out_unpaired_2)
    if not os.path.getsize(out_paired_1) and not os.path.getsize(out_paired_2):
        log_handler.warning("No paired reads found?!")
    return True


def extract_organelle_genome(out_base, spades_output, ignore_kmer_res, slim_out_fg, organelle_prefix,
                             organelle_type, blast_db, read_len_for_log, verbose, log_handler, basic_prefix,
                             expected_maximum_size, expected_minimum_size, do_spades_scaffolding, options):

    from GetOrganelleLib.assembly_parser import ProcessingGraphFailed, Assembly
    random.seed(options.random_seed)
    # import numpy as np
    # np.random.seed(options.random_seed)

    def disentangle_assembly(assembly_obj, fastg_file, tab_file, output, weight_factor, log_dis, time_limit,
                             type_factor=3.,
                             mode="embplant_pt", blast_db_base="embplant_pt", contamination_depth=3.,
                             contamination_similarity=0.95, degenerate=True,
                             degenerate_depth=1.5, degenerate_similarity=0.98,
                             expected_max_size=inf, expected_min_size=0, hard_cov_threshold=10.,
                             min_sigma_factor=0.1, here_only_max_c=True, with_spades_scaffolds=False,
                             here_acyclic_allowed=False,
                             here_verbose=False, timeout_flag_str="'--disentangle-time-limit'", temp_graph=None):
        @set_time_limit(time_limit, flag_str=timeout_flag_str)
        def disentangle_inside(input_graph,
                               fastg_f,
                               tab_f,
                               o_p,
                               w_f,
                               log_in,
                               type_f=3.,
                               mode_in="embplant_pt",
                               in_db_n="embplant_pt", c_d=3., c_s=0.95,
                               deg=True, deg_dep=1.5, deg_sim=0.98, hard_c_t=10., min_s_f=0.1, max_c_in=True,
                               max_s=inf, min_s=0, with_spades_scaffolds_in=False,
                               acyclic_allowed_in=False, verbose_in=False, in_temp_graph=None):
            image_produced = False
            this_K = os.path.split(os.path.split(fastg_f)[0])[-1]
            o_p += "." + this_K
            if with_spades_scaffolds_in:
                log_in.info("Scaffolding disconnected contigs using SPAdes scaffolds ... ")
                log_in.warning("Assembly based on scaffolding may not be as accurate as "
                               "the ones directly exported from the assembly graph.")
            if acyclic_allowed_in:
                log_in.info("Disentangling " + fastg_f + " as a/an " + in_db_n + "-insufficient graph ... ")
            else:
                log_in.info("Disentangling " + fastg_f + " as a circular genome ... ")
            # input_graph = Assembly(fastg_f)
            if with_spades_scaffolds_in:
                if not input_graph.add_gap_nodes_with_spades_res(os.path.join(spades_output, "scaffolds.fasta"),
                                                                 os.path.join(spades_output, "scaffolds.paths"),
                                                                 # min_cov=options.min_depth, max_cov=options.max_depth,
                                                                 log_handler=log_handler):
                    raise ProcessingGraphFailed("No new connections.")
                else:
                    log_handler.info("Re-loading labels along " + slim_out_fg)
                    input_graph.parse_tab_file(
                        tab_f,
                        database_name=in_db_n,
                        type_factor=type_f,
                        max_gene_gap=250,
                        max_cov_diff=hard_c_t,  # contamination_depth?
                        verbose=verbose,
                        log_handler=log_handler,
                        random_obj=random,
                        np_rd_obj=np.random)
                    if in_temp_graph:
                        if in_temp_graph.endswith(".gfa"):
                            this_tmp_graph = in_temp_graph[:-4] + ".scaffolds.gfa"
                        else:
                            this_tmp_graph = in_temp_graph + ".scaffolds.gfa"
                        input_graph.write_to_gfa(this_tmp_graph)
            target_results = input_graph.find_target_graph(  # tab_f,
                                                           mode=mode_in,
                                                           db_name=in_db_n,
                                                           # type_factor=type_f,
                                                           hard_cov_threshold=hard_c_t,
                                                           contamination_depth=c_d,
                                                           contamination_similarity=c_s,
                                                           degenerate=deg, degenerate_depth=deg_dep,
                                                           degenerate_similarity=deg_sim,
                                                           expected_max_size=max_s, expected_min_size=min_s,
                                                           only_keep_max_cov=max_c_in,
                                                           min_sigma_factor=min_s_f,
                                                           weight_factor=w_f,
                                                           broken_graph_allowed=acyclic_allowed_in,
                                                           read_len_for_log=read_len_for_log,
                                                           kmer_for_log=int(this_K[1:]),
                                                           log_handler=log_in, verbose=verbose_in,
                                                           temp_graph=in_temp_graph,
                                                           selected_graph=o_p + ".graph.selected_graph.gfa",
                                                           random_obj=random)
            if not target_results:
                raise ProcessingGraphFailed("No target graph detected!")
            if len(target_results) > 1:
                log_in.warning(str(len(target_results)) + " sets of graph detected!")
            # log_in.info("Slimming and disentangling graph finished!")

            log_in.info("Writing output ...")
            ambiguous_base_used = False
            if acyclic_allowed_in:
                contig_num = set()
                still_complete = []
                for go_res, res in enumerate(target_results):
                    go_res += 1
                    broken_graph = res["graph"]
                    count_path = 0
                    # use options.max_paths_num + 1 to trigger the warning
                    these_paths = broken_graph.get_all_paths(mode=mode_in, log_handler=log_in,
                                                             max_paths_num=options.max_paths_num + 1)
                    # reducing paths
                    if len(these_paths) > options.max_paths_num:
                        log_in.warning("Only exporting " + str(options.max_paths_num) + " out of all " +
                                       str(options.max_paths_num) +
                                       "+ possible paths. (see '--max-paths-num' to change it.)")
                        these_paths = these_paths[:options.max_paths_num]
                    # exporting paths, reporting results
                    for this_paths, other_tag in these_paths:
                        count_path += 1
                        all_contig_str = []
                        contig_num.add(len(this_paths))
                        contigs_are_circular = []
                        for go_contig, this_p_part in enumerate(this_paths):
                            this_contig = broken_graph.export_path(this_p_part)
                            if DEGENERATE_BASES & set(this_contig.seq):
                                ambiguous_base_used = True
                            if this_contig.label.endswith("(circular)"):
                                contigs_are_circular.append(True)
                            else:
                                contigs_are_circular.append(False)
                            if len(this_paths) == 1 and contigs_are_circular[-1]:
                                all_contig_str.append(this_contig.fasta_str())
                            else:
                                all_contig_str.append(">scaffold_" + str(go_contig + 1) + "--" + this_contig.label +
                                                          "\n" + this_contig.seq + "\n")

                        if len(all_contig_str) == 1 and set(contigs_are_circular) == {True}:
                            if "GAP" in all_contig_str:
                                still_complete.append("nearly-complete")
                            else:
                                still_complete.append("complete")
                                # print ir stat
                                if count_path == 1 and in_db_n == "embplant_pt":
                                    detect_seq = broken_graph.export_path(this_paths[0]).seq
                                    ir_stats = detect_plastome_architecture(detect_seq, 1000)
                                    log_in.info("Detecting large repeats (>1000 bp) in PATH1 with " + ir_stats[-1] +
                                                ", Total:LSC:SSC:Repeat(bp) = " + str(len(detect_seq)) + ":" +
                                                ":".join([str(len_val) for len_val in ir_stats[:3]]))
                        else:
                            still_complete.append("incomplete")

                        if still_complete[-1] == "complete":
                            out_n = o_p + ".complete.graph" + str(go_res) + "." + \
                                    str(count_path) + other_tag + ".path_sequence.fasta"
                            log_in.info("Writing PATH" + str(count_path) + " of complete " + mode_in + " to " + out_n)
                        elif still_complete[-1] == "nearly-complete":
                            out_n = o_p + ".nearly-complete.graph" + str(go_res) + "." + \
                                    str(count_path) + other_tag + ".path_sequence.fasta"
                            log_in.info(
                                "Writing PATH" + str(count_path) + " of nearly-complete " + mode_in + " to " + out_n)
                        else:
                            out_n = o_p + ".scaffolds.graph" + str(go_res) + other_tag + "." + \
                                    str(count_path) + ".path_sequence.fasta"
                            log_in.info(
                                "Writing PATH" + str(count_path) + " of " + mode_in + " scaffold(s) to " + out_n)
                        open(out_n, "w").write("\n".join(all_contig_str))

                    if set(still_complete[-len(these_paths):]) == {"complete"}:
                        this_out_base = o_p + ".complete.graph" + str(go_res) + ".path_sequence."
                        log_in.info("Writing GRAPH to " + this_out_base + "gfa")
                        broken_graph.write_to_gfa(this_out_base + "gfa")
                        image_produced = draw_assembly_graph_using_bandage(
                            input_graph_file=this_out_base + "gfa",
                            output_image_file=this_out_base + "png",
                            assembly_graph_ob=broken_graph,
                            log_handler=log_handler, verbose_log=verbose_in, which_bandage=options.which_bandage)
                    elif set(still_complete[-len(these_paths):]) == {"nearly-complete"}:
                        this_out_base = o_p + ".nearly-complete.graph" + str(go_res) + ".path_sequence."
                        log_in.info("Writing GRAPH to " + this_out_base + "gfa")
                        broken_graph.write_to_gfa(this_out_base + "gfa")
                        image_produced = draw_assembly_graph_using_bandage(
                            input_graph_file=this_out_base + "gfa",
                            output_image_file=this_out_base + "png",
                            assembly_graph_ob=broken_graph,
                            log_handler=log_handler, verbose_log=verbose_in, which_bandage=options.which_bandage)
                    else:
                        this_out_base = o_p + ".contigs.graph" + str(go_res) + ".path_sequence."
                        log_in.info("Writing GRAPH to " + this_out_base + "gfa")
                        broken_graph.write_to_gfa(this_out_base + "gfa")
                        # image_produced = draw_assembly_graph_using_bandage(
                        #     input_graph_file=this_out_base + "gfa",
                        #     output_image_file=this_out_base + "png",
                        #     assembly_graph_ob=broken_graph,
                        #     log_handler=log_handler, verbose_log=verbose_in, which_bandage=options.which_bandage)
                if set(still_complete) == {"complete"}:
                    log_in.info("Result status of " + mode_in + ": circular genome")
                elif set(still_complete) == {"nearly-complete"}:
                    log_in.info("Result status of " + mode_in + ": circular genome with gaps")
                else:
                    log_in.info("Result status of " + mode_in + ": " +
                                ",".join(sorted([str(c_n) for c_n in contig_num])) + " scaffold(s)")
            else:
                status_str = "complete"
                for go_res, res in enumerate(target_results):
                    go_res += 1
                    idealized_graph = res["graph"]
                    count_path = 0
                    # use options.max_paths_num + 1 to trigger the warning
                    these_paths = idealized_graph.get_all_circular_paths(
                        mode=mode_in, log_handler=log_in, reverse_start_direction_for_pt=options.reverse_lsc,
                        max_paths_num=options.max_paths_num + 1)
                    # reducing paths
                    if len(these_paths) > options.max_paths_num:
                        log_in.warning("Only exporting " + str(options.max_paths_num) + " out of all " +
                                       str(options.max_paths_num) +
                                       "+ possible paths. (see '--max-paths-num' to change it.)")
                        these_paths = these_paths[:options.max_paths_num]

                    # exporting paths, reporting results
                    for this_path, other_tag in these_paths:
                        count_path += 1
                        this_seq_obj = idealized_graph.export_path(this_path)
                        if DEGENERATE_BASES & set(this_seq_obj.seq):
                            ambiguous_base_used = True
                            status_str = "nearly-complete"
                        out_n = o_p + "." + status_str + ".graph" + str(go_res) + "." + str(
                            count_path) + other_tag + ".path_sequence.fasta"
                        open(out_n, "w").write(this_seq_obj.fasta_str())

                        # print ir stat
                        if count_path == 1 and in_db_n == "embplant_pt" and not ambiguous_base_used:
                            detect_seq = this_seq_obj.seq
                            ir_stats = detect_plastome_architecture(detect_seq, 1000)
                            log_in.info("Detecting large repeats (>1000 bp) in PATH1 with " + ir_stats[-1] +
                                        ", Total:LSC:SSC:Repeat(bp) = " + str(len(detect_seq)) + ":" +
                                        ":".join([str(len_val) for len_val in ir_stats[:3]]))
                        log_in.info(
                            "Writing PATH" + str(count_path) + " of " + status_str + " " + mode_in + " to " + out_n)
                    temp_base_out = o_p + "." + status_str + ".graph" + str(go_res) + ".path_sequence."
                    log_in.info("Writing GRAPH to " + temp_base_out + "gfa")
                    idealized_graph.write_to_gfa(temp_base_out + "gfa")
                    image_produced = draw_assembly_graph_using_bandage(
                        input_graph_file=temp_base_out + "gfa", output_image_file=temp_base_out + "png",
                        assembly_graph_ob=idealized_graph, log_handler=log_handler, verbose_log=verbose_in,
                        which_bandage=options.which_bandage)
                if ambiguous_base_used:
                    log_in.info("Result status of " + mode_in + ": circular genome with gaps")
                else:
                    log_in.info("Result status of " + mode_in + ": circular genome")
            if ambiguous_base_used:
                log_in.warning("Ambiguous base(s) used!")
            o_p_extended = os.path.join(os.path.split(o_p)[0], basic_prefix + "extended_" + this_K + ".")
            os.system("cp " + os.path.join(os.path.split(fastg_f)[0], "assembly_graph.fastg") + " " +
                      o_p_extended + "assembly_graph.fastg")
            os.system("cp " + fastg_f + " " + o_p_extended + os.path.basename(fastg_f))
            os.system("cp " + tab_f + " " + o_p_extended + os.path.basename(tab_f))
            if not acyclic_allowed_in:
                if image_produced:
                    log_in.info("Please check the produced assembly image"
                                " or manually visualize " + o_p_extended + os.path.basename(fastg_f) +
                                " using Bandage to confirm the final result.")
                else:
                    log_in.info("Please visualize " + o_p_extended + os.path.basename(fastg_f) +
                                " using Bandage to confirm the final result.")
            log_in.info("Writing output finished.")

        disentangle_inside(input_graph=deepcopy(assembly_obj),
                           fastg_f=fastg_file, tab_f=tab_file, o_p=output, w_f=weight_factor, log_in=log_dis,
                           type_f=type_factor,
                           mode_in=mode, in_db_n=blast_db_base,
                           c_d=contamination_depth, c_s=contamination_similarity,
                           deg=degenerate, deg_dep=degenerate_depth, deg_sim=degenerate_similarity,
                           hard_c_t=hard_cov_threshold, min_s_f=min_sigma_factor, max_c_in=here_only_max_c,
                           max_s=expected_max_size, min_s=expected_min_size,
                           with_spades_scaffolds_in=with_spades_scaffolds,
                           acyclic_allowed_in=here_acyclic_allowed, verbose_in=here_verbose, in_temp_graph=temp_graph)

    # start
    kmer_values = sorted([int(kmer_d[1:])
                          for kmer_d in os.listdir(spades_output)
                          if os.path.isdir(os.path.join(spades_output, kmer_d))
                          and kmer_d.startswith("K")
                          and os.path.exists(os.path.join(spades_output, kmer_d, "assembly_graph.fastg"))],
                         reverse=True)
    kmer_values = [kmer_val for kmer_val in kmer_values if kmer_val > ignore_kmer_res]
    kmer_dirs = [os.path.join(spades_output, "K" + str(kmer_val)) for kmer_val in kmer_values]
    timeout_flag = "'--disentangle-time-limit'"
    export_succeeded = False
    path_prefix = os.path.join(out_base, organelle_prefix)
    graph_temp_file1 = path_prefix + "R1.temp.gfa" if options.keep_temp_files else None
    file_to_assembly_obj = {}
    for go_k, kmer_dir in enumerate(kmer_dirs):
        out_fastg = slim_out_fg[go_k]
        if out_fastg and os.path.getsize(out_fastg):
            try:
                """disentangle"""
                out_csv = out_fastg[:-5] + "csv"
                # if it is the first round (the largest kmer), copy the slimmed result to the main spades output
                # if go_k == 0:
                #     main_spades_folder = os.path.split(kmer_dir)[0]
                #     os.system("cp " + out_fastg + " " + main_spades_folder)
                #     os.system("cp " + out_csv + " " + main_spades_folder)
                log_handler.info("Parsing " + out_fastg)
                assembly_graph_obj = Assembly(out_fastg)
                log_handler.info("Loading and cleaning labels along " + out_fastg)
                assembly_graph_obj.parse_tab_file(
                    out_csv,
                    database_name=blast_db,
                    type_factor=options.disentangle_type_factor,
                    max_gene_gap=250,
                    max_cov_diff=options.disentangle_depth_factor,  # contamination_depth?
                    verbose=verbose,
                    log_handler=log_handler,
                    random_obj=random)
                file_to_assembly_obj[out_fastg] = assembly_graph_obj
                disentangle_assembly(assembly_obj=assembly_graph_obj,
                                     fastg_file=out_fastg,
                                     blast_db_base=blast_db,
                                     mode=organelle_type,
                                     tab_file=out_csv,
                                     output=path_prefix,
                                     weight_factor=100,
                                     type_factor=options.disentangle_type_factor,
                                     hard_cov_threshold=options.disentangle_depth_factor,
                                     contamination_depth=options.contamination_depth,
                                     contamination_similarity=options.contamination_similarity,
                                     degenerate=options.degenerate,
                                     degenerate_depth=options.degenerate_depth,
                                     degenerate_similarity=options.degenerate_similarity,
                                     expected_max_size=expected_maximum_size,
                                     expected_min_size=expected_minimum_size,
                                     here_only_max_c=True,
                                     here_acyclic_allowed=False, here_verbose=verbose, log_dis=log_handler,
                                     time_limit=options.disentangle_time_limit, timeout_flag_str=timeout_flag,
                                     temp_graph=graph_temp_file1)
            except ImportError as e:
                log_handler.error("Disentangling failed: " + str(e))
                return False
            except AttributeError as e:
                if verbose:
                    raise e
            except RuntimeError as e:
                if verbose:
                    log_handler.exception("")
                log_handler.info("Disentangling failed: RuntimeError: " + str(e).strip())
            except TimeoutError:
                log_handler.info("Disentangling timeout. (see " + timeout_flag + " for more)")
            except ProcessingGraphFailed as e:
                log_handler.info("Disentangling failed: " + str(e).strip())
            except Exception as e:
                log_handler.exception("")
                sys.exit()
            else:
                export_succeeded = True
                break

    if not export_succeeded and do_spades_scaffolding:
        graph_temp_file1s = path_prefix + "R1S.temp.gfa" if options.keep_temp_files else None
        largest_k_graph_f_exist = bool(slim_out_fg[0])
        if kmer_dirs and largest_k_graph_f_exist:
            out_fastg = slim_out_fg[0]
            if out_fastg and os.path.getsize(out_fastg):
                try:
                    """disentangle"""
                    out_csv = out_fastg[:-5] + "csv"
                    disentangle_assembly(assembly_obj=file_to_assembly_obj[out_fastg],
                                         fastg_file=out_fastg,
                                         blast_db_base=blast_db,
                                         mode=organelle_type,
                                         tab_file=out_csv,
                                         output=path_prefix,
                                         weight_factor=100,
                                         type_factor=options.disentangle_type_factor,
                                         hard_cov_threshold=options.disentangle_depth_factor,
                                         contamination_depth=options.contamination_depth,
                                         contamination_similarity=options.contamination_similarity,
                                         degenerate=options.degenerate, degenerate_depth=options.degenerate_depth,
                                         degenerate_similarity=options.degenerate_similarity,
                                         expected_max_size=expected_maximum_size,
                                         expected_min_size=expected_minimum_size,
                                         here_only_max_c=True, with_spades_scaffolds=True,
                                         here_acyclic_allowed=False, here_verbose=verbose, log_dis=log_handler,
                                         time_limit=options.disentangle_time_limit, timeout_flag_str=timeout_flag,
                                         temp_graph=graph_temp_file1s)
                except FileNotFoundError:
                    log_handler.warning("scaffolds.fasta and/or scaffolds.paths not found!")
                except ImportError as e:
                    log_handler.error("Disentangling failed: " + str(e))
                    return False
                except AttributeError as e:
                    if verbose:
                        raise e
                except RuntimeError as e:
                    if verbose:
                        log_handler.exception("")
                    log_handler.info("Disentangling failed: RuntimeError: " + str(e).strip())
                except TimeoutError:
                    log_handler.info("Disentangling timeout. (see " + timeout_flag + " for more)")
                except ProcessingGraphFailed as e:
                    log_handler.info("Disentangling failed: " + str(e).strip())
                except Exception as e:
                    log_handler.exception("")
                    sys.exit()
                else:
                    export_succeeded = True

    if not export_succeeded:
        graph_temp_file2 = path_prefix + "R2.temp.gfa" if options.keep_temp_files else None
        largest_k_graph_f_exist = bool(slim_out_fg[0])
        if kmer_dirs and largest_k_graph_f_exist:
            for go_k, kmer_dir in enumerate(kmer_dirs):
                out_fastg = slim_out_fg[go_k]
                if out_fastg and os.path.getsize(out_fastg):
                    try:
                        """disentangle the graph as scaffold(s)"""
                        out_fastg_list = sorted([os.path.join(kmer_dir, x)
                                                 for x in os.listdir(kmer_dir) if x.count(".fastg") == 2])
                        if out_fastg_list:
                            out_fastg = out_fastg_list[0]
                            out_csv = out_fastg[:-5] + "csv"
                            disentangle_assembly(assembly_obj=file_to_assembly_obj[out_fastg],
                                                 fastg_file=out_fastg,
                                                 blast_db_base=blast_db,
                                                 mode=organelle_type,
                                                 tab_file=out_csv,
                                                 output=path_prefix,
                                                 weight_factor=100,
                                                 type_factor=options.disentangle_type_factor,
                                                 here_verbose=verbose,
                                                 log_dis=log_handler,
                                                 hard_cov_threshold=options.disentangle_depth_factor * 0.6,
                                                 # TODO the adjustment should be changed if it's RNA data
                                                 contamination_depth=options.contamination_depth,
                                                 contamination_similarity=options.contamination_similarity,
                                                 degenerate=options.degenerate,
                                                 degenerate_depth=options.degenerate_depth,
                                                 degenerate_similarity=options.degenerate_similarity,
                                                 expected_max_size=expected_maximum_size,
                                                 expected_min_size=expected_minimum_size,
                                                 here_only_max_c=True, here_acyclic_allowed=True,
                                                 time_limit=3600, timeout_flag_str=timeout_flag,
                                                 temp_graph=graph_temp_file2)
                    except (ImportError, AttributeError) as e:
                        log_handler.error("Disentangling failed: " + str(e))
                        break
                    except RuntimeError as e:
                        if verbose:
                            log_handler.exception("")
                        log_handler.info("Disentangling failed: RuntimeError: " + str(e).strip())
                    except TimeoutError:
                        log_handler.info("Disentangling timeout. (see " + timeout_flag + " for more)")
                    except ProcessingGraphFailed as e:
                        log_handler.info("Disentangling failed: " + str(e).strip())
                    except Exception as e:
                        raise e
                    else:
                        export_succeeded = True
                        out_csv = out_fastg[:-5] + "csv"
                        log_handler.info("Please ...")
                        log_handler.info("load the graph file '" + os.path.basename(out_fastg) +
                                         "' in " + ",".join(["K" + str(k_val) for k_val in kmer_values]))
                        log_handler.info("load the CSV file '" + os.path.basename(out_csv) +
                                         "' in " + ",".join(["K" + str(k_val) for k_val in kmer_values]))
                        log_handler.info("visualize and confirm the incomplete result in Bandage.")
                        # log.info("-------------------------------------------------------")
                        log_handler.info("If the result is nearly complete, ")
                        log_handler.info("you can also adjust the arguments according to "
                                         "https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#what-should-i-do-with-incomplete-resultbroken-assembly-graph")
                        log_handler.info("If you have questions for us, "
                                         "please provide us with the get_org.log.txt file "
                                         "and the post-slimming graph in the format you like!")
                        # log.info("-------------------------------------------------------")
                        break
            if not export_succeeded:
                out_fastg = slim_out_fg[0]
                out_csv = out_fastg[:-5] + "csv"
                log_handler.info("Please ...")
                log_handler.info("load the graph file '" + os.path.basename(out_fastg) + ",assembly_graph.fastg" +
                                 "' in " + ",".join(["K" + str(k_val) for k_val in kmer_values]))
                log_handler.info("load the CSV file '" + os.path.basename(out_csv) +
                                 "' in " + ",".join(["K" + str(k_val) for k_val in kmer_values]))
                log_handler.info("visualize and export your result in Bandage.")
                log_handler.info("If you have questions for us, please provide us with the get_org.log.txt file "
                                 "and the post-slimming graph in the format you like!")
        else:
            # slim failed with unknown error
            log_handler.info("Please ...")
            log_handler.info("load the graph file: " + os.path.join(spades_output, 'assembly_graph.fastg'))
            log_handler.info("visualize and export your result in Bandage.")
            log_handler.info("If you have questions for us, please provide us with the get_org.log.txt file "
                             "and the post-slimming graph in the format you like!")
    return export_succeeded


def main():
    time0 = time.time()
    from GetOrganelleLib.versions import get_versions
    title = "GetOrganelle v" + str(get_versions()) + \
            "\n" \
            "\nget_organelle_from_reads.py assembles organelle genomes from genome skimming data." \
            "\nFind updates in https://github.com/Kinggerm/GetOrganelle and see README.md for more information." \
            "\n"
    options, log_handler, previous_attributes, run_slim, run_disentangle = \
        get_options(description=title, version=get_versions())

    resume = options.script_resume
    verb_log = options.verbose_log
    out_base = options.output_base
    echo_step = options.echo_step
    reads_files_to_drop = []
    # global word_size
    word_size = None
    mean_read_len = None
    mean_error_rate = None
    # all_bases = None
    low_quality_pattern = None
    max_read_len = None
    # max_extending_lens
    max_extending_lens = {inf}
    slim_extending_len = None
    phred_offset = options.phred_offset
    try:
        if options.fq_file_1 and options.fq_file_2:
            reads_paired = {'input': True, 'pair_out': bool}
            original_fq_files = [options.fq_file_1, options.fq_file_2] + \
                                [fastq_file for fastq_file in options.unpaired_fq_files]
            direction_according_to_user_input = [1, 2] + [1] * len(options.unpaired_fq_files)
        else:
            reads_paired = {'input': False, 'pair_out': False}
            original_fq_files = [fastq_file for fastq_file in options.unpaired_fq_files]
            direction_according_to_user_input = [1] * len(options.unpaired_fq_files)
        all_read_nums = [options.maximum_n_reads for foo in original_fq_files]
        other_spd_options = options.other_spades_options.split(' ')
        if '-o' in other_spd_options:
            which_out = other_spd_options.index('-o')
            spades_output = other_spd_options[which_out + 1]
            del other_spd_options[which_out: which_out + 2]
        else:
            spades_output = os.path.join(out_base, options.prefix + "extended_spades")
        if "--phred-offset" in other_spd_options:
            log_handler.warning("--spades-options '--phred-offset' was deprecated in GetOrganelle. ")
            which_po = other_spd_options.index("--phred-offset")
            del other_spd_options[which_po: which_po + 2]
        other_spd_options = ' '.join(other_spd_options)

        """ get reads """
        extended_files_exist = max(
            min([os.path.exists(
                str(os.path.join(out_base, options.prefix + "extended")) + '_' + str(i + 1) + '_unpaired.fq')
                 for i in range(2)] +
                [os.path.exists(str(os.path.join(out_base, options.prefix + "extended")) + '_' + str(i + 1) + '.fq')
                 for i in range(2, len(original_fq_files))]),
            min([os.path.exists(str(os.path.join(out_base, options.prefix + "extended")) + '_' + str(i + 1) + '.fq')
                 for i in range(len(original_fq_files))]))
        extended_fq_gz_exist = max(
            min([os.path.exists(
                str(os.path.join(out_base, options.prefix + "extended")) + '_' + str(i + 1) + '_unpaired.fq.tar.gz')
                 for i in range(2)] +
                [os.path.exists(str(os.path.join(out_base, options.prefix + "extended")) + '_' + str(i + 1) + '.fq.tar.gz')
                 for i in range(2, len(original_fq_files))]),
            min([os.path.exists(str(os.path.join(out_base, options.prefix + "extended")) + '_' + str(i + 1) + '.fq.tar.gz')
                 for i in range(len(original_fq_files))]))

        if resume:
            if "max_read_len" in previous_attributes and "mean_read_len" in previous_attributes and \
                    "phred_offset" in previous_attributes:
                try:
                    max_read_len = int(previous_attributes["max_read_len"])
                    mean_read_len = float(previous_attributes["mean_read_len"])
                    phred_offset = int(previous_attributes["phred_offset"])
                except ValueError:
                    resume = False
            else:
                resume = False
            if not resume and verb_log:
                log_handler.info("Previous attributes: max/mean read lengths/phred offset not found. "
                                 "Restart a new run.\n")
            try:
                word_size = int(previous_attributes["w"])
            except (KeyError, ValueError):
                if extended_files_exist or extended_fq_gz_exist:
                    if verb_log:
                        log_handler.info("Previous attributes: word size not found. Restart a new run.\n")
                    resume = False
                else:
                    pass

        if not (resume and (extended_files_exist or (extended_fq_gz_exist and phred_offset != -1))):
            anti_seed = options.anti_seed
            pre_grp = options.pre_grouped
            in_memory = options.index_in_memory
            log_handler.info("Pre-reading fastq ...")
            # using mapping to estimate maximum_n_reads when options.reduce_reads_for_cov != inf.
            all_read_nums = None
            if resume:
                try:
                    all_read_nums = [int(sub_num) for sub_num in previous_attributes["num_reads_1"].split("+")]
                except (KeyError, ValueError):
                    resume = False
                else:
                    try:
                        low_quality_pattern = "[" + previous_attributes["trim_chars"] + "]"
                        mean_error_rate = float(previous_attributes["mean_error_rate"])
                    except (KeyError, ValueError):
                        low_quality_pattern = "[]"
                        mean_error_rate = None
                    # all_bases = mean_read_len * sum(all_read_nums)
            if all_read_nums is None:
                if options.reduce_reads_for_cov != inf:
                    log_handler.info(
                        "Estimating reads to use ... "
                        "(to use all reads, set '--reduce-reads-for-coverage inf --max-reads inf')")
                    all_read_nums = estimate_maximum_n_reads_using_mapping(
                        twice_max_coverage=options.reduce_reads_for_cov * 2, check_dir=os.path.join(out_base, "check"),
                        original_fq_list=original_fq_files, reads_paired=reads_paired["input"],
                        maximum_n_reads_hard_bound=options.maximum_n_reads,
                        seed_files=options.seed_file, organelle_types=options.organelle_type,
                        in_customs=options.genes_fasta, ex_customs=options.exclude_genes,
                        target_genome_sizes=options.target_genome_size,
                        keep_temp=options.keep_temp_files, resume=options.script_resume,
                        other_spades_opts=other_spd_options,
                        which_blast=options.which_blast, which_spades=options.which_spades,
                        which_bowtie2=options.which_bowtie2, threads=options.threads,
                        random_seed=options.random_seed, verbose_log=options.verbose_log, log_handler=log_handler)
                    log_handler.info("Estimating reads to use finished.")
                else:
                    all_read_nums = [options.maximum_n_reads] * len(original_fq_files)

            if original_fq_files:
                for file_id, read_file in enumerate(original_fq_files):
                    # unzip fq files if needed
                    if read_file.endswith(".gz") or read_file.endswith(".zip"):
                        target_fq = os.path.join(out_base, str(file_id + 1) + "-" +
                                                 os.path.basename(read_file)) + ".fastq"
                        if not (os.path.exists(target_fq) and resume):
                            unzip(read_file, target_fq, 4 * all_read_nums[file_id],
                                  options.verbose_log, log_handler)
                    else:
                        target_fq = os.path.join(out_base, str(file_id + 1) + "-" +
                                                 os.path.basename(read_file))
                        if os.path.exists(target_fq) and os.path.islink(target_fq):
                            if os.path.realpath(target_fq) != os.path.realpath(os.path.join(os.getcwd(), read_file)):
                                log_handler.error("Existed symlink (%s -> %s) does not link to the input file (%s)!" %
                                                  (target_fq,
                                                   os.path.realpath(target_fq),
                                                   os.path.realpath(os.path.join(os.getcwd(), read_file))))
                                exit()
                        elif os.path.realpath(target_fq) == os.path.realpath(os.path.join(os.getcwd(), read_file)):
                            log_handler.error("Do not put original reads file(s) in the output directory!")
                            exit()
                        if not (os.path.exists(target_fq) and resume):
                            if all_read_nums[file_id] > READ_LINE_TO_INF:
                                # os.system("cp " + read_file + " " + target_fq + ".Temp")
                                # os.system("mv " + target_fq + ".Temp " + target_fq)
                                if os.path.exists(target_fq):
                                    os.remove(target_fq)
                                os.system("ln -s " + os.path.abspath(read_file) + " " + target_fq)
                            else:
                                os.system("head -n " + str(int(4 * all_read_nums[file_id])) + " " +
                                          read_file + " > " + target_fq + ".Temp")
                                os.system("mv " + target_fq + ".Temp " + target_fq)
                    if os.path.getsize(target_fq) == 0:
                        raise ValueError("Empty file " + target_fq)
                    original_fq_files[file_id] = target_fq
                    reads_files_to_drop.append(target_fq)
            if not resume:
                sampling_reads_for_quality = 10000
                # pre-reading fastq
                log_handler.info("Counting read qualities ...")
                low_quality_pattern, mean_error_rate, phred_offset = \
                    get_read_quality_info(original_fq_files, sampling_reads_for_quality, options.min_quality_score,
                                          log_handler, maximum_ignore_percent=options.maximum_ignore_percent)
                log_handler.info("Counting read lengths ...")
                mean_read_len, max_read_len, all_read_nums = get_read_len_mean_max_count(
                    original_fq_files, options.maximum_n_reads, n_process=1)
                    # original_fq_files, options.maximum_n_reads, n_process=options.threads)
                log_handler.info("Mean = " + str(round(mean_read_len, 1)) + " bp, maximum = " +
                                 str(max_read_len) + " bp.")
                log_handler.info("Reads used = " + "+".join([str(sub_num) for sub_num in all_read_nums]))
                log_handler.info("Pre-reading fastq finished.\n")
            else:
                log_handler.info("Pre-reading fastq skipped.\n")

            # reading seeds
            log_handler.info("Making seed reads ...")
            seed_fq_files = []
            seed_sam_files = []
            seed_fs_files = []
            for go_t, seed_f in enumerate(options.seed_file):
                seed_fq, seed_sam, new_seed_f = making_seed_reads_using_mapping(
                    seed_file=seed_f,
                    original_fq_files=original_fq_files,
                    out_base=out_base, resume=resume, verbose_log=verb_log, threads=options.threads,
                    random_seed=options.random_seed, organelle_type=options.organelle_type[go_t],
                    prefix=options.prefix, keep_temp=options.keep_temp_files,
                    bowtie2_other_options=options.bowtie2_options, which_bowtie2=options.which_bowtie2,
                    log_handler=log_handler)
                seed_fq_files.append(seed_fq)
                seed_sam_files.append(seed_sam)
                seed_fs_files.append(new_seed_f)
            anti_lines = get_anti_lines_using_mapping(
                anti_seed=anti_seed, seed_sam_files=seed_sam_files,
                original_fq_files=original_fq_files, out_base=out_base, resume=resume,
                verbose_log=verb_log, threads=options.threads,
                random_seed=options.random_seed, prefix=options.prefix,
                keep_temp=options.keep_temp_files, bowtie2_other_options=options.bowtie2_options,
                which_bowtie2=options.which_bowtie2, log_handler=log_handler)
            log_handler.info("Making seed reads finished.\n")

            log_handler.info("Checking seed reads and parameters ...")
            if not resume or options.word_size:
                word_size = options.word_size
            word_size, keep_seq_parts, mean_base_cov_values, max_extending_lens, all_read_limits = \
                check_parameters(word_size=word_size,
                                 original_fq_files=original_fq_files,
                                 seed_fs_files=seed_fs_files,
                                 seed_fq_files=seed_fq_files, seed_sam_files=seed_sam_files,
                                 organelle_types=options.organelle_type,
                                 in_custom_list=options.genes_fasta,
                                 ex_custom_list=options.exclude_genes,
                                 mean_error_rate=mean_error_rate,
                                 target_genome_sizes=options.target_genome_size,
                                 max_extending_len=options.max_extending_len, mean_read_len=mean_read_len,
                                 max_read_len=max_read_len, low_quality_pattern=low_quality_pattern,
                                 all_read_nums=all_read_nums, reduce_reads_for_cov=options.reduce_reads_for_cov,
                                 log_handler=log_handler,
                                 other_spades_opts=other_spd_options,
                                 which_spades=options.which_spades,
                                 which_blast=options.which_blast, which_bowtie2=options.which_bowtie2,
                                 wc_bc_ratio_constant=0.35, larger_auto_ws=options.larger_auto_ws,
                                 threads=options.threads, random_seed=options.random_seed,
                                 resume=resume, verbose_log=verb_log, zip_files=options.zip_files)
            log_handler.info("Checking seed reads and parameters finished.\n")

            # make read index
            log_handler.info("Making read index ...")
            fq_info_in_memory = make_read_index(original_fq_files, direction_according_to_user_input,
                                                all_read_limits, options.rm_duplicates, out_base, word_size,
                                                anti_lines, pre_grp, in_memory, anti_seed,
                                                keep_seq_parts=keep_seq_parts, low_quality=low_quality_pattern,
                                                resume=resume, echo_step=echo_step, log_handler=log_handler)
            len_indices = fq_info_in_memory[2]
            keep_seq_parts = fq_info_in_memory[3]
            if keep_seq_parts:
                log_handler.info("Reads are stored as fragments.")
            # pre-grouping if asked
            if pre_grp:
                preg_word_size = word_size if not options.pregroup_word_size else options.pregroup_word_size
                groups_of_lines, lines_with_dup, group_id_to_read_counts = \
                    pre_grouping(fastq_indices_in_memory=fq_info_in_memory, dupli_threshold=pre_grp, out_base=out_base,
                                 preg_word_size=preg_word_size, index_in_memory=in_memory, log_handler=log_handler)
            else:
                groups_of_lines = lines_with_dup = group_id_to_read_counts = None
            if not in_memory:
                fq_info_in_memory = None
            log_handler.info("Making read index finished.\n")

            # extending process
            log_handler.info("Extending ...")
            if set(max_extending_lens) == {inf}:
                accepted_rd_id = extending_no_lim(word_size=word_size, seed_file=seed_fq_files,
                                                  original_fq_files=original_fq_files, len_indices=len_indices,
                                                  pre_grouped=pre_grp, groups_of_duplicate_lines=groups_of_lines,
                                                  lines_with_duplicates=lines_with_dup,
                                                  fq_info_in_memory=fq_info_in_memory, output_base=out_base,
                                                  max_rounds=options.max_rounds,
                                                  min_rounds=1, fg_out_per_round=options.fg_out_per_round,
                                                  jump_step=options.jump_step, mesh_size=options.mesh_size,
                                                  verbose=verb_log, resume=resume,
                                                  all_read_limits=all_read_limits,
                                                  maximum_n_words=options.maximum_n_words,
                                                  keep_seq_parts=keep_seq_parts, low_qual_pattern=low_quality_pattern,
                                                  echo_step=echo_step, log_handler=log_handler)
            else:
                accepted_rd_id = extending_with_lim(word_size=word_size, seed_file=seed_fq_files,
                                                    original_fq_files=original_fq_files, len_indices=len_indices,
                                                    pre_grouped=pre_grp, groups_of_duplicate_lines=groups_of_lines,
                                                    lines_with_duplicates=lines_with_dup,
                                                    group_id_to_read_counts=group_id_to_read_counts,
                                                    fq_info_in_memory=fq_info_in_memory, output_base=out_base,
                                                    max_rounds=options.max_rounds,
                                                    extending_dist_limit=max_extending_lens,
                                                    min_rounds=1, fg_out_per_round=options.fg_out_per_round,
                                                    jump_step=options.jump_step, mesh_size=options.mesh_size,
                                                    verbose=verb_log, resume=resume,
                                                    all_read_limits=all_read_limits,
                                                    maximum_n_words=options.maximum_n_words,
                                                    keep_seq_parts=keep_seq_parts,
                                                    low_qual_pattern=low_quality_pattern,
                                                    mean_read_len=mean_read_len,
                                                    mean_base_cov=min([cov_v[0] for cov_v in mean_base_cov_values]),
                                                    echo_step=echo_step, log_handler=log_handler)
            mapped_read_ids = set()
            write_fq_results(original_fq_files, accepted_rd_id,
                             os.path.join(out_base, options.prefix + "extended"),
                             os.path.join(out_base, 'temp.indices.2'),
                             fq_info_in_memory, all_read_limits,
                             echo_step, verb_log, in_memory, log_handler, mapped_read_ids)
            del accepted_rd_id, fq_info_in_memory, groups_of_lines, \
                anti_lines, lines_with_dup

            if not options.keep_temp_files:
                try:
                    os.remove(os.path.join(out_base, 'temp.indices.1'))
                    os.remove(os.path.join(out_base, 'temp.indices.2'))
                except OSError:
                    pass

            log_handler.info("Extending finished.\n")
        else:
            log_handler.info("Extending ... skipped.\n")
        if reads_files_to_drop and not options.keep_temp_files:
            for rm_read_file in reads_files_to_drop:
                os.remove(rm_read_file)

        if reads_paired['input']:
            if not (resume and (min([os.path.exists(x) for x in
                                     (os.path.join(out_base, options.prefix + "extended_" + y + "_" + z + "paired.fq")
                                      for y in ('1', '2') for z in ('', 'un'))]) or extended_fq_gz_exist)):
                resume = False
                reads_paired['pair_out'] = separate_fq_by_pair(out_base, options.prefix, verb_log, log_handler)
                if reads_paired['pair_out'] and not options.keep_temp_files:
                    os.remove(os.path.join(out_base, options.prefix + "extended_1.fq"))
                    os.remove(os.path.join(out_base, options.prefix + "extended_2.fq"))
            else:
                log_handler.info("Separating extended fastq file ... skipped.\n")

        """ assembly """
        is_assembled = False
        if options.run_spades:
            if not (resume and os.path.exists(os.path.join(spades_output, 'assembly_graph.fastg'))):
                if extended_fq_gz_exist and not extended_files_exist:
                    files_to_unzip = [os.path.join(out_base, candidate)
                                      for candidate in os.listdir(out_base) if candidate.endswith(".fq.tar.gz")]
                    for file_to_u in files_to_unzip:
                        unzip(source=file_to_u, target=file_to_u[:-7], line_limit=inf)
                options.spades_kmer = check_kmers(options.spades_kmer, word_size, max_read_len, log_handler)
                log_handler.info("Assembling using SPAdes ...")
                if not executable("pigz -h"):
                    log_handler.warning("Compression after read correction will be skipped for lack of 'pigz'")
                    if "--disable-gzip-output" not in other_spd_options:
                        other_spd_options += " --disable-gzip-output"
                if phred_offset in (33, 64):
                    other_spd_options += " --phred-offset %i" % phred_offset
                is_assembled = assembly_with_spades(options.spades_kmer, spades_output, other_spd_options, out_base,
                                                    options.prefix, original_fq_files, reads_paired,
                                                    which_spades=options.which_spades, verbose_log=options.verbose_log,
                                                    resume=resume, threads=options.threads, log_handler=log_handler)
            else:
                is_assembled = True
                log_handler.info("Assembling using SPAdes ... skipped.\n")
        if options.zip_files:
            files_to_zip = [os.path.join(out_base, candidate)
                            for candidate in os.listdir(out_base) if candidate.endswith(".fq")]
            files_to_zip.extend([os.path.join(out_base, "seed", candidate)
                                 for candidate in os.listdir(os.path.join(out_base, "seed"))
                                 if candidate.endswith(".fq") or candidate.endswith(".sam")])
            if files_to_zip:
                log_handler.info("Compressing files ...")
                for file_to_z in files_to_zip:
                    zip_file(source=file_to_z, target=file_to_z + ".tar.gz", remove_source=True)
                log_handler.info("Compressing files finished.\n")

        """ export organelle """
        if is_assembled and run_slim:
            slim_stat_list, ignore_k = slim_spades_result(
                organelle_types=options.organelle_type, in_custom=options.genes_fasta, ex_custom=options.exclude_genes,
                spades_output=spades_output, ignore_kmer_res=options.ignore_kmer_res,
                max_slim_extending_len=slim_extending_len,
                verbose_log=options.verbose_log, log_handler=log_handler, threads=options.threads,
                which_blast=options.which_blast, resume=options.script_resume, keep_temp=options.keep_temp_files)
            slim_stat_codes = [s_code for s_code, fastg_out in slim_stat_list]
            slim_fastg_file = [fastg_out for s_code, fastg_out in slim_stat_list]
            options.ignore_kmer_res = ignore_k
            if set(slim_stat_codes) == {2}:
                log_handler.warning("No sequence hit our LabelDatabase!")
                log_handler.warning("This might due to unreasonable seed/parameter choices or a bug.")
                log_handler.info("Please open an issue at https://github.com/Kinggerm/GetOrganelle/issues "
                                 "with the get_org.log.txt file.\n")
            elif 0 in slim_stat_codes:
                log_handler.info("Slimming assembly graphs finished.\n")
                if run_disentangle:
                    organelle_type_prefix = []
                    duplicated_o_types = {o_type: 1
                                          for o_type in options.organelle_type
                                          if options.organelle_type.count(o_type) > 1}
                    for here_type in options.organelle_type:
                        if here_type in duplicated_o_types:
                            organelle_type_prefix.append(here_type + "-" + str(duplicated_o_types[here_type]))
                            duplicated_o_types[here_type] += 1
                        else:
                            organelle_type_prefix.append(here_type)
                    for go_t, sub_organelle_type in enumerate(options.organelle_type):
                        og_prefix = options.prefix + organelle_type_prefix[go_t]
                        graph_existed = bool([gfa_f for gfa_f in os.listdir(out_base)
                                              if gfa_f.startswith(og_prefix) and gfa_f.endswith(".path_sequence.gfa")])
                        fasta_existed = bool([fas_f for fas_f in os.listdir(out_base)
                                              if fas_f.startswith(og_prefix) and fas_f.endswith(".path_sequence.fasta")])
                        if resume and graph_existed and fasta_existed:
                            log_handler.info("Extracting " + sub_organelle_type + " from the assemblies ... skipped.\n")
                        else:
                            # log_handler.info("Parsing assembly graph and outputting ...")
                            log_handler.info("Extracting " + sub_organelle_type + " from the assemblies ...")
                            if options.genes_fasta:
                                db_base_name = remove_db_postfix(os.path.basename(options.genes_fasta[go_t]))
                            else:
                                db_base_name = sub_organelle_type
                            ext_res = extract_organelle_genome(out_base=out_base, spades_output=spades_output,
                                                               ignore_kmer_res=options.ignore_kmer_res,
                                                               slim_out_fg=slim_fastg_file, organelle_prefix=og_prefix,
                                                               organelle_type=sub_organelle_type,
                                                               blast_db=db_base_name,
                                                               read_len_for_log=mean_read_len,
                                                               verbose=options.verbose_log, log_handler=log_handler,
                                                               basic_prefix=options.prefix,
                                                               expected_minimum_size=options.expected_min_size[go_t],
                                                               expected_maximum_size=options.expected_max_size[go_t],
                                                               options=options,
                                                               do_spades_scaffolding=reads_paired["input"])
                            if ext_res:
                                log_handler.info("Extracting " + sub_organelle_type + " from the assemblies finished.\n")
                            else:
                                log_handler.info("Extracting " + sub_organelle_type + " from the assemblies failed.\n")
            else:
                log_handler.error("No valid assembly graph found!")
                log_handler.warning("This might due to a damaged dependency, "
                                    "to unreasonable seed/parameter choices, or to a bug.")
                log_handler.info("Please first search similar issues at "
                                 "https://github.com/Kinggerm/GetOrganelle/issues, "
                                 "then leave your message following the same issue, "
                                 "or open an issue at https://github.com/Kinggerm/GetOrganelle/issues if it is new, "
                                 "Please always attach the get_org.log.txt file.\n")
        log_handler = simple_log(log_handler, out_base, prefix=options.prefix + "get_org.")
        log_handler.info("\nTotal cost " + "%.2f" % (time.time() - time0) + " s")
        log_handler.info("Thank you!")
    except SystemExit:
        final_error_summary_log(log_handler, out_base, options.prefix, time0)
    except:
        log_handler.exception("")
        final_error_summary_log(log_handler, out_base, options.prefix, time0)
    logging.shutdown()


if __name__ == '__main__':
    main()

"""Copyright 2016 Jianjun Jin"""
