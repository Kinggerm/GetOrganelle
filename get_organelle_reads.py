#!/usr/bin/env python

import datetime
import sys
import os
from copy import deepcopy
try:
    from math import inf
except ImportError:
    inf = float("inf")
from optparse import OptionParser, OptionGroup
from VERSIONS import get_versions
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(PATH_OF_THIS_SCRIPT, ".."))
from Library.pipe_control_func import *
from Library.seq_parser import *
from Library.sam_parser import *
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
NOT_REF_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "Library", "NotationReference")
SEQ_REF_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "Library", "SeqReference")
import time
import random


MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    PYTHON_VERSION = "2.7+"
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    PYTHON_VERSION = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)

if PYTHON_VERSION == "2.7+":
    from commands import getstatusoutput
else:
    from subprocess import getstatusoutput
import subprocess

DEAD_CODE = {"2.7+": 32512, "3.5+": 127}[PYTHON_VERSION]
MAX_RATIO_RL_WS = 0.75
AUTO_MIN_WS = 49
GLOBAL_MIN_WS = 29


def get_options(descriptions, version):
    version = version
    usage = "\n###  Plant Plastome, Normal, 2*(1G raw data, 150 bp) reads\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -o plastome_output " \
            " -R 15 -k 21,45,65,85,105 -F plant_cp\n" \
            "###  Plant Mitochondria\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s mt_reference.fasta -o mitochondria_output " \
            " -R 30 -k 21,45,65,85,105 -F plant_mt\n" \
            "###  Plant Nuclear Ribosomal RNA (18S-ITS1-5.8S-ITS2-26S)\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -o nr_output " \
            " -R 7 -k 35,85,115 -P 0 -F plant_nr"
    parser = OptionParser(usage=usage, version=version, description=descriptions, add_help_option=False)
    # group 1
    group_inout = OptionGroup(parser, "IN-OUT OPTIONS", "Options on inputs and outputs")
    group_inout.add_option("-1", dest="fq_file_1",
                           help="Input file with forward paired-end reads to be added to the pool.")
    group_inout.add_option("-2", dest="fq_file_2",
                           help="Input file with reverse paired-end reads to be added to the pool.")
    group_inout.add_option("-u", dest="unpaired_fq_files",
                           help="Input file(s) with unpaired (single-end) reads to be added to the pool. "
                                "files could be comma-separated lists such as 'seq1,seq2'.")
    group_inout.add_option("-o", dest="output_base", help="Output directory. Overwriting files if directory exists.")
    group_inout.add_option("-s", dest="seed_file", default=os.path.join(SEQ_REF_PATH, "*.fasta"),
                           help="Seed sequence(s)/reference(s). Input fasta format file as initial seed. "
                                "A seed sequence in GetOrganelle is only used for identifying initial "
                                "organelle reads. The assembly process is purely de novo."
                                "Default: '%default' (* depends on the value followed with flag '-F')")
    group_inout.add_option("-a", dest="anti_seed",
                           help="Anti-seed(s). Not suggested unless what you really know what you are doing. "
                                "Input fasta format file as anti-seed, where the extension process "
                                "stop. Typically serves as excluding chloroplast reads when extending mitochondrial "
                                "reads, or the other way around. You should be cautious about using this option, "
                                "because if the anti-seed includes some word in the target but not in the seed, "
                                "the result would have gaps. For example, use the plant_mt and plant_cp from the same "
                                "plant-species as seed and anti-seed.")
    group_inout.add_option("--max-reads", dest="maximum_n_reads", type=float, default=1.5E7,
                           help="Maximum number of reads to be used per file. "
                                "Default: 1.5E7 (-F plant_cp/plant_nr/animal_mt/fungus_mt) "
                                "or 7.5E7 (-F plant_mt/anonym)")
    group_inout.add_option("--max-ignore-percent", dest="maximum_ignore_percent", type=float, default=0.01,
                           help="The maximum percent of bases to be ignore in extension, due to low quality. "
                                "Default: %default")
    group_inout.add_option("--min-quality-score", dest="min_quality_score", type=int, default=1,
                           help="Minimum quality score in extension. This value would be automatically increased "
                                "to prevent ignoring too much raw data (see --max-ignore-percent)."
                                "Default: %default ('\"' in Phred+33; 'A' in Phred+64/Solexa+64)")
    group_inout.add_option("--prefix", dest="prefix", default="",
                           help="Add extra prefix to resulting files under the output directory.")
    group_inout.add_option("--out-per-round", dest="fg_out_per_round", action="store_true", default=False,
                           help="Enable output per round. Choose to save memory but cost more time per round.")
    group_inout.add_option("--keep-temp", dest="keep_temp_files", action="store_true", default=False,
                           help="Choose to keep the running temp/index files.")
    # group 2
    group_scheme = OptionGroup(parser, "SCHEME OPTIONS", "Options on running schemes.")
    group_scheme.add_option("-F", dest="organelle_type",
                            help="This flag should be followed with plant_cp (chloroplast), plant_mt "
                                 "(plant mitochondrion), plant_nr (plant nuclear ribosomal RNA), animal_mt "
                                 "(animal mitochondrion), fungus_mt (fungus mitochondrion), "
                                 "or anonym (uncertain organelle genome type, with customized gene database "
                                 "('--genes'), which is suggested only when the above database is genetically distant "
                                 "from your sample).")
    group_scheme.add_option("--safe", dest="safe_strategy", default=False, action="store_true",
                            help="=\"-R 200 --max-reads 2E8 -J 1 -M 1 --min-quality-score -5 --max-ignore-percent 0 "
                                 "--max-n-words 4E9 -k 21,35,45,55,65,75,85,95,105,115,127\" "
                                 "This option is suggested for lowly-covered or nonhomogeneously-covered data (poor "
                                 "data). You can overwrite the value of a specific option listed above by adding "
                                 "that option along with the \"--safe\" flag. "
                                 "Intensive resources (time and data) would be utilized for a complete result. "
                                 "If you can already get the circular/complete result from other option combinations, "
                                 "you do not need to use this option, because the circular result by GetOrganelle is "
                                 "generally consistent under different options. If you still can not get a complete "
                                 "(say a circular plastome) under this option, you need to adjust '-w' and '-k' "
                                 "carefully to confirm whether it was the problem of the dataset. ")
    group_scheme.add_option("--fast", dest="fast_strategy", default=False, action="store_true",
                            help="=\"-R 10 --max-reads 5E6 -t 4 -J 5 -M 7 --max-n-words 3E7 --larger-auto-ws\" "
                                 "This option is suggested for homogeneously and highly covered data (very fine data). "
                                 "You can overwrite the value of a specific option listed above by adding "
                                 "that option along with the \"--fast\" flag. "
                                 "You could try GetOrganelle with this option for a list of samples and run a second "
                                 "time without this option for the rest with incomplete results. ")
    group_scheme.add_option("--memory-save", dest="memory_save", default=False, action="store_true",
                            help="=\"--out-per-round -P 0 --remove-duplicates 0\" "
                                 "You can overwrite the value of a specific option listed above by adding "
                                 "that option along with the \"--memory-save\" flag. A larger '-R' value is suggested "
                                 "when \"--memory-save\" is chosen.")
    group_scheme.add_option("--memory-unlimited", dest="memory_unlimited", default=False, action="store_true",
                            help="=\"-P 1E7 --index-in-memory --remove-duplicates 2E8 "
                                 "--min-quality-score -5 --max-ignore-percent 0\" "
                                 "You can overwrite the value of a specific option listed above by adding "
                                 "that option along with the \"--memory-unlimited\" flag. ")
    # group 3
    group_extending = OptionGroup(parser, "EXTENDING OPTIONS", "Options on the performance of extending process")
    group_extending.add_option("-w", dest="word_size", type=float,
                               help="Word size (W) for pre-grouping (if not assigned by '--pre-w') and extending "
                                    "process. This script would try to guess (auto-estimate) a proper W "
                                    "using an empirical function based on average read length, reads quality, "
                                    "target genome coverage, and other variables that might influence the extending "
                                    "process. You could assign the ratio (1>input>0) of W to "
                                    "read_length, based on which this script would estimate the W for you; "
                                    "or assign an absolute W value (read length>input>=35). Default: auto-estimated.")
    group_extending.add_option("--pre-w", dest="pregroup_word_size", type=float,
                               help="Word size (W) for pre-grouping. Used to reproduce result when word size is "
                                    "a certain value during pregrouping process and later changed during reads "
                                    "extending process. Similar to word size. Default: auto-estimated.")
    group_extending.add_option("-R", "--max-rounds", dest="max_rounds", type=int, default=1000,
                               help="Maximum number of running rounds (suggested: >=2). Default: %default.")
    group_extending.add_option("-r", "--min-rounds", dest="min_rounds", type=int, default=5,
                               help="Minimum number of running rounds (>=1). (NOT suggested) "
                                    "If 'auto_word_size_step' is enabled "
                                    "(see '--auto-wss' for more) and extending stopped before finishing the designed "
                                    "minimum rounds, automatically restart extending with a smaller word size. "
                                    "Default: %default")
    group_extending.add_option("--max-n-words", dest="maximum_n_words", type=float, default=4E8,
                               help="Maximum number of words to be used in total."
                                    "Default: 4E8 (-F plant_cp), 8E7 (-F plant_nr/fungus_mt), 4E7 (-F animal_mt), "
                                    "2E9 (-F plant_mt)")
    group_extending.add_option("-J", dest="jump_step", type=int, default=3,
                               help="The length of step for checking words in reads during extending process "
                                    "(integer >= 1). When you have reads of high quality, the larger the number is, "
                                    "the faster the extension will be, "
                                    "the more risk of missing reads in low coverage area. "
                                    "Choose 1 to choose the slowest but safest extension strategy. Default: %default")
    group_extending.add_option("-M", dest="mesh_size", type=int, default=2,
                               help="(Beta parameter) "
                                    "The length of step for building words from seeds during extending process "
                                    "(integer >= 1). When you have reads of high quality, the larger the number is, "
                                    "the faster the extension will be, "
                                    "the more risk of missing reads in low coverage area. "
                                    "Another usage of this mesh size is to choose a larger mesh size coupled with a "
                                    "smaller word size, which makes smaller word size feasible when memory is limited."
                                    "Choose 1 to choose the slowest but safest extension strategy. Default: %default")
    group_extending.add_option("--no-bowtie2", dest="utilize_mapping", action="store_false", default=True,
                               help="Choose to disable mapping process."
                                    "By default, this script would map original reads to pre-seed (or bowtie2 index) "
                                    "to acquire the initial seed. This requires bowtie2 to be installed "
                                    "(and thus Linux/Unix only). It is suggested to keep mapping as default "
                                    "when the seed is too diverse from target.")
    group_extending.add_option("--bowtie2-options", dest="bowtie2_options", default="--very-fast -t",
                               help="Bowtie2 options, such as '--ma 3 --mp 5,2 --very-fast -t'. Default: %default.")
    group_extending.add_option("--rough-pre-reading", dest="pre_reading", default=True, action="store_false",
                               help="Roughly pre-reading the read characteristics, in which case a user-defined "
                                    "word size should be given. Default: Fa")
    group_extending.add_option("--auto-wss", dest="auto_word_size_step", default=0, type=int,
                               help="(NOT suggested) The step of word size adjustment during extending process."
                                    "Use 0 to disable this automatic adjustment. Default: %default.")
    group_extending.add_option("--larger-auto-ws", dest="larger_auto_ws", default=False, action="store_true",
                               help="By using this flag, the empirical function for estimating W would tend to "
                                    "produce a relative larger W, which would speed up the matching in extending, "
                                    "reduce the memory cost in extending, but increase the risk of broken final "
                                    "graph. Suggested when the data is good with high and homogenous coverage.")
    group_extending.add_option("--soft-max-words", dest="soft_max_words", type=float, default=8E7,
                               help="Maximum number of words to be used to test the fitness of a word size value."
                                    "Default: 8E7 (-F plant_cp), 1.6E7 (-F plant_nr/animal_mt/fungus_mt), "
                                    "4E8 (-F plant_mt)")
    group_extending.add_option("--target-genome-size", dest="target_genome_size", default=130000, type=int,
                               help="Hypothetical value of target genome size. This is only used for estimating "
                                    "word size when no '-w word_size' is given. No need to be precise."
                                    "Default: 130000 (-F plant_cp) or 390000 (-F plant_mt/anonym) or "
                                    "65000 (-F fungus_mt) or 13000 (-F plant_nr/animal_mt) ")
    # group 4
    group_assembly = OptionGroup(parser, "ASSEMBLY OPTIONS", "These options are about the assembly and "
                                                             "graph disentangling")
    group_assembly.add_option("-k", dest="spades_kmer", default="21,55,85,115",
                              help="SPAdes kmer settings. Use the same format as in SPAdes. illegal kmer values "
                                   "would be automatically discarded by GetOrganelle. "
                                   "Default: %default")
    group_assembly.add_option("--spades-options", dest="other_spades_options", default="",
                              help="Other SPAdes options. Use double quotation marks to include all the arguments "
                                   "and parameters.")
    group_assembly.add_option("--no-spades", dest="run_spades", action="store_false", default=True,
                              help="Disable SPAdes.")
    group_assembly.add_option("--no-gradient-k", dest="auto_gradient_k", action="store_false", default=True,
                              help="Disable adding an extra set of kmer values according to word size.")
    group_assembly.add_option("--genes", dest="genes_fasta",
                              help="Customized database (fasta format) containing ONE set of protein coding genes "
                                   "and ribosomal RNAs from ONE reference genome. "
                                   "This database is only used for identifying/tagging target organelle contigs, "
                                   "while the gene order of output, meaning the scaffolding/disentangling process "
                                   "is purely decided by the topology of de Bruijn graph and the coverage info. "
                                   "Used when '-F anonym'.")
    group_assembly.add_option("--disentangle-df", dest="disentangle_depth_factor", default=10.0, type=float,
                              help="Depth factor for differentiate genome type of contigs. "
                                   "The genome type of contigs are determined by blast. "
                                   "You can also make the blast index by your self and add those index to '" +
                                   NOT_REF_PATH + "'. Default: %default")
    group_assembly.add_option("--contamination-depth", dest="contamination_depth", default=3., type=float,
                              help="Depth factor for confirming contamination in parallel contigs. Default: %default")
    group_assembly.add_option("--contamination-similarity", dest="contamination_similarity", default=0.9, type=float,
                              help="Similarity threshold for confirming contaminating contigs. Default: %default")
    group_assembly.add_option("--no-degenerate", dest="degenerate", default=True, action="store_false",
                              help="Disable making consensus from parallel contig based on nucleotide degenerate table.")
    group_assembly.add_option("--degenerate-depth", dest="degenerate_depth", default=1.5, type=float,
                              help="Depth factor for confirming parallel contigs. Default: %default")
    group_assembly.add_option("--degenerate-similarity", dest="degenerate_similarity", default=0.98, type=float,
                              help="Similarity threshold for confirming parallel contigs. Default: %default")
    group_assembly.add_option("--disentangle-time-limit", dest="disentangle_time_limit", default=600, type=int,
                              help="Time limit (second) for each try of disentangling a graph file as a circular "
                                   "genome. Disentangling a graph as contigs is not limited. Default: %default")
    group_assembly.add_option("--expected-max-size", dest="expected_max_size", default=200000, type=int,
                              help="Expected maximum target genome size. Default: 200000 (-F plant_cp/fungus_mt), "
                                   "50000 (-F plant_nr/animal_mt/fungus_mt), 600000 (-F plant_mt)")
    group_assembly.add_option("--expected-min-size", dest="expected_min_size", default=10000, type=int,
                              help="Expected mininum target genome size. Default: %default")
    # group 5
    group_computational = OptionGroup(parser, "ADDITIONAL OPTIONS", "")
    group_computational.add_option("-t", dest="threads", type=int, default=1,
                                   help="Maximum threads to use.")
    group_computational.add_option("-P", dest="pre_grouped", type=float, default=2E5,
                                   help="The maximum number (integer) of high-covered reads to be pre-grouped "
                                        "before extending process. pre_grouping is suggested when the whole genome "
                                        "coverage is shallow but the organ genome coverage is deep. "
                                        "The default value is 2E5. "
                                        "For personal computer with 8G memory, we suggest no more than 3E5. "
                                        "A larger number (ex. 6E5) would run faster but exhaust memory "
                                        "in the first few minutes. Choose 0 to disable this process.")
    group_computational.add_option("--continue", dest="script_resume", default=False, action="store_true",
                                   help="Several check points based on files produced, rather than log, "
                                        "so keep in mind that this script will not detect the difference "
                                        "between this input parameters and the previous ones.")
    group_computational.add_option("--index-in-memory", dest="index_in_memory", action="store_true", default=False,
                                   help="Keep index in memory. Choose save index in memory than in disk.")
    group_computational.add_option("--remove-duplicates", dest="rm_duplicates", default=1E7, type=float,
                                   help="By default this script use unique reads to extend. Choose the number of "
                                        "duplicates (integer) to be saved in memory. A larger number (ex. 2E7) would "
                                        "run faster but exhaust memory in the first few minutes. "
                                        "Choose 0 to disable this process. "
                                        "Note that whether choose or not will not disable "
                                        "the calling of replicate reads. Default: %default.")
    group_computational.add_option("--flush-frequency", dest="echo_frequency", default=54321, type=int,
                                   help="Flush frequency for presenting progress. Default: %default")
    group_computational.add_option("--random-seed", dest="random_seed", default=12345, type=int,
                                   help="Random seed (only for disentangling at this moment). Default: %default")
    group_computational.add_option("--verbose", dest="verbose_log", action="store_true", default=False,
                                   help="Verbose output. Choose to enable verbose running log.")

    parser.add_option("-h", dest="simple_help", default=False, action="store_true",
                      help="print brief introduction for frequently-used options.")
    parser.add_option("--help", dest="verbose_help", default=False, action="store_true",
                      help="print verbose introduction for all options.")
    if "--help" in sys.argv:
        parser.add_option_group(group_inout)
        parser.add_option_group(group_scheme)
        parser.add_option_group(group_extending)
        parser.add_option_group(group_assembly)
        parser.add_option_group(group_computational)
        parser.print_help()
        exit()
    elif "-h" in sys.argv:
        for not_often_used in ("-a", "--max-reads", "--max-ignore-percent",
                               "--min-quality-score", "--prefix", "--out-per-round", "--keep-temp",
                               "--memory-save", "--memory-unlimited", "--pre-w", "-r", "--max-n-words",
                               "-J", "-M", "--no-bowtie2", "--bowtie2-options", "--rough-pre-reading",
                               "--auto-wss", "--larger-auto-ws",
                               "--soft-max-words", "--target-genome-size", "--spades-options", "--no-spades",
                               "--no-gradient-k", "--genes", "--disentangle-df", "--contamination-depth",
                               "--contamination-similarity", "--no-degenerate",
                               "--degenerate-depth", "--degenerate-similarity", "--disentangle-time-limit",
                               "--expected-max-size", "--expected-min-size", "--continue", "--index-in-memory",
                               "--remove-duplicates", "--flush-frequency", "--verbose"):
            parser.remove_option(not_often_used)
        parser.remove_option("-1")
        parser.add_option("-1", dest="fq_file_1", help="Input file with forward paired-end reads.")
        parser.remove_option("-2")
        parser.add_option("-2", dest="fq_file_2", help="Input file with reverse paired-end reads.")
        parser.remove_option("-u")
        parser.add_option("-u", dest="unpaired_fq_files", help="Input file(s) with unpaired (single-end) reads. ")
        parser.remove_option("-o")
        parser.add_option("-o", dest="output_base", help="Output directory.")
        parser.remove_option("-s")
        parser.add_option("-s", dest="seed_file", help="Reference. Input fasta format file as initial seed. "
                                                       "Default: " + os.path.join(SEQ_REF_PATH, "*.fasta"))
        parser.remove_option("-w")
        parser.add_option("-w", dest="word_size", help="Word size (W) for extension. Default: auto-estimated")
        parser.remove_option("-R")
        parser.add_option("-R", dest="max_rounds", help="Maximum running rounds (suggested: >=2). Default: unlimited")
        parser.remove_option("-F")
        parser.add_option("-F", dest="organelle_type", default="plant_cp",
                          help="Target organelle genome type: "
                               "plant_cp/plant_mt/plant_nr/animal_mt/fungus_mt/anonym. Default: %default")
        parser.remove_option("--safe")
        parser.add_option("--safe", dest="safe_strategy",
                          help="=\"-R 200 --max-reads 2E8 -J 1 -M 1 --min-quality-score -5 --max-ignore-percent 0 "
                               "--max-n-words 4E9 -k 21,35,45,55,65,75,85,95,105,115,127 "
                               "--disentangle-time-limit 3600\"")
        parser.remove_option("--fast")
        parser.add_option("--fast", dest="fast_strategy",
                          help="=\"-R 10 --max-reads 5E6 -t 4 -J 5 -M 7 --max-words 3E7 --larger-auto-ws "
                               "--disentangle-time-limit 180 --no-gradient-k\"")
        parser.remove_option("-k")
        parser.add_option("-k", dest="spades_kmer", default="21,55,85,115", help="SPAdes kmer settings. Default: %default")
        parser.remove_option("-t")
        parser.add_option("-t", dest="threads", type=int, default=1, help="Maximum threads to use. Default: %default")
        parser.remove_option("-P")
        parser.add_option("-P", dest="pre_grouped", default=int(2E5), help="Pre-grouping value. Default: %default")
        parser.print_help()
        sys.stdout.write("\n")
        exit()
    else:
        parser.add_option_group(group_inout)
        parser.add_option_group(group_scheme)
        parser.add_option_group(group_extending)
        parser.add_option_group(group_assembly)
        parser.add_option_group(group_computational)

    try:
        (options, args) = parser.parse_args()
    except Exception as e:
        sys.stdout.write("\n############################################################################" + str(e))
        sys.stdout.write("\n\"-h\" for more usage")
        exit()
    else:
        if not (options.seed_file # or options.bowtie2_seed)
                and ((options.fq_file_1 and options.fq_file_2) or options.unpaired_fq_files)
                and options.output_base and options.organelle_type):
            sys.stdout.write("\n############################################################################"
                             "\nERROR: Insufficient arguments!\nUsage:")
            sys.stdout.write(usage + "\n\n")
            exit()
        if int(bool(options.fq_file_1)) + int(bool(options.fq_file_2)) == 1:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: unbalanced paired reads!\n\n")
            exit()
        if options.organelle_type not in {"plant_cp", "plant_mt", "plant_nr", "animal_mt", "fungus_mt", "anonym"}:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: \"-F\" MUST be one of 'plant_cp', 'plant_mt', 'plant_nr', "
                             "'animal_mt', 'fungus_mt', 'anonym'!\n\n")
            exit()
        if "*" in options.seed_file:
            if options.organelle_type == "anonym":
                sys.stdout.write("\n############################################################################"
                                 "\nERROR: \"-s\" and \"--genes\" must be specified when \"-F anonym\"!\n\n")
                exit()
            elif options.organelle_type == "animal_mt":
                sys.stdout.write("\n############################################################################"
                                 "\nERROR: Currently \"-s\" and \"--genes\" must be specified when "
                                 "\"-F animal_mt\"!\n\n")
                exit()
            else:
                options.seed_file = options.seed_file.replace("*", options.organelle_type)
        if options.organelle_type in {"anonym", "animal_mt"} and not options.genes_fasta:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: \"-s\" and \"--genes\" must be specified when \"-F anonym\"!\n\n")
            exit()
        for check_file in (options.fq_file_1, options.fq_file_2, options.seed_file, options.anti_seed):
            if check_file:
                if not os.path.exists(check_file):
                    sys.stdout.write("\n############################################################################"
                                     "\nERROR: " + check_file + " not found!\n\n")
                    exit()
        if options.unpaired_fq_files:
            options.unpaired_fq_files = options.unpaired_fq_files.split(",")
            for fastq_file in options.unpaired_fq_files:
                if not os.path.exists(fastq_file):
                    sys.stdout.write("\n############################################################################"
                                     "\nERROR: " + fastq_file + " not found!\n\n")
                    exit()
        else:
            options.unpaired_fq_files = []
        if options.jump_step < 1:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: Jump step MUST be an integer that >= 1")
            exit()
        if options.mesh_size < 1:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: Mesh size MUST be an integer that >= 1")
            exit()
        if options.fq_file_1 == options.fq_file_2 and options.fq_file_1:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: 1st fastq file is the same with 2nd fastq file!")
            exit()
        if options.safe_strategy and options.fast_strategy:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: \"--safe\" and \"--fast\" are not compatible!")
        if options.memory_save and options.memory_unlimited:
            sys.stdout.write("\n############################################################################"
                             "\nERROR: \"--memory-save\" and \"--memory-unlimited\" are not compatible!")
        options.prefix = os.path.basename(options.prefix)
        if options.script_resume and os.path.isdir(options.output_base):
            previous_attributes = LogInfo(options.output_base, options.prefix)
        else:
            previous_attributes = None
            options.script_resume = False
        if not os.path.isdir(options.output_base):
            os.mkdir(options.output_base)
            options.script_resume = False
        log = simple_log(logging.getLogger(), options.output_base, options.prefix + "get_org.")
        log.info("")
        log.info(descriptions)
        log.info("Python " + sys.version)
        log.info(" ".join(sys.argv) + "\n")
        log = timed_log(log, options.output_base, options.prefix + "get_org.")
        if options.word_size is None:
            # if options.auto_word_size_step:
            #     pass
            # else:
            #     log.error("Do not turn off auto_word_size_step if no input initial word size!")
            #     exit()
            pass
        elif 0 < options.word_size < 1:
            pass
        elif options.word_size >= GLOBAL_MIN_WS:
            options.word_size = int(options.word_size)
        else:
            log.error("Illegal word size (\"-w\") value!")
            exit()
        
        if options.pregroup_word_size:
            if 0 < options.pregroup_word_size < 1:
                pass
            elif options.pregroup_word_size >= GLOBAL_MIN_WS:
                options.pregroup_word_size = int(options.pregroup_word_size)
            else:
                log.error("Illegal word size (\"--pre-w\") value!")
                exit()

        if options.safe_strategy:
            if "-R" not in sys.argv and "--max-rounds" not in sys.argv:
                options.max_rounds = 200
            if "--max-reads" not in sys.argv:
                options.maximum_n_reads = 2E8
            if "-J" not in sys.argv:
                options.jump_step = 1
            if "-M" not in sys.argv:
                options.mesh_size = 1
            if "--min-quality-score" not in sys.argv:
                options.min_quality_score = -5
            if "--max-ignore-percent" not in sys.argv:
                options.maximum_ignore_percent = 0
            # if "--auto-wss" not in sys.argv:
            #     options.auto_word_size_step = 2
            if "--max-n-words" not in sys.argv:
                options.maximum_n_words = 4E9
            if "-k" not in sys.argv:
                options.spades_kmer = "21,35,45,55,65,75,85,95,105,115,127"
            if "--disentangle-time-limit" not in sys.argv:
                options.disentangle_time_limit = 3600

        if options.fast_strategy:
            if "-R" not in sys.argv and "--max-rounds" not in sys.argv:
                options.max_rounds = 10
            if "--max-reads" not in sys.argv:
                options.maximum_n_reads = 5E6
            if "-t" not in sys.argv:
                options.threads = 4
            if "-J" not in sys.argv:
                options.jump_step = 5
            if "-M" not in sys.argv:
                options.mesh_size = 7
            if "--max-n-words" not in sys.argv:
                options.maximum_n_words = 3E7
            options.larger_auto_ws = True
            options.auto_gradient_k = False
            if "--disentangle-time-limit" not in sys.argv:
                options.disentangle_time_limit = 180

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
        
        options.auto_word_size_step = abs(options.auto_word_size_step)

        if options.organelle_type in ("animal_mt", "fungus_mt"):
            log.warning("Currently GetOrganelle had been tested on limited animal & fungus samples!")
            log.warning("The default seed/references for animal & fungus are for temporary usage!")
            log.warning("No guarantee for success rate in animal & fungus samples!\n")
            if options.organelle_type == "fungus_mt":
                log.warning("Too improve this, a customized close-related reference (-s reference.fasta) with its "
                            "protein coding & ribosomal genes extracted (--genes reference.genes.fasta) is "
                            "highly suggested for your own animal & fungus samples!\n")

        # using the default
        if "--max-reads" not in sys.argv:
            if options.organelle_type == ("plant_mt", "anonym"):
                options.maximum_n_reads *= 5
        if "--max-n-words" not in sys.argv:
            if options.organelle_type in ("plant_mt", "anonym"):
                options.maximum_n_words *= 5
            elif options.organelle_type in ("plant_nr", "fungus_mt"):
                options.maximum_n_words /= 5
            elif options.organelle_type == "animal_mt":
                options.maximum_n_words /= 10
        if "--soft-max-words" not in sys.argv:
            if options.organelle_type in ("plant_mt", "anonym"):
                options.soft_max_words *= 5
            elif options.organelle_type in ("plant_nr", "animal_mt", "fungus_mt"):
                options.soft_max_words /= 5
        if "--target-genome-size" not in sys.argv:
            if options.organelle_type == "plant_mt":
                options.target_genome_size *= 3
            elif options.organelle_type == "fungus_mt":
                options.target_genome_size /= 2.
            elif options.organelle_type in ("plant_nr", "animal_mt"):
                options.target_genome_size /= 10.
            elif options.organelle_type == "anonym":
                ref_seqs = read_fasta(options.genes_fasta)[1]
                options.target_genome_size = 2 * sum([len(this_seq) for this_seq in ref_seqs])
                log.info("Setting '--target-genome-size " + str(options.target_genome_size) +
                         "' for estimating the word size value.")
        if "--expected-max-size" not in sys.argv:
            if options.organelle_type == "plant_mt":
                options.expected_max_size *= 3
            elif options.organelle_type in ("plant_nr", "animal_mt", "fungus_mt"):
                options.expected_max_size /= 4

        if options.organelle_type in ("fungus_mt", "animal_mt", "anonym"):
            global MAX_RATIO_RL_WS
            MAX_RATIO_RL_WS = 0.8

        if not options.pre_reading and not options.word_size:
            log.error("When '--rough-pre-reading' was chosen, a user defined word size value ('-w') must be given!")
            exit()
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
            if not executable("bowtie2-build"):
                options.utilize_mapping = False
                if options.seed_file:
                    log.warning('bowtie2-build not in the path! Take the seed file as initial seed.')
                else:
                    log.error('bowtie2-build not in the path!')
                    exit()
                if options.anti_seed:
                    log.warning('bowtie2-build not in the path! Take the anti-seed file as initial anti-seed.')
                else:
                    log.warning('bowtie2-build not in the path! Anti-seed disabled!')
            if not executable("bowtie2-build-l"):
                options.utilize_mapping = False
                if options.seed_file:
                    log.warning('bowtie2-build-l not in the path! Take the seed file as initial seed.')
                else:
                    log.error('bowtie2-build-l not in the path!')
                    exit()
                if options.anti_seed:
                    log.warning('bowtie2-build-l not in the path! Take the anti-seed file as initial anti-seed.')
                else:
                    log.warning('bowtie2-build-l not in the path! Anti-seed disabled!')
        if options.run_spades:
            if not executable("spades.py -h"):
                log.warning("spades.py not found in the path. Only get the reads and skip assembly.")
                options.run_spades = False
            if not executable("blastn"):
                log.warning("blastn not found in the path!")
            if options.organelle_type == "anonym" and not executable("makeblastdb"):
                log.warning("makeblastdb not found in the path!")
        options.rm_duplicates = int(options.rm_duplicates)
        options.pre_grouped = int(options.pre_grouped)
        if not options.rm_duplicates and options.pre_grouped:
            log.warning("removing duplicates was inactive, so that the pre-grouping was disabled.")
            options.pre_grouped = False
        if options.min_rounds < 1:
            log.warning("illegal minimum rounds! Set to default: 1.")
            options.min_rounds = 1
        if options.max_rounds and options.max_rounds < 1:
            log.warning("illegal maximum rounds! Set to infinite")
            options.max_rounds = inf
        random.seed(options.random_seed)
        try:
            import numpy as np
        except ImportError:
            pass
        else:
            np.random.seed(options.random_seed)
        return options, log, previous_attributes


def combination_res(all_choices_num, chosen_num):
    res = 1
    for ch_n in range(chosen_num, 0, -1):
        res *= all_choices_num - ch_n + 1
        res /= ch_n
    return res


def trans_word_cov(word_cov, base_cov, mean_base_error_rate, read_length):
    wrong_words_percent = 0
    for error_site_num in range(1, int(min(read_length * mean_base_error_rate * 10, read_length))):
        prob_of_err_site_num = combination_res(read_length, error_site_num) \
                               * (mean_base_error_rate ** error_site_num) \
                               * ((1 - mean_base_error_rate) ** (read_length - error_site_num))
        wrong_words_percent += (1 - 2 ** (-error_site_num)) * prob_of_err_site_num
    # if word size < read_len/2, wrong words percent decreases
    increase_word_cov = word_cov / (1 - wrong_words_percent) - word_cov
    if word_cov > 0.5 * base_cov:
        word_cov += increase_word_cov ** 0.34
    elif word_cov + increase_word_cov > 0.5 * base_cov:
        word_cov = 0.5 * base_cov + (word_cov + increase_word_cov - 0.5 * base_cov) ** 0.34
    else:
        word_cov += increase_word_cov
    return word_cov


def estimate_word_size(base_cov_values, read_length, target_size, mean_error_rate=0.015, log_handler=None,
                       max_discontinuous_prob=0.01, min_word_size=AUTO_MIN_WS, max_effective_word_cov=60,
                       wc_bc_ratio_constant=0.35):
    echo_problem = False
    all_word_sizes = []
    for base_cov in base_cov_values:
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
        #
        wc_bc_ratio_max = 1 - (min_word_size - 1) / read_length
        if base_cov * wc_bc_ratio_max < min_word_cov:
            min_word_cov = base_cov * wc_bc_ratio_max
            echo_problem = True
        word_cov = max(min_word_cov, word_cov)
        word_cov = trans_word_cov(word_cov, base_cov, mean_error_rate / 2., read_length)
        # 1. relationship between kmer coverage and base coverage, k_cov = base_cov * (read_len - k_len + 1) / read_len
        estimated_word_size = int(read_length * (1 - word_cov / base_cov)) + 1
        estimated_word_size = min(int(read_length * MAX_RATIO_RL_WS), max(min_word_size, estimated_word_size))
        all_word_sizes.append(int(round(estimated_word_size, 0)))
    if echo_problem:
        if log_handler:
            log_handler.warning("Guessing that you are using too few data for assembly!")
            log_handler.warning("GetOrganelle is still trying ...")
        else:
            sys.stdout.write("Guessing that you are using too few data for assembly!\n")
            sys.stdout.write("GetOrganelle is still trying ...\n")
    return all_word_sizes


def calculate_word_size_according_to_ratio(word_size_ratio, mean_read_len, log):
    if word_size_ratio < 1:
        new_word_size = int(round(word_size_ratio * mean_read_len, 0))
        if new_word_size < GLOBAL_MIN_WS:
            new_word_size = GLOBAL_MIN_WS
            log.warning("Too small ratio " + str(new_word_size) + ", setting '-w " + str(GLOBAL_MIN_WS) + "'")
        else:
            log.info("Setting '-w " + str(new_word_size) + "'")
        return new_word_size
    else:
        max_ws = int(round(mean_read_len * 0.9))
        if word_size_ratio > max_ws:
            word_size_ratio = max_ws
            log.warning("Too large word size for mean read length " + str(mean_read_len) +
                        ", setting '-w " + str(word_size_ratio) + "'")
        return word_size_ratio


def check_parameters(word_size, out_base, utilize_mapping, maximum_n_reads, original_fq_files,
                     organelle_type, all_bases, auto_word_size_step, mean_error_rate,
                     target_genome_size, mean_read_len, max_read_len, low_quality_pattern,
                     log, wc_bc_ratio_constant=0.35, larger_auto_ws=False):
    if utilize_mapping:
        coverage_info = get_coverage_from_sam_fast(os.path.join(out_base, "seed_bowtie.sam"))
        all_coverages = [coverage_info[ref][pos] for ref in coverage_info for pos in coverage_info[ref]
                         if coverage_info[ref][pos] > 2]
        if not all_coverages:
            all_coverages = [coverage_info[ref][pos] for ref in coverage_info for pos in coverage_info[ref]]
            if not all_coverages:
                if log:
                    log.error("No seed reads found!")
                    log.error("Please check your raw data or change your reference!")
                exit()
        base_cov_values = get_cover_range(all_coverages, guessing_percent=0.07)  # top 0.07 from mapped reads
        log.info("Estimated " + organelle_type + " base-coverage = " + str(base_cov_values[1]))
        if base_cov_values[0] < 200 and organelle_type in {"animal_mt", "fungus_mt", "anonym"}:
            log.info("Re-estimating " + organelle_type + " base-coverage using word frequency counting ...")
            counting_word_size = min(49, 2 * int(mean_read_len * 0.35) - 1)
            smp_percent = min(1., max(0.2, 2E7/all_bases))  # at least 20M raw data should be used
            words_dict = chop_seqs_as_empty_dict(
                fq_seq_simple_generator(os.path.join(out_base, "Initial.mapped.fq")), word_size=counting_word_size)
            counting_words(fq_seq_simple_generator(original_fq_files,
                                                   max_n_reads=all_bases/mean_read_len * smp_percent),
                           words_initial_dict=words_dict, word_size=counting_word_size)
            all_word_coverages = []
            for this_word in list(words_dict):
                if this_word in words_dict:
                    reverse_word = complementary_seq(this_word)
                    if reverse_word == this_word:
                        all_word_coverages.append(words_dict.pop(this_word))
                    else:
                        all_word_coverages.append(words_dict.pop(this_word) + words_dict.pop(reverse_word))
            word_cov_values = get_cover_range(all_word_coverages, guessing_percent=0.5)  # mean from kmer counting
            base_cov_values = [round(this_word_cov * mean_read_len / (mean_read_len - counting_word_size + 1)
                                     / smp_percent, 2)
                               for this_word_cov in word_cov_values]
            log.info("Estimated " + organelle_type + " base-coverage = " + str(base_cov_values[1]))
    else:
        organelle_base_percent = 0.05
        guessing_base_cov = organelle_base_percent * all_bases / target_genome_size
        base_cov_values = [guessing_base_cov * 0.8, guessing_base_cov, guessing_base_cov * 1.2]
        log.info("Guessing " + organelle_type + " base-coverage = " + str(base_cov_values[1]))
    if word_size is None:
        if larger_auto_ws:
            min_word_size, word_size, max_word_size = estimate_word_size(
                base_cov_values=base_cov_values, read_length=mean_read_len, target_size=target_genome_size,
                max_discontinuous_prob=0.05, min_word_size=70,
                mean_error_rate=mean_error_rate, log_handler=log, wc_bc_ratio_constant=wc_bc_ratio_constant - 0.03)
        else:
            # standard dev
            min_word_size, word_size, max_word_size = estimate_word_size(
                base_cov_values=base_cov_values, read_length=mean_read_len, target_size=target_genome_size,
                max_discontinuous_prob=0.01, min_word_size=AUTO_MIN_WS,
                mean_error_rate=mean_error_rate, log_handler=log, wc_bc_ratio_constant=wc_bc_ratio_constant)
        log.info("Setting '-w " + str(word_size) + "'")
        log.info("The automatically-estimated word size does not ensure the best choice.")
        log.info("If the result graph is not a circular organelle genome, ")
        log.info("you could adjust the word size for another new run.")
    elif word_size < 1:
        new_word_size = int(round(word_size * mean_read_len, 0))
        if new_word_size < GLOBAL_MIN_WS:
            word_size = GLOBAL_MIN_WS
            log.warning("Too small ratio " + str(word_size) + ", setting '-w " + str(GLOBAL_MIN_WS) + "'")
        else:
            word_size = new_word_size
            log.info("Setting '-w " + str(word_size) + "'")
        min_word_size = max(GLOBAL_MIN_WS, word_size - auto_word_size_step * 4)
        max_word_size = min(word_size + auto_word_size_step * 4, int(mean_read_len * MAX_RATIO_RL_WS))
    else:
        min_word_size = max(GLOBAL_MIN_WS, word_size - auto_word_size_step * 4)
        max_word_size = min(word_size + auto_word_size_step * 4, int(mean_read_len * MAX_RATIO_RL_WS))
    # arbitrarily adjust word size range
    if not auto_word_size_step:
        min_word_size = max_word_size = word_size
    else:
        max_word_size = max(int(mean_read_len * (1 - wc_bc_ratio_constant)), max_word_size)

    if float(word_size) / max_read_len <= 0.5 and len(low_quality_pattern) > 2:
        keep_seq_parts = True
    else:
        keep_seq_parts = False
    return min_word_size, word_size, max_word_size, keep_seq_parts


def check_kmers(kmer_str, auto_gradient_k, word_s, max_r_len, log):
    if kmer_str:
        # delete illegal kmer
        kmer_values = [int(kmer_v) for kmer_v in kmer_str.split(",")]

        if auto_gradient_k:
            # add kmer around estimated word_size
            k_tens_digits = [int(kmer_v / 10) for kmer_v in kmer_values]
            w_tens_digit = int(word_s / 10)
            if w_tens_digit not in k_tens_digits:
                # kmer_values.append(int(w_tens_digit * 10 + 5))
                avail_tens_digits = k_tens_digits + [w_tens_digit]
                for add_tens_digit in range(min(avail_tens_digits), max(avail_tens_digits)):
                    if add_tens_digit not in k_tens_digits:
                        kmer_values.append(int(add_tens_digit * 10 + 5))
        # delete illegal kmer
        kmer_values = [kmer_v for kmer_v in kmer_values if 21 <= kmer_v <= min(max_r_len, 127)]

        spades_kmer = ",".join([str(kmer_v) for kmer_v in sorted(kmer_values)])
        log.info("Setting '-k " + spades_kmer + "'")
        return spades_kmer
    else:
        return None


# test whether an external binary is executable
def executable(test_this):
    return True if os.access(test_this, os.X_OK) or getstatusoutput(test_this)[0] != DEAD_CODE else False


try:
    import psutil
except ImportError:
    this_process = None
else:
    this_process = psutil.Process(os.getpid())


def write_fq_results(original_fq_files, accepted_contig_id, out_file_name, temp2_clusters_dir, fq_info_in_memory,
                     maximum_n_reads, verbose, index_in_memory, log, extra_accepted_lines=set()):
    if verbose:
        sys.stdout.write(' ' * 100 + '\b' * 100)
        sys.stdout.flush()
        log.info("Producing output ...")
        log.info("reading indices ...")
    accepted_lines = []
    if index_in_memory:
        # read cluster indices
        for this_index in accepted_contig_id:
            accepted_lines += fq_info_in_memory[1][this_index]
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
    # add initial mapped read ids
    for line_id in extra_accepted_lines:
        accepted_lines.add(line_id)

    # write by line
    if verbose:
        log.info("writing fastq lines ...")
    post_reading = [open(fq_file, 'rU') for fq_file in original_fq_files]
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
            if count_r >= maximum_n_reads:
                break
        files_out[i].close()
        post_reading[i].close()
    del accepted_lines
    for i in range(len(original_fq_files)):
        os.rename(out_file_name + '_' + str(i + 1) + '.temp', out_file_name + '_' + str(i + 1) + '.fq')
    if verbose:
        log.info("writing fastq lines finished.")


def make_read_index(original_fq_files, direction_according_to_user_input, maximum_n_reads, rm_duplicates, output_base,
                    word_size, anti_lines, pre_grouped, index_in_memory, anti_seed, keep_seq_parts,
                    low_quality, echo_frequency, resume, log):
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
            log.info("Reading existed indices for fastq ...")
            #
            if keep_seq_parts:
                forward_reverse_reads = [x.strip().split("\t") for x in open(temp1_contig_dir[1], 'rU')]
                cancel_seq_parts = True if max([len(x) for x in forward_reverse_reads]) == 1 else False
            else:
                forward_reverse_reads = [x.strip() for x in open(temp1_contig_dir[1], 'rU')]
            #
            line_clusters = [[int(x) for x in y.split('\t')] for y in open(temp2_clusters_dir[1], 'rU')]
            if rm_duplicates:
                line_count = sum([len(x) for x in line_clusters]) * 4
            # log
            len_indices = len(line_clusters)
            if this_process:
                memory_usage = "Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3)) + " G, "
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
        lengths = []
        use_user_direction = False
        for id_file, file_name in enumerate(original_fq_files):
            file_in = open(file_name, "rU")
            count_this_read_n = 0
            line = file_in.readline()
            # if anti seed input, name & direction should be recognized
            if anti_seed:
                while line and count_this_read_n < maximum_n_reads:
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
                                else:
                                    log.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                    log.info("Using user-defined directions. ")
                                    use_user_direction = True
                                    this_name = line[1:].strip()
                                    direction = direction_according_to_user_input[id_file]
                            except (ValueError, IndexError):
                                log.info('Unrecognized head: ' + file_name + ': ' + str(line.strip()))
                                log.info("Using user-defined directions. ")
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
                                lengths.extend([len(seq_part) for seq_part in this_seq])
                            else:
                                this_seq = this_seq[0]
                                this_c_seq = complementary_seq(this_seq)
                                lengths.append(len(this_seq))
                        else:
                            this_c_seq = complementary_seq(this_seq)
                            lengths.append(len(this_seq))
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
                        log.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if line_count % echo_frequency == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            else:
                while line and count_this_read_n < maximum_n_reads:
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
                                lengths.extend([len(seq_part) for seq_part in this_seq])
                            else:
                                this_seq = this_seq[0]
                                this_c_seq = complementary_seq(this_seq)
                                lengths.append(len(this_seq))
                        else:
                            this_c_seq = complementary_seq(this_seq)
                            lengths.append(len(this_seq))
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
                        log.error("Illegal fq format in line " + str(line_count) + ' ' + str(line))
                        exit()
                    if line_count % echo_frequency == 0:
                        to_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: " + str(
                            (line_count + 4) // 4) + " reads"
                        sys.stdout.write(to_print + '\b' * len(to_print))
                        sys.stdout.flush()
                    line_count += 4
                    line = file_in.readline()
            line = file_in.readline()
            file_in.close()
            if line:
                log.warning("Number of reads exceeded " + str(int(maximum_n_reads)) + " in " + file_name + ", only "
                            "top " + str(int(maximum_n_reads)) + " reads are used in downstream analysis.")
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
                log.error("No qualified reads found!")
                max_read_len = max(lengths)
                if max_read_len < word_size:
                    log.error("Word size (" + str(word_size) + ") CANNOT be larger than your post-trimmed maximum read "
                              "length (" + str(max_read_len) + ")!")
                exit()
            log.info(memory_usage + str(len_indices) + " candidates in all " + str(line_count // 4) + " reads")
        else:
            del lengths
            log.info(memory_usage + str(len_indices) + " reads")
    if keep_seq_parts and cancel_seq_parts:
        keep_seq_parts = False
        for go_to, all_seq_parts in enumerate(forward_reverse_reads):
            forward_reverse_reads[go_to] = all_seq_parts[0]
    return forward_reverse_reads, line_clusters, len_indices, keep_seq_parts


def pre_grouping(fastq_indices_in_memory, dupli_threshold, out_base, index_in_memory, preg_word_size, log):
    forward_and_reverse_reads, line_clusters, len_indices, keep_seq_parts = fastq_indices_in_memory
    log.info("Pre-grouping reads ...")
    log.info("Setting '--pre-w " + str(preg_word_size) + "'")
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
    log.info(memory_usage + str(len(lines_with_duplicates)) + "/" + str(count_dupli) + " used/duplicated")

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
    log.info(memory_usage + str(len(groups_of_duplicate_lines)) + " groups made.")
    return groups_of_duplicate_lines, lines_with_duplicates


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


class WordsSoftLimitException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NoMoreReads(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def extending_reads(word_size, seed_file, seed_is_fq, original_fq_files, len_indices, pre_grouped,
                    groups_of_duplicate_lines, lines_with_duplicates, fq_info_in_memory, output_base,
                    max_rounds, min_rounds, fg_out_per_round, jump_step, mesh_size, verbose, resume,
                    maximum_n_reads, maximum_n_words,
                    auto_word_size_step, soft_max_n_words, max_word_size, min_word_size,
                    keep_seq_parts, low_quality_pattern, echo_frequency, log):
    used_word_sizes = {word_size}
    temp_max_n_words = soft_max_n_words
    not_reasonable_word_size = True

    check_times = 1000
    check_frequency = int(len_indices / check_times)

    while not_reasonable_word_size:
        # adding initial word
        log.info("Adding initial words ...")
        if seed_is_fq:
            if keep_seq_parts:
                accepted_words = chop_seq_list(
                    fq_seq_simple_generator(seed_file, split_pattern=low_quality_pattern, min_sub_seq=word_size),
                    word_size)
            else:
                accepted_words = chop_seqs(fq_seq_simple_generator(seed_file), word_size)
        else:
            accepted_words = chop_seqs(read_fasta(seed_file)[1], word_size)
        log.info("AW " + str(len(accepted_words)))

        accepted_rd_id = set()
        accepted_rd_id_this_round = set()
        g_duplicate_lines = deepcopy(groups_of_duplicate_lines)
        l_with_duplicates = deepcopy(lines_with_duplicates)
        line_to_accept = set()
        round_count = 1
        initial_aw_count = len(accepted_words)
        prev_aw_count = initial_aw_count
        accumulated_num_words = initial_aw_count
        check_points = [(int(percent * check_times) * check_frequency, percent) for percent in (0.25, 0.25, 0.5, 0.75)]
        if fg_out_per_round:
            round_dir = os.path.join(output_base, "Reads_per_round")
            if not os.path.exists(round_dir):
                os.mkdir(round_dir)
        if this_process and verbose:
            log.warning("Package psutil is not installed, so that memory usage will not be logged\n"
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
                                     os.path.join(output_base, 'temp.indices.2'), fq_info_in_memory, maximum_n_reads,
                                     verbose, bool(fq_info_in_memory), log)
                    # clear former accepted words from memory
                    del acc_words
                    # then add new accepted words into memory
                    if keep_seq_parts:
                        acc_words = chop_seq_list(
                            fq_seq_simple_generator(
                                [os.path.join(round_dir, "Round." + str(r_count) + '_' + str(x + 1) + '.fq') for x in
                                 range(len(original_fq_files))],
                                split_pattern=low_quality_pattern, min_sub_seq=word_size),
                            word_size, mesh_size)
                    else:
                        acc_words = chop_seqs(
                            fq_seq_simple_generator(
                                [os.path.join(round_dir, "Round." + str(r_count) + '_' + str(x + 1) + '.fq') for x in
                                 range(len(original_fq_files))]),
                            word_size, mesh_size)
                    acc_contig_id_this_round = set()
                log.info("Round " + str(r_count) + ': ' + str(unique_id + 1) + '/' + str(len_indices) + " AI " + str(
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

            def check_words_limit(inside_max_n_words):
                if accumulated_num_words + len(accepted_words) - prev_aw_count > inside_max_n_words:
                    if inside_max_n_words == maximum_n_words:
                        if this_process:
                            inside_memory_usage = " Mem " + str(
                                round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
                        else:
                            inside_memory_usage = ''
                        log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                 str(len_indices) + " AI " + str(len(accepted_rd_id_this_round)) +
                                 " AW " + str(len(accepted_words)) + inside_memory_usage)
                        raise WordsLimitException("")
                    elif word_size + auto_word_size_step < max_word_size and round_count <= min_rounds and \
                            (word_size + auto_word_size_step not in used_word_sizes):
                        if this_process:
                            inside_memory_usage = " Mem " + str(
                                round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
                        else:
                            inside_memory_usage = ''
                        log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                 str(len_indices) + " AI " + str(len(accepted_rd_id_this_round)) +
                                 " AW " + str(len(accepted_words)) + inside_memory_usage)
                        raise WordsSoftLimitException("")
                    else:
                        return True
                else:
                    return False

            def echo_to_screen():
                inside_this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round " \
                                    + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                    " AI " + str(len(accepted_rd_id_this_round)) + " AW " + str(len(accepted_words))
                sys.stdout.write(inside_this_print + '\b' * len(inside_this_print))
                sys.stdout.flush()

            # core extending code
            # here efficiency is more important than code conciseness,
            # so there are four similar structure with minor differences
            while True:
                if verbose:
                    log.info("Round " + str(round_count) + ": Start ...")

                if fq_info_in_memory:
                    reads_generator = (this_read for this_read in fq_info_in_memory[0])
                else:
                    if keep_seq_parts:
                        reads_generator = (this_read.strip().split("\t") for this_read in
                                           open(os.path.join(output_base, 'temp.indices.1'), 'rU'))
                    else:
                        reads_generator = (this_read.strip() for this_read in
                                           open(os.path.join(output_base, 'temp.indices.1'), 'rU'))
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
                            if unique_read_id % echo_frequency == 0:
                                echo_to_screen()
                            if unique_read_id % check_frequency == 0:
                                if temp_max_n_words != maximum_n_words and \
                                        check_points and unique_read_id % check_points[0][0] == 0:
                                    if check_words_limit(initial_aw_count +
                                                         check_points[0][1] * (temp_max_n_words - initial_aw_count)):
                                        temp_max_n_words = maximum_n_words
                                    check_points.pop(0)
                                else:
                                    if check_words_limit(temp_max_n_words):
                                        temp_max_n_words = maximum_n_words
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
                            if unique_read_id % echo_frequency == 0:
                                echo_to_screen()
                            if unique_read_id % check_frequency == 0:
                                if temp_max_n_words != maximum_n_words and \
                                        check_points and unique_read_id % check_points[0][0] == 0:
                                    if check_words_limit(initial_aw_count +
                                                         check_points[0][1] * (temp_max_n_words - initial_aw_count)):
                                        temp_max_n_words = maximum_n_words
                                    check_points.pop(0)
                                else:
                                    if check_words_limit(temp_max_n_words):
                                        temp_max_n_words = maximum_n_words
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
                            if unique_read_id % echo_frequency == 0:
                                echo_to_screen()
                            if unique_read_id % check_frequency == 0:
                                if temp_max_n_words != maximum_n_words and \
                                        check_points and unique_read_id % check_points[0][0] == 0:
                                    if check_words_limit(initial_aw_count +
                                                         check_points[0][1] * (temp_max_n_words - initial_aw_count)):
                                        temp_max_n_words = maximum_n_words
                                    check_points.pop(0)
                                else:
                                    if check_words_limit(temp_max_n_words):
                                        temp_max_n_words = maximum_n_words
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
                            if unique_read_id % echo_frequency == 0:
                                echo_to_screen()
                            if unique_read_id % check_frequency == 0:
                                if temp_max_n_words != maximum_n_words and \
                                        check_points and unique_read_id % check_points[0][0] == 0:
                                    if check_words_limit(initial_aw_count +
                                                         check_points[0][1] * (temp_max_n_words - initial_aw_count)):
                                        temp_max_n_words = maximum_n_words
                                    check_points.pop(0)
                                else:
                                    if check_words_limit(temp_max_n_words):
                                        temp_max_n_words = maximum_n_words
                        accepted_words, accepted_rd_id_this_round, prev_aw_count, round_count, accumulated_num_words \
                            = summarise_round(accepted_words, accepted_rd_id_this_round, prev_aw_count,
                                              round_count,
                                              accumulated_num_words, unique_read_id)
                reads_generator.close()
        except KeyboardInterrupt:
            reads_generator.close()
            sys.stdout.write(' ' * 100 + '\b' * 100)
            sys.stdout.flush()
            log.info(
                "Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + " AI " + str(
                    len(accepted_rd_id_this_round)) + " AW " + str(len(accepted_words)))
            log.info("KeyboardInterrupt")
        except NoMoreReads:
            reads_generator.close()
            if round_count < min_rounds:
                if auto_word_size_step:
                    word_size -= auto_word_size_step
                    if word_size > min_word_size and word_size not in used_word_sizes:
                        log.info("No more reads found and terminated ...")
                        log.info("Setting '-w " + str(word_size) + "'")
                        used_word_sizes.add(word_size)
                        continue
                    else:
                        log.info("No more reads found and terminated ...")
                        log.warning("Terminated at an insufficient number of rounds while "
                                    "'-w " + str(word_size + auto_word_size_step) + "'")
                else:
                    log.info("No more reads found and terminated ...")
                    log.warning("Terminated at an insufficient number of rounds")  # see '--auto-wss' and '-r'.")
            else:
                log.info("No more reads found and terminated ...")
        except WordsLimitException:
            reads_generator.close()
            if round_count <= min_rounds:
                if auto_word_size_step:
                    log.info("Hit the words limit and terminated ...")
                    log.warning("Terminated at an insufficient number of rounds while '-w " + str(word_size) + "'")
                else:
                    log.info("Hit the words limit and terminated ...")
                    log.warning("Terminated at an insufficient number of rounds, see "
                                "'--max-n-words' for more.")
            else:
                log.info("Hit the words limit and terminated ...")
        except WordsSoftLimitException:
            reads_generator.close()
            word_size += auto_word_size_step
            # log.info("Hit the words soft limit and terminated ...")
            log.info("Setting '-w " + str(word_size) + "'")
            used_word_sizes.add(word_size)
            continue
        except RoundLimitException as r_lim:
            reads_generator.close()
            log.info("Hit the round limit " + str(r_lim) + " and terminated ...")
        del reads_generator
        not_reasonable_word_size = False
        accepted_words = set()
        accepted_rd_id_this_round = set()
    del l_with_duplicates
    return accepted_rd_id


def get_anti_with_fas(word_size, anti_words, anti_input, original_fq_files, log):
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
            log.error('Unrecognized fq format in ' + str(line_count))
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
                    for i in range(0, temp_length + 1):
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


def mapping_with_bowtie2(seed_file, anti_seed, max_num_reads, original_fq_files, original_fq_beyond_read_limit,
                         out_base, resume, verbose_log, threads, random_seed,
                         prefix, keep_temp, bowtie2_other_options, log):
    if os.path.exists(seed_file + '.index.1.bt2l'):
        log.info("Bowtie2 index existed!")
    else:
        log.info("Making seed - bowtie2 index ...")
        build_seed_index = subprocess.Popen("bowtie2-build --large-index " + seed_file + " " + seed_file + '.index',
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log:
            log.info("bowtie2-build --large-index " + seed_file + " " + seed_file + '.index')
        output, err = build_seed_index.communicate()
        if "unrecognized option" in output.decode("utf8"):
            build_seed_index = subprocess.Popen("bowtie2-build " + seed_file + " " + seed_file + '.index',
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output, err = build_seed_index.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            log.error('\n' + output.decode("utf8"))
            exit()
        if verbose_log:
            log.info("\n" + output.decode("utf8").strip())
        log.info("Making seed - bowtie2 index finished.")
    seed_index_base = seed_file + '.index'

    query_fq_files = []
    for count_fq, ori_fq in enumerate(original_fq_files):
        if original_fq_beyond_read_limit[count_fq]:
            query_fq = ori_fq.strip() + ".Reduced"
            os.system("head -n " + str(int(max_num_reads * 4)) + " " + ori_fq + " > " + query_fq)
            query_fq_files.append(query_fq)
        else:
            query_fq_files.append(ori_fq)

    total_seed_fq = [os.path.join(out_base, x + prefix + "Initial.mapped.fq") for x in ("temp.", "")]
    total_seed_sam = [os.path.join(out_base, x + prefix + "seed_bowtie.sam") for x in ("temp.", "")]
    if resume and os.path.exists(total_seed_fq[1]):
        log.info("Initial seeds existed!")
    else:
        log.info("Mapping reads to seed - bowtie2 index ...")
        this_command = "bowtie2 --seed " + str(random_seed) + " --mm -p " + str(threads) + " " \
                       + bowtie2_other_options + " --al " + total_seed_fq[0] + \
                       " -x " + seed_index_base + " -U " + ",".join(query_fq_files) + " -S " + total_seed_sam[0] + \
                       " --no-unal --no-hd --no-sq --omit-sec-seq"
        make_seed_bowtie2 = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log:
            log.info(this_command)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            log.error('\n' + output.decode("utf8"))
            exit()
        if verbose_log:
            log.info("\n" + output.decode("utf8").strip())
        if os.path.exists(total_seed_sam[0]):
            os.rename(total_seed_sam[0], total_seed_sam[1])
            os.rename(total_seed_fq[0], total_seed_fq[1])
            log.info("Mapping finished.")
        else:
            log.error("Cannot find bowtie2 result!")
            exit()
    if anti_seed:
        if os.path.exists(anti_seed + '.index.1.bt2l'):
            log.info("anti-seed - bowtie2 index existed!")
        else:
            log.info("Making anti-seed - bowtie2 index ...")
            build_anti_index = subprocess.Popen(
                "bowtie2-build --large-index " + anti_seed + " " + anti_seed + '.index',
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            if verbose_log:
                log.info("bowtie2-build --large-index " + anti_seed + " " + anti_seed + '.index')
            output, err = build_anti_index.communicate()
            if "unrecognized option" in output.decode("utf8"):
                build_anti_index = subprocess.Popen("bowtie2-build " + anti_seed + " " + anti_seed + '.index',
                                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                if verbose_log:
                    log.info("bowtie2-build " + anti_seed + " " + anti_seed + '.index')
                output, err = build_anti_index.communicate()
            if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
                log.error('\n' + output.decode("utf8"))
                exit()
            if verbose_log:
                log.info("\n" + output.decode("utf8").strip())
            log.info("Making anti-seed - bowtie2 index finished.")
        anti_index_base = anti_seed + '.index'

        anti_seed_sam = [os.path.join(out_base, x + prefix + "anti_seed_bowtie.sam") for x in ("temp.", "")]
        if resume and os.path.exists(anti_seed_sam[1]):
            log.info("Anti-seed mapping information existed!")
        else:
            log.info("Mapping reads to anti-seed - bowtie2 index ...")
            this_command = "bowtie2 --seed " + str(random_seed) + " --mm -p " + str(threads) + " " + \
                           bowtie2_other_options + " -x " + anti_index_base + \
                           " -U " + ",".join(query_fq_files) + " -S " + anti_seed_sam[0] + \
                           " --no-unal --no-hd --no-sq --omit-sec-seq"
            make_anti_seed_bowtie2 = subprocess.Popen(this_command,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            if verbose_log:
                log.info(this_command)
            output, err = make_anti_seed_bowtie2.communicate()
            if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
                log.error('\n' + output.decode("utf8"))
                exit()
            if verbose_log:
                log.info("\n" + output.decode("utf8").strip())
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

    if not keep_temp:
        for count_fq, query_fq in enumerate(query_fq_files):
            if original_fq_beyond_read_limit[count_fq]:
                os.remove(query_fq)
    seed_fq_size = os.path.getsize(total_seed_fq[1])
    if not seed_fq_size:
        if log:
            log.error("No seed reads found!")
            log.error("Please check your raw data or change your reference!")
        exit()
    if seed_fq_size < 1024:
        log.info("Seed reads made: " + total_seed_fq[1] + " (" + str(int(seed_fq_size)) + " bytes)")
    elif 1024 * 1024 > seed_fq_size >= 1024:
        log.info("Seed reads made: " + total_seed_fq[1] + " (" + "%.2f" % round(seed_fq_size / 1024, 2) + " K)")
    elif seed_fq_size >= 1024 * 1024:
        log.info("Seed reads made: " + total_seed_fq[1] + " (" + "%.2f" % round(seed_fq_size / 1024 / 1024, 2) + " M)")
    return total_seed_fq[1], anti_lines


def assembly_with_spades(spades_kmer, spades_out_put, parameters, out_base, prefix, original_fq_files, reads_paired,
                         verbose_log, resume, threads, log):
    if '-k' in parameters or not spades_kmer:
        kmer = ''
    else:
        kmer = '-k ' + spades_kmer
    if resume and os.path.exists(spades_out_put):
        spades_command = 'spades.py --continue -o ' + spades_out_put
    else:
        spades_out_command = '-o ' + spades_out_put
        if reads_paired['input'] and reads_paired['pair_out']:
            all_unpaired = []
            # spades does not accept empty files
            if os.path.getsize(os.path.join(out_base, prefix + "filtered_1_unpaired.fq")):
                all_unpaired.append(os.path.join(out_base, prefix + "filtered_1_unpaired.fq"))
            if os.path.getsize(os.path.join(out_base, prefix + "filtered_2_unpaired.fq")):
                all_unpaired.append(os.path.join(out_base, prefix + "filtered_2_unpaired.fq"))
            for iter_unpaired in range(len(original_fq_files) - 2):
                if os.path.getsize(str(os.path.join(out_base, prefix + "filtered_" + str(iter_unpaired + 3) + ".fq"))):
                    all_unpaired.append(str(os.path.join(out_base, prefix + "filtered_" + str(iter_unpaired + 3) + ".fq")))
            if os.path.getsize(os.path.join(out_base, prefix + "filtered_1_paired.fq")):
                spades_command = ' '.join(
                    ['spades.py', '-t', str(threads), parameters, '-1',
                     os.path.join(out_base, prefix + "filtered_1_paired.fq"), '-2',
                     os.path.join(out_base, prefix + "filtered_2_paired.fq")] +
                    ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                    [kmer, spades_out_command]).strip()
            else:
                log.warning("No paired reads found for the target!?")
                spades_command = ' '.join(
                    ['spades.py', '-t', str(threads), parameters] +
                    ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                    [kmer, spades_out_command]).strip()
        else:
            all_unpaired = []
            for iter_unpaired in range(len(original_fq_files)):
                if os.path.getsize(str(os.path.join(out_base, prefix + "filtered_" + str(iter_unpaired + 1) + ".fq"))):
                    all_unpaired.append(str(os.path.join(out_base, prefix + "filtered_" + str(iter_unpaired + 1) + ".fq")))
            spades_command = ' '.join(
                ['spades.py', '-t', str(threads), parameters] +
                ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                [kmer, spades_out_command]).strip()
    spades_running = subprocess.Popen(spades_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log:
        log.info(spades_command)
    output, err = spades_running.communicate()
    if "not recognized" in output.decode("utf8") or not os.path.exists(spades_out_put):
        if verbose_log:
            log.error('Problem with running SPAdes:')
            log.error(output.decode("utf8"))
        log.error("WAINING in SPAdes: unrecognized option in " + parameters)
        log.error('Assembling failed.')
        return False
    elif "== Error ==" in output.decode("utf8"):
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
                log.warning("SPAdes failed for '-k " + failed_at_k + "'!")
                log.warning("If you need result based on kmer=" + failed_at_k + " urgently, "
                            "please check " + os.path.join(spades_out_put, "spades.log"))
                del real_kmer_values[real_kmer_values.index(failed_at_k)]
                log.warning("GetOrganelle would continue to process results based on "
                            "kmer=" + ",".join(real_kmer_values) + ".")
                # os.system("cp " + os.path.join(spades_out_put, "K" + real_kmer_values[-1], "assembly_graph.fastg")
                #           + " " + spades_out_put)
                log.info('Assembling finished with warnings.\n')
                return True
            else:
                log.warning("SPAdes failed with unknown errors!")
                log.warning("If you need to know more details, please check " +
                            os.path.join(spades_out_put, "spades.log") + " and contact SPAdes developers.")
                log.warning("GetOrganelle would continue to process results based on "
                            "kmer=" + ",".join(real_kmer_values) + ".")
                # os.system("cp " + os.path.join(spades_out_put, "K" + real_kmer_values[-1], "assembly_graph.fastg")
                #           + " " + spades_out_put)
                log.info('Assembling finished with warnings.\n')
                return True
        else:
            log.error("Error in SPAdes: \n== Error ==" + output.decode("utf8").split("== Error ==")[-1].split("In case you")[0])
            log.error('Assembling failed.')
            return False
    elif not os.path.exists(os.path.join(spades_out_put, "assembly_graph.fastg")):
        if verbose_log:
            log.info(output.decode("utf8"))
        log.warning("Assembling exited halfway.\n")
        return True
    else:
        spades_log = output.decode("utf8").split("\n")
        if verbose_log:
            log.info(output.decode("utf8"))
        for line in spades_log:
            line = line.strip()
            if line.count(":") > 2 and "Insert size = " in line and \
                    line.split()[0].replace(":", "").replace(".", "").isdigit():
                try:
                    log.info(line.split("   ")[-1].split(", read length =")[0].strip())
                except IndexError:
                    pass
        log.info('Assembling finished.\n')
        return True


def slim_spades_result(scheme, custom_db, spades_output, verbose_log, log, threads, depth_threshold=0, resume=False):
    if not executable("blastn"):
        if log:
            log.warning('blastn not in the path! Skip slimming assembly result ...')
        return
    if not executable("makeblastdb"):
        if log:
            log.warning('makeblastdb not in the path! Skip slimming assembly result ...')
        return
    graph_file = os.path.join(spades_output, "assembly_graph.fastg")
    if resume:
        for existed_file in os.listdir(spades_output):
            if existed_file.count(".fastg") == 2:
                if log:
                    log.info("Slimming      " + graph_file + " ... skipped.")
                return 0
    scheme_translation = {
        'plant_cp': ' --include-priority ' +
                    os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'plant_cp') + ' --exclude ' +
                    os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'plant_mt'),
        'plant_mt': ' --include-priority ' +
                    os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'plant_mt') + ' --exclude ' +
                    os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'plant_cp'),
        'plant_nr': ' --include-priority ' +
                    os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'plant_nr'),
        'animal_mt': ' --include-priority ' +
                     os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'animal_mt'),
        'fungus_mt': ' --include-priority ' +
                     os.path.join(PATH_OF_THIS_SCRIPT, 'Library', 'NotationReference', 'fungus_mt'),
        'anonym': ' -F anonym --genes ' + custom_db if custom_db else ""
    }
    if scheme in scheme_translation:
        run_command = scheme_translation[scheme]
    else:
        run_command = scheme
    run_command = os.path.join(PATH_OF_THIS_SCRIPT, 'Utilities', 'slim_fastg.py') + ' -t ' + str(threads) + ' ' \
                  + graph_file + run_command  #\
                  # + ' -o ' + out_base + (' --prefix ' + prefix if prefix else "")
    slim_spades = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log:
        log.info(run_command)
    output, err = slim_spades.communicate()
    output_file_list = [os.path.join(spades_output, x) for x in os.listdir(spades_output) if x.count(".fastg") == 2]
    if "not recognized" in output.decode("utf8") or "command not found" in output.decode("utf8"):
        if log:
            if verbose_log:
                log.warning(os.path.join(PATH_OF_THIS_SCRIPT, "Utilities", "slim_spades_fastg_by_blast.py") + ' not found!')
                log.warning(output.decode("utf8"))
            log.warning("Slimming      " + graph_file + " failed.")
        return 1
    elif "failed" in output.decode("utf8") or "error" in output.decode("utf8") or "Error" in output.decode("utf8"):
        if log:
            if verbose_log:
                log.error(output.decode("utf8"))
            log.warning("Slimming      " + graph_file + " failed.")
        return 1
    elif output_file_list and os.path.getsize(output_file_list[0]) == 0:
        if log:
            log.warning("Slimming      " + graph_file + " finished with no target organelle contigs found!")
        return 2
    else:
        if log:
            if verbose_log:
                log.info(output.decode("utf8"))
            log.info("Slimming      " + graph_file + " finished!")
        return 0


def separate_fq_by_pair(out_base, prefix, verbose_log, log):
    log.info("Separating filtered fastq file ... ")
    run_command = os.path.join(PATH_OF_THIS_SCRIPT, "Utilities", "get_pair_reads.py") \
                  + ' ' + os.path.join(out_base, prefix + "filtered_1.fq") \
                  + ' ' + os.path.join(out_base, prefix + "filtered_2.fq")
    separating = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if verbose_log:
        log.info(run_command)
    output, err = separating.communicate()
    if "not recognized" in output.decode("utf8") or "command not found" in output.decode("utf8"):
        if verbose_log:
            log.warning(os.path.join(PATH_OF_THIS_SCRIPT, "Utilities", "get_pair_reads.py") + "not found!")
            log.warning(output.decode("utf8"))
        log.warning("Separating filtered fastq file failed.")
        return False
    elif not os.path.exists(os.path.join(out_base, prefix + "filtered_2_paired.fq")):
        log.warning("Separating filtered fastq file failed with unknown error: " + run_command)
        return False
    else:
        if verbose_log:
            log.info(output.decode("utf8"))
        log.info("Separating filtered fastq file finished!")
        return True


def unzip(source, target, line_limit, verbose_log, log):
    target_temp = target + ".Temp"
    try_commands = ["gunzip -c " + source + " | head -n " + str(line_limit) + " > " + target_temp,
                    "tar -x -f " + source + " -O | head -n " + str(line_limit) + " > " + target_temp]
    log.info("Unzipping reads file: " + source)
    success = False
    output = b""
    for run_command in try_commands:
        unzipping = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log:
            log.info(run_command)
        output, err = unzipping.communicate()
        if "Unrecognized" not in output.decode("utf8") and \
                "Error" not in output.decode("utf8") and \
                "error" not in output.decode("utf8"):
            success = True
            break
    if success:
        if verbose_log:
            log.info(output.decode("utf8"))
        os.rename(target_temp, target)
    else:
        if verbose_log:
            log.error("\n" + output.decode("utf8"))
        try:
            os.remove(target_temp)
        except:
            pass
        raise NotImplementedError("unzipping failed!")


def extract_organelle_genome(out_base, spades_output, prefix, organelle_type, ref_seq, read_len_for_log,
                             verbose, log_in, threads, options):
    def disentangle_assembly(fastg_file, tab_file, output, weight_factor, log_dis, time_limit, type_factor=3.,
                             mode="plant_cp", contamination_depth=3., contamination_similarity=0.95, degenerate=True,
                             degenerate_depth=1.5, degenerate_similarity=0.98,
                             expected_max_size=inf, expected_min_size=0, hard_cov_threshold=10.,
                             min_sigma_factor=0.1, here_only_max_c=True, here_acyclic_allowed=False,
                             here_verbose=False, timeout_flag_str="'--disentangle-time-limit'"):
        @set_time_limit(time_limit, flag_str=timeout_flag_str)
        def disentangle_inside(fastg_f, tab_f, o_p, w_f, log_in, type_f=3., mode_in="plant_cp", c_d=3., c_s=0.95,
                               deg=True, deg_dep=1.5, deg_sim=0.98, hard_c_t=10., min_s_f=0.1, max_c_in=True,
                               max_s=inf, min_s=0, acyclic_allowed_in=False, verbose_in=False):
            from Library.assembly_parser import Assembly
            this_K = os.path.split(os.path.split(fastg_f)[0])[-1]
            o_p += "." + this_K
            if acyclic_allowed_in:
                log_in.info("Disentangling " + fastg_f + " as contig(s) ... ")
            else:
                log_in.info("Disentangling " + fastg_f + " as a circular genome ... ")
            input_graph = Assembly(fastg_f)
            target_results = input_graph.find_target_graph(tab_f,
                                                           mode=mode_in, type_factor=type_f,
                                                           log_hard_cov_threshold=hard_c_t,
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
                                                           log_handler=log_in, verbose=verbose_in)
            if len(target_results) > 1:
                log_in.warning(str(len(target_results)) + " sets of graph detected!")
            log_in.info("Slimming and disentangling graph finished!\n")

            log_in.info("Writing output ...")
            degenerate_base_used = False
            if acyclic_allowed_in:
                contig_num = set()
                for go_res, res in enumerate(target_results):
                    broken_graph = res["graph"]
                    count_path = 0
                    test_paths = broken_graph.get_all_paths(mode=mode_in, log_handler=log_in)
                    for this_paths, other_tag in test_paths:
                        count_path += 1
                        all_contig_str = []
                        contig_num.add(len(this_paths))
                        for go_contig, this_p_part in enumerate(this_paths):
                            this_contig = broken_graph.export_path(this_p_part)
                            if DEGENERATE_BASES & set(this_contig.seq):
                                degenerate_base_used = True
                            all_contig_str.append(">contig_" + str(go_contig + 1) + "--" + this_contig.label + "\n" +
                                                  this_contig.seq + "\n")
                        open(o_p + ".contigs.graph" + str(go_res + 1) + other_tag + "." + str(count_path) +
                             ".path_sequence.fasta", "w").write("\n".join(all_contig_str))
                        log_in.info("Writing PATH" + str(count_path) + " of " + mode_in + " contig(s) to " + o_p +
                                    ".contigs.graph" + str(go_res + 1) + other_tag + "." + str(count_path) +
                                    ".path_sequence.fasta")
                    log_in.info("Writing GRAPH to " + o_p + ".contigs.graph" + str(go_res + 1) + ".selected_graph.gfa")
                    broken_graph.write_to_gfa(o_p + ".contigs.graph" + str(go_res + 1) + ".selected_graph.gfa")
                log_in.info("Result status: " + ",".join(sorted([str(c_n) for c_n in contig_num])) + " contig(s)")
            else:
                for go_res, res in enumerate(target_results):
                    go_res += 1
                    idealized_graph = res["graph"]
                    # average_kmer_cov = res["cov"]
                    # log.info("Detecting target graph" + str(go_res) +
                    #          " finished with average kmer coverage: " + str(round(average_kmer_cov, 4)))
                    count_path = 0
                    for this_path, other_tag in idealized_graph.get_all_circular_paths(mode=mode_in, log_handler=log_in):
                        count_path += 1
                        out_n = o_p + ".complete.graph" + str(go_res) + "." + str(count_path) + other_tag + \
                                ".path_sequence.fasta"
                        this_seq_obj = idealized_graph.export_path(this_path)
                        if DEGENERATE_BASES & set(this_seq_obj.seq):
                            degenerate_base_used = True
                        open(out_n, "w").write(this_seq_obj.fasta_str())
                        log_in.info("Writing PATH" + str(count_path) + " of complete genome to " + out_n)
                    log_in.info("Writing GRAPH to " + o_p + ".complete.graph" + str(go_res) + ".selected_graph.gfa")
                    idealized_graph.write_to_gfa(o_p + ".complete.graph" + str(go_res) + ".selected_graph.gfa")
                log_in.info("Result status: circular genome")
            if degenerate_base_used:
                log_in.warning("Degenerate base(s) used!")
            os.system("cp " + os.path.join(os.path.split(fastg_f)[0], "assembly_graph.fastg") + " " +
                      o_p + ".assembly_graph.fastg")
            os.system("cp " + fastg_f + " " + o_p + "." + os.path.basename(fastg_f))
            os.system("cp " + tab_f + " " + o_p + "." + os.path.basename(tab_f))
            if not acyclic_allowed_in:
                log_in.info("Please visualize " + o_p + "." + os.path.basename(fastg_f) + " to confirm the final result.")
            log_in.info("Writing output finished.")

        disentangle_inside(fastg_f=fastg_file, tab_f=tab_file, o_p=output, w_f=weight_factor, log_in=log_dis,
                           type_f=type_factor, mode_in=mode, c_d=contamination_depth, c_s=contamination_similarity,
                           deg=degenerate, deg_dep=degenerate_depth, deg_sim=degenerate_similarity,
                           hard_c_t=hard_cov_threshold, min_s_f=min_sigma_factor, max_c_in=here_only_max_c,
                           max_s=expected_max_size, min_s=expected_min_size,
                           acyclic_allowed_in=here_acyclic_allowed, verbose_in=here_verbose)

    # start
    export_succeeded = False
    kmer_vals = sorted([int(kmer_d[1:])
                        for kmer_d in os.listdir(spades_output)
                        if os.path.isdir(os.path.join(spades_output, kmer_d))
                        and kmer_d.startswith("K")
                        and os.path.exists(os.path.join(spades_output, kmer_d, "assembly_graph.fastg"))],
                       reverse=True)
    kmer_dirs = [os.path.join(spades_output, "K" + str(kmer_val)) for kmer_val in kmer_vals]
    run_stat_set = set()
    timeout_flag = "'--disentangle-time-limit'"
    for go_k, kmer_dir in enumerate(kmer_dirs):
        try:
            running_stat = slim_spades_result(organelle_type, ref_seq, kmer_dir, verbose, log_in, threads=threads,
                                              resume=options.script_resume)
            run_stat_set.add(running_stat)
            """disentangle"""
            if running_stat == 0:
                out_fastg = sorted([os.path.join(kmer_dir, x)
                                    for x in os.listdir(kmer_dir) if x.count(".fastg") == 2])[0]
                out_csv = out_fastg[:-5] + "csv"
                # if it is the first round (the largest kmer), copy the slimmed result to the main spades output folder
                # if go_k == 0:
                #     main_spades_folder = os.path.split(kmer_dir)[0]
                #     os.system("cp " + out_fastg + " " + main_spades_folder)
                #     os.system("cp " + out_csv + " " + main_spades_folder)
                path_prefix = os.path.join(out_base, prefix + organelle_type)
                disentangle_assembly(fastg_file=out_fastg, mode=organelle_type, tab_file=out_csv, output=path_prefix,
                                     weight_factor=100, hard_cov_threshold=options.disentangle_depth_factor,
                                     contamination_depth=options.contamination_depth,
                                     contamination_similarity=options.contamination_similarity,
                                     degenerate=options.degenerate, degenerate_depth=options.degenerate_depth,
                                     degenerate_similarity=options.degenerate_similarity,
                                     expected_max_size=options.expected_max_size,
                                     expected_min_size=options.expected_min_size,
                                     here_only_max_c=True,
                                     here_acyclic_allowed=False, here_verbose=verbose, log_dis=log_in,
                                     time_limit=options.disentangle_time_limit, timeout_flag_str=timeout_flag)
                # currently time is not limited for exporting contigs
        except ImportError:
            log_in.warning("Disentangling failed: numpy/scipy/sympy not installed!")
            break
        except RuntimeError:
            log_in.info("Disentangling timeout. (see " + timeout_flag + " for more)")
        except Exception as e:
            log_in.info("Disentangling failed: " + str(e).strip())
        else:
            if running_stat == 0:
                export_succeeded = True
                break

    if not export_succeeded:
        if run_stat_set == {2}:
            log_in.warning("No target organelle contigs found!")
            log_in.warning("This might due to a bug or unreasonable reference/parameter choices")
            log_in.info("Please email phylojianjun@163.com or jinjianjun@mail.kib.ac.cn with the get_org.log.txt file.")
            return export_succeeded

        largest_k_graph_f_exist = sorted([x for x in os.listdir(kmer_dirs[0]) if x.count(".fastg") == 2])
        if kmer_dirs and largest_k_graph_f_exist:
            for go_k, kmer_dir in enumerate(kmer_dirs):
                try:
                    """disentangle the graph as contig(s)"""
                    out_fastg_list = sorted([os.path.join(kmer_dir, x)
                                             for x in os.listdir(kmer_dir) if x.count(".fastg") == 2])
                    if out_fastg_list:
                        out_fastg = out_fastg_list[0]
                        out_csv = out_fastg[:-5] + "csv"
                        path_prefix = os.path.join(out_base, prefix + organelle_type)
                        disentangle_assembly(fastg_file=out_fastg, mode=organelle_type, tab_file=out_csv,
                                             output=path_prefix, weight_factor=100, here_verbose=verbose,
                                             log_dis=log_in, hard_cov_threshold=options.disentangle_depth_factor * 0.8,
                                             contamination_depth=options.contamination_depth,
                                             contamination_similarity=options.contamination_similarity,
                                             degenerate=options.degenerate, degenerate_depth=options.degenerate_depth,
                                             degenerate_similarity=options.degenerate_similarity,
                                             expected_max_size=options.expected_max_size,
                                             expected_min_size=options.expected_min_size,
                                             here_only_max_c=True, here_acyclic_allowed=True,
                                             time_limit=3600, timeout_flag_str=timeout_flag)
                except ImportError:
                    break
                except RuntimeError:
                    log_in.info("Disentangling timeout. (see " + timeout_flag + " for more)")
                except Exception as e:
                    log_in.info("Disentangling failed: " + str(e).strip())
                else:
                    export_succeeded = True
                    out_fastg = largest_k_graph_f_exist[0]
                    out_csv = out_fastg[:-5] + "csv"
                    log_in.info("Please ...")
                    log_in.info("load the graph file '" + out_fastg +
                                "' in " + ",".join(["K" + str(k_val) for k_val in kmer_vals]))
                    log_in.info("load the CSV file '" + out_csv +
                                "' in " + ",".join(["K" + str(k_val) for k_val in kmer_vals]))
                    log_in.info("visualize and confirm the incomplete result in Bandage.")
                    # log.info("-------------------------------------------------------")
                    log_in.info("If the result is nearly complete, ")
                    log_in.info("you can also adjust the arguments to try again.")
                    log_in.info("If you have questions for us, "
                                "please provide us with the get_org.log.txt file and the graph in the format you like!")
                    # log.info("-------------------------------------------------------")
                    break
            if not export_succeeded:
                out_fastg = largest_k_graph_f_exist[0]
                out_csv = out_fastg[:-5] + "csv"
                log_in.info("Please ...")
                log_in.info("load the graph file '" + out_fastg + "/assembly_graph.fastg" +
                            "' in " + ",".join(["K" + str(k_val) for k_val in kmer_vals]))
                log_in.info("load the CSV file '" + out_csv +
                            "' in " + ",".join(["K" + str(k_val) for k_val in kmer_vals]))
                log_in.info("visualize and export your result in Bandage.")
                log_in.info("If you have questions for us, "
                            "please provide us with the get_org.log.txt file and the graph in the format you like!")
        else:
            # slim failed with unknown error
            log_in.info("Please ...")
            log_in.info("load the graph file: " + os.path.join(spades_output, 'assembly_graph.fastg'))
            log_in.info("visualize and export your result in Bandage.")
            log_in.info("If you have questions for us, "
                        "please provide us with the get_org.log.txt file and the graph in the format you like!")
    return export_succeeded


def main():
    time0 = time.time()
    title = "GetOrganelle v" + str(get_versions()) + \
            "\n" \
            "\ngets organelle reads and genomes from genome skimming data by extending." \
            "\nFind updates in https://github.com/Kinggerm/GetOrganelle and see README.md for more information." \
            "\n"
    options, log, previous_attributes = get_options(descriptions=title, version=get_versions())
    resume = options.script_resume
    verb_log = options.verbose_log
    out_base = options.output_base
    echo_frequency = options.echo_frequency
    reads_files_to_drop = []
    # global word_size
    word_size = None
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
        other_options = options.other_spades_options.split(' ')
        if '-o' in other_options:
            which_out = other_options.index('-o')
            spades_output = other_options[which_out + 1]
            del other_options[which_out: which_out + 2]
        else:
            spades_output = os.path.join(out_base, options.prefix + "filtered_spades")
        other_options = ' '.join(other_options)

        """ get reads """
        filtered_files_exist = max(
            min([os.path.exists(str(os.path.join(out_base, options.prefix + "filtered")) + '_' + str(i + 1) + '_unpaired.fq')
                 for i in range(2)] +
                [os.path.exists(str(os.path.join(out_base, options.prefix + "filtered")) + '_' + str(i + 1) + '.fq')
                 for i in range(2, len(original_fq_files))]),
            min([os.path.exists(str(os.path.join(out_base, options.prefix + "filtered")) + '_' + str(i + 1) + '.fq')
                 for i in range(len(original_fq_files))]))

        if resume:
            try:
                max_read_len = int(previous_attributes.max_read_len)
                mean_read_len = float(previous_attributes.mean_read_len)
            except (AttributeError, ValueError):
                resume = False
            try:
                word_size = int(previous_attributes.w)
            except (AttributeError, ValueError):
                if filtered_files_exist:
                    resume = False
                else:
                    pass

        if not (resume and filtered_files_exist):
            # bowt_seed = options.bowtie2_seed
            anti_seed = options.anti_seed
            # b_at_seed = options.bowtie2_anti_seed
            pre_grp = options.pre_grouped
            in_memory = options.index_in_memory

            if original_fq_files:
                for file_id, read_file in enumerate(original_fq_files):
                    # unzip fq files if needed
                    if read_file.endswith(".gz") or read_file.endswith(".zip"):
                        target_fq = read_file + ".fastq"
                        if not (os.path.exists(target_fq) and resume):
                            unzip(read_file, target_fq, int(4 * options.maximum_n_reads), options.verbose_log, log)
                        if os.path.getsize(target_fq) == 0:
                            raise ValueError("Empty file " + target_fq)
                        original_fq_files[file_id] = target_fq
                        reads_files_to_drop.append(target_fq)

            log.info("Pre-reading fastq ...")
            if resume:
                try:
                    original_fq_beyond_read_limit = [previous_attributes.record_fq_beyond_limit
                                                     for j in range(len(original_fq_files))]
                    all_read_num = int(previous_attributes.num_reads)
                except (AttributeError, ValueError):
                    resume = False
                else:
                    try:
                        low_quality_pattern = "[" + previous_attributes.trim_chars + "]"
                        mean_error_rate = float(previous_attributes.mean_error_rate)
                    except (AttributeError, ValueError):
                        low_quality_pattern = "[]"
                        mean_error_rate = None
                    all_bases = mean_read_len * all_read_num
            if not resume:
                sampling_percent = 0.1
                if options.pre_reading:
                    # pre-reading fastq
                    log.info("Counting read qualities ...")
                    low_quality_pattern, mean_error_rate, original_fq_beyond_read_limit = \
                        get_read_quality_info(original_fq_files, options.maximum_n_reads, options.min_quality_score,
                                              log, maximum_ignore_percent=options.maximum_ignore_percent,
                                              sampling_percent=sampling_percent)
                    log.info("Counting read lengths ...")
                    mean_read_len, max_read_len, all_read_num = get_read_len_mean_max_count(original_fq_files,
                                                                                            options.maximum_n_reads)
                    all_bases = mean_read_len * all_read_num
                    log.info("Mean = " + str(round(mean_read_len, 1)) + " bp, maximum = " + str(max_read_len) + " bp.")
                    log.info("Reads used = " + str(all_read_num))
                else:
                    low_quality_pattern, mean_error_rate, original_fq_beyond_read_limit = \
                        "[]", None, [True for j in range(len(original_fq_files))]
                    mean_read_len, max_read_len, all_read_num = get_read_len_mean_max_count(original_fq_files,
                                                                                            options.maximum_n_reads,
                                                                                            sampling_percent)
                    log.info("Mean = " + str(round(mean_read_len, 1)) + " bp, maximum = " + str(max_read_len) + " bp.")
                    log.info("Reads used = " + str(all_read_num))
                    all_bases = None
                log.info("Pre-reading fastq finished.\n")
            else:
                log.info("Pre-reading fastq skipped.\n")

            # reading seeds
            log.info("Making seed reads ...")
            seed_file = check_fasta_seq_names(options.seed_file, log)
            if options.utilize_mapping:
                seed_file, anti_lines = mapping_with_bowtie2(seed_file, anti_seed, options.maximum_n_reads,
                                                             original_fq_files,
                                                             original_fq_beyond_read_limit, out_base, resume,
                                                             verb_log, options.threads, options.random_seed,
                                                             options.prefix, options.keep_temp_files,
                                                             bowtie2_other_options=options.bowtie2_options,
                                                             log=log)
            log.info("Making seed reads finished.\n")

            log.info("Checking seed reads and parameters ...")
            if not resume or options.word_size:
                word_size = options.word_size
            if options.pre_reading:
                min_word_size, word_size, max_word_size, keep_seq_parts = \
                    check_parameters(word_size=word_size, out_base=out_base, utilize_mapping=options.utilize_mapping,
                                     maximum_n_reads=options.maximum_n_reads, original_fq_files=original_fq_files,
                                     organelle_type=options.organelle_type, all_bases=all_bases,
                                     auto_word_size_step=options.auto_word_size_step, mean_error_rate=mean_error_rate,
                                     target_genome_size=options.target_genome_size, mean_read_len=mean_read_len,
                                     max_read_len=max_read_len, low_quality_pattern=low_quality_pattern, log=log,
                                     wc_bc_ratio_constant=0.35, larger_auto_ws=options.larger_auto_ws)
            else:
                word_size = calculate_word_size_according_to_ratio(word_size, mean_read_len, log)
                min_word_size, max_word_size, keep_seq_parts = word_size, word_size, False
            log.info("Checking seed reads and parameters finished.\n")

            # make read index
            log.info("Making read index ...")
            if not options.utilize_mapping:
                anti_lines = get_anti_with_fas(word_size, chop_seqs(read_fasta(anti_seed)[1], word_size),
                                               bool(anti_seed),
                                               original_fq_files, log) if anti_seed else set()
            fq_info_in_memory = make_read_index(original_fq_files, direction_according_to_user_input,
                                                options.maximum_n_reads, options.rm_duplicates, out_base, word_size,
                                                anti_lines, pre_grp, in_memory, anti_seed,
                                                keep_seq_parts=keep_seq_parts, low_quality=low_quality_pattern,
                                                resume=resume, echo_frequency=echo_frequency, log=log)
            len_indices = fq_info_in_memory[2]
            keep_seq_parts = fq_info_in_memory[3]
            if keep_seq_parts:
                log.info("Reads are stored as fragments.")
            # pre-grouping if asked
            if pre_grp:
                preg_word_size = word_size if not options.pregroup_word_size else options.pregroup_word_size
                groups_of_lines, lines_with_dup = pre_grouping(fq_info_in_memory, pre_grp, out_base, in_memory,
                                                               preg_word_size=preg_word_size, log=log)
            else:
                groups_of_lines = None
                lines_with_dup = None
            if not in_memory:
                fq_info_in_memory = None
            log.info("Making read index finished.\n")

            # extending process
            log.info("Extending ...")
            # accepted_ids = set()
            if options.auto_word_size_step:
                if options.maximum_n_words <= options.soft_max_words:
                    log.info("Setting '--soft-max-words " + str(int(options.maximum_n_words)) + "'")
                    options.soft_max_words = options.maximum_n_words
            accepted_rd_id = extending_reads(word_size=word_size, seed_file=seed_file,
                                             seed_is_fq=options.utilize_mapping,
                                             original_fq_files=original_fq_files, len_indices=len_indices,
                                             pre_grouped=pre_grp, groups_of_duplicate_lines=groups_of_lines,
                                             lines_with_duplicates=lines_with_dup,
                                             fq_info_in_memory=fq_info_in_memory, output_base=out_base,
                                             max_rounds=options.max_rounds, min_rounds=options.min_rounds,
                                             fg_out_per_round=options.fg_out_per_round,
                                             jump_step=options.jump_step, mesh_size=options.mesh_size,
                                             verbose=verb_log, resume=resume,
                                             maximum_n_reads=options.maximum_n_reads,
                                             maximum_n_words=options.maximum_n_words,
                                             auto_word_size_step=options.auto_word_size_step,
                                             soft_max_n_words=options.soft_max_words,
                                             min_word_size=min_word_size, max_word_size=max_word_size,
                                             keep_seq_parts=keep_seq_parts, low_quality_pattern=low_quality_pattern,
                                             echo_frequency=echo_frequency, log=log)
            mapped_read_ids = set()
            write_fq_results(original_fq_files, accepted_rd_id,
                             os.path.join(out_base, options.prefix + "filtered"),
                             os.path.join(out_base, 'temp.indices.2'),
                             fq_info_in_memory, options.maximum_n_reads, verb_log, in_memory, log, mapped_read_ids)
            del accepted_rd_id, fq_info_in_memory, groups_of_lines, \
                anti_lines, lines_with_dup

            if not options.keep_temp_files:
                try:
                    os.remove(os.path.join(out_base, 'temp.indices.1'))
                    os.remove(os.path.join(out_base, 'temp.indices.2'))
                except OSError:
                    pass

            log.info("Extending finished.\n")
        else:
            log.info("Extending ... skipped.\n")
        if reads_files_to_drop and not options.keep_temp_files:
            for rm_read_file in reads_files_to_drop:
                os.remove(rm_read_file)

        if reads_paired['input']:
            if not (resume and min([os.path.exists(x) for x in
                                    (os.path.join(out_base, options.prefix + "filtered_" + y + "_" + z + "paired.fq")
                                     for y in ('1', '2') for z in ('', 'un'))])):
                resume = False
                reads_paired['pair_out'] = separate_fq_by_pair(out_base, options.prefix, verb_log, log)
                if reads_paired['pair_out'] and not options.keep_temp_files:
                    os.remove(os.path.join(out_base, options.prefix + "filtered_1.fq"))
                    os.remove(os.path.join(out_base, options.prefix + "filtered_2.fq"))
            else:
                log.info("Separating filtered fastq file ... skipped.\n")

        """ assembly """
        is_assembled = False
        if options.run_spades:
            if not (resume and os.path.exists(os.path.join(spades_output, 'assembly_graph.fastg'))):
                options.spades_kmer = check_kmers(options.spades_kmer, options.auto_gradient_k, word_size,
                                                  max_read_len, log)
                log.info('Assembling using SPAdes ...')
                is_assembled = assembly_with_spades(options.spades_kmer, spades_output, other_options, out_base,
                                                    options.prefix, original_fq_files, reads_paired,
                                                    options.verbose_log, resume, options.threads, log)
            else:
                is_assembled = True
                log.info('Assembling using SPAdes ... skipped.\n')

        """ export organelle """
        if is_assembled and options.organelle_type != '0':
            graph_existed = bool([gfa_f for gfa_f in os.listdir(out_base) if gfa_f.endswith(".selected_graph.gfa")])
            fasta_existed = bool([fas_f for fas_f in os.listdir(out_base) if fas_f.endswith(".path_sequence.fasta")])
            if resume and graph_existed and fasta_existed:
                log.info("Slimming and disentangling graph ... skipped.")
            else:
                log.info("Slimming and disentangling graph ...")
                extract_organelle_genome(out_base=out_base, spades_output=spades_output, prefix=options.prefix,
                                         organelle_type=options.organelle_type,
                                         ref_seq=options.genes_fasta,
                                         read_len_for_log=mean_read_len,
                                         verbose=options.verbose_log, log_in=log,
                                         threads=options.threads, options=options)

        log = simple_log(log, out_base, prefix=options.prefix + "get_org.")
        log.info("\nTotal cost " + "%.2f"%(time.time() - time0) + " s")
        log.info("Thank you!")
    except:
        log.exception("")
        log = simple_log(log, out_base, prefix=options.prefix + "get_org.")
        log.info("\nTotal cost " + "%.2f"%(time.time() - time0) + " s")
        log.info("Please email phylojianjun@163.com or jinjianjun@mail.kib.ac.cn if you find bugs!")
        log.info("Please provide me with the get_org.log.txt file!")
    logging.shutdown()


if __name__ == '__main__':
    main()

"""Copyright 2016 Jianjun Jin"""
