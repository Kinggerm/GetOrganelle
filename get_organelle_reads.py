#!/usr/bin/env python

import datetime
import sys
import os
import re
from optparse import OptionParser, OptionGroup
from VERSIONS import get_versions
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
from Library.seq_parser import *
from Library.pipe_control_func import *
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
import time

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

word_size = 0


def require_commands(print_title, version):
    version = version
    usage = "\n###  Chloroplast, Normal, 2G raw data, 150 bp reads\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -o chloroplast_output " \
            " -R 10 -k 75,85,95,105\n" \
            "###  Mitochondria\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s mt_reference.fasta -w 0.6 -o mitochondria_output " \
            " -R 30 -k 75,85,95,105 -F mt\n" \
            "###  Nuclear Ribosomal RNA (18S-ITS1-5.8S-ITS2-26S)\n" + str(os.path.basename(__file__)) + \
            " -1 sample_1.fq -2 sample_2.fq -s nr_reference.fasta -w 0.85 -o nr_output " \
            " -R 7 -k 95,105,115 -P 0 -F nr\n"
    description = print_title
    parser = OptionParser(usage=usage, version=version, description=description)
    # group1
    group_need = OptionGroup(parser, "COMMON OPTIONS", "All these arguments are required unless alternations provided")
    group_need.add_option('-1', dest='fastq_file_1', help='Input file with forward paired-end reads as pool.')
    group_need.add_option('-2', dest='fastq_file_2', help='Input file with reverse paired-end reads as pool.')
    group_need.add_option('-s', dest='seed_file', help='Reference. Input fasta format file as initial seed '
                                                       'or input bowtie index base name as pre-seed (see flag "--bs")')
    group_need.add_option('-w', dest='word_size', type=float, default=0.7,
                          help='Word size (W) for extension. You could assign the ratio (1>input>0) of W to '
                               'read_length, based on which this script would estimate the W for you; '
                               'or assign an absolute W value (read length-1>input>=21). Default:0.7')
    group_need.add_option('-o', dest='output_base', help='Out put directory. Overwriting files if directory exists.')
    # group2
    group_result = OptionGroup(parser, "INFLUENTIAI OPTIONS", "These option will affect the final results"
                                                              " or serve as alternations of the required options")
    group_result.add_option('-R', dest='round_limit', type=int,
                            help='Limit running rounds (>=2).')
    group_result.add_option('-u', dest='unpaired_fastq_files',
                            help='Input file(s) with unpaired reads as pool. '
                                 'files could be comma-separated lists such as "seq1,seq2".')
    group_result.add_option('--max-reads', dest='maximum_n_reads', type=float, default=1E7,
                            help="Maximum number of reads to be used per file. Default: 1E7 (-F cp,nr) or 5E7 (-F mt)")
    group_result.add_option('--max-words', dest='maximum_n_words', type=float, default=5E7,
                            help="Maximum number of words to be used in total."
                                 "Default: 5E7 (-F cp), 1E7 (-F nr) or 3E8 (-F mt)")
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
    group_result.add_option('-F', dest='organelle_type', default='cp',
                            help='This flag should be followed with cp (if you want get chloroplast), '
                                 'mt (mitochondria), nr (nuclear ribosomal RNA), 0 (disable this). Default: cp. '
                                 'You can also make the index by your self and add those index to ' +
                                 os.path.join(path_of_this_script, 'Library', '/NotationReference') + '')
    group_result.add_option("--disentangle-df", dest="disentangle_depth_factor", default=10.0, type=float,
                            help="Depth factor for differentiate genome type of contigs. Default:%default")
    group_result.add_option("--contamination-depth", dest="contamination_depth", default=5., type=float,
                            help="Depth factor for confirming contamination in parallel contigs. Default:%default")
    group_result.add_option("--contamination-similarity", dest="contamination_similarity", default=0.9, type=float,
                            help="Similarity threshold for confirming contaminating contigs. Default:%default")
    group_result.add_option("--no-degenerate", dest="degenerate", default=True, action="store_false",
                            help="Disable making consensus from parallel contig based on nucleotide degenerate table.")
    group_result.add_option("--degenerate-depth", dest="degenerate_depth", default=1.5, type=float,
                            help="Depth factor for confirming parallel contigs. Default:%default")
    group_result.add_option("--degenerate-similarity", dest="degenerate_similarity", default=0.95, type=float,
                            help="Similarity threshold for confirming parallel contigs. Default:%default")
    # group_result.add_option("--trim", dest='trim_values',
    #                         help='Assign the number of bases in the ends to trim in extending process. '
    #                              'This function will not change the length of the out put reads. '
    #                              'Input format: int,int (Example: 4,4). Default: 0,0')
    group_result.add_option("--min-quality-score", dest="min_quality_score", type=int, default=15,
                            help="Minimum quality score in extending extension. "
                                 "Default:%default ('+' in Phred+33; 'J' in Phred+64/Solexa+64)")
    group_result.add_option('-k', dest='spades_kmer', default='65,75,85',
                            help='SPAdes kmer settings. Use the same format as in SPAdes. Default=65,75,85')
    group_result.add_option('--spades-options', dest='other_spades_options', default='',
                            help='Other SPAdes options. Use double quotation marks to include all the arguments '
                                 'and parameters, such as "--careful -o test"')
    group_result.add_option('--no-spades', dest='run_spades', action="store_false", default=True,
                            help='Disable SPAdes.')
    group_result.add_option('--no-bowtie2', dest='utilize_mapping', action="store_false", default=True,
                            help='Choose to disable mapping process.'
                                 'By default, this script would map original reads to pre-seed (or bowtie2 index) '
                                 'to acquire the initial seed. This requires bowtie2 to be installed '
                                 '(and thus Linux/Unix only). It is suggested to keep mapping as default '
                                 'when the seed is too diverse from target.')
    # group3
    group_computational = OptionGroup(parser, "ADDITIONAL OPTIONS", "These options will NOT affect the final results. "
                                                                    "Take easy to pick some according your computer's flavour")
    group_computational.add_option('-t', dest='threads', type=int, default=4,
                                   help="Maximum threads to use.")
    group_computational.add_option('-P', dest='pre_grouped', type=float, default=2E5,
                                   help='The maximum number (integer) of high-covered reads to be pre-grouped '
                                        'before extension. '
                                        'pre_grouping is suggested when the whole genome coverage is shallow but '
                                        'the organ genome coverage is deep.'
                                        'The default value is 2E5. '
                                        'For personal computer with 8G memory, we suggest no more than 3E5. '
                                        'A larger number (ex. 6E5) would run faster but exhaust memory '
                                        'in the first few minutes. Choose 0 to disable this process.')
    group_computational.add_option('--continue', dest='script_resume', default=False, action="store_true",
                                   help='Several check point based on files produced, rather than log, '
                                        'so keep in mind that this script will not detect the difference '
                                        'between this input parameters and the previous ones.')
    group_computational.add_option('--index-in-memory', dest='index_in_memory', action="store_true", default=False,
                                   help="Keep index in memory. Choose save index in memory than in disk.")
    group_computational.add_option('--out-per-round', dest='fg_out_per_round', action="store_true", default=False,
                                   help='Enable output per round. Choose to save memory but cost more time per round.')
    group_computational.add_option('--remove-duplicates', dest='rm_duplicates', default=1E7, type=float,
                                   help='By default this script use unique reads to extend. Choose the number of '
                                        'duplicates (integer) to be saved in memory. A larger number (ex. 2E7) would '
                                        'run faster but exhaust memory in the first few minutes. '
                                        'Choose 0 to disable this process. '
                                        'Note that whether choose or not will not disable '
                                        'the calling of replicate reads. Default: %default.')
    group_computational.add_option("--prefix", dest="prefix", default="",
                                   help="Add extra prefix to resulting files under the output directory.")
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
        sys.stdout.write('\n############################################################################' + str(e))
        sys.stdout.write('\n"-h" for more usage')
        exit()
    else:
        if not ((options.seed_file or options.bowtie2_seed) and ((
                                                                         options.fastq_file_1 and options.fastq_file_2) or options.unpaired_fastq_files) and options.word_size and options.output_base):
            parser.print_help()
            sys.stdout.write('\n############################################################################'
                             '\nERROR: Insufficient arguments!\nUsage:')
            sys.stdout.write(usage + '\n\n')
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
        options.prefix = os.path.basename(options.prefix)
        log = simple_log(logging.getLogger(), options.output_base, options.prefix + "get_org.")
        log.info(print_title)
        log.info(' '.join(sys.argv) + '\n')
        log = timed_log(log, options.output_base, options.prefix + "get_org.")
        if 0 < options.word_size < 1:
            pass
        elif options.word_size >= 21:
            options.word_size = int(options.word_size)
        else:
            log.error("Illegal '-w' value!")
            exit()
        if "--max-reads" not in sys.argv:
            if options.organelle_type == "mt":
                options.maximum_n_reads *= 5
                # log.info("--max-reads " + str(options.maximum_n_reads) + " (mt)")
        if "--max-words" not in sys.argv:
            if options.organelle_type == "mt":
                options.maximum_n_words *= 6
            elif options.organelle_type == "nr":
                options.maximum_n_words /= 5
                # log.info()
        if options.seed_file and options.bowtie2_seed:
            log.error('Simultaneously using "-s" and "--bs" is not allowed!')
            exit()
        if options.anti_seed and options.bowtie2_anti_seed:
            log.error('Simultaneously using "-a" and "--as" is not allowed!')
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
        options.rm_duplicates = int(options.rm_duplicates)
        options.pre_grouped = int(options.pre_grouped)
        if not options.rm_duplicates and options.pre_grouped:
            log.warning("removing duplicates was inactive, so that the pre-grouping was disabled.")
            options.pre_grouped = False
        if options.round_limit and options.round_limit < 2:
            log.warning("illegal limit for rounds! Been set to default: unlimited.")
            options.round_limit = None
        return options, log


# test whether an external binary is executable
def executable(test_this):
    return True if os.access(test_this, os.X_OK) or getstatusoutput(test_this)[0] != dead_code else False


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
            for i in range(0, this_seq_len - word_size + 1):
                forward = seed[i:i + word_size]
                return_words.add(forward)
                reverse = cpt_seed[i:i + word_size]
                return_words.add(reverse)
    return return_words


def chop_seq_list(seq_generator_or_list):
    return_words = set()
    for seed in seq_generator_or_list:
        for seq_part in seed:
            this_seq_len = len(seq_part)
            if this_seq_len >= word_size:
                cpt_seed = complementary_seq(seq_part)
                for i in range(0, this_seq_len - word_size + 1):
                    forward = seq_part[i:i + word_size]
                    return_words.add(forward)
                    reverse = cpt_seed[i:i + word_size]
                    return_words.add(reverse)
    return return_words


def get_read_len_mean_max(fq_files, maximum_n_reads):
    read_lengths = []
    for fq_f in fq_files:
        count_r = 0
        for seq in fq_seq_simple_generator(fq_f):
            count_r += 1
            read_lengths.append(len(seq.strip("N")))
            if count_r >= maximum_n_reads:
                break
    return sum(read_lengths)/len(read_lengths), max(read_lengths)


def get_low_quality_char_pattern(fq_files, maximum_n_reads, digit_threshold, log, sampling_percent=0.1):
    max_quality = "!"
    min_quality = "~"
    sampling_percent = int(1/sampling_percent)
    for fq_f in fq_files:
        count_r = 0
        for quality_str in fq_seq_simple_generator(fq_f, go_to_line=3):
            if count_r % sampling_percent == 0:
                max_quality = max(max_quality, max(quality_str))
                min_quality = min(max_quality, min(quality_str))
            count_r += 1
            if count_r >= maximum_n_reads:
                break
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
    the_name, the_c_min, the_c_max, the_s_min, the_s_max, deviation = sorted(decision_making, key=lambda x: x[-1])[0]
    log.info("Identified quality encoding format: " + the_name)
    if max_quality > the_c_max:
        log.warning("Max quality score " + repr(chr(max_quality)) +
                    "(" + str(max_quality) + ":" + str(max_quality - the_c_min + the_s_min) +
                    ") in your fastq file exceeds the usual boundary " + str((the_c_min, the_c_max)))
    if min_quality < the_c_min:
        log.warning("Min quality score " + repr(chr(min_quality)) +
                    "(" + str(min_quality) + ":" + str(min_quality - the_c_min + the_s_min) +
                    ") in your fastq file is under the usual lower boundary " + str((the_c_min, the_c_max)))
    low_quality_chars = [chr(low_quality_score)
                         for low_quality_score in range(min_quality, the_c_min + digit_threshold - the_s_min)]
    return "[" + "".join(low_quality_chars) + "]"


def write_fq_results(original_fq_files, accepted_contig_id, out_file_name, temp2_clusters_dir, fq_info_in_memory,
                     maximum_n_reads, verbose, index_in_memory, log):
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
                    anti_lines, pre_grouped, index_in_memory, bowtie2_anti_seed, anti_seed, keep_seq_parts,
                    low_quality, echo_frequency, resume, log):
    # if trim_values:
    #     trim1, trim2 = [int(trim_value) for trim_value in trim_values.split(',')]
    # else:
    #     trim1, trim2 = 0, 0

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
        lengths = []
        use_user_direction = False
        for id_file, file_name in enumerate(original_fq_files):
            file_in = open(file_name, "rU")
            count_this_read_n = 0
            line = file_in.readline()
            # if anti seed input, name & direction should be recognized
            if bowtie2_anti_seed or anti_seed:
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
                        # if trim_values:
                        #     this_seq = this_seq[trim1:(len(this_seq) - trim2)].strip("N")
                        # else:
                        #     this_seq = this_seq.strip("N")
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

                        # if trim_values:
                        #     this_seq = this_seq[trim1:(len(this_seq) - trim2)].strip("N")
                        # else:
                        #     this_seq = this_seq.strip("N")

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
                log.warning("Number of reads exceeded " + str(int(maximum_n_reads)) + " in " + file_name + ", only top "
                            + str(int(maximum_n_reads)) + " reads are used in downstream analysis (suggested).")
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


def pre_grouping(fastq_indices_in_memory, dupli_threshold, out_base, index_in_memory, log):
    global word_size
    forward_and_reverse_reads, line_clusters, len_indices, keep_seq_parts = fastq_indices_in_memory
    log.info("Pre-grouping reads...")
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
                # skip those reads that are not unique/represented by others
                while here_go_to[0] < 2 * here_unique_id:
                    forward_seq_line = temp_seq_file.readline()
                    here_go_to[0] += 1
                reverse_seq_line = temp_seq_file.readline()
                here_go_to[0] += 1
                return forward_seq_line.strip().split("\t"), reverse_seq_line.strip().split("\t")
        else:
            def generate_forward_and_reverse(here_unique_id):
                forward_seq_line = temp_seq_file.readline()
                # skip those reads that are not unique/represented by others
                while here_go_to[0] < 2 * here_unique_id:
                    forward_seq_line = temp_seq_file.readline()
                    here_go_to[0] += 1
                reverse_seq_line = temp_seq_file.readline()
                here_go_to[0] += 1
                return forward_seq_line.strip(), reverse_seq_line.strip()

    for this_unique_read_id in list(lines_with_duplicates):
        this_seq, this_c_seq = generate_forward_and_reverse(this_unique_read_id)
        these_group_id = set()
        this_words = []
        if keep_seq_parts:
            for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
                seq_len = len(this_seq_part)
                temp_length = seq_len - word_size
                for i in range(0, temp_length + 1):
                    forward = this_seq_part[i:i + word_size]
                    reverse = this_c_seq_part[temp_length - i:seq_len - i]
                    if forward in these_words:
                        these_group_id.add(these_words[forward])
                    else:
                        this_words.append(forward)
                        this_words.append(reverse)
        else:
            seq_len = len(this_seq)
            temp_length = seq_len - word_size
            for i in range(0, temp_length + 1):
                forward = this_seq[i:i + word_size]
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


class NoMoreReads(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def extending_reads(accepted_words, accepted_contig_id, original_fq_dir, len_indices, pre_grouped,
                    groups_of_duplicate_lines, lines_with_duplicates, fq_info_in_memory, output_base,
                    round_limit, fg_out_per_round, jump_step, mesh_size, verbose, resume,
                    maximum_n_reads, maximum_n_words, keep_seq_parts, low_quality_pattern, echo_frequency, log):
    global word_size
    accumulated_num_words = 0
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
        def summarise_round(acc_words, acc_contig_id_this_round, pre_aw, r_count, acc_num_words, unique_id):
            len_aw = len(acc_words)
            len_al = len(acc_contig_id_this_round)
            acc_num_words += len_aw - pre_aw
            if this_process:
                memory_usage = " Mem " + str(round(this_process.memory_info().rss / 1024.0 / 1024 / 1024, 3))
            else:
                memory_usage = ''
            if fg_out_per_round:
                write_fq_results(original_fq_dir, acc_contig_id_this_round,
                                 os.path.join(round_dir, "Round." + str(r_count)),
                                 os.path.join(output_base, 'temp.indices.2'), fq_info_in_memory, maximum_n_reads,
                                 verbose, bool(fq_info_in_memory), log)
                # clear former accepted words from memory
                del acc_words
                # then add new accepted words into memory
                acc_words = chop_seqs(fq_seq_simple_generator(
                    [os.path.join(round_dir, "Round." + str(r_count) + '_' + str(x + 1) + '.fq') for x in
                     range(len(original_fq_dir))], split_pattern=low_quality_pattern, min_sub_seq=word_size))
                acc_contig_id_this_round = set()
            log.info("Round " + str(r_count) + ': ' + str(unique_id + 1) + '/' + str(len_indices) + " AI " + str(
                len_al) + " AW " + str(len_aw) + memory_usage)
            #
            if len_aw == pre_aw:
                raise NoMoreReads('')
            pre_aw = len(acc_words)
            #
            if r_count == round_limit:
                raise RoundLimitException(r_count)
            r_count += 1
            return acc_words, acc_contig_id_this_round, pre_aw, r_count, acc_num_words

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
                if pre_grouped and groups_of_duplicate_lines:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_contig_id:
                            if unique_read_id in line_to_accept:
                                accepted_contig_id.add(unique_read_id)
                                accepted_contig_id_this_round.add(unique_read_id)
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
                                for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
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
                                    accepted_contig_id.add(unique_read_id)
                                    accepted_contig_id_this_round.add(unique_read_id)
                                    if unique_read_id in lines_with_duplicates:
                                        which_group = lines_with_duplicates[unique_read_id]
                                        for id_to_accept in groups_of_duplicate_lines[which_group]:
                                            line_to_accept.add(id_to_accept)
                                            del lines_with_duplicates[id_to_accept]
                                        line_to_accept.remove(unique_read_id)
                                        del groups_of_duplicate_lines[which_group]
                        if unique_read_id % echo_frequency == 0:
                            this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round "\
                                         + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                         " AI " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words))
                            sys.stdout.write(this_print + '\b' * len(this_print))
                            sys.stdout.flush()
                            if accumulated_num_words + len(accepted_words) - previous_aw_count > maximum_n_words:
                                log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                         str(len_indices) + " AI " + str(len(accepted_contig_id_this_round)) +
                                         " AW " + str(len(accepted_words)))
                                raise WordsLimitException("")
                    accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count, accumulated_num_words\
                        = summarise_round(accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count,
                                          accumulated_num_words, unique_read_id)
                else:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_contig_id:
                            accepted = False
                            for this_seq_part, this_c_seq_part in zip(this_seq, this_c_seq):
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
                                accepted_contig_id.add(unique_read_id)
                                accepted_contig_id_this_round.add(unique_read_id)
                        if unique_read_id % echo_frequency == 0:
                            this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round " \
                                         + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                         " AI " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words))
                            sys.stdout.write(this_print + '\b' * len(this_print))
                            sys.stdout.flush()
                            if accumulated_num_words + len(accepted_words) - previous_aw_count > maximum_n_words:
                                log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                         str(len_indices) + " AI " + str(len(accepted_contig_id_this_round)) + " AW "
                                         + str(len(accepted_words)))
                                raise WordsLimitException("")
                    accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count, accumulated_num_words\
                        = summarise_round(accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count,
                                          accumulated_num_words, unique_read_id)
            else:
                if pre_grouped and groups_of_duplicate_lines:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_contig_id:
                            seq_len = len(this_seq)
                            temp_length = seq_len - word_size
                            if unique_read_id in line_to_accept:
                                accepted_contig_id.add(unique_read_id)
                                accepted_contig_id_this_round.add(unique_read_id)
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
                                    accepted_contig_id.add(unique_read_id)
                                    accepted_contig_id_this_round.add(unique_read_id)
                                    if unique_read_id in lines_with_duplicates:
                                        which_group = lines_with_duplicates[unique_read_id]
                                        for id_to_accept in groups_of_duplicate_lines[which_group]:
                                            line_to_accept.add(id_to_accept)
                                            del lines_with_duplicates[id_to_accept]
                                        line_to_accept.remove(unique_read_id)
                                        del groups_of_duplicate_lines[which_group]
                        if unique_read_id % echo_frequency == 0:
                            this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round "\
                                         + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                         " AI " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words))
                            sys.stdout.write(this_print + '\b' * len(this_print))
                            sys.stdout.flush()
                            if accumulated_num_words + len(accepted_words) - previous_aw_count > maximum_n_words:
                                log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                         str(len_indices) + " AI " + str(len(accepted_contig_id_this_round)) +
                                         " AW " + str(len(accepted_words)))
                                raise WordsLimitException("")
                    accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count, accumulated_num_words\
                        = summarise_round(accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count,
                                          accumulated_num_words, unique_read_id)
                else:
                    for unique_read_id in range(len_indices):
                        this_seq = next(reads_generator)
                        this_c_seq = next(reads_generator)
                        if unique_read_id not in accepted_contig_id:
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
                                accepted_contig_id.add(unique_read_id)
                                accepted_contig_id_this_round.add(unique_read_id)
                        if unique_read_id % echo_frequency == 0:
                            this_print = str("%s" % datetime.datetime.now())[:23].replace('.', ',') + " - INFO: Round " \
                                         + str(round_count) + ': ' + str(unique_read_id + 1) + '/' + str(len_indices) + \
                                         " AI " + str(len(accepted_contig_id_this_round)) + " AW " + str(len(accepted_words))
                            sys.stdout.write(this_print + '\b' * len(this_print))
                            sys.stdout.flush()
                            if accumulated_num_words + len(accepted_words) - previous_aw_count > maximum_n_words:
                                log.info("Round " + str(round_count) + ': ' + str(unique_read_id + 1) + '/' +
                                         str(len_indices) + " AI " + str(len(accepted_contig_id_this_round)) + " AW "
                                         + str(len(accepted_words)))
                                raise WordsLimitException("")
                    accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count, accumulated_num_words\
                        = summarise_round(accepted_words, accepted_contig_id_this_round, previous_aw_count, round_count,
                                          accumulated_num_words, unique_read_id)
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
        log.info("Hit the round limit " + str(r_lim) + " and terminated ...")
    except WordsLimitException:
        reads_generator.close()
        log.info("Hit the words limit and terminated ...")
    del reads_generator

    del accepted_words
    del accepted_contig_id_this_round
    del lines_with_duplicates
    return accepted_contig_id


def get_anti_with_fas(anti_words, anti_input, original_fq_files, log):
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


def mapping_with_bowtie2(seed_file, bowtie2_seed, anti_seed, bowtie2_anti_seed, original_fq_files, out_base, resume,
                         verbose_log, threads, prefix, log):
    if seed_file:
        if os.path.exists(seed_file + '.index.1.bt2l'):
            log.info("Bowtie2 index existed!")
        else:
            log.info("Making seed - bowtie2 index ...")
            build_seed_index = subprocess.Popen("bowtie2-build --large-index " + seed_file + " " + seed_file + '.index',
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
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
    else:
        seed_index_base = bowtie2_seed
    total_seed_file = [os.path.join(out_base, x + prefix + "Initial.mapped.fq") for x in ("temp.", "")]
    total_seed_sam = [os.path.join(out_base, x + prefix + "seed_bowtie.sam") for x in ("temp.", "")]
    if resume and os.path.exists(total_seed_file[1]):
        log.info("Initial seeds existed!")
    else:
        log.info("Mapping reads to seed - bowtie2 index ...")
        make_seed_bowtie2 = subprocess.Popen(
            "bowtie2 -p " + str(threads) + " --very-fast-local --al " + total_seed_file[
                0] + " -x " + seed_index_base + " -U " +
            ",".join(original_fq_files) + " -S " + total_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            log.error('\n' + output.decode("utf8"))
            exit()
        if verbose_log:
            log.info("\n" + output.decode("utf8").strip())
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
                log.info("anti-seed - bowtie2 index existed!")
            else:
                log.info("Making anti-seed - bowtie2 index ...")
                build_anti_index = subprocess.Popen(
                    "bowtie2-build --large-index " + anti_seed + " " + anti_seed + '.index',
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                output, err = build_anti_index.communicate()
                if "unrecognized option" in output.decode("utf8"):
                    build_anti_index = subprocess.Popen("bowtie2-build " + anti_seed + " " + anti_seed + '.index',
                                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                    output, err = build_anti_index.communicate()
                if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
                    log.error('\n' + output.decode("utf8"))
                    exit()
                if verbose_log:
                    log.info("\n" + output.decode("utf8").strip())
                log.info("Making anti-seed - bowtie2 index finished.")
            anti_index_base = anti_seed + '.index'
        else:
            anti_index_base = bowtie2_anti_seed

        anti_seed_sam = [os.path.join(out_base, x + prefix + "anti_seed_bowtie.sam") for x in ("temp.", "")]
        if resume and os.path.exists(anti_seed_sam[1]):
            log.info("Anti-seed mapping information existed!")
        else:
            log.info("Mapping reads to anti-seed - bowtie2 index ...")
            make_anti_seed_bowtie2 = subprocess.Popen(
                "bowtie2 -p " + str(threads) + " --very-fast-local -x " + anti_index_base + " -U " +
                ",".join(original_fq_files) + " -S " +
                anti_seed_sam[0] + " --no-unal --no-hd --no-sq -t",
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
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
    return total_seed_file[1], anti_lines


def assembly_with_spades(spades_kmer, spades_out_put, parameters, out_base, prefix, original_fq_files, reads_paired,
                         verbose_log, resume, threads, log):
    if '-k' in parameters:
        kmer = ''
    else:
        kmer = '-k ' + spades_kmer
    if resume and os.path.exists(spades_out_put):
        continue_command = '--continue'
    else:
        continue_command = ''
    spades_out_put = '-o ' + spades_out_put
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
                ['spades.py', '-t', str(threads), continue_command, parameters, '-1',
                 os.path.join(out_base, prefix + "filtered_1_paired.fq"), '-2',
                 os.path.join(out_base, prefix + "filtered_2_paired.fq")] +
                ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                [kmer, spades_out_put]).strip()
        else:
            log.warning("No paired reads found for the target!?")
            spades_command = ' '.join(
                ['spades.py', '-t', str(threads), continue_command, parameters] +
                ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
                [kmer, spades_out_put]).strip()
    else:
        all_unpaired = []
        for iter_unpaired in range(len(original_fq_files)):
            if os.path.getsize(str(os.path.join(out_base, prefix + "filtered_" + str(iter_unpaired + 1) + ".fq"))):
                all_unpaired.append(str(os.path.join(out_base, prefix + "filtered_" + str(iter_unpaired + 1) + ".fq")))
        spades_command = ' '.join(
            ['spades.py', '-t', str(threads), continue_command, parameters] +
            ['--s' + str(i + 1) + ' ' + out_f for i, out_f in enumerate(all_unpaired)] +
            [kmer, spades_out_put]).strip()
    spades_running = subprocess.Popen(spades_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = spades_running.communicate()
    if "not recognized" in output.decode("utf8"):
        if verbose_log:
            log.error('Problem with running SPAdes:')
            log.error(output.decode("utf8"))
        log.error("WAINING in SPAdes: unrecognized option in " + parameters)
        log.error('Assembling failed.')
        return False
    elif "== Error ==" in output.decode("utf8"):
        log.error("Error in SPAdes: \n== Error ==" + output.decode("utf8").split("== Error ==")[-1].split("In case you")[0])
        log.error('Assembling failed.')
        return False
    else:
        if verbose_log:
            log.info(output.decode("utf8"))
        log.info('Assembling finished.\n')
        return True


def slim_spades_result(scheme, spades_output, verbose_log, log, threads, depth_threshold=0, resume=False):
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
        'cp': ' --include-priority ' + os.path.join(path_of_this_script, 'Library', 'NotationReference',
                                                    'cp') + ' --exclude ' + os.path.join(path_of_this_script, 'Library',
                                                                                         'NotationReference', 'mt'),
        'mt': ' --include-priority ' + os.path.join(path_of_this_script, 'Library', 'NotationReference',
                                                    'mt') + ' --exclude ' + os.path.join(path_of_this_script, 'Library',
                                                                                         'NotationReference', 'cp'),
        'nr': ' --include-priority ' + os.path.join(path_of_this_script, 'Library', 'NotationReference', 'nr')}
    if scheme in scheme_translation:
        run_command = scheme_translation[scheme]
    else:
        run_command = scheme
    run_command = os.path.join(path_of_this_script, 'Utilities', 'slim_fastg.py') + ' -t ' + str(threads) + ' ' \
                  + graph_file + run_command + ' --depth-threshold ' + str(depth_threshold)  #\
                  # + ' -o ' + out_base + (' --prefix ' + prefix if prefix else "")
    slim_spades = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = slim_spades.communicate()
    if "not recognized" in output.decode("utf8") or "command not found" in output.decode("utf8"):
        if log:
            if verbose_log:
                log.warning(os.path.join(path_of_this_script, "Utilities", "slim_spades_fastg_by_blast.py") + ' not found!')
                log.warning(output.decode("utf8"))
            log.warning("Slimming      " + graph_file + " failed.")
        return 1
    elif "failed" in output.decode("utf8") or "error" in output.decode("utf8"):
        if log:
            if verbose_log:
                log.error(output.decode("utf8"))
            log.warning("Slimming      " + graph_file + " failed.")
        return 1
    else:
        if log:
            if verbose_log:
                log.info(output.decode("utf8"))
            log.info("Slimming      " + graph_file + " finished!")
        return 0


def separate_fq_by_pair(out_base, prefix, verbose_log, log):
    log.info("Separating filtered fastq file ... ")
    run_command = os.path.join(path_of_this_script, "Utilities", "get_pair_reads.py") \
                  + ' ' + os.path.join(out_base, prefix + "filtered_1.fq") \
                  + ' ' + os.path.join(out_base, prefix + "filtered_2.fq")
    separating = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = separating.communicate()
    if "not recognized" in output.decode("utf8") or "command not found" in output.decode("utf8"):
        if verbose_log:
            log.warning(os.path.join(path_of_this_script, "Utilities", "get_pair_reads.py") + "not found!")
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


def unzip(source, target, verbose_log, log):
    target_temp = target + ".Temp"
    try_commands = ["gunzip -c " + source + " > " + target_temp, "tar -x -f " + source + " -O > " + target_temp]
    log.info("Unzipping reads file: " + source)
    success = False
    output = b""
    for run_command in try_commands:
        unzipping = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
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


@set_time_limit(300)
def disentangle_circular_assembly(fastg_file, tab_file, output, weight_factor, log, type_factor=3.,
                                  mode="cp", verbose=False,
                                  contamination_depth=5., contamination_similarity=5.,
                                  degenerate=True, degenerate_depth=1.5, degenerate_similarity=1.5,
                                  hard_cov_threshold=10.,
                                  min_sigma_factor=0.1):
    from Library.assembly_parser import Assembly
    log.info("Disentangling " + fastg_file + " ... ")
    input_graph = Assembly(fastg_file)
    target_results = input_graph.find_target_graph(tab_file,
                                                   mode=mode, type_factor=type_factor,
                                                   log_hard_cov_threshold=hard_cov_threshold,
                                                   contamination_depth=contamination_depth,
                                                   contamination_similarity=contamination_similarity,
                                                   degenerate=degenerate, degenerate_depth=degenerate_depth,
                                                   degenerate_similarity=degenerate_similarity,
                                                   min_sigma_factor=min_sigma_factor,
                                                   weight_factor=weight_factor,
                                                   log_handler=log, verbose=verbose)
    if len(target_results) > 1:
        log.warning(str(len(target_results)) + " sets of graph detected!")
    for go_res, res in enumerate(target_results):
        go_res += 1
        idealized_graph = res["graph"]
        # average_kmer_cov = res["cov"]
        # log.info("Detecting target graph" + str(go_res) +
        #          " finished with average kmer coverage: " + str(round(average_kmer_cov, 4)))

        # should add making one-step-inversion pairs for paths,
        # which would be used to identify existence of a certain isomer using mapping information
        count_path = 0
        for this_path, other_tag in idealized_graph.get_all_circular_paths(mode=mode, log_handler=log):
            count_path += 1
            out_n = output + ".graph" + str(go_res) + "." + str(count_path) + other_tag + ".path_sequence.fasta"
            open(out_n, "w").write(idealized_graph.export_path(this_path).fasta_str())
            log.info("Writing PATH" + str(count_path) + " to " + out_n)
        log.info("Writing GRAPH to " + output + ".graph" + str(go_res) + ".selected_graph.gfa")
        idealized_graph.write_to_gfa(output + ".graph" + str(go_res) + ".selected_graph.gfa")

    log.info("Solving and unfolding graph finished!")
    log.info("Please visualize " + fastg_file + " to confirm the final result.")


def extract_organelle_genome(out_base, spades_output, go_round,
                             prefix, organelle_type, verbose, log, threads, options):
    running_stat = slim_spades_result(organelle_type, spades_output, verbose, log, threads=threads,
                                      resume=options.script_resume)
    """disentangle"""
    if running_stat == 0:
        out_fastg = sorted([os.path.join(spades_output, x)
                            for x in os.listdir(spades_output) if x.count(".fastg") == 2])[0]
        out_csv = out_fastg[:-5] + "csv"
        # if it is the first round (the largest kmer), copy the slimmed result to the main spades output folder
        if go_round == 0:
            main_spades_folder = os.path.split(spades_output)[0]
            os.system("cp " + out_fastg + " " + main_spades_folder)
            os.system("cp " + out_csv + " " + main_spades_folder)

        path_prefix = os.path.join(out_base, prefix + organelle_type)
        disentangle_circular_assembly(fastg_file=out_fastg, mode=organelle_type,
                                      tab_file=out_csv, output=path_prefix, weight_factor=100,
                                      hard_cov_threshold=options.disentangle_depth_factor,
                                      contamination_depth=options.contamination_depth,
                                      contamination_similarity=options.contamination_similarity,
                                      degenerate=options.degenerate,
                                      degenerate_depth=options.degenerate_depth,
                                      degenerate_similarity=options.degenerate_similarity,
                                      verbose=verbose, log=log)
        return 0
    else:
        return 1


def main():
    time0 = time.time()
    title = "\nGetOrganelle v" + str(get_versions()) + \
            "\n" \
            "\nThis pipeline get organelle reads and genomes from genome skimming data by extending." \
            "\nFind updates in https://github.com/Kinggerm/GetOrganelle and see README.md for more information." \
            "\n"
    options, log = require_commands(print_title=title, version=get_versions())
    resume = options.script_resume
    verb_out = options.verbose_log
    out_base = options.output_base
    reads_files_to_drop = []
    try:
        """ initialization """
        global word_size
        word_size = options.word_size
        if options.fastq_file_1 and options.fastq_file_2:
            reads_paired = {'input': True, 'pair_out': bool}
            original_fq_files = [options.fastq_file_1, options.fastq_file_2] + \
                                [fastq_file for fastq_file in options.unpaired_fastq_files]
            direction_according_to_user_input = [1, 2] + [1] * len(options.unpaired_fastq_files)
        else:
            reads_paired = {'input': False, 'pair_out': False}
            original_fq_files = [fastq_file for fastq_file in options.unpaired_fastq_files]
            direction_according_to_user_input = [1] * len(options.unpaired_fastq_files)
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
        if not (resume and filtered_files_exist):

            seed_file = options.seed_file
            bowt_seed = options.bowtie2_seed
            anti_seed = options.anti_seed
            b_at_seed = options.bowtie2_anti_seed
            pre_grp = options.pre_grouped
            # trim_ends = options.trim_values
            in_memory = options.index_in_memory

            if original_fq_files:
                for file_id, read_file in enumerate(original_fq_files):
                    # unzip fq files if needed
                    if read_file.endswith(".gz") or read_file.endswith(".zip"):
                        # log.info("Unzipping " + read_file + " ...")
                        target_fq = read_file + ".fastq"
                        if not (os.path.exists(target_fq) and resume):
                            unzip(read_file, target_fq, options.verbose_log, log)
                        original_fq_files[file_id] = target_fq
                        reads_files_to_drop.append(target_fq)
                        # log.info("Unzipping " + read_file + " finished.\n")

            # pre-reading fastq
            log.info("Pre-reading fastq ...")
            log.info("Counting read qualities ...")
            low_quality_pattern = get_low_quality_char_pattern(original_fq_files,
                                                               options.maximum_n_reads,
                                                               options.min_quality_score, log, sampling_percent=0.1)
            log.info("Counting read lengths ...")
            mean_read_len, max_read_len = get_read_len_mean_max(original_fq_files, options.maximum_n_reads)
            if word_size < 1:
                new_word_size = int(word_size * mean_read_len)
                if new_word_size < 21:
                    word_size = 21
                    log.warning("Too small ratio " + str(word_size) + ", setting word_size = 21.")
                else:
                    word_size = min(new_word_size, {"nr": 141, "cp": 121, "mt": 101}[options.organelle_type])
                    log.info("Setting word_size = " + str(word_size))
            if float(word_size) / max_read_len <= 0.5 and len(low_quality_pattern) > 2:
                keep_seq_parts = True
            else:
                keep_seq_parts = False
            log.info("Pre-reading fastq finished.\n")

            # reading seeds
            log.info("Making seed reads ...")
            if not options.utilize_mapping:
                anti_lines = get_anti_with_fas(chop_seqs(read_fasta(anti_seed)[1]),
                                               (anti_seed or b_at_seed),
                                               original_fq_files, log) if anti_seed else set()
            else:
                seed_fastq, anti_lines = mapping_with_bowtie2(seed_file, bowt_seed, anti_seed,
                                                              b_at_seed, original_fq_files,
                                                              out_base, resume,
                                                              verb_out, options.threads, options.prefix, log)
            log.info("Making seed reads finished.\n")

            # make read index
            log.info("Making read index ...")
            echo_frequency = 54321
            fq_info_in_memory = make_read_index(original_fq_files, direction_according_to_user_input,
                                                options.maximum_n_reads, options.rm_duplicates, out_base,
                                                anti_lines, pre_grp, in_memory, b_at_seed, anti_seed,
                                                keep_seq_parts=keep_seq_parts, low_quality=low_quality_pattern,
                                                resume=resume, echo_frequency=echo_frequency, log=log)
            len_indices = fq_info_in_memory[2]
            keep_seq_parts = fq_info_in_memory[3]
            if keep_seq_parts:
                log.info("Reads are stored as fragments.")
            # pre-grouping if asked
            if pre_grp:
                groups_of_lines, lines_in_a_group = pre_grouping(fq_info_in_memory, pre_grp, out_base, in_memory, log)
            else:
                groups_of_lines = None
                lines_in_a_group = None
            if not in_memory:
                fq_info_in_memory = None
            log.info("Making read index finished.\n")

            # adding initial word
            log.info("Adding initial words ...")
            if not options.utilize_mapping:
                initial_accepted_words = chop_seqs(read_fasta(seed_file)[1])
            else:
                if keep_seq_parts:
                    initial_accepted_words = chop_seq_list(
                        fq_seq_simple_generator(seed_fastq, split_pattern=low_quality_pattern, min_sub_seq=word_size))
                else:
                    initial_accepted_words = chop_seqs(fq_seq_simple_generator(seed_fastq))
            log.info("Adding initial words finished.\n")

            # extending process
            log.info("Extending ...")
            accepted_ids = set()
            echo_frequency = 1234321//(mean_read_len - word_size + 1)//2
            accepted_contig_id = extending_reads(initial_accepted_words, accepted_ids, original_fq_files, len_indices,
                                                 pre_grp, groups_of_lines, lines_in_a_group,
                                                 fq_info_in_memory, out_base, options.round_limit,
                                                 options.fg_out_per_round,
                                                 options.jump_step,
                                                 options.mesh_size, verb_out, resume,
                                                 options.maximum_n_reads, options.maximum_n_words,
                                                 keep_seq_parts=keep_seq_parts,
                                                 low_quality_pattern=low_quality_pattern,
                                                 echo_frequency=echo_frequency,
                                                 log=log)
            write_fq_results(original_fq_files, accepted_contig_id,
                             os.path.join(out_base, options.prefix + "filtered"),
                             os.path.join(out_base, 'temp.indices.2'),
                             fq_info_in_memory, options.maximum_n_reads, verb_out, in_memory, log)
            del accepted_contig_id, fq_info_in_memory, groups_of_lines, \
                anti_lines, initial_accepted_words, lines_in_a_group

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
                reads_paired['pair_out'] = separate_fq_by_pair(out_base, options.prefix, verb_out, log)
                if reads_paired['pair_out'] and not options.keep_temp_files:
                    os.remove(os.path.join(out_base, options.prefix + "filtered_1.fq"))
                    os.remove(os.path.join(out_base, options.prefix + "filtered_2.fq"))
            else:
                log.info("Separating filtered fastq file ... skipped.\n")

        """ assembly """
        if options.run_spades:
            if not (resume and os.path.exists(os.path.join(spades_output, 'assembly_graph.fastg'))):
                # resume = False
                log.info('Assembling using SPAdes ...')
                assembly_with_spades(options.spades_kmer, spades_output, other_options, out_base, options.prefix,
                                     original_fq_files, reads_paired, options.verbose_log, resume, options.threads, log)
            else:
                log.info('Assembling using SPAdes ... skipped.\n')

        """ export organelle """
        if os.path.exists(os.path.join(spades_output, 'assembly_graph.fastg')) \
                and options.organelle_type != '0':
            export_succeeded = False
            kmer_vals = sorted([int(kmer_val) for kmer_val in options.spades_kmer.split(",")], reverse=True)
            kmer_dirs = [os.path.join(spades_output, "K" + str(kmer_val)) for kmer_val in kmer_vals]
            for go_k, kmer_dir in enumerate(kmer_dirs):
                try:
                    run_stat = extract_organelle_genome(out_base=out_base, spades_output=kmer_dir, go_round=go_k,
                                                        prefix=options.prefix, organelle_type=options.organelle_type,
                                                        verbose=options.verbose_log, log=log, threads=options.threads,
                                                        options=options)
                except ImportError:
                    log.warning("Disentangling failed: numpy/scipy/sympy not installed!")
                    break
                except:
                    log.info("Disentangling failed.")
                else:
                    if run_stat == 0:
                        export_succeeded = True
                        break
            if not export_succeeded:
                out_fastg = sorted([os.path.join(spades_output, x)
                                    for x in os.listdir(spades_output) if x.count(".fastg") == 2])[0]
                out_csv = out_fastg[:-5] + "csv"
                log.info("Please ...")
                log.info("load the graph file: " + out_fastg)
                log.info("load the CSV file: " + out_csv)
                log.info("visualize and export your result in Bandage.\n")

        log = simple_log(log, out_base, prefix=options.prefix + "get_org.")
        log.info("\nTotal Calc-cost " + str(time.time() - time0))
        log.info("Thank you!")
    except:
        log.exception("")
        log = simple_log(log, out_base, prefix=options.prefix + "get_org.")
        log.info("\nTotal cost " + str(time.time() - time0))
        log.info("Please email jinjianjun@mail.kib.ac.cn if you find bugs!")
    logging.shutdown()


if __name__ == '__main__':
    main()

"""Copyright 2016 Jianjun Jin"""
