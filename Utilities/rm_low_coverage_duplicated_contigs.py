#!/usr/bin/env python
import sys
from optparse import OptionParser
import subprocess
try:
    import commands
except:
    pass
import os
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(path_of_this_script, ".."))
import GetOrganelleLib
from GetOrganelleLib.seq_parser import *
from GetOrganelleLib.pipe_control_func import executable, make_blast_db, execute_blast
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
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


def check_db(reference_fa_base, which_blast=""):
    in_index = reference_fa_base + '.index'
    if min([os.path.exists(in_index+postfix) for postfix in ('.nhr', '.nin', '.nsq')]):
        pass
    elif reference_fa_base:
        print('Making BLAST db ... ')
        make_blast_db(input_file=reference_fa_base, output_base=in_index, which_blast=which_blast)
        print('Making BLAST db finished.')
    else:
        print('ERROR: No reference input!')
        exit()
    return in_index


def require_options():
    usage = "Usage: rm_low_coverage_duplicated_contigs.py *.fastg"
    parser = OptionParser(usage=usage)
    parser.add_option('--cov-t', dest='coverage_threshold', default=0.12,
                      help='With ratio (coverage of query/coverage of subject) below which, '
                           'the query would be exposed to discarded. Default: 0.12')
    parser.add_option('--len-t', dest='length_threshold', default=0.9,
                      help='With overlap (length of hit of query/ length of query) above which, '
                           'the query would be exposed to discarded. Default: 0.9')
    parser.add_option('--blur', dest='blur_bases', default=False, action='store_true',
                      help='Replace hit low-coverage bases with N.')
    parser.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true',
                      help='Keep temp blast files.')
    parser.add_option("--which-blast", dest="which_blast", default="",
                      help="Assign the path to BLAST binary files if not added to the path.")
    parser.add_option('-o', dest='output_dir',
                      help='Output directory. Default: along with the original file')
    parser.add_option('-t', '--threads', dest="threads", default=4, type=int,
                      help="Threads of blastn.")
    options, args = parser.parse_args()
    if not args:
        parser.print_help()
        sys.stdout.write('\n######################################\nERROR: Insufficient REQUIRED arguments!\n\n')
        exit()
    if not options.which_blast:
        try_this_bin = os.path.join(GO_DEP_PATH, "ncbi-blast", "blastn")
        if os.path.isfile(try_this_bin) and executable(try_this_bin):
            output, err = subprocess.Popen(
                try_this_bin + " -version", stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, shell=True).communicate()
            if "not found" in output.decode("utf8"):
                sys.stdout.write(output.decode("utf8") + "\n")
            else:
                options.which_blast = os.path.split(try_this_bin)[0]
    if not executable(os.path.join(options.which_blast, "blastn")):
        sys.stdout.write(os.path.join(options.which_blast, "blastn") + " not accessible!")
        exit()
    if not executable(os.path.join(options.which_blast, "makeblastdb")):
        sys.stdout.write(os.path.join(options.which_blast, "makeblastdb") + " not accessible!")
        exit()
    if options.treat_no_hits not in ["ex_no_con", "ex_no_hit", "keep_all"]:
        sys.stdout.write('\n\nOption Error: you should choose assign one of "ex_no_con", "ex_no_hit"'
                         ' and "keep_all" to variable treat_no_hits\n')
        exit()
    return options, args


def purify_fastg(fastg_file, cov_threshold, len_threshold, blur, keep_temp, output, threads, which_blast):
    index_base = check_db(fastg_file, which_blast=which_blast)
    execute_blast(fastg_file, index_base, fastg_file + '.blast', threads=threads, outfmt=6,
                  e_value="1E-30", word_size=10, which_blast=which_blast)
    blast_result = open(fastg_file + '.blast', 'r')
    suspicious_nodes = {}
    for line in blast_result:
        if not line.startswith("#") and line.strip():
            records = line.strip().split('\t')
            query_id = records[0]
            subject_id = records[1]
            if query_id != subject_id:
                query_cov = float(query_id.split(':')[0].split(';')[0].rstrip('\'').split('cov_')[1])
                subject_cov = float(subject_id.split(':')[0].split(';')[0].rstrip('\'').split('cov_')[1])
                node_name = query_id.split(':')[0]
                if subject_cov and query_cov/subject_cov < cov_threshold:
                    if node_name in suspicious_nodes:
                        for base in range(int(records[6]), int(records[7]) + 1):
                            suspicious_nodes[node_name].add(base)
                    else:
                        suspicious_nodes[node_name] = set(range(int(records[6]), int(records[7]) + 1))
    fastg_matrix = read_fasta_as_list(fastg_file)
    i = 0
    while i < len(fastg_matrix[0]):
        node_name = fastg_matrix[0][i].split(':')[0].strip(';').strip('-')
        if node_name in suspicious_nodes:
            if len(suspicious_nodes[node_name])/float(len(fastg_matrix[1][i])) > len_threshold:
                del fastg_matrix[0][i]
                del fastg_matrix[1][i]
                continue
            elif blur:
                for base_to_blur in suspicious_nodes[node_name]:
                    fastg_matrix[1][i][base_to_blur-1] = 'N'
        i += 1
    if output:
        output_f = os.path.join(output, os.path.basename(fastg_file+'.purified.fastg'))
    else:
        output_f = fastg_file + '.purified.fastg'
    write_fasta_with_list(output_f, fastg_matrix, True)
    blast_result.close()
    if not keep_temp:
        os.remove(fastg_file + '.blast')
        for postfix in ('.nhr', '.nin', '.nsq'):
            os.remove(fastg_file + '.index' + postfix)


def main():
    options, args = require_options()
    for fastg_file in args:
        purify_fastg(fastg_file, options.coverage_threshold, options.length_threshold,
                     options.blur_bases, options.keep_temp, options.output, options.threads, options.which_blast)


if __name__ == '__main__':
    main()


"""Copyright 2016 Jianjun Jin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License."""