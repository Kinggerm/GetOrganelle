#!/usr/bin/env python
import os
import sys
import platform
from argparse import ArgumentParser as Option
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
import GetOrganelleLib
from GetOrganelleLib.pipe_control_func import executable
from GetOrganelleLib.versions import get_versions
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
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

sys.stdout.write("\nThis is a script for making batch file for "
                 "iteratively mapping with bowtie2 and assemblying with SPAdes.\n\n")
sys.stdout.flush()

usage = "\n" + str(os.path.basename(__file__)) + \
        " -1 original_1.fq -2 original_2.fq -s seed.fasta -R 5 -k 21,45,65,85,105 -o output_base"
parser = Option(usage=usage)
parser.add_argument('-1', dest='fastq_file_1', help='Input 1st fastq format file as pool')
parser.add_argument('-2', dest='fastq_file_2', help='Input 2nd fastq format file as pool')
parser.add_argument('-s', dest='seed_dir', help='Input fasta format file as initial seed')
parser.add_argument('-R', dest='rounds', default=3, type=int,
                  help='How many iterations would you like to have? Default=3')
parser.add_argument('-t', dest="threads", default=1, type=int,
                  help="theads used for bowtie2 and SPAdes. Default=1")
parser.add_argument('-k', dest='spades_kmer', default='21,45,65,85,105',
                  help='SPAdes k-mer settings. Use the same format as in SPAdes. Default=21,45,65,85,105')
parser.add_argument('-o', dest='output_sh_file',
                  help='Executable output batch file.')
parser.add_argument('--un', dest='unpaired', default=False, action='store_true',
                  help='Try to map and assembly without paired information.')
parser.add_argument('--random-seed', dest="random_seed", type=int, default=12345,
                  help="seed for random generator for bowtie2. Default: %(default)s")
parser.add_argument("--which-bowtie2", dest="which_bowtie2", default="",
                  help="Assign the path to Bowtie2 binary files if not added to the path. "
                       "Default: try GetOrganelleDep/" + SYSTEM_NAME + "/bowtie2 first, then $PATH")
parser.add_argument("--which-spades", dest="which_spades", default="",
                  help="Assign the path to SPAdes binary files if not added to the path. "
                       "Default: try GetOrganelleDep/" + SYSTEM_NAME + "/SPAdes first, then $PATH")
parser.add_argument("-v", "--version", action="version",
                        version="GetOrganelle v{version}".format(version=get_versions()))
options = parser.parse_args()
if not (options.seed_dir and options.fastq_file_1 and options.fastq_file_2 and options.output_sh_file):
    parser.print_help()
    sys.stdout.write('\nERROR: Insufficient arguments!\n')
    exit()
if options.fastq_file_1 == options.fastq_file_2:
    raise IOError("1st fastq file should NOT be the same with 2nd fastq file!")
if not options.which_bowtie2:
    try_this_bin = os.path.join(GO_DEP_PATH, "bowtie2", "bowtie2")
    if os.path.isfile(try_this_bin) and executable(try_this_bin):
        options.which_bowtie2 = os.path.split(try_this_bin)[0]
if not options.which_spades:
    try_this_bin = os.path.join(GO_DEP_PATH, "SPAdes", "bin", "spades.py")
    if os.path.isfile(try_this_bin) and executable(try_this_bin):
        options.which_spades = os.path.split(try_this_bin)[0]
if not executable(os.path.join(options.which_bowtie2, "bowtie2")):
    sys.stdout.write("Warning: " + os.path.join(options.which_bowtie2, "bowtie2") + " not accessible!")
if not executable(os.path.join(options.which_spades, "spades.py")):
    sys.stdout.write("Warning: " + os.path.join(options.which_spades, "spades.py") + " not accessible!")

out_f_h = open(options.output_sh_file+'.sh', 'w')
if not os.path.exists(options.output_sh_file + '.mapped.RUN1'):
    os.mkdir(options.output_sh_file + '.mapped.RUN1')
else:
    sys.stdout.write('\nError: file exists!\n')
    exit()
if options.unpaired:
    out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2-build') + ' --seed ' + str(options.random_seed) + ' --large-index ' + options.seed_dir + ' ' + options.seed_dir + '.index\n')
    out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2') + ' --seed ' + str(options.random_seed) + ' -I 0 -X 1000 -p ' + str(options.threads) + ' --very-fast-local --al ' + options.output_sh_file + '.mapped.RUN1/mapped.fq -q -x ' + options.seed_dir + '.index -U ' + options.fastq_file_1 + ',' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam\n')
    out_f_h.write(os.path.join(options.which_spades, 'spades.py') + ' -t ' + str(options.threads) + ' -k ' + options.spades_kmer + ' -s ' + options.output_sh_file + '.mapped.RUN1/mapped.fq -o ' + options.output_sh_file + '.mapped.RUN1/spades\n')
    for i in range(2, options.rounds + 1):
        os.mkdir(options.output_sh_file + '.mapped.RUN' + str(i))
        out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2-build') + ' --seed ' + str(options.random_seed) + ' --large-index ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg.index\n')
        out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2') + ' --seed ' + str(options.random_seed) + ' -I 0 -X 1000 -p ' + str(options.threads) + ' --very-fast-local --al ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/mapped.fq -q -x ' + options.seed_dir + '.index -U ' + options.fastq_file_1 + ',' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam\n')
        out_f_h.write(os.path.join(options.which_spades, 'spades.py') + ' -t ' + str(options.threads) + ' -k ' + options.spades_kmer + ' -s ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/mapped.fq -o ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/spades\n')
else:
    out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2-build') + ' --seed ' + str(options.random_seed) + ' --large-index ' + options.seed_dir + ' ' + options.seed_dir + '.index\n')
    out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2') + ' --seed ' + str(options.random_seed) + ' -I 0 -X 1000 -p ' + str(options.threads) + ' --very-fast-local --al-conc ' + options.output_sh_file + '.mapped.RUN1/%.mapped.fq -q -x ' + options.seed_dir + '.index -1 ' + options.fastq_file_1 + ' -2 ' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam\n')
    out_f_h.write(os.path.join(options.which_spades, 'spades.py') + ' -t ' + str(options.threads) + ' -k ' + options.spades_kmer + ' -1 ' + options.output_sh_file + '.mapped.RUN1/1.mapped.fq -2 ' + options.output_sh_file + '.mapped.RUN1/2.mapped.fq -o ' + options.output_sh_file + '.mapped.RUN1/spades\n')
    for i in range(2, options.rounds+1):
        os.mkdir(options.output_sh_file + '.mapped.RUN' + str(i))
        out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2-build') + ' --seed ' + str(options.random_seed) + ' --large-index ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg.index\n')
        out_f_h.write(os.path.join(options.which_bowtie2, 'bowtie2') + ' --seed ' + str(options.random_seed) + ' -I 0 -X 1000 -p ' + str(options.threads) + ' --very-fast-local --al-conc ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/%.mapped.fq -q -x ' + options.seed_dir + '.index -1 ' + options.fastq_file_1 + ' -2 ' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam\n')
        out_f_h.write(os.path.join(options.which_spades, 'spades.py') + ' -t ' + str(options.threads) + ' -k ' + options.spades_kmer + ' -1 ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/1.mapped.fq -2 ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/2.mapped.fq -o ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/spades\n')
out_f_h.close()
os.system('chmod 777 ' + options.output_sh_file+'.sh')


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