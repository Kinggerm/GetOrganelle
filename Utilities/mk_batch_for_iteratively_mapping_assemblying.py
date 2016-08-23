#!/usr/bin/env python
import os
import sys
from optparse import OptionParser as Option

sys.stdout.write("\nThis is a script for making batch file for "
                 "iteratively mapping with bowtie2 and assemblying with SPAdes.\n\n")
sys.stdout.flush()

usage = "\n" + str(os.path.basename(__file__)) + \
        " -1 original_1.fq -2 original_2.fq -s seed.fasta -R 5 -k 75,85,95 -o output_base"
parser = Option(usage=usage)
parser.add_option('-1', dest='fastq_file_1', help='Input 1st fastq format file as pool')
parser.add_option('-2', dest='fastq_file_2', help='Input 2nd fastq format file as pool')
parser.add_option('-s', dest='seed_dir', help='Input fasta format file as initial seed')
parser.add_option('-R', dest='rounds', default=3, type=int,
                  help='How many iterations would you like to have? Default=3')
parser.add_option('-k', dest='spades_kmer', default='65,75,85',
                  help='SPAdes k-mer settings. Use the same format as in SPAdes. Default=65,75,85')
parser.add_option('-o', dest='output_sh_file',
                  help='Executable output batch file.')
parser.add_option('--un', dest='unpaired', default=False, action='store_true',
                  help='Try to map and assembly without paired information.')
options, args = parser.parse_args()
if not (options.seed_dir and options.fastq_file_1 and options.fastq_file_2 and options.output_sh_file):
    parser.print_help()
    sys.stdout.write('\nERROR: Insufficient arguments!\n')
    exit()

out_f_h = open(options.output_sh_file+'.sh', 'w')
if not os.path.exists(options.output_sh_file + '.mapped.RUN1'):
    os.mkdir(options.output_sh_file + '.mapped.RUN1')
else:
    sys.stdout.write('\nError: file exists!\n')
    exit()
if options.unpaired:
    out_f_h.write('bowtie2-build --large-index ' + options.seed_dir + ' ' + options.seed_dir + '.index\n')
    out_f_h.write('bowtie2 -I 0 -X 800 -p 5 --very-fast-local --al ' + options.output_sh_file + '.mapped.RUN1/mapped.fq -q -x ' + options.seed_dir + '.index -U ' + options.fastq_file_1 + ',' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam\n')
    out_f_h.write('spades.py --careful -k ' + options.spades_kmer + ' -s ' + options.output_sh_file + '.mapped.RUN1/mapped.fq -o ' + options.output_sh_file + '.mapped.RUN1/spades\n')
    for i in range(2, options.rounds + 1):
        os.mkdir(options.output_sh_file + '.mapped.RUN' + str(i))
        out_f_h.write('bowtie2-build --large-index ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg.index\n')
        out_f_h.write('bowtie2 -I 0 -X 800 -p 5 --very-fast-local --al ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/mapped.fq -q -x ' + options.seed_dir + '.index -U ' + options.fastq_file_1 + ',' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam\n')
        out_f_h.write('spades.py --careful -k ' + options.spades_kmer + ' -s ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/mapped.fq -o ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/spades\n')
else:
    out_f_h.write('bowtie2-build --large-index ' + options.seed_dir + ' ' + options.seed_dir + '.index\n')
    out_f_h.write('bowtie2 -I 0 -X 800 -p 5 --very-fast-local --al-conc ' + options.output_sh_file + '.mapped.RUN1/%.mapped.fq -q -x ' + options.seed_dir + '.index -1 ' + options.fastq_file_1 + ' -2 ' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN1/bowtie2.bam\n')
    out_f_h.write('spades.py --careful -k ' + options.spades_kmer + ' -1 ' + options.output_sh_file + '.mapped.RUN1/1.mapped.fq -2 ' + options.output_sh_file + '.mapped.RUN1/2.mapped.fq -o ' + options.output_sh_file + '.mapped.RUN1/spades\n')
    for i in range(2, options.rounds+1):
        os.mkdir(options.output_sh_file + '.mapped.RUN' + str(i))
        out_f_h.write('bowtie2-build --large-index ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg ' + options.output_sh_file + '.mapped.RUN' + str(i - 1) + '/spades/assembly_graph.fastg.index\n')
        out_f_h.write('bowtie2 -I 0 -X 800 -p 5 --very-fast-local --al-conc ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/%.mapped.fq -q -x ' + options.seed_dir + '.index -1 ' + options.fastq_file_1 + ' -2 ' + options.fastq_file_2 + ' -S ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam --no-unal -t && rm ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/bowtie2.bam\n')
        out_f_h.write('spades.py --careful -k ' + options.spades_kmer + ' -1 ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/1.mapped.fq -2 ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/2.mapped.fq -o ' + options.output_sh_file + '.mapped.RUN' + str(i) + '/spades\n')
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