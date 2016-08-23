#!/usr/bin/env python

import os
from optparse import OptionParser

usage = "mk_get_organelle.py -o basename -p \"*\" " \
        "[a_list_of_directories_containing_balanced_fastq_file, default='./*']"
parser = OptionParser(usage=usage)
parser.add_option('-o', dest='output_base',
                  help='output base name for each sample')
parser.add_option('-p', dest='other_arguments',
                  help='Other arguments that get_organelle_reads.py would take.'
                       'Use double quotation marks to include all the arguments'
                       'Example: "-s chloroplast.fasta -a mitochondrial.fasta -m cp -w 105"')
parser.add_option('--all', dest='skip_done', default=True, action='store_false',
                  help='Choose to make command for all samples including samples with results.'
                       'Default: skip those with results.')
parser.add_option('--strict', dest='strict_name', default=False, action='store_true',
                  help='Choose to only search for the fastq with the same base name with the directory '
                       '(*/*_1.fq). Default: relaxed searching.')
options, args = parser.parse_args()
if not (options.output_base and options.other_arguments):
    print("\nUsage: "+usage+'\n')
    exit()

if args:
    dirs = args
else:
    dirs = [x for x in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), x))]

dirs.sort()
lines = []
for dire in dirs:
    if options.strict_name:
        these_files = [os.path.join(dire, str(os.path.split(dire)[1])+'_'+str(y)+'.fq') for y in (1, 2)]
        for this_file in these_files:
            if not os.path.exists(this_file):
                print('Warning: '+this_file+' not found! Omitted!')
                continue
    else:
        these_files = [os.path.join(dire, y) for y in os.listdir(dire) if y.endswith('.fq') or y.endswith('.fastq')]
        if these_files:
            if len(these_files) != 2:
                print('Warning: Undetermined fq file in '+dire+'. Omitted!')
                continue
    if options.skip_done and os.path.isdir(os.path.join(dire, options.output_base)):
        print('Warning: '+os.path.join(dire, options.output_base)+' already exists. Annotated!')
        lines.append('# get_organelle_reads.py -1 '+these_files[0]+' -2 '+these_files[1]+' -o ' +
                     os.path.join(dire, options.output_base) + ' ' + options.other_arguments + '\n')
    else:
        lines.append('get_organelle_reads.py -1 ' + these_files[0] + ' -2 ' + these_files[1] + ' -o ' +
                     os.path.join(dire, options.output_base) + ' ' + options.other_arguments + '\n')

out_file = './get_organelle'
while os.path.exists(out_file+'.sh'):
    out_file += '_'
out_file += '.sh'
open(out_file, 'w').writelines(lines)
os.system('chmod 777 '+out_file)
