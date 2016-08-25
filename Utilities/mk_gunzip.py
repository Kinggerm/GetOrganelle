#!/usr/bin/env python

import os
from optparse import OptionParser

usage = "mkgunzip.py [a_list_of_directories_containing_balanced_fastq_file, default='./*']"
parser = OptionParser(usage=usage)
parser.add_option('-e', dest='execute', action='store_true', default=False,
                        help='Execute the gunzipping process directly in the script.'
                             'Default: False. Mainly create the batch file and wait for user\'s command.')
parser.add_option('-o', dest='batch_file', default=os.path.join(os.getcwd(), 'gunzip.bash.sh'),
                        help='Output file.')
parser.add_option('-n', dest='first_reads', type=int,
                        help='Only keep the first n (n >= 1) reads if the file is too large. '
                             'For get_organelle_reads.py and angiosperm cp genome, 7500000 is typically enough.')
parser.add_option('-f', dest='overwrite', default=False, action='store_true',
                        help='Overwrite even the fastq file exists.')
options, args = parser.parse_args()

if args:
    dirs = args
else:
    dirs = [x for x in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), x))]
lines = []
for dire in dirs:
    file_1 = []
    file_2 = []
    these_files = [os.path.join(dire, y) for y in os.listdir(dire) if y.endswith('.gz')]
    if these_files:
        these_files.sort()
        try:
            for i in range(0, len(these_files), 2):
                file_1.append(these_files[i])
                file_2.append(these_files[i+1])
        except IndexError:
            print('Error: Unbalanced file in '+dire+'. Omitted!')
            continue
        if options.overwrite and os.path.exists(os.path.join(dire, os.path.split(dire)[1]+'_1.fq')) and os.path.exists(os.path.join(dire, os.path.split(dire)[1]+'_2.fq')):
            print('Warning: '+os.path.join(dire, dire+'_*.fq')+' already exists. Omitted!')
        else:
            new_files = [os.path.join(dire, os.path.split(dire)[1]+'_'+str(k)+'.fq') for k in (1, 2)]
            
            if options.execute:
                lines.append('echo "gunzip -c '+' '.join(file_1)+' > '+new_files[0]+' ...\c"\n')
            else:
                lines.append('echo -n "gunzip -c '+' '.join(file_1)+' > '+new_files[0]+' ..."\n')
            lines.append('gunzip -c '+' '.join(file_1)+' > '+new_files[0]+'\n')
            if options.first_reads:
                lines.append('head -n '+str(4*options.first_reads)+' '+new_files[0]+' > '+new_files[0]+'.temp\n')
                lines.append('mv '+new_files[0]+'.temp '+new_files[0]+'\n')
            lines.append('echo " finished!"\n')
            
            if options.execute:
                lines.append('echo "gunzip -c '+' '.join(file_2)+' > '+new_files[1]+' ...\c"\n')
            else:
                lines.append('echo -n "gunzip -c '+' '.join(file_2)+' > '+new_files[1]+' ..."\n')
            lines.append('gunzip -c '+' '.join(file_2)+' > '+new_files[1]+'\n')
            if options.first_reads:
                lines.append('head -n '+str(4*options.first_reads)+' '+new_files[1]+' > '+new_files[1]+'.temp\n')
                lines.append('mv '+new_files[1]+'.temp '+new_files[1]+'\n')
            lines.append('echo " finished!"\n')

if options.execute:
    for line in lines:
        os.system(line)
else:
    open(options.batch_file, 'w').writelines(lines)
    os.system('chmod 777 '+options.batch_file)
