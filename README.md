# GetOrganelle


This pipeline assemblies organelle genome from genomic skimming data with an input reference, which should be the same organelle genome sequence but not necessarily a relate species within the focal clade.

Similar pipelines (I will post the differences latter.):</p>
1.<a href='https://github.com/chrishah/MITObim'>MITObim</a> is where the first realization of this idea.</p>
2.<a href='http://metabarcoding.org/org-asm'>ORG.asm</a> is what I suggest be the best choice for normal samples.</p>
3.<a href='http://metabarcoding.org/obitools'>OBItools</a> is no more suggested to do the same thing.</p>
4.<a href='http://ibest.github.io/ARC'>ARC</a></p>
5.<a href='https://github.com/holmrenser/IOGA'>IOGA</a></p>
6.<a href='https://github.com/quxiaojian/PERR'>PERR</a> is my senior of the same group, Xiaojian Qu's perl version work, which implements similar strategy but takes a different recruiting method and usually have a constant memory cost.</p>

Many thanks to Chaonan Fu, Dr Wenbin Yu, Hantao Qin and Shuo Wang!

==========================================================================
# Installation

My script was written in python 3.5.1, but compatible with 2.7.11.

Execute following simple git commands to download and keep update:

    git clone "https://github.com/Kinggerm/GetOrganelle"

and keep update:
    
    cd GetOrganelle

    git pull

It would be convenient to use my script if you add */GetOrganelle and */GetOrganelle/Utilities to the path.

You could run the main script (get_organelle_reads.py) to get organelle reads (*.fastq) successfully, without any third-party libraries or software.

However, to get a complete organ genome (such as a chloroplast genome) rather than organ reads, other files in GetOrganelle are needed in the original relative path. Also, the following software are suggested to be installed and configured in the path, since they could be called automatically by my script:

<a href='http://bioinf.spbau.ru/spades'>SPAdes</a>

<a href='http://bowtie-bio.sourceforge.net/bowtie2/index.shtml'>bowtie2</a>

<a href='http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews'>BLAST</a>

Besides, <a href='http://rrwick.github.io/Bandage/'>Bandage</a> is suggested to view the final contig graph (*.fastg). And if you install python library psutil, the memory cost of get_organelle_reads.py will be automatically logged.

==========================================================================
# HowTo

1. Preparing Data: Cut raw data into certain size (<2G per-end is enough for chloroplast for most normal angiosperm samples) if it is too large dataset. You could use the Linux or Mac OS build-in command to easily get a reduced file. Currently, this script was written for illumina pair-end data (fastq).

2. Filtering and Assembly: Take your input reference (fasta or bowtie index) as probe, the script would recruit target reads in successive rounds (iterations). The value word size (followed with "-w"), like the kmer in assembly, is crucial to the feasibility and efficiency of this process. The best word size changes from data to data and will be affected by read length, read quality, base coverage, organ DNA percent and other factors. After recruitment, this script will automatically call SPAdes to assembly the target reads produced by the former step. The best kmer depends on a wide variety of factors too.

3. Producing result: By default, the main script would call Utilities/slim_spades_fastg_by_blast.py to modify the assembly_graph.fastg file and produce a new fastg file along with a csv file. View the new fastg file and load the csv file in Bandage and choose the best path as the final result. 

==========================================================================
# Example

For 2G raw data, 150 bp reads, to assembly chloroplast, typically I use:

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 10 -k 75,85,95,105 -P 300000

or in a fast but memory-consuming way:

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 5 -k 75,85,95,105 -P 1000000 -a mitochondria.fasta -J 3 -M 5

or in a slow and memory-economic way:

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 10 -k 75,85,95,105 -P 0 --out-per-round --no-remove-duplicates

For 2G raw data, 150 bp reads, to assembly mitochondria

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s mt_reference.fasta -w 93 -o mitochondria_output -R 30 -k 65,75,85,95 -P 1000000 -F mt 
    
For 2G raw data, 150 bp reads, to assembly nuclear ribosomal RNA (18S-ITS1-5.8S-ITS2-26S)

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s nr_reference.fasta -w 115 -o nr_output -R 7 -k 95,105,115 -P 0 -F nr
