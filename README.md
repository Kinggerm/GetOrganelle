# GetOrganelle


This pipeline assemblies organelle genome from genomic skimming data with an input reference, which should be the same organelle genome sequence but not necessarily a relate species within the focal clade.

Citation: [Jian-Jun Jin*, Wen-Bin Yu*, Jun-Bo Yang, Yu Song, Ting-Shuang Yi, De-Zhu Li (2018). GetOrganelle: a simple and fast pipeline for de novo assembly of a complete circular chloroplast genome using genome skimming data. bioRxiv, 256479. http://doi.org/10.1101/256479](https://www.biorxiv.org/content/early/2018/03/12/256479)

License: GPL https://www.gnu.org/licenses/gpl-3.0.html

Bug&Usage contact: [jinjianjun@mail.kib.ac.cn](jinjianjun@mail.kib.ac.cn); [yuwenbin@xtbg.ac.cn](yuwenbin@xtbg.ac.cn)

Please cite the dependencies if they are used:

SPAdes: [Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Nikolenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M. A. Alekseyev and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology 19: 455-477.](https://www.liebertpub.com/doi/abs/10.1089/cmb.2012.0021)

Bowtie2: [Langmead, B. and S. L. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. Nat Meth 9: 357-359.](https://www.nature.com/articles/nmeth.1923)

BLAST+: [Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, K. Bealer and T. L. Madden. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10: 421.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421)

Bandage: [Wick, R. R., M. B. Schultz, J. Zobel and K. E. Holt. 2015. Bandage: interactive visualization of de novo genome assemblies. Bioinformatics 31: 3350-3352.](https://academic.oup.com/bioinformatics/article/31/20/3350/196114)


## Installation

This pipeline was written in python 3.5.1, but compatible with 2.7.11.

Execute following simple git commands to download:

    git clone "https://github.com/Kinggerm/GetOrganelle"

then add */GetOrganelle and */GetOrganelle/Utilities to the path:

    # for MacOS
    
    echo "{where_you_clone_GetOrganelle}/GetOrganelle:$PATH" >> ~/.bash_profile
    
    echo "{where_you_clone_GetOrganelle}/GetOrganelle/Utilities:$PATH" >> ~/.bash_profile
    
    # for Linux
    
    echo "{where_you_clone_GetOrganelle}/GetOrganelle:$PATH" >> ~/.bashrc
    
    echo "{where_you_clone_GetOrganelle}/GetOrganelle/Utilities:$PATH" >> ~/.bashrc


and keep update when necessary:
    
    cd {where_you_clone_GetOrganelle}/GetOrganelle

    git pull

You could run the main script (get_organelle_reads.py) to get organelle reads (*.fastq) successfully, without any third-party libraries or software.

However, to get a complete organelle genome (such as a chloroplast genome) rather than organelle reads, other files in GetOrganelle are needed in the original relative path. Also, the following software are needed to be installed and configured in the path, since they could be called automatically:

<a href='http://bioinf.spbau.ru/spades'>SPAdes</a> is the assembler

<a href='http://bowtie-bio.sourceforge.net/bowtie2/index.shtml'>Bowtie2</a> is used to speed up initial recruitment of target-like reads

<a href='http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews'>BLAST+</a> is used to filter target-like contigs and simplify the final assembly graph

<a href='http://rrwick.github.io/Bandage/'>Bandage</a> is suggested to view the final contig graph (*.fastg).

Besides, if you installed python library psutil (pip install psutil), the memory cost of get_organelle_reads.py will be automatically logged.


## HowTo

1. Preparing Data: Cut raw data into certain size (ca. 2G per-end is enough for chloroplast for most normal angiosperm samples) if it is too large dataset. You could use the Linux or Mac OS build-in command to easily get a reduced file. Currently, this script was written for illumina pair-end data (fastq).

2. Filtering and Assembly: Take your input reference (fasta or bowtie index) as probe, the script would recruit target reads in successive rounds (iterations). The value word size (followed with "-w"), like the kmer in assembly, is crucial to the feasibility and efficiency of this process. The best word size changes from data to data and will be affected by read length, read quality, base coverage, organ DNA percent and other factors. After recruitment, this script will automatically call SPAdes to assembly the target reads produced by the former step. The best kmer depends on a wide variety of factors too.

3. Producing result: By default, the main script would call Utilities/slim_spades_fastg_by_blast.py to modify the assembly_graph.fastg file and produce a new fastg file along with a csv file. View the new fastg file and load the csv file in Bandage and choose the best path as the final result. 


## Example

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


## Similar pipelines 

I will post the differences latter:</p>
1.<a href='https://github.com/chrishah/MITObim'>MITObim</a> is where the first realization of this idea.</p>
2.<a href='http://metabarcoding.org/org-asm'>ORG.asm</a> is what I suggest be the best choice for normal samples.</p>
3.<a href='http://metabarcoding.org/obitools'>OBItools</a> is no more suggested to do the same thing.</p>
4.<a href='http://ibest.github.io/ARC'>ARC</a></p>
5.<a href='https://github.com/holmrenser/IOGA'>IOGA</a></p>
6.<a href='https://github.com/quxiaojian/PERR'>PERR</a></p>

## Acknowledgement

Chao-Nan Fu, Han-Tao Qin, Shuo Wang, Rong Zhang, Xiao-Jian Qu