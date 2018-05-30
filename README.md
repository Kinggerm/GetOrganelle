# GetOrganelle


This pipeline assemblies organelle genome from genomic skimming data with an input reference, which should be the same organelle genome sequence but not necessarily a relate species within the focal clade.

<div id="citation"></div>

<b>Citation:</b> Jian-Jun Jin*, Wen-Bin Yu*, Jun-Bo Yang, Yu Song, Ting-Shuang Yi, De-Zhu Li. 2018. GetOrganelle: a simple and fast pipeline for de novo assembly of a complete circular chloroplast genome using genome skimming data. bioRxiv, 256479. [http://doi.org/10.1101/256479](https://www.biorxiv.org/content/early/2018/03/14/256479)

<b>License:</b> GPL https://www.gnu.org/licenses/gpl-3.0.html

<b>Bug&Usage contact:</b> [jinjianjun@mail.kib.ac.cn](mailto:jinjianjun@mail.kib.ac.cn); [yuwenbin@xtbg.ac.cn](mailto:yuwenbin@xtbg.ac.cn)

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
    
    echo "export PATH" >> ~/.bash_profile
    
    # for Linux
    
    echo "{where_you_clone_GetOrganelle}/GetOrganelle:$PATH" >> ~/.bashrc
    
    echo "{where_you_clone_GetOrganelle}/GetOrganelle/Utilities:$PATH" >> ~/.bashrc
    
    echo "export PATH" >> ~/.bashrc


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


## How To

<b>Preparing Data</b>
Cut raw data into certain size (ca. 2G per-end is enough for chloroplast for most normal angiosperm samples) if it is too large dataset. You could use the Linux or Mac OS build-in command to easily get a reduced file. Currently, this script was written for illumina pair-end data (fastq).

<b>Filtering and Assembly</b>
Take your input reference (fasta or bowtie index) as probe, the script would recruit target reads in successive rounds (iterations). You could also using the references in `Library/SeqReference`, but a more related reference is safer if the sequence quality is bad (say, degraded DNA samples). The value word size (followed with "-w"), like the kmer in assembly, is crucial to the feasibility and efficiency of this process. The best word size changes from data to data and will be affected by read length, read quality, base coverage, organ DNA percent and other factors. After recruitment, this script will automatically call SPAdes to assembly the target reads produced by the former step. The best kmer depends on a wide variety of factors too.

<b>Producing Result</b>
By default, the main script would call Utilities/slim_fastg.py to modify the `assembly_graph.fastg` file and produce a new fastg file (would be `assembly_graph.fastg.extend+cp-mt.fastg` if -F cp been called) along with a csv file (`assembly_graph.fastg.extend+cp-mt.csv`). View `assembly_graph.fastg.extend+cp-mt.fastg` and load the `assembly_graph.fastg.extend+cp-mt.csv` in Bandage, choose the best path as the final result. 
[Here](http://player.youku.com/embed/XMzUxODc3MDQyOA) (or [here](https://youtu.be/NqOIi-fBma4)) is a short video showing a standard way to extract the plastome from the assembly graph with Bandage. See [here](https://v.qq.com/x/page/g0602unrcsf.html) or [here](https://www.youtube.com/watch?v=cXUV7k-F26w) for more examples with more complicated (do not miss `3m01s - 5m53s`) situations.


## Example

For 2G raw data, 150 bp reads, to assembly chloroplast, typically I use:

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 10 -k 75,85,95,105 -P 300000

or in a fast but memory-consuming way:

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 5 -k 75,85,95,105 -P 1000000 -a mitochondria.fasta -J 3 -M 5

or in a slow and memory-economic way:

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 10 -k 75,85,95,105 -P 0 --out-per-round --remove-duplicates 0

For 2G raw data, 150 bp reads, to assembly mitochondria

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s mt_reference.fasta -w 93 -o mitochondria_output -R 30 -k 65,75,85,95 -P 1000000 -F mt 
    
For 2G raw data, 150 bp reads, to assembly nuclear ribosomal RNA (18S-ITS1-5.8S-ITS2-26S)

    get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s nr_reference.fasta -w 115 -o nr_output -R 7 -k 95,105,115 -P 0 -F nr


## Published Works Using GetOrganelle

It was previously cited as GetOrganelle (https://github.com/Kinggerm/GetOrganelle), but now we have a report paper (<a href="#citation">see above</a>) to cite.

Yu Song, Wen-Bin Yu, Yun-Bong Tan, Bing Liu, Xin Yao, Jian-Jun Jin, Michael Padmanaba, Jun-Bo Yang, Richard T. Corlett. 2017. Evolutionary comparisons of the chloroplast genome in Lauraceae and insights into loss events in the Magnoliids. Genome biology and evolution. 9(9): 2354-64. doi: [https://doi.org/10.1093/gbe/evx180](https://doi.org/10.1093/gbe/evx180)

Twyford AD, Ness RW. 2017. Strategies for complete plastid genome sequencing. Molecular Ecology Resources. 17(5):858-68. doi: [https://doi.org/10.1111/1755-0998.12626](https://doi.org/10.1111/1755-0998.12626)

Guan-Song Yang, Yin-Huan Wang, Yue-Hua Wang, Shi-Kang Shen. 2017. The complete chloroplast genome of a vulnerable species Champereia manillana (Opiliaceae). Conservation Genetics Resources. 9(3): 415-418. doi: [https://doi.org/10.1007/s12686-017-0697-1](https://doi.org/10.1007/s12686-017-0697-1)

[See here for more (10+)](http://www.wbyu.net/getorganelle.html)

## Acknowledgement

Thanks to Chao-Nan Fu, Han-Tao Qin, Xiao-Jian Qu, Shuo Wang, and Rong Zhang for giving tests or suggestions.