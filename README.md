# GetOrganelle


This pipeline assemblies organelle genome from genomic skimming data with an input reference, which should be the same organelle genome sequence but not necessarily a relate species within the focal clade.

At first, I was trying to follow <a href='http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12246/abstract'>Pierre-Jean G. Malé et al., 2014</a> to achieve chloroplast genome but found extractread in <a href='http://metabarcoding.org/obitools'>OBItools</a> hard to use. Inspired by the idea in that paper, I tried to write a script to express my understanding of getting organelle reads and practise my python programming skill. I didn't have any data that time and sent it to my friend Chaonan Fu to test it. While testing my script with her data, Chaonan contacted Dr Eric Coissac with the problem of extractread in IBOtools and found it no more maintenance. But nice Dr Eric Coissac suggested their another cool project <a href='http://metabarcoding.org/org-asm'>ORG.asm</a>, which has a sexy name. Although ORG.sam was found to be quite efficient and convenient with normal data (qualified total genomic pair-end reads) in my latter comparing tests, ORG.sam did not help my friend with her special samples (a kind of saprobe, I will put species name here when she publish her paper someday). However, to my surprise, my poorly written script did! The script filter_fq_with_probes.py found reads covering most of the plastid genome, which was subsequently checked with PCR, while many other pipelines found her no more than half of the plastid genome.

That's the beginning of the story. Latter, I revised the script for several times. As I tested my pipeline (bowtie2+script+SPAdes) with more data (120 samples including Fabaceae, Elaeagnaceae, Poaceae and Orobanchaceae species with fresh leave or specimen materials), I knew more about the associated priciples and the pros and cons of my pipeline. Now I post it here for fun even it's still memory exhausted with large data (Maybe I will optimize it someday, maybe not).

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

    python get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s reference.fasta -w 101 -o chloroplast

    python get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s reference.fasta -w 101 -o chloroplast -a mitochondrial.fasta

    python get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq --bs reference.fasta.index -w 103 -o chloroplast -P 500000 -R 20 -k 75,85,95,105
