# GetOrganelle

[![Anaconda-Server Badge](https://anaconda.org/bioconda/getorganelle/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/getorganelle/badges/version.svg)](https://anaconda.org/bioconda/getorganelle)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/getorganelle/badges/latest_release_date.svg)](https://anaconda.org/bioconda/getorganelle)

[![GitHub release](https://img.shields.io/github/release/Kinggerm/GetOrganelle.svg)](https://GitHub.com/Kinggerm/GetOrganelle/releases/)
[![GitHub version](https://img.shields.io/github/commits-since/Kinggerm/GetOrganelle/1.7.3.3.svg)](https://github.com/Kinggerm/GetOrganelle/commit/master)

This toolkit assemblies organelle genome from genomic skimming data. 

It achieved the best performance overall both on simulated and real data and was recommended as the default for chloroplast genome assemblies in a third-party comparison paper ([Freudenthal et al. 2020. Genome Biology](https://doi.org/10.1186/s13059-020-02153-6)).

<div id="citation"></div>

Please denote the version of GetOrganelle as well as the dependencies in your manuscript for reproducible science.

<b>Citation:</b> Jian-Jun Jin*, Wen-Bin Yu*, Jun-Bo Yang, Yu Song, Claude W. dePamphilis, Ting-Shuang Yi, De-Zhu Li. <b>GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes.</b> <i>Genome Biology</i> <b>21</b>, 241 (2020). [https://doi.org/10.1186/s13059-020-02154-5](https://doi.org/10.1186/s13059-020-02154-5)

<b>License:</b> GPL https://www.gnu.org/licenses/gpl-3.0.html

Please also cite the dependencies if used:

SPAdes: [Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Nikolenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M. A. Alekseyev and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology 19: 455-477.](https://www.liebertpub.com/doi/abs/10.1089/cmb.2012.0021)

Bowtie2: [Langmead, B. and S. L. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. Nature Methods 9: 357-359.](https://www.nature.com/articles/nmeth.1923)

BLAST+: [Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, K. Bealer and T. L. Madden. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10: 421.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421)

Bandage: [Wick, R. R., M. B. Schultz, J. Zobel and K. E. Holt. 2015. Bandage: interactive visualization of de novo genome assemblies. Bioinformatics 31: 3350-3352.](https://academic.oup.com/bioinformatics/article/31/20/3350/196114)


## Installation & Initialization

GetOrganelle is currently maintained under Python 3.7.0, but designed to be compatible with versions higher than 3.5.1 and 2.7.11. It was built for Linux and macOS. Windows Subsystem Linux is currently not supported, we are working on this.

* The easiest way to install GetOrganelle and its [dependencies](https://github.com/Kinggerm/GetOrganelle/wiki/Installation#requirement--dependencies) is using conda:
       
       
      conda install -c bioconda getorganelle

  You have to install [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) before using the above command. If you don't like conda, or want to follow the latest versions, you can find [more installation options here](https://github.com/Kinggerm/GetOrganelle/wiki/Installation#installation).

* After installation of GetOrganelle v1.7.0+, please download and initialize the database of your preferred organelle type (embplant_pt, embplant_mt, embplant_nr, fungus_mt, fungus_nr, animal_mt, and/or other_pt). Supposing you are assembling chloroplast genomes:

    
      get_organelle_config.py --add embplant_pt,embplant_mt
    
  If connection keeps failing, please manually download the latest database from [GetOrganelleDB](https://github.com/Kinggerm/GetOrganelleDB) and [initialization from local files](https://github.com/Kinggerm/GetOrganelle/wiki/Initialization#option-2-initialization-from-local-files).
  
  The database will be located at `~/.GetOrganelle` by default, which can be changed via the command line parameter `--config-dir`, or via the shell environment variable `GETORG_PATH` (see more [here](https://github.com/Kinggerm/GetOrganelle/wiki/Initialization)).
    

## Test

Download [a simulated _Arabidopsis thaliana_ WGS dataset](https://github.com/Kinggerm/GetOrganelleGallery/tree/master/Test/reads):

    wget https://github.com/Kinggerm/GetOrganelleGallery/raw/master/Test/reads/Arabidopsis_simulated.1.fq.gz
    wget https://github.com/Kinggerm/GetOrganelleGallery/raw/master/Test/reads/Arabidopsis_simulated.2.fq.gz

then verify the integrity of downloaded files using `md5sum`:

    md5sum Arabidopsis_simulated.*.fq.gz
    # 935589bc609397f1bfc9c40f571f0f19  Arabidopsis_simulated.1.fq.gz
    # d0f62eed78d2d2c6bed5f5aeaf4a2c11  Arabidopsis_simulated.2.fq.gz
    # Please re-download the reads if your md5 values unmatched above

then do the fast plastome assembly (memory: ~600MB, CPU time: ~60s):

    get_organelle_from_reads.py -1 Arabidopsis_simulated.1.fq.gz -2 Arabidopsis_simulated.2.fq.gz -t 1 -o Arabidopsis_simulated.plastome -F embplant_pt -R 10

You are going to get a similar running log as [here](https://github.com/Kinggerm/GetOrganelle/wiki/Example-1#running-log) and the same result as [here](https://github.com/Kinggerm/GetOrganelleGallery/tree/master/Test/results/Arabidopsis_simulated.plastome).

Find more real data examples at [GetOrganelle wiki](https://github.com/Kinggerm/GetOrganelle/wiki/Examples), [GetOrganelleGallery](https://github.com/Kinggerm/GetOrganelleGallery) and [GetOrganelleComparison](https://github.com/Kinggerm/GetOrganelleComparison).


## Instruction

<b>What you actually need to do is just typing in one simple command as suggested in <a href="#recipes">Recipes</a ></b>. But you are still recommended to read the following introductions:

<b>Preparing Data</b>

Currently, `get_organelle_from_reads.py` was written for illumina pair-end/single-end data (fastq or fastq.gz). Usually, >1G per end is enough for plastome for most normal angiosperm samples, and >5G per end is enough for mitochondria genome assembly. Since v1.6.2, `get_organelle_from_reads.py` will automatically estimate the read data it needs, without user assignment nor data reducing (see flags `--reduce-reads-for-coverage` and `--max-reads` for more options). 

<b>Filtering and Assembly</b>

Take your input seed (fasta; the default is `GetOrganelleLib/SeedDatabase/*.fasta`) as probe, the script would recruit target reads in successive rounds (extending process). The default seed works for most samples, but using a complete organelle genome sequence of a related species as the seed would help rescue data of bad sequence quality (e.g. degraded DNA samples). The value word size (followed with "-w"), like the kmer in assembly, is crucial to the feasibility and efficiency of this process. The best word size changes upon data and will be affected by read length, read quality, base coverage, organ DNA percent and other factors. Since version 1.4.0, if there is no user assigned word size value, GetOrganelle would automatically estimate a proper word size based on the data characters. Although the automatically-estimated word size value does not ensure the best performance nor the best result, you do not need to adjust the value if a complete/circular organelle result is produced, because the circular result by GetOrganelle is highly consistent under different options and seeds. After extending, this script will automatically call SPAdes to assembly the target reads produced by the former step. The best kmer depends on a wide variety of factors too.

<b>Producing Result</b>

By default, `SPAdes` is automatically called to produce the assembly graph file `filtered_spades/assembly_graph.fastg`. Then, `Utilities/slim_graph.py` is called to modify the `filtered_spades/assembly_graph.fastg` file and produce a new fastg file (would be `assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` if "-F embplant_pt" been used by `get_organelle_from_reads.py`) along with a tab-format annotation file (`assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv`). 

The `assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` file along with the `assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv` file would be further parsed by `disentangle_organelle_assembly.py`, and your target sequence file(s) `*complete*path_sequence.fasta` would be produced as the <b>final result</b>, if disentangle_organelle_assembly.py successfully solve the path. 

Otherwise, if GetOrganelle failed to solve the path (produce `*scaffolds*path_sequence.fasta`), you could use the incomplete sequence to conduct downstream analysis or manually view `assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` and load the `assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv` in [Bandage](http://rrwick.github.io/Bandage/), choose the best path(s) as the <b>final result</b>. You could execute `slim_graph.py -F embplant_pt -E embplant_mt assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` to further remove mitogenome contigs for this easier visualization and manual completion.
[Here](http://player.youku.com/embed/XMzUxODc3MDQyOA) (or [here](https://youtu.be/NqOIi-fBma4)) is a short video showing a standard way to manually extract the plastome from the assembly graph with Bandage. See [here](https://v.qq.com/x/page/g0602unrcsf.html) or [here](https://www.youtube.com/watch?v=cXUV7k-F26w) for more examples.


<b>GetOrganelle flowchart</b>

![flowchart](https://user-images.githubusercontent.com/8598031/83836465-85afa080-a6c1-11ea-8b08-b08623d974f4.png)

## Recipes

To assembly Embryophyta plant plastome (e.g. using 2G raw data of 150 bp paired reads), typically I use:

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o plastome_output -R 15 -k 21,45,65,85,105 -F embplant_pt

or in a draft way:

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o plastome_output --fast -k 21,65,105 -w 0.68 -F embplant_pt

or in a slow and memory-economic way:

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o plastome_output -R 30 -k 21,45,65,85,105  -F embplant_pt --memory-save

To assembly Embryophyta plant mitochondria (usually you need more than 5G raw data):

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o mitochondria_output -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
    
To assembly Embryophyta plant nuclear ribosomal RNA (18S-ITS1-5.8S-ITS2-26S):

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o nr_output -R 10 -k 35,85,115 -F embplant_nr

To assembly fungus mitochondria:

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -R 10 -k 21,45,65,85,105 -F fungus_mt -o fungus_mt_out  # if you fails with the default database, use your own seed database and label database with "-s" and "--genes" 

To assembly fungus nuclear ribosomal RNA (18S-ITS1-5.8S-ITS2-28S):

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -R 10 -k 21,45,65,85,105 -F fungus_nr -o fungus_nr_out  # if you fails with the default database, use your own seed database and label database with "-s" and "--genes" 

To assembly animal mitochondria:

    get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -R 10 -k 21,45,65,85,105 -F animal_mt -o animal_mt_out   # if you fails with the default database, use your own seed database and label database with "-s" and "--genes"

See a brief illustrations of those arguments by typing in:

    get_organelle_from_reads.py -h
    
or see the detailed illustrations:
    
    get_organelle_from_reads.py --help
    
To extract the plastid genome from an existing assembly graph (`*.fastg`/`*.gfa`; e.g. from long-read sequencing assemblies):

    get_organelle_from_assembly.py -F embplant_pt -g ONT_assembly_graph.gfa


## Contact

If your question is running specific, please attach the `get_org.log.txt` file and the post-slimming assembly graph (`assembly_graph.fastg.extend_*.fastg`, could be Bandage-visualized *.png format to protect your data privacy).

* Report bugs & Open issues [here](https://github.com/Kinggerm/GetOrganelle/issues).

* Send email to us ([jianjun.jin@columbia.edu](mailto:jianjun.jin@columbia.edu)&[jinjianjun@mail.kib.ac.cn](mailto:jinjianjun@mail.kib.ac.cn), [yuwenbin@xtbg.ac.cn](mailto:yuwenbin@xtbg.ac.cn))

* Join the QQ group (ID: 859158590)
