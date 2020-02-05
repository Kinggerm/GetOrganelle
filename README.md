# GetOrganelle


This toolkit assemblies organelle genome from genomic skimming data.

<div id="citation"></div>

Please denote the versions of GetOrganelle as well as the dependencies in your manuscript for reproducible science.

<b>Citation:</b> Jian-Jun Jin*, Wen-Bin Yu*, Jun-Bo Yang, Yu Song, Ting-Shuang Yi, De-Zhu Li. 2018. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. bioRxiv, 256479. [http://doi.org/10.1101/256479](https://www.biorxiv.org/content/early/2018/03/14/256479)

<b>License:</b> GPL https://www.gnu.org/licenses/gpl-3.0.html

<b>Bug&Usage contact:</b> [jinjianjun@mail.kib.ac.cn](mailto:jinjianjun@mail.kib.ac.cn) or [phylojin@163.com](mailto:phylojin@163.com); [yuwenbin@xtbg.ac.cn](mailto:yuwenbin@xtbg.ac.cn)

Please also cite the dependencies if used:

SPAdes: [Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Nikolenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M. A. Alekseyev and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology 19: 455-477.](https://www.liebertpub.com/doi/abs/10.1089/cmb.2012.0021)

Bowtie2: [Langmead, B. and S. L. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. Nature Methods 9: 357-359.](https://www.nature.com/articles/nmeth.1923)

BLAST+: [Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, K. Bealer and T. L. Madden. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10: 421.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421)

Bandage: [Wick, R. R., M. B. Schultz, J. Zobel and K. E. Holt. 2015. Bandage: interactive visualization of de novo genome assemblies. Bioinformatics 31: 3350-3352.](https://academic.oup.com/bioinformatics/article/31/20/3350/196114)


## Installation

This toolkit is currently maintained under Python 3.7.0, but designed to be compatible with versions higher than 3.5.1 and 2.7.11. GetOrganelle is generally more efficient under Python 3.*.

There are generally two ways to install GetOrganelle: 1) `Using the setup.py` is the way with GetOrganelleLib installed in the $PYTHONPATH ; 2) `In situ configuration` is the classic and heavy way, but easier to keep updated.

### Using the setup.py

Execute following curl commands to download suitable version (see more versions [here](https://github.com/Kinggerm/GetOrganelle/releases)). You can also use [git](https://www.atlassian.com/git/tutorials/install-git) to download as explained latter in the `In situ configuration`, but without the need of cloning into the installation directory.

    # To dowload GetOrganelle using curl and decompress it. 
    # Supposing your system is linux, otherwise change the 'linux' into 'macOS'; If you do not need the attached dependency, change 'linux' into 'light';
    # Supposing you download GetOrganelle to ~/Downloads
    cd ~/Downloads
    curl -L https://github.com/Kinggerm/GetOrganelle/releases/download/v1.6.2e/v1.6.2e-linux.tar.gz | tar zx
    
install [pip](https://pip.pypa.io/en/stable/installing/) and then install downloaded GetOrganelle with pip.

    # install pip, NOT neccessary if pip was already available
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    python get-pip.py
    
    # install GetOrganelle; add the option "--user" after "install" to install for current user
    pip install ./GetOrganelle
    
Alternatively, if you have Python library setuptools installed (`sudo apt install -y python-setuptools` or `sudo yum install -y python-setuptools`), you can install GetOrganelle with `cd GetOrganelle && python setup.py install`. 
    
For some fresh linux systems, after above commands you still cannot execute `get_organelle_from_reads.py` in a new terminal directly, meaning `~/.local/bin` was not added to the $PATH, you have to manually add `~/.local/bin` by executing `echo "PATH=~/.local/bin:\$PATH" >> ~/.bashrc`. For fresh MacOS environment with similar situation, for example, if you installed GetOrganelle with `Python 3.6` and find scripts not in the $PATH, please execute `echo "PATH=/Library/Frameworks/Python.framework/Versions/3.6/bin:\$PATH" >> ~/.bash_profile`.
    
### In situ configuration

You could use curl as explained above, however [git](https://www.atlassian.com/git/tutorials/install-git) would be more suggested for update and version control.

    # Supposing you are going to install it at ~/Applications/bin
    mkdir ~/Applications && mkdir ~/Applications/bin  # create directories if not existed
    GetOrganellePATH=~/Applications/bin
    cd $GetOrganellePATH
    git clone git://github.com/Kinggerm/GetOrganelle
    
use following commands to make GetOrganelle scripts executable; and make blast-databases and bowtie2 indices for default seeds.
    
    cd GetOrganelle
    python setup.py --in-situ

add GetOrganelle and GetOrganelle/Utilities to the $PATH.
    
    # for Linux
    echo "PATH=$GetOrganellePATH/GetOrganelle:\$PATH" >> ~/.bashrc
    echo "PATH=$GetOrganellePATH/GetOrganelle/Utilities:\$PATH" >> ~/.bashrc
    echo "export PATH" >> ~/.bashrc
    source ~/.bashrc
    
    ## for MacOS
    echo "PATH=$GetOrganellePATH/GetOrganelle:\$PATH" >> ~/.bash_profile
    echo "PATH=$GetOrganellePATH/GetOrganelle/Utilities:\$PATH" >> ~/.bash_profile
    echo "export PATH" >> ~/.bash_profile
    source ~/.bash_profile
    
At last, install python libraries numpy, scipy, and sympy using [pip](https://pip.pypa.io/en/stable/installing/). Alternatively you could install package/environment management systems such as [conda](https://conda.io/en/latest/), which already have those python packages installed. [Pyenv](https://github.com/pyenv/pyenv) is highly suggested to control python versions/environments.
    
    pip install numpy scipy sympy

### Updating GetOrganelle

You are always recommended to use the latest GetOrganelle, although you could find many old versions of GetOrganelle [here](https://github.com/Kinggerm/GetOrganelle/releases).

1. If you follow the way of `Using the setup.py`, to update, just do the same thing as the installation process.

2. If you follow the way of `In situ configuration` with git, go to the directory where you cloned GetOrganelle:
    
        # Supposing you are going to installed it at ~/Applications/bin
        cd ~/Applications/bin/GetOrganelle
        git stash
        git pull
        python setup.py --in-situ

### Requirement & Dependencies

Since v1.6*, GetOrganelle included binary files of all dependencies (SPAdes, Bowtie2, BLAST+) in its repository. Although making GetOrganelle use your own installed dependencies is not suggested for compatibility consideration, but you could still do this. For example, if SPAdes v3.6.2 is already available in the $PATH and you would like GetOrganelle to use the installed SPAdes v3.6.2, you could remove the SPAdes folder before executing `pip install ./GetOrganelle`. If all dependencies were previously installed (using `sudo apt install spades bowtie2 ncbi-blast+` for Ubuntu), you could download the [light version](https://github.com/Kinggerm/GetOrganelle/releases/download/v1.6.2e/v1.6.2e-light.tar.gz) upon installing GetOrganelle. 

Besides, no worries about interference from GetOrganelle's dependencies. Because during the installing process mentioned above, GetOrganelle would add those dependencies (SPAdes, Bowtie2, BLAST+) to the GetOrganelle-*.egg rather than to the $PATH, thereby not influence your own usage. For example, if you already installed SPAdes v3.6.2, after installing GetOrganelle, the spades version for your system would still be v3.6.2, while only GetOrganelle uses the [version in GetOrganelleDep](https://github.com/Kinggerm/GetOrganelle/blob/master/GetOrganelleDep/linux/SPAdes/share/spades/VERSION). 

Python libraries (numpy, scipy, sympy) is covered in the installation part. 

Perl is required for the wrapper of Bowtie2, but we assume that it was builtin in your system.

#### NOT Required Dependencies

<a href=' '>Bandage</a > is a fantastic tool to view the assembly graph (`*.fastg`/`*.gfa`). If you have Bandage correctly configured and add the binary folder of Bandage (which is `Bandage.app/Contents/MacOS` for MacOS) to the $PATH, get_organelle_from_*.py would automatically generate the a png formatted image of the assembly graph. 

If you installed python library psutil (version >= 3.0; pip install psutil), the memory cost of get_organelle_from_reads.py will be automatically logged. If you want to evaluate your results and plot the evaluation with `evaluate_assembly_using_mapping.py` and `round_statistics.py`, you have to further install python library matplotlib (pip install matplotlib).


## How To

<b>What you actually need to do is just typing in one simple command as suggested in <a href="#example">Example</a ></b>. But you are still recommended to read the following introductions:

<b>Preparing Data</b>

Currently, this script was written for illumina pair-end/single-end data (fastq or fastq.gz). 1G per end is enough for plastome for most normal angiosperm samples, and 5G per end is enough for mitochondria data. You could simply assign a maximum number of reads (number of seqs, not number of bases) for `get_organelle_from_reads.py` to use with flag `--max-reads` (Default value: 1.5E7 for "-F embplant_pt/embplant_nr/fungus_mt"; 7.5E7 for "-F embplant_mt/animal_mt/other_pt/anonym"; 3E8 for "-F animal_mt") or manually cut raw data into certain size before running GetOrganelle using the Linux or Mac OS build-in command (eg. `head -n 20000000 large.fq > small.fq`). 

<b>Filtering and Assembly</b>

Take your input seed (fasta; the default is `GetOrganelleLib/SeedDatabase/*.fasta`) as probe, the script would recruit target reads in successive rounds (extending process). You could also use a seed sequence of a related species, which would be safer if the sequence quality is bad (say, degraded DNA samples). The value word size (followed with "-w"), like the kmer in assembly, is crucial to the feasibility and efficiency of this process. The best word size changes upon data and will be affected by read length, read quality, base coverage, organ DNA percent and other factors. Since version 1.4.0, if there is no user assigned word size value, GetOrganelle would automatically estimate a proper word size based on the data characters. Although the automatically-estimated word size value does not ensure the best performance nor the best result, you do not need to adjust the value if a complete/circular organelle result is produced, because the circular result by GetOrganelle is generally consistent under different options. After extending, this script will automatically call SPAdes to assembly the target reads produced by the former step. The best kmer depends on a wide variety of factors too.

<b>Producing Result</b>

By default, SPAdes is automatically called to produce the assembly graph file `filtered_spades/assembly_graph.fastg`. Then, Utilities/slim_fastg.py is called to modify the `filtered_spades/assembly_graph.fastg` file and produce a new fastg file (would be `assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` if "-F embplant_pt" been used by get_organelle_from_reads.py) along with a tab-format annotation file (`assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv`). 

The `assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` file along with the `assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv` file would be further parsed by disentangle_organelle_assembly.py, and your target sequence file(s) `*complete*path_sequence.fasta` would be produced as the <b>final result</b>, if disentangle_organelle_assembly.py successfully solve the path. 

Otherwise, if disentangle_organelle_assembly.py failed to solve the path (produce `*contigs*path_sequence.fasta`), you could use the incomplete sequence to conduct downstream analysis or manually view `assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` and load the `assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv` in Bandage, choose the best path(s) as the <b>final result</b>. You could execute `slim_fastg.py -F embplant_pt -E embplant_mt assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg` to further remove mitogenome contigs for this easier visualization and manual completion.
[Here](http://player.youku.com/embed/XMzUxODc3MDQyOA) (or [here](https://youtu.be/NqOIi-fBma4)) is a short video showing a standard way to manually extract the plastome from the assembly graph with Bandage. See [here](https://v.qq.com/x/page/g0602unrcsf.html) or [here](https://www.youtube.com/watch?v=cXUV7k-F26w) for more examples with more complicated (do not miss `3m01s - 5m53s`) situations.


<b>GetOrganelle flowchart</b>

![flowchart](https://user-images.githubusercontent.com/8598031/65656060-086f5e80-e051-11e9-97a2-fb1d2a79375b.png)

## Example

To assembly Embryophyta plant plastome (e.g. using 2G raw data of 150 bp paired reads), typically I use:

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -o plastome_output -R 15 -k 21,45,65,85,105 -F embplant_pt

or in a draft way:

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -o plastome_output --fast -k 21,65,105 -w 0.68 -F embplant_pt

or in a slow and memory-economic way:

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -o plastome_output -R 30 -k 21,45,65,85,105  -F embplant_pt --memory-save

To assembly Embryophyta plant mitochondria (usually you need more than 5G raw data):

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -o mitochondria_output -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
    
To assembly Embryophyta plant nuclear ribosomal RNA (18S-ITS1-5.8S-ITS2-26S):

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -o nr_output -R 7 -k 35,85,115 -F embplant_nr

To assembly fungus mitochondria:

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -s fungus_mt_seed.fasta --genes fungus_mt_genes.fasta -R 10 -k 21,45,65,85,105 -F fungus_mt

To assembly animal mitochondria:

    get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -s animal_mt_seed.fasta --genes animal_mt_genes.fasta -R 10 -k 21,45,65,85,105 -F animal_mt

See a brief illustrations of those arguments by typing in:

    get_organelle_from_reads.py -h
    
or see the detailed illustrations:
    
    get_organelle_from_reads.py --help

Also see [GetOrganelleComparison](https://github.com/Kinggerm/GetOrganelleComparison) for a benchmark test of `GetOrganelle` and `NOVOPlasty` using 50 online samples.

## Acknowledgement

Thanks to Chao-Nan Fu, Han-Tao Qin, Yang Pan, Xiao-Jian Qu, Shuo Wang, Rong Zhang, Fei Zhao and lots of users for giving tests or suggestions.