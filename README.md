# GetOrganelle


This script, written in python2, get organelle reads by extending with input seeds.

At first, I was trying to follow <a href='http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12246/abstract'>Pierre-Jean G. Malé et al., 2014</a> to achieve chloroplast genome but found extractread in <a href='http://metabarcoding.org/obitools'>OBItools</a> hard to use. Inspired by the idea in that paper, I tried to write a script to express my understanding of getting organelle reads and practise my python programming skill in Dec 3, 2015. I didn't have any data that time and sent it to my friend Chaonan Fu to test it. While testing my script with her data, Chaonan contacted Dr Eric Coissac with the problem of extractread in IBOtools and found it no more maintenance. But nice Dr Eric Coissac suggested their another cool project <a href='http://metabarcoding.org/org-asm'>ORG.asm</a>, which has a sexy name. Although ORG.sam was found to be supernaturally efficient and convenient with normal data (qualified total genomic pair-end reads) in my latter comparing tests, ORG.sam did not help my friend with her special samples (a kind of saprobe, I will put species name here when she publish her paper someday). However, to my surprise, my poorly written script did! The script filter_fq_with_probes.py found reads covering most of the plastid genome, which was subsequently checked with PCR, while many other pipelines found her no more than half of the plastid genome.

That's the beginning of the story. Latter, I revised the script for several times. As I tested my pipeline (bowtie2+script+SPAdes) with more data (120 samples including Fabaceae, Elaeagnaceae, Poaceae and Orobanchaceae species with fresh leave or specimen materials), I knew more about the associated priciples and the pros and cons of my pipeline. Now I post it here for fun even it's still memory exhausted with large data (Maybe I will optimize it someday, maybe not).

Similar pipelines (I will post the differences latter.):</p>
1.<a href='https://github.com/chrishah/MITObim'>MITObim</a> is where the first realization of this idea.</p>
2.<a href='http://metabarcoding.org/org-asm'>ORG.asm</a> is what I suggest be the best choice for normal samples.</p>
3.<a href='http://metabarcoding.org/obitools'>OBItools</a></p>
4.<a href='http://ibest.github.io/ARC'>ARC</a></p>
5.<a href='https://github.com/holmrenser/IOGA'>IOGA</a></p>

Many thanks to Chaonan Fu, Dr Wenbin Yu and Hantao Qin!

==========================================================================

<b>Brief workflow of assemblying plastid genome from total reads</b>

1. Preparing Data: Cut raw data into certain size (<2G per-end).
   
   ref command: cat *_1.fq > sample_1.fq
   
   ref command: head sample_1.fq -n 30000000 > sample_cut_1.fq
   
   software: build-in command under Linux or Mac OS

2. Mapping reads: Map raw data to reference plastid genome (could be distant speceis).
   
   ref command: bowtie2-build --large-index reference.fasta cp.index
   
   ref command: bowtie2 --very-fast-local --al-conc %.mapped.fq -q -x cp.index -1 sample_cut_1.fq -2 sample_cut_2.fq -S bowtie.sam --no-unal -t

   software: bowtie2
   
   note: bowtie.sam could be discard since not used after this process

3. Filtering by Extension in silica: Take the mapped reads (1.mapped.fq and 2.mapped.fq) as probe, and use the custom script filter_fq_with_probes.py to extend more target reads.
   
   ref command: python2 filter_fq_with_probes.py -w 93 -A -1 sample_cut_1.fq -2 sample_cut_2.fq --fastq -3 1.mapped.fq -4 2.mapped.fq
   
   software: filter_fq_with_probes.py
   
   note: the value word size (followed with "-w"), like the kmer in assembly, is crucial to the feasibility and efficiency of this process. However, theoretical best word size changes from data to data and will be affected by read length, read quality, base coverage, organ DNA percent and other factors.

4. Denovo Assembly: Utilize SPAdes to assembly the target reads produced by the former step. View contig graph and choose the path with Bandage.
   
   ref command: spades.py --careful -k 65,75,85,95 -1 sample.filtered_1.fq -2 sample.filtered_2.fq -o sample_spades
   
   software: <a href='http://bioinf.spbau.ru/spades'>SPAdes</a>, <a href='https://github.com/rrwick/Bandage'>Bandage</a>
   
   note: the best kmer depends on a wide variety of factors too.
