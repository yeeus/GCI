# Genome Continuity Index
## About
Genome Continuity Index (GCI) is a program for assessing the T2T genome, which has a high resolution (in theory). After analyzing the alignment files generated by mapping long reads (PacBio HiFi and/or Oxford Nanopore long reads) to the T2T genome, GCI can provide an index to represent the continuity of the genome.

**As for the version 1.0, we can only analyze the case of inputting one paf and one bam and can only parse one type of reads alignment. Also, we can only use one thread..**

![7fe814001fd9da8e1ef94fb3f386129](https://github.com/yeeus/genome_assessment_tool/assets/118142448/75b978b6-a29f-4ade-b9c2-51a1c0ff60b0)



## Contents
- [Requirements](https://github.com/yeeus/GCI#requirements)
- [Parameters](https://github.com/yeeus/GCI#parameters)
- [Usage](https://github.com/yeeus/GCI#usage)
- [Outputs](https://github.com/yeeus/GCI#outputs)
- [Citation](https://github.com/yeeus/GCI#citation)
- [Help](https://github.com/yeeus/GCI#help)
- [To do](https://github.com/yeeus/GCI#to-do)

### Requirements
For the complete pipeline, there are several necessary softwares:

- [canu](https://github.com/marbl/canu) (for trio-binning)
- [seqkit](https://github.com/shenwei356/seqkit) (for fa processing)
- [minimap2](https://github.com/lh3/minimap2)
- [winnowmap](https://github.com/marbl/Winnowmap)
- [paftools.js](https://github.com/lh3/minimap2/blob/master/misc/paftools.js)

As for **GCI**, it requires:
- python3.*
- pysam
- numpy
- karyotypeR (for plotting) ###

### Parameters
```
python GCI.py --help

usage: ../GCI.py [--hifi FILE FILE] [--nano FILE FILE] [-d [PATH]] [-o [STR]] [-t [INT]] [-mq [INT]] [-ip [FLOAT]] [-op [FLOAT]] [-cp [FLOAT]] [-fl [INT]]
                 [-ts [INT]] [-f] [-p] [-h] [-v]

A program for assessing the T2T genome

Input/Output:
  --hifi FILE FILE      PacBio HiFi reads alignment files (.bam and/or .paf) including at least one bam file
  --nano FILE FILE      Oxford Nanopore long reads alignment files (.bam and/or .paf) including at least one bam file
  -d [PATH]             The dictionary of output files [.]
  -o [STR], --output [STR]
                        Prefix of output files [GCI]
  -t [INT], --threads [INT]
                        Number of threads [1]

Filter Options:
  -mq [INT], --map-qual [INT]
                        Minium mapping quality for alignments in both bam and paf files [30]
  -ip [FLOAT], --iden-percent [FLOAT]
                        Minimum identity (num_match_res/len_aln) of the reads in paf files [0.9]
  -op [FLOAT], --ovlp-percent [FLOAT]
                        Minimum overlapping percentage of the reads if inputting more than one file [0.9]
  -cp [FLOAT], --clip-percent [FLOAT]
                        Maximum clipped percentage of the reads in bam files [0.1]
  -fl [INT], --flank-len [INT]
                        The flanking length of the clipped bases [10]
  -ts [INT], --threshold [INT]
                        The threshold of depth in the final bed file [0]

Other Options:
  -f, --force           Force rewriting of existing files
  -p, --plot            Visualize the final result
  -h, --help            Show this help message and exit
  -v, --version         Show program's version number and exit
```

### Usage
1. **(For haplotype-resolved genome)** Prepare parental (specific) reads (if parental sequencing data are available, please skip this step) 
```
# we recommend to use canu for binning
canu -haplotype \
    -p $prefix -d $dictionary \
    genomeSize=3g \
    maxThreads=$threads \
    -haplotypePat $pat \
    -haplotypeMat $mat \
    -pacbio-hifi $hifi \   ## binning ONT reads with -nanopore $ont
    useGrid=false

# because there would be unknown reads which could't be reliably binned, we suggest to combine them with haplotype-specific reads
seqkit shuffle -2 -o ${canu_unknown_shuffle.fa.gz} -j $threads ${canu_unknown.fa.gz}
seqkit split2 ${canu_unknown_shuffle.fa.gz} -p 2 -j $threads
cat ${canu_mat.fa.gz} ${canu_unknown_shuffle.part_001.fa.gz} > ${canu_mat.final.fa.gz}
cat ${canu_pat.fa.gz} ${canu_unknown_shuffle.part_002.fa.gz} > ${canu_pat.final.fa.gz}
```

2. Map HiFi and/or ONT reads to assemblies (using minimap2 and winnowmap)
```
# minimap2 
minimap2 -t $threads -ax map-hifi $mat_asm $mat_hifi > ${mat.minimap2.hifi.sam}   ## mapping ONT reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.minimap2.hifi.sam} | samtools sort -@ $threads -o ${mat.minimap2.hifi.bam}
samtools index ${mat.minimap2.hifi.bam} ## this is necessary!!!
paftools.js sam2paf ${mat.minimap2.hifi.sam} | sort -k6,6V -k8,8n > ${mat.minimap2.hifi.paf} ## please sort the paf file because our program don't automatically sort the file by the targets names!

# winnowmap
meryl count k=15 output $mat_merylDB $mat_asm
meryl print greater-than distinct=0.9998 $mat_merylDB > $mat_repetitive_k15.txt
winnowmap -W $mat_repetitive_k15.txt -ax map-pb $mat_asm $mat_hifi > ${mat.winnowmap.hifi.sam}   ## mapping ONT reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.winnowmap.hifi.sam} | samtools sort -@ $threads -o ${mat.winnowmap.hifi.bam}
samtools index ${mat.minimap2.hifi.bam} ## this is necessary!!!
paftools.js sam2paf ${mat.winnowmap.hifi.sam} | sort -k6,6V -k8,8n > ${mat.winnowmap.hifi.paf} ## please sort the paf file because our program don't automatically sort the file by the targets names!
```

3. Filter the mapping files and get the genome continuity index

We recommend to input only one alignment file per software (minimap2 and winnowmap) using the same set of long reads. **Importantly,** there needs at least one bam file for one type of long reads, which means you'll get errors when providing only paf files.
```
# Before this, make sure you've generated the index file (.bai) for bam files
# we recommend to input one bam and one paf file produced by two softwares (for example, one bam file from minimap2 and one paf file from winnowmap)
python GCI.py --hifi hifi.bam hifi.paf (--nano ont.bam ont.paf) ## v.1.0 can only process one type of reads alignment files 
```

### Outputs
- ${dictionay}
    - ${prefix}.depth ## the whole-genome depth file
    - ${prefix}.${threshold}.depth.bed ## the merged depth file in bed format
    - ${prefix}.gci ## an index file containing the original N50, new N50, genome continuity index

### Citation


### Help
If you get any problems, please raise an [issue](https://github.com/yeeus/GCI/issues) first.
For another helps, please contact quanyu_chen@outlook.com.


### To do
- visualization
- contribution for N50 per chromosome
- bam + bam
- ONT intergration
- multiple threads
- benchmark
  - clip_percent: 0.05, 0.1, 0.2
  - flank_len: 10, 20, 50
  - performance: distribution, index
- edit the function (compute_index)
  - formula
  - merge gaps for calculating the number of contigs
- simulated
