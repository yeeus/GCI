# Genome Continuity Index
## About
Genome Continuity Index (GCI) is a script to assess the assembly continuity for high-quality genomes (e.g. T2T genomes), in a high resolution (base-level). After analyzing the alignment files generated by mapping long reads (PacBio HiFi and/or Oxford Nanopore long reads) to the T2T genome, GCI can provide an index to represent the continuity of the genome.

![7fe814001fd9da8e1ef94fb3f386129](https://github.com/yeeus/genome_assessment_tool/assets/118142448/75b978b6-a29f-4ade-b9c2-51a1c0ff60b0)



## Contents
- [Requirements](https://github.com/yeeus/GCI#requirements)
- [Parameters](https://github.com/yeeus/GCI#parameters)
- [Usage](https://github.com/yeeus/GCI#usage) 
- [Outputs](https://github.com/yeeus/GCI#outputs)
- [Benchmark](https://github.com/yeeus/GCI#benchmark)
- [Utility](https://github.com/yeeus/GCI#utility)
- [Citation](https://github.com/yeeus/GCI#citation)
- [Help](https://github.com/yeeus/GCI#help)
- [To do](https://github.com/yeeus/GCI#to-do)

### Requirements
For the complete pipeline, there are several necessary softwares:

- [canu](https://github.com/marbl/canu) (for trio-binning)
- [minimap2](https://github.com/lh3/minimap2) (for mapping)
- [winnowmap](https://github.com/marbl/Winnowmap) (for mapping)
- [samtools](https://github.com/samtools/samtools) (for sam/bam processing)
- [paftools.js](https://github.com/lh3/minimap2/blob/master/misc/paftools.js) (for converting sam to paf)

As for **GCI**, it requires:
- **python3.x** (tested in python3.10)
- pysam (stable version)
- numpy (stable version)
- matplotlib (stable version)
- [bamsnap](https://github.com/yeeus/bamsnap) (for plotting)


### Parameters
```
python GCI.py --help

usage: GCI.py [--hifi  [...]] [--nano  [...]] [-ts INT] [-dp FLOAT] [-d PATH] [-o STR] [-mq INT] [--mq-cutoff INT] [-ip FLOAT] [-op FLOAT]
              [-cp FLOAT] [-fl INT] [-p] [-dmin FLOAT] [-dmax FLOAT] [-ws FLOAT] [-it STR] [-g] [-f] [-h] [-v]

A program for assessing the T2T genome

Input/Output:
  --hifi  [ ...]        PacBio HiFi reads alignment files (at least one bam file)
  --nano  [ ...]        Oxford Nanopore long reads alignment files (at least one bam file)
  -ts INT, --threshold INT
                        The threshold of depth in the final bed file [0]
  -dp FLOAT, --dist-percent FLOAT
                        The percentage of the distance between the candidate gap intervals in the whole chromosome (contig) [0.005]
  -d PATH               The directory of output files [.]
  -o STR, --output STR  Prefix of output files [GCI]

Filter Options:
  -mq INT, --map-qual INT
                        Minium mapping quality for alignments [30]
  --mq-cutoff INT       The cutoff of mapping quality for keeping the alignment [50]
  -ip FLOAT, --iden-percent FLOAT
                        Minimum identity (num_match_res/len_aln) of the reads [0.9]
  -op FLOAT, --ovlp-percent FLOAT
                        Minimum overlapping percentage of the reads if inputting more than one alignment files [0.9]
  -cp FLOAT, --clip-percent FLOAT
                        Maximum clipped percentage of the reads [0.1]
  -fl INT, --flank-len INT
                        The flanking length of the clipped bases [15]

Plot Options:
  -p, --plot            Visualize the finally filtered whole genome depth
  -dmin FLOAT, --depth-min FLOAT
                        Minimum depth in folds of mean coverage for plotting [0.1]
  -dmax FLOAT, --depth-max FLOAT
                        Maximum depth in folds of mean coverage for plotting [4.0]
  -ws FLOAT, --window-size FLOAT
                        The window size in chromosome units (0-1) when plotting [0.001]
  -it STR, --image-type STR
                        The format of the output images: png or pdf [png]

Other Options:
  -g, --generate        Generate the depth files
  -f, --force           Force rewriting of existing files
  -h, --help            Show this help message and exit
  -v, --version         Show program's version number and exit

Examples:
python GCI.py --hifi hifi.bam hifi.paf ... --nano nano.bam nano.paf ...
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
cat ${canu_mat.fa.gz} ${canu_unknown.fa.gz} > ${canu_mat.final.fa.gz}
cat ${canu_pat.fa.gz} ${canu_unknown.fa.gz} > ${canu_pat.final.fa.gz}
```

2. Map HiFi and/or ONT reads to assemblies (using minimap2 and winnowmap)
```
# minimap2 
minimap2 -t $threads -ax map-hifi $mat_asm $mat_hifi > ${mat.minimap2.hifi.sam}   ## mapping ONT reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.minimap2.hifi.sam} | samtools sort -@ $threads -o ${mat.minimap2.hifi.bam}
samtools index ${mat.minimap2.hifi.bam} ## this is necessary!!!
paftools.js sam2paf (-p) ${mat.minimap2.hifi.sam} | sort -k6,6V -k8,8n > ${mat.minimap2.hifi.paf} ## We recommend to use "-p" to filter the supplementary alignments
                                                                                                  ## please sort the paf file because our program don't automatically sort the file by the targets names!

# winnowmap
meryl count k=15 output $mat_merylDB $mat_asm
meryl print greater-than distinct=0.9998 $mat_merylDB > $mat_repetitive_k15.txt
winnowmap -W $mat_repetitive_k15.txt -ax map-pb $mat_asm $mat_hifi > ${mat.winnowmap.hifi.sam}   ## mapping ONT reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.winnowmap.hifi.sam} | samtools sort -@ $threads -o ${mat.winnowmap.hifi.bam}
samtools index ${mat.minimap2.hifi.bam} ## this is necessary!!!
paftools.js sam2paf (-p) ${mat.winnowmap.hifi.sam} | sort -k6,6V -k8,8n > ${mat.winnowmap.hifi.paf} ## We recommend to use "-p" to filter the supplementary alignments
                                                                                                    ## please sort the paf file because our program don't automatically sort the file by the targets names!
```

3. Filter the mapping files and get the genome continuity index

We recommend to input only one alignment file per software (minimap2 and winnowmap) using the same set of long reads. **Importantly,** there needs at least one bam file for one type of long reads, which means you'll get errors when providing only paf files.
```
# Before this, make sure you've generated the index file (.bai) for bam files
# we recommend to input one bam and one paf file produced by two softwares (for example, one bam file from winnowmap and one paf file from minimap2)
python GCI.py --hifi hifi.bam hifi.paf --nano ont.bam ont.paf -d mat -o mat -p ...
```

### Outputs
#### only providing one type of reads
- ${dictionay}/
    - (${prefix}.depth) ## the whole-genome depth file (only generated with the parameter **-g**)
    - ${prefix}.${threshold}.depth.bed ## the merged depth file in bed format
    - ${prefix}.gci ## an index file containing expected N50, observed N50, expected number of contigs, observed number of contigs, genome continuity index of each chromosome and reads type
    - (images/) ## the depth plots across the whole chromosome (if providing the parameter **-p**)
      -  ${prefix}.${target}.${image_type}
      > ![https://github.com/yeeus/GCI/images/MH63.chr01.png](https://github.com/yeeus/GCI/blob/main/images/MH63.chr01.png)
  
#### providing two types
- ${dictionay}/
  - (${prefix}_hifi.depth) ## the whole-genome depth file generated by the hifi alignment file
  - (${prefix}_nano.depth) ## the whole-genome depth file generated by the ont alignment file
  - (${prefix}_two_type.depth) ## the whole-genome depth file generated by the two types of reads alignment file
  - ${prefix}_hifi.${threshold}.depth.bed ## the merged depth file in bed format generated by the hifi alignment file
  - ${prefix}_nano.${threshold}.depth.bed ## the merged depth file in bed format generated by the ont alignment file
  - ${prefix}_two_type.${threshold}.depth.bed ## the merged depth file in bed format generated by the two types of reads alignment file
  - ${prefix}.gci
  - (images/)
    - ${prefix}.${target}.${image_type}
    > ![https://github.com/yeeus/GCI/images/chm13.chr19.png](https://github.com/yeeus/GCI/blob/main/images/chm13.chr19.png)

### Benchmark
We benchmarked GCI in many genomes:
| Type of reads                |  CHM13.v.2.0   | CN1.mat.v0.9  | CN1.pat.v0.9  | HG002.mat.cur.20211005 | HG002.pat.cur.20211005 | GGswu          | Col-CEN.v1.2   | MH63RS3       |
| :--------------------------: | :----------:   | :----------:  | :----------:  | :--------------------: | :--------------------: | :---:          | :----------:   | :-----:       |
| HiFi (depth; GCI)            | ~58x; 41.8259  | ~44x; 22.8391 | ~44x; 22.4743 | ~83x; 7.9651           | ~83x; 14.9464          | ~51x; 7.9901   | ~90x; 30.7545  | ~39x; 49.8945 | 
| Nano (depth; GCI)            | ~134x; 87.0425 | ~39x; 51.5398 | ~39x; 63.0391 | ~257x; 34.9807         | ~257x; 35.4896         | ~103x; 30.5893 | ~480x; 99.9999 |               | 
| HiFi + Nano                  | 87.0425        | 66.7940       | 77.8956       | 36.6257                | 41.6960                | 29.5146        | 99.9999        |               | 

Note: all the results are computed using one bam file from winnowmap and one paf file from minimap2

### Utility
- filter_bam.py
  - Usage

    This script is used to get the filtered bam files (and get the final plots).
    - Example1

      After getting the ${prefix}.${threshold}.depth.bed file, we'd like to get the detailed filtering information. So, we can just extract the alignments in the bed file:
      ```
      samtools view -@ $threads -Sb -L ${prefix}.${threshold}.depth.bed ${hifi.bam} > test.bam
      samtools index test.bam
      ```
      Then we can input the bam file with other paf file(s):
      ```
      python filter_bam.py test.bam test.paf -d test -o test ## if no prefix provided, the output file would be `test.filter.bam`
      ```
      Finally we would get the filtered bam file `test.bam`. Next, we can visualize the raw and filtered bam files in [IGV](https://github.com/igvteam/igv):
      > ![https://github.com/yeeus/GCI/images/igv_test.png](https://github.com/yeeus/GCI/blob/main/images/igv_test.png)
      > The top track is the raw bam file and the bottom is the filtered.
      >
      > Don't be confused about the results, which are generated by filtering both test.bam and test.paf files, while only the former is plotted here. If there is any contradiction with the final bed file, please check the paf file to figure out the conflict.
      
    - Example2

      We can immediately visualize the alignments after getting the filtered file without using IGV. We provide a convenient method for visualizing via [bamsnap](https://github.com/yeeus/bamsnap).
      ```
      # first install the prerequisites and bamsnap
      # then we can visualize the alignments in one command
      python filter_bam.py test.bam test.paf -d test -o test -p -ref test.fasta -r chrxxx:xxx-xxx  ## Be cautious, in this case, the filtered file is `test.filter.bam`
      ```
   
      The result:
      > ![https://github.com/yeeus/GCI/images/bamsnap_test.png](https://github.com/yeeus/GCI/blob/main/images/bamsnap_test.png)
      
### Citation
writing...

### Help
If you get any problems, please raise an [issue](https://github.com/yeeus/GCI/issues) first.
For another helps, please contact quanyu_chen@outlook.com.


### To do
- multiple threads
