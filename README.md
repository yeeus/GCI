# Genome Continuity Inspector
## About
Genome Continuity Inspector (GCI) is an assembly assessment tool for high-quality genomes (e.g. T2T genomes), in base resolution. After stringently filtering the alignments generated by mapping long reads (PacBio HiFi and/or Oxford Nanopore long reads) back to the genome assembly, GCI will report potential assembly issues and also a score to quantify the continuity of assembly.

![https://github.com/yeeus/GCI/images/pipeline.png](https://github.com/yeeus/GCI/blob/main/images/pipeline.png)


## Contents
- [Requirements](https://github.com/yeeus/GCI#requirements)
- [Parameters](https://github.com/yeeus/GCI#parameters)
- [Usage](https://github.com/yeeus/GCI#usage)
- [Test data](https://github.com/yeeus/GCI#Test-data)  
- [Outputs](https://github.com/yeeus/GCI#outputs)
- [Benchmark](https://github.com/yeeus/GCI#benchmark)
- [Utility](https://github.com/yeeus/GCI#utility)
- [FAQ](https://github.com/yeeus/GCI#Frequently-asked-questions)
- [Citation](https://github.com/yeeus/GCI#citation)
- [Help](https://github.com/yeeus/GCI#help)
- [To do](https://github.com/yeeus/GCI#to-do)

### Requirements

- [canu](https://github.com/marbl/canu) (for trio-binning)
- [minimap2](https://github.com/lh3/minimap2) (for mapping)
- [winnowmap](https://github.com/marbl/Winnowmap) (for mapping)
- [veritymap](https://github.com/ablab/VerityMap) (for mapping)
- [samtools](https://github.com/samtools/samtools) (for sam/bam processing)
- [paftools.js](https://github.com/lh3/minimap2/blob/master/misc/paftools.js) (for converting sam to paf)

As for **GCI**, it requires:
- **python3.x** (tested in python3.10)
- pysam (stable version)
- biopython (stable version)
- numpy (stable version)
- matplotlib (stable version)
- [bamsnap](https://github.com/yeeus/bamsnap) (for plotting in utility `filter_bam.py`)


### Parameters
```
python GCI.py --help

usage: GCI.py [-r FILE] [--hifi  [...]] [--nano  [...]] [--chrs] [-R FILE] [-ts INT] [-dp FLOAT] [-d PATH] [-o STR] [-mq INT] [--mq-cutoff INT] [-ip FLOAT] [-op FLOAT] [-cp FLOAT] [-fl INT] [-p]
              [-dmin FLOAT] [-dmax FLOAT] [-ws INT] [-it STR] [-f] [-h] [-v]

A program for assessing the T2T genome

Input/Output:
  -r FILE, --reference FILE
                        The reference file
  --hifi  [ ...]        PacBio HiFi reads alignment files (at least one bam file)
  --nano  [ ...]        Oxford Nanopore long reads alignment files (at least one bam file)
  --chrs                A list of chromosomes separated by comma
  -R FILE, --regions FILE
                        Bed file containing regions
                        Be cautious! If both specify `--chrs` and `--regions`, chromosomes in regions bed file should be included in the chromosomes list
  -ts INT, --threshold INT
                        The threshold of depth to be reported as issues [0]
  -dp FLOAT, --dist-percent FLOAT
                        The distance between the candidate gap intervals for combining in chromosome units [0.005]
  -d PATH               The directory of output files [.]
  -o STR, --output STR  Prefix of output files [GCI]

Filter Options:
  -mq INT, --map-qual INT
                        Minium mapping quality for alignments [30]
  --mq-cutoff INT       The cutoff of mapping quality for keeping the alignment [50]
                        (only used when inputting more than one alignment files)
  -ip FLOAT, --iden-percent FLOAT
                        Minimum identity (num_match_res/len_aln) of alignments [0.9]
  -op FLOAT, --ovlp-percent FLOAT
                        Minimum overlapping percentage of the same read alignment if inputting more than one alignment files [0.9]
  -cp FLOAT, --clip-percent FLOAT
                        Maximum clipped percentage of the alignment [0.1]
  -fl INT, --flank-len INT
                        The flanking length of the clipped bases [15]

Plot Options:
  -p, --plot            Visualize the finally filtered whole genome (and regions if providing the option `-R`) depth [False]
  -dmin FLOAT, --depth-min FLOAT
                        Minimum depth in folds of mean coverage for plotting [0.1]
  -dmax FLOAT, --depth-max FLOAT
                        Maximum depth in folds of mean coverage for plotting [4.0]
  -ws INT, --window-size INT
                        The window size when plotting [50000]
  -it STR, --image-type STR
                        The format of the output images: png or pdf [png]

Other Options:
  -f, --force           Force rewriting of existing files [False]
  -h, --help            Show this help message and exit
  -v, --version         Show program's version number and exit

Examples:
python GCI.py -r ref.fa --hifi hifi.bam hifi.paf ... --nano nano.bam nano.paf ...
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

## If providing paf file, we recommend using paftools to convert sam to paf
paftools.js sam2paf ${mat.minimap2.hifi.sam} | sort -k6,6V -k8,8n > ${mat.minimap2.hifi.paf} ## please sort the paf file because our program don't automatically sort the file by the targets names!

# winnowmap
meryl count k=15 output $mat_merylDB $mat_asm
meryl print greater-than distinct=0.9998 $mat_merylDB > $mat_repetitive_k15.txt
winnowmap -W $mat_repetitive_k15.txt -ax map-pb $mat_asm $mat_hifi > ${mat.winnowmap.hifi.sam}
samtools view -@ $threads -Sb ${mat.winnowmap.hifi.sam} | samtools sort -@ $threads -o ${mat.winnowmap.hifi.bam}
samtools index ${mat.minimap2.hifi.bam}
paftools.js sam2paf ${mat.winnowmap.hifi.sam} | sort -k6,6V -k8,8n > ${mat.winnowmap.hifi.paf}
```

3. Filter the mapping files and get the genome continuity index

We recommend to input only one alignment file per software (minimap2 and winnowmap) using the same set of long reads. 
**Importantly,** GCI needs at least one bam file for one type of long reads, which means you'll get errors when providing only paf files.
```
# Before this, make sure you've generated the index file (.bai) for bam files
# We recommend to input one bam and one paf file produced by two softwares (for example, one bam file from winnowmap and one paf file from minimap2)
# which would be theoretically faster and consume less memory but at the cost of lower sensitivity relative to all bams 
# PDF is recommended because PNG file may lose some details though GCI will output png files by default
python GCI.py -r ref.fa --hifi hifi.bam hifi.paf --nano ont.bam ont.paf -d mat -o mat -p -it pdf ...
```

### Test data

You can first download the test files from [zenodo](https://zenodo.org/records/12748594)
```
tar zxf example.tar.gz

# Then you can run GCI on the test data (it will take some seconds)
python GCI.py -r example/MH63.fasta \
       --hifi example/MH63_winnowmap_hifi.subsample.bam example/MH63.minimap2_hifi.subsample.paf \
       -d my_test -o MH63 -f

# And you will get output files my_test/MH63.depth.gz, my_test/MH63.0.depth.bed, and my_test/MH63.gci
# which will be the same as files provided in folder example
```

### Outputs
#### only providing one type of reads
- ${dictionay}/
    - (${prefix}.gaps.bed) ## the positions of Ns in the assembly
    - ${prefix}.depth.gz ## the gzipped whole-genome depth file
    - ${prefix}.${threshold}.depth.bed ## the merged depth file in bed format
    - ${prefix}.gci ## an index file containing expected N50, observed N50, expected number of contigs, observed number of contigs, genome continuity index of each chromosome and reads type
    - (${prefix}.regins.gci) ## gci file for regions (if providing the parameter **-R**)
    - (images/) ## the depth plots across the whole chromosome (if providing the parameter **-p**)
      -  ${prefix}.${target}.${image_type}
      > ![https://github.com/yeeus/GCI/images/MH63.chr01.png](https://github.com/yeeus/GCI/blob/main/images/MH63.chr01.png)

      - (${prefix}.${target}:${start}-${end}.${image_type}) (if providing the parameter **-R**)
      > ![https://github.com/yeeus/GCI/images/MH63.Chr02_MH63:13000000-15000000.png](https://github.com/yeeus/GCI/blob/main/images/MH63.Chr02_MH63_13000000-15000000.png)
  
#### providing two types
- ${dictionay}/
  - (${prefix}.gaps.bed)
  - ${prefix}_hifi.depth.gz ## the gzipped whole-genome depth file generated by the hifi alignment file
  - ${prefix}_nano.depth.gz ## the gzipped whole-genome depth file generated by the ont alignment file
  - ${prefix}_two_type.depth.gz ## the gzipped whole-genome depth file generated by the two types of reads alignment file
  - ${prefix}_hifi.${threshold}.depth.bed ## the merged depth file in bed format generated by the hifi alignment file
  - ${prefix}_nano.${threshold}.depth.bed ## the merged depth file in bed format generated by the ont alignment file
  - ${prefix}_two_type.${threshold}.depth.bed ## the merged depth file in bed format generated by the two types of reads alignment file
  - ${prefix}.gci
  - (${prefix}.regins.gci)
  - (images/)
    - ${prefix}.${target}.${image_type}
    > ![https://github.com/yeeus/GCI/images/chm13.chr19.png](https://github.com/yeeus/GCI/blob/main/images/chm13.chr19.png)

    - (${prefix}.${target}:${start}-${end}.${image_type})
    > ![https://github.com/yeeus/GCI/images/arabidopsis.Chr2:0-500000.png](https://github.com/yeeus/GCI/blob/main/images/arabidopsis.Chr2_0-500000.png)


### Benchmark
We benchmarked GCI in many genomes (details in folder [benchmark](https://github.com/yeeus/GCI/tree/main/benchmark) and [citation](https://github.com/yeeus/GCI#citation)):
| Type of reads                |  CHM13.v.2.0   | CN1.mat.v0.9  | CN1.pat.v0.9  | HG002.mat.cur.20211005 | HG002.pat.cur.20211005 | GGswu          | Col-CEN.v1.2   | MH63RS3       |
| :--------------------------: | :----------:   | :----------:  | :----------:  | :--------------------: | :--------------------: | :---:          | :----------:   | :-----:       |
| HiFi (depth; GCI)            | ~58x; 41.8259  | ~44x; 22.8391 | ~44x; 22.4743 | ~83x; 7.2645           | ~83x; 11.9397          | ~51x; 7.9901   | ~90x; 30.7545  | ~39x; 49.8945 | 
| ONT (depth; GCI)             | ~134x; 87.0425 | ~39x; 51.5398 | ~39x; 63.0391 | ~257x; 18.3920         | ~257x; 27.1588         | ~103x; 30.4181 | ~480x; 99.9999 |       NA      | 
| HiFi + ONT                   | 87.0425        | 66.7940       | 77.8956       | 18.7177                | 27.7796                | 29.3659        | 99.9999        |       NA      | 

*Note: all the results are computed using one bam file from winnowmap and one paf file from minimap2 which would be sightly higher than all bams*

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

- GCI_score.py
  - Usage

    > Now users can provide bed files generated by GCI.py (but can't compute score for regions in this case).
  
    This script is used for computing GCI score using depth files generated by GCI.py instead of running it again ([#5](https://github.com/yeeus/GCI/issues/5))
    - Example

    ```
    python GCI_score.py -r ref.fa --hifi hifi.depth.gz --nano nano.depth.gz --two-type two_type.depth.gz  ## `--two-type` is recommended which is used to compute the final score (HiFi + Nano) in GCI.py
    ```
    > Feel free to compute score for regions ([#5](https://github.com/yeeus/GCI/issues/5)) or specific chromosomes ([#3](https://github.com/yeeus/GCI/issues/3)) as in GCI.py
  
- plot_depth.py
  - Usage
    
    This script is used to generate the depth plot as a stand-alone function ([#1](https://github.com/yeeus/GCI/issues/1)) and pdf is recommended.
    - Example
    
    ```
    python plot_depth.py -r ref.fasta --hifi hifi.depth.gz --nano nano.depth.gz -it pdf
    ```
    > Unlike GCI.py, this script will output plots only for regions ([#5](https://github.com/yeeus/GCI/issues/5)).

- convert_samtools_depth.py
  - Usage

    This script is used for converting depth file generated by `samtools depth` into the format compatible with `plot_depth.py` ([#6](https://github.com/yeeus/GCI/issues/6))
    - Example

    ```
    python convert_samtools_depth.py $samtools_depth_file output_refix
    ```
    > This will generate $output_prefix.depth.gz which can be used for `plot_depth.py`





### Frequently asked questions


**1. Why is the numerator based on N50, a single point on the Nx step plot, rather than on a more stable cumulative auN/E-size?**

Contig N50 is a well-established and widely recognized metric for assessing genome assembly continuity. Although auN/E-size provides valuable insight into assembly continuity, our analysis has shown that auN/E-size is highly correlated with N50, as demonstrated by both real and simulated data (Fig. 2 in paper). We calculated the GCI scores of several genome assemblies using both N50 and auN/E-size, and observed that the results were similar (sheet 1 of [benchmark/supplementary_tables.xlsx](https://github.com/yeeus/GCI/blob/main/benchmark/supplementary_tables.xlsx)).



**2. How to select two aligners? Compared to minimap2 and Winnowmap2, does VerityMap work ?**

VerityMap (abbreviated as VM) was designed for mapping long reads to assemblies with extra-long tandem repeats (Mikheenko et al., 2020, Bioinformatics), and similarly Winnowmap2 (WM2) was specially optimized for mapping long reads to repetitive reference sequences (Jain et al., 2020, Nat Methods). Compared to minimap2 (MM2), VerityMap and Winnowmap2 are both specially developed for complex regions. We tested and compared the performance of aligner VerityMap compared to minimap2 and Winnomap2, using rice assembly MH63 as instance. We mapped HiFi reads against the assembly and observed that VerityMap (4.5h) took more running time than Winnowmap2 (3.07h, including the k-mer library building using meryl) and minimap2 (0.17h). Using alignments from any two of the three tools, we ran the GCI workflow. WM2+MM2 and VM+MM2 yielded similar potential assembly issues and GCI scores, while WM2+VM detected fewer issues with a higher GCI score. Therefore, the combination between WM2 and MM2 is recommended. See details in [benchmark/comparing_alignment_tools.pdf](https://github.com/yeeus/GCI/blob/main/benchmark/comparing_alignment_tools.pdf).



**3. How about the performance of the GCI pipeline (RAM & time requirements).**

See the computing requirements of GCI for three model genomes (human, Arabidopsis and rice) in sheet 3 of [benchmark/supplementary_tables.xlsx](https://github.com/yeeus/GCI/blob/main/benchmark/supplementary_tables.xlsx).
GCI was given one thread in all runs (which will be improved in the future). Note that the runtime only includes running GCI and does not include the computational cost of read mapping.




### Citation
[Chen, Quanyu, et al. "GCI: a continuity inspector for complete genome assembly." bioRxiv (2024): 2024-04.](https://doi.org/10.1101/2024.04.06.588431)

### Help
If you get any problems, please raise an [issue](https://github.com/yeeus/GCI/issues) first.
For another helps, please contact quanyu_chen@outlook.com.


### To do
- for paf files: add cigar operations, as well as add query to high_qual_querys in the end instead of at reading paf files
- speed up and reduce memory usage
