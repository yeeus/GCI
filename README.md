# genome_assessment_tool
![7fe814001fd9da8e1ef94fb3f386129](https://github.com/yeeus/genome_assessment_tool/assets/118142448/75b978b6-a29f-4ade-b9c2-51a1c0ff60b0)



## Contents
- [Requirements](https://github.com/yeeus/genome_assessment_tool#requirements)
- [Usage](https://github.com/yeeus/genome_assessment_tool#usage)


### Requirements
For the complete pipeline, there are several necessary softwares:

- [canu](https://github.com/marbl/canu)
- seqkit
- [minimap2](https://github.com/lh3/minimap2)
- [winnowmap](https://github.com/marbl/Winnowmap)
- paftools.js

As for .., it just needs:
- python3.*

### Usage
1. Prepare parental (specific) reads (if parental sequencing data are accessible, please skip this step) 
```
# we suggest to use canu for binning
canu -haplotype \
    -p $prefix -d $dictionary \
    genomeSize=3g \
    maxThreads=$threads \
    -haplotypePat $pat \
    -haplotypeMat $mat \
    -pacbio-hifi $hifi \   ## binning ul reads with -nanopore $ul
    useGrid=false

# because there would be unknown reads which could't be reliably binned, we suggest to combine them with haplotype-specific reads
seqkit shuffle -2 -o ${canu_unknown_shuffle.fa.gz} -j $threads ${canu_unknown.fa.gz}
seqkit split2 ${canu_unknown_shuffle.fa.gz} -p 2 -j $threads
cat ${canu_mat.fa.gz} ${canu_unknown_shuffle.part_001.fa.gz} > ${canu_mat.final.fa.gz}
cat ${canu_pat.fa.gz} ${canu_unknown_shuffle.part_002.fa.gz} > ${canu_pat.final.fa.gz}
```

2. Map HiFi and/or UL reads to assemblies (using minimap2 and winnowmap)
```
# minimap2 
minimap2 -t $threads -ax map-hifi $mat_asm $mat_hifi > ${mat.minimap2.hifi.sam}   ## mapping UL reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.minimap2.hifi.sam} | samtools sort -@ $threads -o ${mat.minimap2.hifi.bam}
paftools.js sam2paf ${mat.minimap2.hifi.sam} > ${mat.minimap2.hifi.paf}


# winnowmap
meryl count k=15 output $mat_merylDB $mat_asm
meryl print greater-than distinct=0.9998 $mat_merylDB > $mat_repetitive_k15.txt
winnowmap -W $mat_repetitive_k15.txt -ax map-pb $mat_asm $mat_hifi > ${mat.winnowmap.hifi.sam}   ## mapping UL reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.winnowmap.hifi.sam} | samtools sort -@ $threads -o ${mat.winnowmap.hifi.bam}
paftools.js sam2paf ${mat.winnowmap.hifi.sam} > ${mat.winnowmap.hifi.paf}
```

3. Filter the mapping files
