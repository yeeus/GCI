# genome_assessment_tool
![7fe814001fd9da8e1ef94fb3f386129](https://github.com/yeeus/genome_assessment_tool/assets/118142448/75b978b6-a29f-4ade-b9c2-51a1c0ff60b0)



## Contents
- [Requirements](https://github.com/yeeus/genome_assessment_tool#requirements)
- [Usage](https://github.com/yeeus/genome_assessment_tool#usage)


### Requirements
For the complete pipeline, there are several necessary softwares:

- minimap2
- winnowmap
- paftools.js

As for .., it just needs:
- python3.*

### Usage
1. Prepare HiFi and/or UL mapping files (using [minimap2](https://github.com/lh3/minimap2) and [winnowmap](https://github.com/marbl/Winnowmap)) for two haplotypes (if don't have the parental sequencing data, we suggest to use [canu](https://github.com/marbl/canu) for binning)
```
# mapping maternal-specific HiFi and/or UL reads with minimap2 
minimap2 -t $threads -ax map-hifi $mat_asm $mat_hifi > ${mat.minimap2.hifi.sam}   ## mapping UL reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.minimap2.hifi.sam} | samtools sort -@ $threads -o ${mat.minimap2.hifi.bam}
paftools.js sam2paf ${mat.minimap2.hifi.sam} > ${mat.minimap2.hifi.paf}


# mapping maternal-specific HiFi and/or UL reads with winnowmap
meryl count k=15 output $mat_merylDB $mat_asm
meryl print greater-than distinct=0.9998 $mat_merylDB > $mat_repetitive_k15.txt
winnowmap -W $mat_repetitive_k15.txt -ax map-pb $mat_asm $mat_hifi > ${mat.winnowmap.hifi.sam}   ## mapping UL reads with -ax map-ont

samtools view -@ $threads -Sb ${mat.winnowmap.hifi.sam} | samtools sort -@ $threads -o ${mat.winnowmap.hifi.bam}
paftools.js sam2paf ${mat.winnowmap.hifi.sam} > ${mat.winnowmap.hifi.paf}
```

2. Filter the mapping files
