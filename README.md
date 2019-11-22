# YAGCloser

**Y**et **A**nother **G**ap **Closer** based on spanning long reads


## Introduction

A program for correction of genome assemblies focused on gap closing and filling, starting from long-reads, and based on overelapping reads, lenght agreement and a lenient consensus.

## The input alignment

We recommend the input alignment to be generated following these instructions. 

1. Map the long reads against the genome assembly using minimap2 (Li, 2018) and tag secondary alignments (option --secondary=yes). 
2. Convert the alignment into its binary form (BAM file), sort it by coordinates and index it with samtools(Li et al., 2009) [Version tested: 1.7 (using htslib  1.8)]. 

```
minimap2 --secondary=yes -ax [map-pb | map-ont] reference.fasta reads.fastq | samtools view -hSb - | samtools sort - > aln.s.bam
```

3. [Optional | Will help with processing speed later] Filter a subset of primary alignments with MAPQ > 20 to identify the closable gaps. 

```
minimap2 --secondary=yes -ax [map-pb | map-ont] reference.fasta reads.fastq | samtools view --hSb -q 20 - | samtools sort - > aln.s.bam
```


## TODO: 
yagloser

test data
simulation 
assembly
make arbitrary gaps of arbitrary lengths
different coverages
how well they are closed
how the consensus sequences are

simulated
bioinformatics application notes

