# PLAEC
**P**roximity **L**igation post-processing for **A**ssembly **E**rror-**C**orrection
---

## Introduction

A set of tools that allows correct assemblies using proximity ligation data.

## Goals
 - Scaffold Fixing (breaking)
 - Inversion detection and fixing

## Post-processing pipeline

### Software requirements
- [bwa](https://github.com/lh3/bwa)
- [samtools](https://github.com/samtools/samtools)
- [picard tools](https://github.com/broadinstitute/picard)


### Pipeline

1. Index previously assembled reference, using [bwa](https://github.com/lh3/bwa)
```
bwa index -a bwtsw assembled_reference.fasta
```
2. Map HiC data library to the previously assembled reference, using [bwa](https://github.com/lh3/bwa), compress and sort the results using [samtools](https://github.com/samtools/samtools)

```
bwa mem -M assembled_reference.fasta HiCLibrary.R1.fq HiCLibrary.R2.fq | samtools view -hB - | samtools sort - > HiCLibrary.s.bam
```

3. Index resulting BAM, using [samtools](https://github.com/samtools/samtools)

```
samtools index HiCLibrary.s.bam
```

4. Mark duplicates, using [picard tools](https://github.com/broadinstitute/picard)

```
java -jar /path/to/picard/picard.jar MarkDuplicates \
    I=HiCLibrary.s.bam \
    O=HiCLibrary.s.d.bam \
    METRICS_FILE=HiCLibrary.dup.metrics.txt \
    CREATE_INDEX=true \
    ASSUME_SORTED=true 
```
5. Generation of [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file with the gap information, with [describe-assembly](https://github.com/conchoecia/kmer_scaffolder/blob/master/utilities/describe-assembly.c)

```
describe-assembly -g assembled_reference.fasta -B > assembled_reference.gaps.bed
```

## PLAEC tools: Interaction matrices and scaffold breaking

**Input **

- [Required] Assembled reference ([FASTA](https://en.wikipedia.org/wiki/FASTA_format) file)
- [Required] Gap info on the assembled reference ([BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file)
- [Required] HiC mapped library ([BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file)
- [Required] Flank size (integer)
- [Optional] Mapping quality filter (Default: 20. MAPQ<20 alignments filtered out)



