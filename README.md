<img src="https://img.shields.io/badge/build-Under%20dev.-yellow"/>

# YAGCloser

**Y**et **A**nother **G**ap **Closer** based on spanning long reads

A program for correction of genome assemblies focused on gap closing and filling, based on gap-spanning of long reads.

## Installation 

Requires the following `python` modules:

- `pysam`,`biopython`, `argparse`,`collections`,`copy`,`csv`,`datetime`,`filetype`, `glob`,`gzip`,`logging`,`os`,`random`,`string`,`sys`, `numpy`

Easy installation with `pip`:

```
pip3 install --user --upgrade pip
pip3 install --user pysam biopython numpy 
```
  

## The input alignment

We recommend the input alignment to be generated following these instructions. 

1. Map the long reads against the genome assembly using minimap2 (Li, 2018) and tag secondary alignments (option `--secondary=yes`). 
2. Convert the alignment into its binary form (BAM file), sort it by coordinates and index it with `samtools` (Li et al., 2009) [Version tested: 1.7 (using htslib  1.8)]. 

```
minimap2 --secondary=yes -ax [map-pb | map-ont] reference.fasta reads.fastq | samtools view -hSb - | samtools sort - > aln.s.bam
```

3. [Optional | Will help with processing speed later] Filter a subset of primary alignments with MAPQ > 20 to identify the closable gaps. 

```
minimap2 --secondary=yes -ax [map-pb | map-ont] reference.fasta reads.fastq | samtools view --hSb -q 20 - | samtools sort - > aln.s.bam
```

4. Index BAM file

```
samtools index aln.s.bam
```

5. Get a description of the gaps present in the reference file, you can use the `detgaps` too from [`asset`](https://github.com/dfguan/asset) to do so. Check on the corresponding repository how to install it.

```
detgaps reference.fasta > .gaps.bed &
```

# Usage

```
usage: yagcloser [-h] -g FASTA_FILE_PATH -b BED_FILE_PATH -a BAM_FILE_PATH 
                  -o FOLDER_PATH -s STR [-mins INT] [-f INT]
                 [-mapq <MAPQ_threshold>] [-mcc INT] [-prt FLOAT] [-eft FLOAT]
                 [-l <log_level>] [-v]
```

## Required arguments:
 - `-g FASTA FILE PATH, --genome FASTA FILE PATH`: Filepath of the reference genome file to be used. Accepts compressed files (GZIP) (default: None)
 - `-b BED FILE PATH, --bed BED FILE PATH`: Filepath of the bed file describing the gaps of the genome. Accepts compressed files (GZIP) (default:None)
- `-a BAM FILE PATH, --aln BAM FILE PATH`: Filepath of the alignment of reads to the reference genome in BAM format. This file needs to be indexed before running. (default: None)
- `-o FOLDER_PATH, --output FOLDER_PATH`: Output path folder. (default: None)
- `-s STR, --samplename STR`:  Short sample name that will be used for naming OUTPUT files. (default: None)

## Optional arguments:
- `-mins INT, --min-support INT`: Minimum number of reads needed spanning a gap to be considered for closing or filling. (default: 5)
- `-f INT, --flanksize INT`: Flank size to be used to select the reads that are in the surroundings of the gap and determine whether there are reads that span the gap or not. (default: 20)
- `-mapq <MAPQ_threshold>, --mapping-guality-threshold <MAPQ_threshold>`: MAPQ value used to filter alignments for posterior processing.Discarding alingments where: alignment_mapq < MAPQ_threshold. (default: 20)
- `-mcc INT, --min-coverage-consensus INT`: Require that more than INT sequences to be part of an alignment to put in the consensus. (default: 2)
- `-prt FLOAT, --percent-reads-threshold FLOAT`: Require that more than INT sequences to remain after the length agreement to be considered for consensus. (default: 0.5)
- `-eft FLOAT, --empty-flanks-threshold FLOAT Percentage of empty flanks required to skip an ambiguous decision on a gap. (default: 0.2)
- `-l <log_level>, --log <log_level>`:  Verbosity levels. (default: INFO)

## Information arguments:
- `-v, --version`: Show program's version number and exit
- `-h, --help`: show this help message and exit






## Output

| File/directory | Name                                                 |             Description                            | 
|----------------|------------------------------------------------------|----------------------------------------------------|
|      D         | `<samplename>.consensus`                             | Folder with FASTA files that are used in the consensus phase. |                  
|      D         | `<samplename>.flanks`                                | Folder with FASTA files with the sequences corresponding to the flank regions. |                  
|      D         | `<samplename>.fullsupport`                           | Folder with FASTA files. Content are the sequences that include the flanks and the sequencing between the flanks. |                  
|      D         | `<samplename>.log`                                   | |                  
|      D         | `<samplename>.msa`                                   | Output of the alignment of the subreads. |                  
|      D         | `<samplename>.pre`                                   | |                  
|      D         | `<samplename>.reads`                                 | Read subsequences that are aligned to the flank and sequences between the flanks. |                  
|      D         | `<samplename>.support`                               | Folder with FASTA files. Content are the sequences between the flanks. | 
|      D         | `<samplename>.consensus.log.txt`                     | |                  
|      F         | `<samplename>.alignment.err`                         | Output of all the alignments generated by MAFFT  |                  
|      F         | `<samplename>.ambiguous.txt`                         | |                  
|      F         | `<samplename>.edits.txt`                             | |                  
|      F         | `<samplename>.gaps.closed.original_coordinates.txt`  | |                  
|      F         | `<samplename>.no.length.agreement.txt`               | |                  
|      F         | `<samplename>.no.support.gaps.txt`                   | |                  
|      F         | `<samplename>.potential.fillable.gaps.txt`           | |                  


# Examples

1. Run to identify potential gaps/edits that will be done to he reference
```
 python yagcloser -g reference.fasta \
    -a reference.sorted.bam \
    -b gaps.bed \
    -o yagcloser_output \
    -f 40 - mins 2 -s reference_data
 ```

2. Edit the reference file

```
python update_assembly_edits_and_breaks.py \
    -i reference.fasta \
    -o reference.v2.fasta \
    -e yagcloser_output/yagcloser_output.edits.txt
```



