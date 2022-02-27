Custom analysis software for DISCOVER-seq+
====

## Software requirements
- [Anaconda Python 3.7](https://www.anaconda.com/distribution/) (Anaconda's python distribution comes with the required numpy and scipy libraries)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- Ensure that both `samtools` and `bowtie2` are added to path and can be called directly from bash

## Data requirements
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA801688

## Installation
1. Download sequencing reads in FASTQ format from SRA
2. Download the prebuilt bowtie2 indices for human hg19 and hg38 genome assemblies
    - [Human hg38](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip)
    - [Human hg19](https://genome-idx.s3.amazonaws.com/bt/hg19.zip)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/` and `hg19_bowtie2/`
3. Download two human hg19 and hg38 genome assemblies in FASTA format
    - [hg38.fa](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - [hg19.fa](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/` and `hg19_bowtie2/`
4. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
    - `samtools faidx hg19_bowtie2/hg19.fa`

## Usage
Bash scripts are used to automate the processing of sequencing data.

`subset_reads.sh` is first used to subset the sequencing reads such that each FASTQ file of a set 
has the same # of reads

Shell file of type `blender_*_xxxx_cc.sh` are used to perform 
[BLENDER](https://github.com/staciawyman/blender). `*` indicates sequencing file dataset,
`xxxx` indicates genome (`hg19` or `mm10`), and cc indicates either `c2` or `c3` BLENDER parameter.

A Python script `script_1.py` is used to perform analysis of various data featured in the manuscript.
