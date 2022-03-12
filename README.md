Custom analysis software for DISCOVER-Seq+
====

## System requirements
- [Anaconda Python 3.7](https://www.anaconda.com/distribution/) (Anaconda's python distribution comes with the required numpy and scipy libraries)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
- Ensure that both `samtools` and `bowtie2` are added to path and can be called directly from bash

## Data requirements
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA801688

## Installation guide (est. 30 min)
1. Download the prebuilt bowtie2 indices for various genome assemblies
    - [Human hg38](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip)
    - [Human hg19](https://genome-idx.s3.amazonaws.com/bt/hg19.zip)
    - [Mouse mm10](https://genome-idx.s3.amazonaws.com/bt/mm10.zip)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/`,`hg19_bowtie2/`, or `mm10_bowtie2/`
2. Download two human hg19 and hg38 genome assemblies in FASTA format
    - [hg38.fa](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - [hg19.fa](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
    - [mm10.fa](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/`,`hg19_bowtie2/`, or `mm10_bowtie2/`
3. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
    - `samtools faidx hg19_bowtie2/hg19.fa`
    - `samtools faidx mm10_bowtie2/mm10.fa`
4. Download [BLENDER](https://github.com/staciawyman/blender) (Wienert et al., 2020)

## Demo
This demo will perform DISCOVER-Seq+ on C57BL/6J mouse liver, 24h after induction with adenovirus
expressing Cas9 and gRNA targeting PCSK9.
1. Open `download_reads_demo.sh`, go to "User Entry Section", enter the directory to download the sequencing data. 
2. Run `download_reads_demo.sh` on command line, which will download the relevant FASTQ files from 
the public database and rename the files appropriately.
3. Copy raw FASTQ files into a newly created `demo` folder in the main directory 
4. Run `subset_reads_demo.sh`, which will subset the relevant FASTQ files to ensure equal # of sequencing
reads between DISCOVER-Seq+ and DISCOVER-Seq, for fair comparison.
5. Open `process_reads_demo.sh`, go to "User Entry Section", enter the directory which has the 
mm10 genome previously installed via the Installation Guide.
6. Run `process_reads_demo.sh`, which will use `bowtie2` to align the files to the mouse mm10 genome, 
followed by `samtools` for further processing before conversion to indexed BAM files (est. time 12 hours).
7. Open `blender_c3_demo.sh`, go to "User Entry Section", enter the paths to 
(1) BLENDER software, (2-3) path to new output folders.
8. Run `blender_c3_demo.sh`, which runs BLENDER to output list of off-target sites from 
DISCOVER-Seq+ versus DISCOVER-Seq (est. time 3 days).
9. Expected results are included in this repository:
   - [DISCOVER-Seq+](https://github.com/rogerzou/DSeqPlus/tree/main/peaks/mouse_PCSK9_KU_r3_c3) 
   - [DISCOVER-Seq](https://github.com/rogerzou/DSeqPlus/tree/main/peaks/mouse_PCSK9_nD_r3_c3)

## Usage instructions
Similar to instructions listed in the Demo, but more general and applicable to all data sets
1. Download raw sequencing data from [SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA801688), using the SRA Toolkit.
2. Copy raw FASTQ files into a folder that holds raw FASTQ files.
3. Open `subset_reads.sh` to modify paths
4. Run `subset_reads.sh` to subset the relevant FASTQ files to ensure equal # of sequencing reads
between DISCOVER-Seq+ and DISCOVER-Seq for fair comparison.
5. Open `process_reads.sh` to modify paths
6. Run `process_reads.sh` to use `bowtie2` to align the files to the mouse mm10 genome, followed by 
`samtools` for further processing before conversion to indexed BAM files (est. time 12 hours).
7. Shell file of type `blender_*_xxxx_cc.sh` are used to perform 
[BLENDER](https://github.com/staciawyman/blender). `*` indicates sequencing file dataset,
`xxxx` indicates genome (`hg19` or `mm10`), and `cc` indicates either `c2` or `c3` BLENDER parameter.
8. A Python script `script_1.py` is used to perform downstream analysis, generating the majority
of the data and results presented in the associated manuscript.
