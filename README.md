Custom analysis software for DISCOVER-Seq+
====

## System requirements
- [Anaconda Python 3.7](https://www.anaconda.com/distribution/) (Anaconda's python distribution comes with the required numpy and scipy libraries)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
  - For linux systems, download compiled binaries from [here](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).
  In this example, will use Ubuntu Linux 64 bit architecture, downloading the file `sratoolkit.3.0.0-ubuntu64.tar.gz`.
  - Run `tar -xvzf sratoolkit.3.0.0-ubuntu64.tar.gz` to extract folder.
  - Copy extracted folder `sratoolkit.3.0.0-ubuntu64` to permanent location, e.g. `/home/roger/bioinformatics`.
  - Add relevant functions to PATH variable by adding the following text to a new last line of `~/.bashrc`, e.g.
  `export PATH=$PATH:/home/roger/bioinformatics/sratoolkit.3.0.0-ubuntu64/bin`.

## Data requirements
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA801688

## Installation guide (est. 30 min)
1. Download the prebuilt bowtie2 indices for various genome assemblies
    - [Human hg38](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip)
    - [Human hg19](https://genome-idx.s3.amazonaws.com/bt/hg19.zip)
    - [Mouse mm10](https://genome-idx.s3.amazonaws.com/bt/mm10.zip)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/`,`hg19_bowtie2/`, or `mm10_bowtie2/`
2. Download genome assemblies in FASTA format
    - [hg38.fa](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - [hg19.fa](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
    - [mm10.fa](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/`,`hg19_bowtie2/`, or `mm10_bowtie2/`
3. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
    - `samtools faidx hg19_bowtie2/hg19.fa`
    - `samtools faidx mm10_bowtie2/mm10.fa`
4. Download [BLENDER](https://github.com/staciawyman/blender) (Wienert et al., 2020)

## Demo (please ensure 100 GB of free disk space)
This demo will perform DISCOVER-Seq+ on C57BL/6J mouse liver, 24h after induction with adenovirus
expressing Cas9 and gRNA targeting PCSK9. Test run performed on personal Windows 10 Desktop with
Ubuntu, Intel i7-8700k 3.7 GHz 6 cores, 32GB RAM.
1. On Linux/Mac command line, navigate to the home directory for this repository. 
2. Open `download_reads_demo.sh`, go to "User Entry Section", enter the desired directory to 
download the sequencing data. For the consistency of this demo, I will use `/mnt/c/Users/rzou4/Downloads/`. 
3. Run `download_reads_demo.sh` on command line, which will download the relevant FASTQ files from 
the public database, rename the files appropriately, and move files to newly created `demo` folder 
within `/mnt/c/Users/rzou4/Downloads/` *(est. time <10 min/file)*.
4. Open `subset_reads_demo.sh`, go to "User Entry Section", ensure the paths are correct. 
5. Run `subset_reads_demo.sh`, which will subset the relevant FASTQ files to ensure equal # of 
sequencing reads between DISCOVER-Seq+ and DISCOVER-Seq, for fair comparison *(est. time 1 min/file)*.
6. Open `process_reads_demo.sh`, go to "User Entry Section", enter the directory which has the 
mm10 genome previously installed via the Installation Guide.
7. Run `process_reads_demo.sh`, which will use `bowtie2` to align the files to the mouse mm10 genome, 
followed by `samtools` for further processing before conversion to indexed BAM files *(est. time 12 hours/file)*.
Performs processing on all 3 samples in parallel.
8. Open `blender_c3_demo.sh`, go to "User Entry Section", enter the paths to 
(1) BLENDER software, (2) mm10 genome, and (3-4) path to new output folders.
9. Run `blender_c3_demo.sh`, which runs BLENDER to output list of off-target sites from 
DISCOVER-Seq+ versus DISCOVER-Seq (est. time 3 days/file). Performs processing on both samples in parallel.
10. Expected results are included in this repository:
    - [DISCOVER-Seq+](https://github.com/rogerzou/DSeqPlus/tree/main/peaks/mouse_PCSK9_KU_r3_c3) 
    - [DISCOVER-Seq](https://github.com/rogerzou/DSeqPlus/tree/main/peaks/mouse_PCSK9_nD_r3_c3)

## Usage instructions
Similar to instructions listed in the Demo, but more general and applicable to all data sets.
The list of sequencing datasets is listed [here](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP362082&o=acc_s%3Aa).
#### `download_reads.sh` downloads the desired FASTQ sequencing reads from Sequencing Read Archive (SRA) (est. time 10 min/file), parallelizable
1. Open `download_reads.sh`, which has the list of FASTQ datasets to download, of format 
`SRR********`, such as `SRR18188706`, along with its corresponding name: `mm_mP9_KU24h_r3`. Enter 
the desired directory to hold the downloaded files at `download_path`. Uncomment only the files you 
would like to download. Expect the processing for each sample to take up ~33 GB of disk space.
2. Run `download_reads.sh`, which downloads the desired data from SRA using the SRA Toolkit, then 
renames and copies the downloaded files into a subfolder inside `download_path` called `SRA_download`
#### `subset_reads.sh` subsets the group of sequencing reads to ensure equal # of reads (est. time 1 min/file), parallelizable
3. Open `subset_reads.sh`, make sure `sra_path` points to the path for `SRA_download` folder.
Under `main() {}`, uncomment only the files sets you would like to subset - likely they are the same
ones previously downloaded in `download_reads.sh`. NOTE: negative control samples generally do not 
need to be subsetted b/c they are only used as a control for off-target detection using BLENDER.
4. Run `subset_reads.sh` to subset the relevant FASTQ files to ensure equal # of sequencing reads
between DISCOVER-Seq+ and DISCOVER-Seq for fair comparison.
#### `process_reads.sh` aligns FASTQ to the appropriate genome, then converts alignment to BAM files (est. time 12 hours/file), parallelizable
5. Open `process_reads.sh`, uncomment only the files you would like to align and process - likely they
are the same ones previously downloaded in `download_reads.sh`, some of which were subsetted with `subset_reads.sh`.
Also, enter valid paths to the downloaded genomes
6. Run `process_reads.sh` to use `bowtie2` to align the files to the mouse mm10 genome, followed by 
`samtools` for further processing before conversion to indexed BAM files.
#### `blender_*_xxxx_cc.sh` determines genome-wide off-target sites using [BLENDER](https://github.com/staciawyman/blender) (est. time 1-3 days/file), parallelizable
7. Shell file of type `blender_*_xxxx_cc.sh`. `*` indicates sequencing file dataset,
`xxxx` indicates genome (`hg19` or `mm10`), and `cc` indicates either `c2` or `c3` BLENDER parameter.
8. Modify `blender_*_xxxx_cc.sh` appropriately for compatibility with your desired datasets, and run.
#### `script_1.py` perform downstream analysis, generating the majority of the data and results presented in the associated manuscript
9. Modify `script_1.py` appropriately for compatibility with your desired datasets, and run.
