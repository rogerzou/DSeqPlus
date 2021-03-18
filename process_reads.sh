#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
#  "/mnt/d/210216_Dseq+/A04_sub" \
#  "/mnt/d/210216_Dseq+/A05_sub" \
#  "/mnt/d/210216_Dseq+/A06_sub" \
#  "/mnt/d/210216_Dseq+/A10_sub" \
#  "/mnt/d/210216_Dseq+/A11_sub" \
#  "/mnt/d/210216_Dseq+/A12_sub" \
#  "/mnt/d/210216_Dseq+/A13_sub" \
#  "/mnt/d/210216_Dseq+/A14_sub" \
#  "/mnt/d/210216_Dseq+/A15_sub" \
#  "/mnt/d/210216_Dseq+/A17_sub" \
#  "/mnt/d/210216_Dseq+/A20_sub" \
#  "/mnt/d/210301_Dseq+/A13_sub" \
#  "/mnt/d/210301_Dseq+/A14_sub" \
#  "/mnt/d/210301_Dseq+/A15_sub" \
#  "/mnt/d/210301_Dseq+/A16_sub" \
#  "/mnt/d/210301_Dseq+/A17_sub" \
#  "/mnt/d/210301_Dseq+/A18_sub" \
#  "/mnt/d/210301_Dseq+/A19" \
#  "/mnt/d/210301_Dseq+/A20" \
)

# Enter path to indexed genome
hg38path="/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38"
hg19path="/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19"
mm10path="/mnt/c/Users/Roger/bioinformatics/mm10_bowtie2/mm10"
##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=1                      # default number of parallel processes
  genome=""
  genomepath=""
  while getopts 'p:g:' opt; do       # pass number of threads with -p | which genome with -g
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      g) genome="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads] [-g genome for alignment (hg38, hg19, mm10)]"
      exit 1
        ;;
    esac
  done
  if [ "$genome" = "hg19" ]; then
    genomepath="$hg19path"
  elif [ "$genome" = "mm10" ]; then
    genomepath="$mm10path"
  else
    genomepath="$hg38path"
    genome="hg38"
  fi
  echo "Number of parallel processes: $numthreads"
  echo "Genome file for alignment: $genome | $genomepath"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    align2bam ${genomepath} ${file} ${genome} &
  done
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
align2bam() {

  # Align reads to either hg38 or mm10 using bowtie2
  bowtie2 -p 6 -q --local -X 1000 -x $1 \
  -1 "$2_1.fastq" -2 "$2_2.fastq" -S "$2_$3.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
  samtools view -h -S -b -F 0x08 -q 25 \
  "$2_$3.sam" > "$2_$3_unsorted.bam" ;

  # Add mate score tags to ensure that paired-end reads contain correct information about the mate reads
  samtools fixmate -m \
  "$2_$3_unsorted.bam" "$2_$3_fixmate.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_$3_fixmate.bam" > "$2_$3_sorted.bam" ;

  # Remove potential PCR duplicates
  samtools markdup -r \
  "$2_$3_sorted.bam" "$2_$3_final.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_$3_final.bam" ;

  # Retrieve read count statistics
  samtools flagstat \
  "$2_$3_final.bam" > "$2_$3_flagstats.txt"

}

main "$@"; exit
