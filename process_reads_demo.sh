#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
  "/mnt/c/Users/rzou4/Downloads/demo/mm_mP9_nD24h_r3_sub" \
  "/mnt/c/Users/rzou4/Downloads/demo/mm_mP9_KU24h_r3_sub" \
  "/mnt/c/Users/rzou4/Downloads/demo/mm_negctrl_r1" \
)

# Enter path to indexed genome
mm10path="/home/roger/bioinformatics/mm10/mm10"
##########################################

# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=3                     # default number of parallel processes
  genome="mm10"
  genomepath="$mm10path"
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
