#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
#sra_path=/mnt/c/Users/rzou4/Downloads/SRA_download/
sra_path=/mnt/d/220701_Dseq+/
# sra_path=/mnt/d/220714_Dseq+/
# sra_path=/mnt/d/211213_Dseq_public/

declare -a filelist=(\
#"HEK_WT"
#"HEK_Vs3_nD12h_sub"
#"HEK_Vs3_NU12h_sub"
#"HEK_Vs3_KU12h_sub"
#"HEK_Vs2_nD12h_sub"
#"HEK_Vs2_KU12h_sub"
#"HEK_Hs4_nD04h_sub"
#"HEK_Hs4_nD24h_sub"
#"HEK_Hs4_nD12h_sub"
#"HEK_Hs4_KU04h_sub"
#"HEK_Hs4_KU24h_sub"
#"HEK_Hs4_KU12h_sub"
#"iPSC_Vs2_nD12h_sub"
#"iPSC_Vs2_KU12h_sub"
#"K562_WT"
#"K562_Vs2_nD12h_sub"
#"K562_Vs2_KU12h_sub"
#"K562_Fs2_nD12h_sub"
#"K562_Fs2_KU12h_sub"
#"mm_mP9_nD24h_r4_sub"
#"mm_mP9_nD24h_r3_sub"
#"mm_mP9_nD24h_r2_sub"
#"mm_mP9_nD24h_r1_sub"
#"mm_mP9_nD24h_r0_sub"
#"mm_mP9_KU24h_r4_sub"
#"mm_mP9_KU24h_r3_sub"
#"mm_mP9_KU24h_r2_sub"
#"mm_mP9_KU24h_r1_sub"
#"mm_mP9_KU24h_r0_sub"
#"mm_GFP_KU24h"
#"mm_negctrl_r1"
#"mm_negctrl_r2"
#"Tcell_Cas9_g_KU_L1_sub"
#"Tcell_Cas9_g_KU_L2_sub"
#"Tcell_Cas9_g_nD_L1_sub"
#"Tcell_Cas9_g_nD_L2_sub"
#"Tcell_Cas9_ng_KU_L1_sub"
#"Tcell_Cas9_ng_KU_L2_sub"
#"Tcell_Cas9_ng_nD_L1_sub"
#"Tcell_Cas9_ng_nD_L2_sub"
#"Tcell_Cpf1_g_KU_L1_sub"
#"Tcell_Cpf1_g_KU_L2_sub"
#"Tcell_Cpf1_g_nD_L1_sub"
#"Tcell_Cpf1_g_nD_L2_sub"
#"Tcell_Cpf1_ng_KU_L1_sub"
#"Tcell_Cpf1_ng_KU_L2_sub"
#"Tcell_Cpf1_ng_nD_L1_sub"
#"Tcell_Cpf1_ng_nD_L2_sub"
# "SRR8550683"
# "SRR8550685"
# "SRR8550686"
"HEK_WT_KU_r1_L1"
"HEK_WT_KU_r1_L2"
"HEK_WT_KU_r2_L1"
"HEK_WT_KU_r2_L2"
"HEK_WT_nD_r1_L1"
"HEK_WT_nD_r1_L2"
"HEK_WT_nD_r2_L1"
"HEK_WT_nD_r2_L2"
"K562_WT_KU_r1_L1"
"K562_WT_KU_r1_L2"
"K562_WT_KU_r2_L1"
"K562_WT_KU_r2_L2"
"K562_WT_nD_r1_L1"
"K562_WT_nD_r1_L2"
"K562_WT_nD_r2_L1"
"K562_WT_nD_r2_L2"
"mice_WT_KU_r1_L1"
"mice_WT_KU_r1_L2"
"mice_WT_KU_r2_L1"
"mice_WT_KU_r2_L2"
"mice_WT_KU_r3_L1"
"mice_WT_KU_r3_L2"
"mice_WT_nD_r1_L1"
"mice_WT_nD_r1_L2"
"mice_WT_nD_r2_L1"
"mice_WT_nD_r2_L2"
"mice_WT_nD_r3_L1"
"mice_WT_nD_r3_L2"
#"Tcell_Cas9_ng_nD_r2_sub"
#"Tcell_Cas9_ng_KU_r2_sub"
#"Tcell_Cas9_g_nD_r2_sub"
#"Tcell_Cas9_g_KU_r2_sub"
#"Tcell_Cpf1_ng_nD_r2_sub"
#"Tcell_Cpf1_ng_KU_r2_sub"
#"Tcell_Cpf1_g_nD_r2_sub"
#"Tcell_Cpf1_g_KU_r2_sub"
#"Tcell_Cas9_g_nD_NT-HDRT_sub"
#"Tcell_nD_T-HDRT"
)

# Enter path to indexed genome
hg38path="/mnt/c/users/rzou4/bioinformatics/hg38_bowtie2/hg38"
hg19path="/mnt/c/users/rzou4/bioinformatics/hg19_bowtie2/hg19"
mm10path="/mnt/c/users/rzou4/bioinformatics/mm10_bowtie2/mm10"

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

  cd $sra_path

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
