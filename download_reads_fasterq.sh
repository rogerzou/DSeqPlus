#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
download_path=/mnt/c/Users/rzou4/Downloads/
# List of all SRA samples with their labels. Comment out samples that are not needed.
declare -a filelist=(\
"SRR18188661" "HEK_WT"
"SRR18188662" "HEK_Vs3_nD12h"
"SRR18188646" "HEK_Vs3_NU12h"
"SRR18188647" "HEK_Vs3_KU12h"
#"SRR18188648" "HEK_Vs2_nD12h"
#"SRR18188649" "HEK_Vs2_KU12h"
#"SRR18188650" "HEK_Hs4_nD04h"
#"SRR18188651" "HEK_Hs4_nD24h"
#"SRR18188652" "HEK_Hs4_nD12h"
#"SRR18188653" "HEK_Hs4_KU04h"
#"SRR18188663" "HEK_Hs4_KU24h"
#"SRR18188664" "HEK_Hs4_KU12h"
#"SRR18188654" "iPSC_Vs2_nD12h"
#"SRR18188655" "iPSC_Vs2_KU12h"
#"SRR18188656" "K562_WT"
#"SRR18188657" "K562_Vs2_nD12h"
#"SRR18188658" "K562_Vs2_KU12h"
#"SRR18188659" "K562_Fs2_nD12h"
#"SRR18188660" "K562_Fs2_KU12h"
#"SRR18188709" "mm_mP9_nD24h_r4"
#"SRR18188701" "mm_mP9_nD24h_r3"
#"SRR18188702" "mm_mP9_nD24h_r2"
#"SRR18188703" "mm_mP9_nD24h_r1"
#"SRR18188704" "mm_mP9_nD24h_r0"
#"SRR18188705" "mm_mP9_KU24h_r4"
#"SRR18188706" "mm_mP9_KU24h_r3"
#"SRR18188707" "mm_mP9_KU24h_r2"
#"SRR18188708" "mm_mP9_KU24h_r1"
#"SRR18188710" "mm_mP9_KU24h_r0"
#"SRR18188711" "mm_GFP_KU24h"
)
##########################################

# processing arguments, proceed with bioinformatics pipeline
main() {
  downloadfastq "${filelist[@]}"
}

# pipeline for downloading FASTQ reads from SRA
downloadfastq() {

  cd $download_path
  mkdir SRA_download

  arraylength=${#filelist[@]}
  for (( i=0; i<${arraylength}; i+=2 ));
  do
    var1="${filelist[$i]}"
    var2="${filelist[$i+1]}"
    cd $download_path
    prefetch $var1 -O $var1
    cd $var1/$var1
    fasterq-dump "${var1}.sra"
    mv "${var1}_1.fastq" "${download_path}SRA_download/${var2}_1.fastq"
    mv "${var1}_2.fastq" "${download_path}SRA_download/${var2}_2.fastq"
  done

#  # Download DISCOVER-Seq (Wienert et al) replicate 1 (negative control)
#  cd $download_path
#  this_file=SRR8553804
#  prefetch $this_file -O $this_file
#  cd $this_file/$this_file
#  fasterq-dump $this_file.sra
#  mv "${this_file}_1.fastq" "${download_path}SRA_download/mm_negctrl_r1_1.fastq"
#  mv "${this_file}_2.fastq" "${download_path}SRA_download/mm_negctrl_r1_2.fastq"
#
#  # Download DISCOVER-Seq (Wienert et al) replicate 2 (negative control)
#  cd $download_path
#  this_file=SRR8553806
#  prefetch $this_file -O $this_file
#  cd $this_file/$this_file
#  fasterq-dump $this_file.sra
#  mv "${this_file}_1.fastq" "${download_path}SRA_download/mm_negctrl_r2_1.fastq"
#  mv "${this_file}_2.fastq" "${download_path}SRA_download/mm_negctrl_r2_2.fastq"

}

main "$@"; exit
