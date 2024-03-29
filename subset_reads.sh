#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
#sra_path=/mnt/c/Users/rzou4/Downloads/SRA_download/
#sra_path=/mnt/d/220701_Dseq+/
sra_path=/mnt/d/220714_Dseq+/

declare -a filelist01=(\
  "HEK_Vs3_nD12h" \
  "HEK_Vs3_NU12h" \
  "HEK_Vs3_KU12h" \
)
declare -a filelist02=(\
  "HEK_Hs4_nD04h" \
  "HEK_Hs4_nD24h" \
  "HEK_Hs4_nD12h" \
  "HEK_Hs4_KU04h" \
  "HEK_Hs4_KU24h" \
  "HEK_Hs4_KU12h" \
)
declare -a filelist03=(\
  "HEK_Vs2_nD12h" \
  "HEK_Vs2_KU12h" \
)
declare -a filelist04=(\
  "K562_Fs2_nD12h" \
  "K562_Fs2_KU12h" \
)
declare -a filelist05=(\
  "K562_Vs2_nD12h" \
  "K562_Vs2_KU12h" \
)
declare -a filelist06=(\
  "iPSC_Vs2_nD12h" \
  "iPSC_Vs2_KU12h" \
)
declare -a filelist07=(\
  "mm_mP9_nD24h_r2" \
  "mm_mP9_KU24h_r2" \
)
declare -a filelist08=(\
  "mm_mP9_nD24h_r3" \
  "mm_mP9_KU24h_r3" \
)
declare -a filelist09=(\
  "mm_mP9_nD24h_r4" \
  "mm_mP9_KU24h_r4" \
)
declare -a filelist10=(\
  "mm_mP9_nD24h_r0" \
  "mm_mP9_KU24h_r0" \
)
declare -a filelist11=(\
  "mm_mP9_nD24h_r1" \
  "mm_mP9_KU24h_r1" \
)

declare -a filelist12=(\
  "Tcell_Cas9_g_KU_L1" \
  "Tcell_Cas9_g_KU_L2" \
  "Tcell_Cas9_g_nD_L1" \
  "Tcell_Cas9_g_nD_L2" \
)

declare -a filelist13=(\
  "Tcell_Cas9_ng_KU_L1" \
  "Tcell_Cas9_ng_KU_L2" \
  "Tcell_Cas9_ng_nD_L1" \
  "Tcell_Cas9_ng_nD_L2" \
)

declare -a filelist14=(\
  "Tcell_Cpf1_g_KU_L1" \
  "Tcell_Cpf1_g_KU_L2" \
  "Tcell_Cpf1_g_nD_L1" \
  "Tcell_Cpf1_g_nD_L2" \
)

declare -a filelist15=(\
  "Tcell_Cpf1_ng_KU_L1" \
  "Tcell_Cpf1_ng_KU_L2" \
  "Tcell_Cpf1_ng_nD_L1" \
  "Tcell_Cpf1_ng_nD_L2" \
)

declare -a filelist16=(\
  "Tcell_Cas9_g_KU_r2" \
  "Tcell_Cas9_g_nD_r2" \
  "Tcell_Cas9_g_nD_NT-HDRT" \
)

declare -a filelist17=(\
  "Tcell_Cas9_ng_KU_r2" \
  "Tcell_Cas9_ng_nD_r2" \
)

declare -a filelist18=(\
  "Tcell_Cpf1_g_KU_r2" \
  "Tcell_Cpf1_g_nD_r2" \
)

declare -a filelist19=(\
  "Tcell_Cpf1_ng_KU_r2" \
  "Tcell_Cpf1_ng_nD_r2" \
)

# processing arguments, proceed with bioinformatics pipeline
main() {
#  subsetfastq "${filelist01[@]}"
#  subsetfastq "${filelist02[@]}"
#  subsetfastq "${filelist03[@]}"
#  subsetfastq "${filelist04[@]}"
#  subsetfastq "${filelist05[@]}"
#  subsetfastq "${filelist06[@]}"
#  subsetfastq "${filelist07[@]}"
#  subsetfastq "${filelist08[@]}"
#  subsetfastq "${filelist09[@]}"
#  subsetfastq "${filelist10[@]}"
#  subsetfastq "${filelist11[@]}"
#  subsetfastq "${filelist12[@]}"
#  subsetfastq "${filelist13[@]}"
#  subsetfastq "${filelist14[@]}"
#  subsetfastq "${filelist15[@]}"
  subsetfastq "${filelist16[@]}"
  subsetfastq "${filelist17[@]}"
  subsetfastq "${filelist18[@]}"
  subsetfastq "${filelist19[@]}"
}

##########################################

# pipeline for subsetting FASTQ reads to
subsetfastq() {

  # calculate the mininum number of lines in a FASTQ file for one set of files
  echo "Number of fastq files: $#"
  numlines=1000000000
  for var in "$@"
  do
    numlines_tmp=$(wc -l < "${sra_path}${var}_1.fastq")
    if [ $numlines_tmp -lt $numlines ]
    then
      numlines=${numlines_tmp}
    fi
  done
  echo "Minimum number of lines: $numlines"

  # subset all FASTQ files in set to have the minimum number of lines
  for var in "$@"
  do
    head -n $numlines "${sra_path}${var}_1.fastq" > "${sra_path}${var}_sub_1.fastq"
    head -n $numlines "${sra_path}${var}_2.fastq" > "${sra_path}${var}_sub_2.fastq"
  done

}

main "$@"; exit
