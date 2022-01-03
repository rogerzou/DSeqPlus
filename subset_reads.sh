#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array

#declare -a filelist1=(\
#  "/mnt/d/210216_Dseq+/A04" \
#  "/mnt/d/210216_Dseq+/A05" \
#  "/mnt/d/210216_Dseq+/A06" \
#)
#declare -a filelist2=(\
#  "/mnt/d/210216_Dseq+/A10" \
#  "/mnt/d/210216_Dseq+/A11" \
#  "/mnt/d/210216_Dseq+/A12" \
#  "/mnt/d/210216_Dseq+/A13" \
#  "/mnt/d/210216_Dseq+/A14" \
#  "/mnt/d/210216_Dseq+/A15" \
#)
#declare -a filelist3=(\
#  "/mnt/d/210216_Dseq+/A17" \
#  "/mnt/d/210216_Dseq+/A20" \
#)

#declare -a filelist1=(\
#  "/mnt/d/210301_Dseq+/A13" \
#  "/mnt/d/210301_Dseq+/A14" \
#)
#declare -a filelist2=(\
#  "/mnt/d/210301_Dseq+/A15" \
#  "/mnt/d/210301_Dseq+/A16" \
#)
#declare -a filelist3=(\
#  "/mnt/d/210301_Dseq+/A17" \
#  "/mnt/d/210301_Dseq+/A18" \
#)

#declare -a filelist1=(\
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A08_mPCSK9_nd24h_r2" \
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A11_mPCSK9_KU24h_r2" \
#)
#declare -a filelist2=(\
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A09_mPCSK9_nd24h_r3" \
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A12_mPCSK9_KU24h_r3" \
#)
#declare -a filelist3=(\
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A05_mPCSK9_nd24h_r4" \
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A06_mPCSK9_KU24h_r4" \
#)
#declare -a filelist1=(\
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A02_mPCSK9_nd24h_r0" \
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A03_mPCSK9_KU24h_r0" \
#)
#declare -a filelist2=(\
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A07_mPCSK9_nd24h_r1" \
#  "/mnt/c/Users/rzou4/Desktop/211213_Dseq+/A10_mPCSK9_KU24h_r1" \
#)


##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  subsetfastq "${filelist1[@]}"
  subsetfastq "${filelist2[@]}"
  subsetfastq "${filelist3[@]}"
}

# pipeline for subsetting FASTQ reads to
subsetfastq() {

  # calculate the mininum number of lines in a FASTQ file for one set of files
  echo "Number of fastq files: $#"
  numlines=1000000000
  for var in "$@"
  do
    numlines_tmp=$(wc -l < "${var}_1.fastq")
    if [ $numlines_tmp -lt $numlines ]
    then
      numlines=${numlines_tmp}
    fi
  done
  echo "Minimum number of lines: $numlines"

  # subset all FASTQ files in set to have the minimum number of lines
  for var in "$@"
  do
    head -n $numlines "${var}_1.fastq" > "${var}_sub_1.fastq"
    head -n $numlines "${var}_2.fastq" > "${var}_sub_2.fastq"
  done

}

main "$@"; exit
