#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array

declare -a filelist1=(\
  "demo/mPCSK9_nd24h_r3" \
  "demo/mPCSK9_KU24h_r3" \
)

##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  subsetfastq "${filelist1[@]}"
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
