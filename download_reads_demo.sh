#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
cd /mnt/c/Users/rzou4/Downloads/
##########################################


# Download DISCOVER-Seq+ replicate 3 (Cas9 targeting PCSK9 in mouse liver, with DNA-PKcs inhibitor)
prefetch SRR18188706 -O SRR18188706
cd SRR18188706/SRR18188706
fasterq-dump SRR18188706.sra
mv SRR18188706_1.fastq mPCSK9_KU24h_r3_1_fastq
mv SRR18188706_2.fastq mPCSK9_KU24h_r3_1_fastq

# Download DISCOVER-Seq replicate 3 (Cas9 targeting PCSK9 in mouse liver, without drug)
prefetch SRR18188701 -O SRR18188701
cd SRR18188701/SRR18188701
fasterq-dump SRR18188701.sra
mv SRR18188701_1.fastq mPCSK9_nD24h_r3_1_fastq
mv SRR18188701_2.fastq mPCSK9_nD24h_r3_1_fastq

# Download DISCOVER-Seq (Wienert et al) replicate 1 (negative control)
prefetch SRR8553804 -O SRR8553804
cd SRR8553804/SRR8553804
fasterq-dump SRR8553804.sra
mv SRR8553804_1.fastq negctrl_r1_1.fastq
mv SRR8553804_2.fastq negctrl_r1_2.fastq
