#!/usr/bin/env bash
cd /Volumes/Lab-Home/rzou4/NGS_data/4_damage/210720_ampNGS/
mkdir -p merged

flash -M 300 -o merged/Fs2_OFF1_ct  raw/Dseq-Ctrl1_S1_L001_R1_001.fastq  raw/Dseq-Ctrl1_S1_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ct  raw/Dseq-Ctrl2_S2_L001_R1_001.fastq  raw/Dseq-Ctrl2_S2_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ct  raw/Dseq-Ctrl3_S3_L001_R1_001.fastq  raw/Dseq-Ctrl3_S3_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ct  raw/Dseq-Ctrl4_S4_L001_R1_001.fastq  raw/Dseq-Ctrl4_S4_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ct  raw/Dseq-CtrlON_S5_L001_R1_001.fastq  raw/Dseq-CtrlON_S5_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF1_ep  raw/Dseq-EP1_S6_L001_R1_001.fastq  raw/Dseq-EP1_S6_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ep  raw/Dseq-EP2_S7_L001_R1_001.fastq  raw/Dseq-EP2_S7_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ep  raw/Dseq-EP3_S8_L001_R1_001.fastq  raw/Dseq-EP3_S8_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ep  raw/Dseq-EP4_S9_L001_R1_001.fastq  raw/Dseq-EP4_S9_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ep  raw/Dseq-EPON_S10_L001_R1_001.fastq  raw/Dseq-EPON_S10_L001_R2_001.fastq  | tee -a flash-stats.txt
