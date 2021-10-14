#!/usr/bin/env bash

cd /Volumes/Lab-Home/rzou4/NGS_data/4_damage/210720_ampNGS/
mkdir -p merged

flash -M 300 -o merged/Fs2_OFF1_ct1  raw/Dseq-Ctrl1_S1_L001_R1_001.fastq  raw/Dseq-Ctrl1_S1_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ct1  raw/Dseq-Ctrl2_S2_L001_R1_001.fastq  raw/Dseq-Ctrl2_S2_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ct1  raw/Dseq-Ctrl3_S3_L001_R1_001.fastq  raw/Dseq-Ctrl3_S3_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ct1  raw/Dseq-Ctrl4_S4_L001_R1_001.fastq  raw/Dseq-Ctrl4_S4_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ct1  raw/Dseq-CtrlON_S5_L001_R1_001.fastq  raw/Dseq-CtrlON_S5_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF1_ep1  raw/Dseq-EP1_S6_L001_R1_001.fastq  raw/Dseq-EP1_S6_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ep1  raw/Dseq-EP2_S7_L001_R1_001.fastq  raw/Dseq-EP2_S7_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ep1  raw/Dseq-EP3_S8_L001_R1_001.fastq  raw/Dseq-EP3_S8_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ep1  raw/Dseq-EP4_S9_L001_R1_001.fastq  raw/Dseq-EP4_S9_L001_R2_001.fastq  | tee -a flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ep1  raw/Dseq-EPON_S10_L001_R1_001.fastq  raw/Dseq-EPON_S10_L001_R2_001.fastq  | tee -a flash-stats.txt


cd /Volumes/Lab-Home/rzou4/NGS_data/4_damage/210912_ampNGS/
mkdir -p merged

flash -M 300 -o merged/Fs2_OFF1_ct2  raw/C2-1_S1_L001_R1_001.fastq  raw/C2-1_S1_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ct2  raw/C2-2_S2_L001_R1_001.fastq  raw/C2-2_S2_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ct2  raw/C2-3_S3_L001_R1_001.fastq  raw/C2-3_S3_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ct2  raw/C2-4_S4_L001_R1_001.fastq  raw/C2-4_S4_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ct2  raw/C2-ON_S5_L001_R1_001.fastq  raw/C2-ON_S5_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF1_ep2  raw/E2-1_S6_L001_R1_001.fastq  raw/E2-1_S6_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ep2  raw/E2-2_S7_L001_R1_001.fastq  raw/E2-2_S7_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ep2  raw/E2-3_S8_L001_R1_001.fastq  raw/E2-3_S8_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ep2  raw/E2-4_S9_L001_R1_001.fastq  raw/E2-4_S9_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ep2  raw/E2-ON_S10_L001_R1_001.fastq  raw/E2-ON_S10_L001_R2_001.fastq  | tee flash-stats.txt

flash -M 300 -o merged/Fs2_OFF1_ct3  raw/C3-1_S11_L001_R1_001.fastq  raw/C3-1_S11_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ct3  raw/C3-2_S12_L001_R1_001.fastq  raw/C3-2_S12_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ct3  raw/C3-3_S13_L001_R1_001.fastq  raw/C3-3_S13_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ct3  raw/C3-4_S14_L001_R1_001.fastq  raw/C3-4_S14_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ct3  raw/C3-ON_S15_L001_R1_001.fastq  raw/C3-ON_S15_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF1_ep3  raw/E3-1_S16_L001_R1_001.fastq  raw/E3-1_S16_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF2_ep3  raw/E3-2_S17_L001_R1_001.fastq  raw/E3-2_S17_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF3_ep3  raw/E3-3_S18_L001_R1_001.fastq  raw/E3-3_S18_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF4_ep3  raw/E3-4_S19_L001_R1_001.fastq  raw/E3-4_S19_L001_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_ON_ep3  raw/E3-ON_S20_L001_R1_001.fastq  raw/E3-ON_S20_L001_R2_001.fastq  | tee flash-stats.txt
