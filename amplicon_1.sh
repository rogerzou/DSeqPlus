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


cd /Volumes/Lab-Home/rzou4/NGS_data/4_damage/220118_ampNGS/
mkdir -p merged

flash -M 300 -o merged/Fs2_OFF0_ct1  raw/C1_0_S1_R1_001.fastq  raw/C1_0_S1_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF5_ct1  raw/C1_5_S2_R1_001.fastq  raw/C1_5_S2_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF7_ct1  raw/C1_7_S3_R1_001.fastq  raw/C1_7_S3_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF8_ct1  raw/C1_8_S4_R1_001.fastq  raw/C1_8_S4_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF9_ct1  raw/C1_9_S5_R1_001.fastq  raw/C1_9_S5_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF10_ct1  raw/C1_10_S6_R1_001.fastq  raw/C1_10_S6_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF11_ct1  raw/C1_11_S7_R1_001.fastq  raw/C1_11_S7_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF12_ct1  raw/C1_12_S8_R1_001.fastq  raw/C1_12_S8_R2_001.fastq  | tee flash-stats.txt

flash -M 300 -o merged/Fs2_OFF0_ep1  raw/EP1_0_S9_R1_001.fastq  raw/EP1_0_S9_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF5_ep1  raw/EP1_5_S10_R1_001.fastq  raw/EP1_5_S10_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF7_ep1  raw/EP1_7_S11_R1_001.fastq  raw/EP1_7_S11_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF8_ep1  raw/EP1_8_S12_R1_001.fastq  raw/EP1_8_S12_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF9_ep1  raw/EP1_9_S13_R1_001.fastq  raw/EP1_9_S13_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF10_ep1  raw/EP1_10_S14_R1_001.fastq  raw/EP1_10_S14_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF11_ep1  raw/EP1_11_S15_R1_001.fastq  raw/EP1_11_S15_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF12_ep1  raw/EP1_12_S16_R1_001.fastq  raw/EP1_12_S16_R2_001.fastq  | tee flash-stats.txt

flash -M 300 -o merged/Fs2_OFF0_ct2  raw/C2_0_S17_R1_001.fastq  raw/C2_0_S17_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF5_ct2  raw/C2_5_S18_R1_001.fastq  raw/C2_5_S18_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF7_ct2  raw/C2_7_S19_R1_001.fastq  raw/C2_7_S19_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF8_ct2  raw/C2_8_S20_R1_001.fastq  raw/C2_8_S20_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF9_ct2  raw/C2_9_S21_R1_001.fastq  raw/C2_9_S21_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF10_ct2  raw/C2_10_S22_R1_001.fastq  raw/C2_10_S22_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF11_ct2  raw/C2_11_S23_R1_001.fastq  raw/C2_11_S23_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF12_ct2  raw/C2_12_S24_R1_001.fastq  raw/C2_12_S24_R2_001.fastq  | tee flash-stats.txt

flash -M 300 -o merged/Fs2_OFF0_ep2  raw/EP2_0_S25_R1_001.fastq  raw/EP2_0_S25_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF5_ep2  raw/EP2_5_S26_R1_001.fastq  raw/EP2_5_S26_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF7_ep2  raw/EP2_7_S27_R1_001.fastq  raw/EP2_7_S27_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF8_ep2  raw/EP2_8_S28_R1_001.fastq  raw/EP2_8_S28_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF9_ep2  raw/EP2_9_S29_R1_001.fastq  raw/EP2_9_S29_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF10_ep2  raw/EP2_10_S30_R1_001.fastq  raw/EP2_10_S30_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF11_ep2  raw/EP2_11_S31_R1_001.fastq  raw/EP2_11_S31_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF12_ep2  raw/EP2_12_S32_R1_001.fastq  raw/EP2_12_S32_R2_001.fastq  | tee flash-stats.txt

flash -M 300 -o merged/Fs2_OFF0_ct3  raw/C3_0_S33_R1_001.fastq  raw/C3_0_S33_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF5_ct3  raw/C3_5_S34_R1_001.fastq  raw/C3_5_S34_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF7_ct3  raw/C3_7_S35_R1_001.fastq  raw/C3_7_S35_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF8_ct3  raw/C3_8_S36_R1_001.fastq  raw/C3_8_S36_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF9_ct3  raw/C3_9_S37_R1_001.fastq  raw/C3_9_S37_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF10_ct3  raw/C3_10_S38_R1_001.fastq  raw/C3_10_S38_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF11_ct3  raw/C3_11_S39_R1_001.fastq  raw/C3_11_S39_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF12_ct3  raw/C3_12_S40_R1_001.fastq  raw/C3_12_S40_R2_001.fastq  | tee flash-stats.txt

flash -M 300 -o merged/Fs2_OFF0_ep3  raw/EP3_0_S41_R1_001.fastq  raw/EP3_0_S41_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF5_ep3  raw/EP3_5_S42_R1_001.fastq  raw/EP3_5_S42_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF7_ep3  raw/EP3_7_S43_R1_001.fastq  raw/EP3_7_S43_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF8_ep3  raw/EP3_8_S44_R1_001.fastq  raw/EP3_8_S44_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF9_ep3  raw/EP3_9_S45_R1_001.fastq  raw/EP3_9_S45_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF10_ep3  raw/EP3_10_S46_R1_001.fastq  raw/EP3_10_S46_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF11_ep3  raw/EP3_11_S47_R1_001.fastq  raw/EP3_11_S47_R2_001.fastq  | tee flash-stats.txt
flash -M 300 -o merged/Fs2_OFF12_ep3  raw/EP3_12_S48_R1_001.fastq  raw/EP3_12_S48_R2_001.fastq  | tee flash-stats.txt