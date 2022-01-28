#!/usr/bin/env bash

cd /mnt/d/211213_Dseq+
# Merge mouse PCSK9 no drug 5 replicates
samtools merge \
"mPCSK9_nD24h_merged.bam" \
"A02_mPCSK9_nd24h_r0_sub_mm10_final.bam" \
"A07_mPCSK9_nd24h_r1_sub_mm10_final.bam" \
"A08_mPCSK9_nd24h_r2_sub_mm10_final.bam" \
"A09_mPCSK9_nd24h_r3_sub_mm10_final.bam" \
"A05_mPCSK9_nd24h_r4_sub_mm10_final.bam" ;
samtools index "mPCSK9_nD24h_merged.bam" ;

# Merge mouse PCSK9 KU-60648 5 replicates
samtools merge \
"mPCSK9_KU24h_merged.bam" \
"A03_mPCSK9_KU24h_r0_sub_mm10_final.bam" \
"A10_mPCSK9_KU24h_r1_sub_mm10_final.bam" \
"A11_mPCSK9_KU24h_r2_sub_mm10_final.bam" \
"A12_mPCSK9_KU24h_r3_sub_mm10_final.bam" \
"A06_mPCSK9_KU24h_r4_sub_mm10_final.bam" ;
samtools index "mPCSK9_KU24h_merged.bam" ;
