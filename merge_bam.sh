#!/usr/bin/env bash

#cd /mnt/d/211213_Dseq+
## Merge mouse PCSK9 no drug 5 replicates
#samtools merge \
#"mPCSK9_nD24h_merged.bam" \
#"A02_mPCSK9_nd24h_r0_sub_mm10_final.bam" \
#"A07_mPCSK9_nd24h_r1_sub_mm10_final.bam" \
#"A08_mPCSK9_nd24h_r2_sub_mm10_final.bam" \
#"A09_mPCSK9_nd24h_r3_sub_mm10_final.bam" \
#"A05_mPCSK9_nd24h_r4_sub_mm10_final.bam" ;
#samtools index "mPCSK9_nD24h_merged.bam" ;
#
## Merge mouse PCSK9 KU-60648 5 replicates
#samtools merge \
#"mPCSK9_KU24h_merged.bam" \
#"A03_mPCSK9_KU24h_r0_sub_mm10_final.bam" \
#"A10_mPCSK9_KU24h_r1_sub_mm10_final.bam" \
#"A11_mPCSK9_KU24h_r2_sub_mm10_final.bam" \
#"A12_mPCSK9_KU24h_r3_sub_mm10_final.bam" \
#"A06_mPCSK9_KU24h_r4_sub_mm10_final.bam" ;
#samtools index "mPCSK9_KU24h_merged.bam" ;


cd /mnt/d/220701_Dseq+
## Merge Tcell TRA
#samtools merge \
#  "Tcell_Cas9_g_KU_sub_hg19_merged.bam" \
#  "Tcell_Cas9_g_KU_L1_sub_hg19_final.bam" \
#  "Tcell_Cas9_g_KU_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cas9_g_KU_sub_hg19_merged.bam"
#samtools merge \
#  "Tcell_Cas9_g_nD_sub_hg19_merged.bam" \
#  "Tcell_Cas9_g_nD_L1_sub_hg19_final.bam" \
#  "Tcell_Cas9_g_nD_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cas9_g_nD_sub_hg19_merged.bam"
#samtools merge \
#  "Tcell_Cas9_ng_KU_sub_hg19_merged.bam" \
#  "Tcell_Cas9_ng_KU_L1_sub_hg19_final.bam" \
#  "Tcell_Cas9_ng_KU_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cas9_ng_KU_sub_hg19_merged.bam"
#samtools merge \
#  "Tcell_Cas9_ng_nD_sub_hg19_merged.bam" \
#  "Tcell_Cas9_ng_nD_L1_sub_hg19_final.bam" \
#  "Tcell_Cas9_ng_nD_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cas9_ng_nD_sub_hg19_merged.bam"
#
#samtools merge \
#  "Tcell_Cpf1_g_KU_sub_hg19_merged.bam" \
#  "Tcell_Cpf1_g_KU_L1_sub_hg19_final.bam" \
#  "Tcell_Cpf1_g_KU_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cpf1_g_KU_sub_hg19_merged.bam"
#samtools merge \
#  "Tcell_Cpf1_g_nD_sub_hg19_merged.bam" \
#  "Tcell_Cpf1_g_nD_L1_sub_hg19_final.bam" \
#  "Tcell_Cpf1_g_nD_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cpf1_g_nD_sub_hg19_merged.bam"
#samtools merge \
#  "Tcell_Cpf1_ng_KU_sub_hg19_merged.bam" \
#  "Tcell_Cpf1_ng_KU_L1_sub_hg19_final.bam" \
#  "Tcell_Cpf1_ng_KU_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cpf1_ng_KU_sub_hg19_merged.bam"
#samtools merge \
#  "Tcell_Cpf1_ng_nD_sub_hg19_merged.bam" \
#  "Tcell_Cpf1_ng_nD_L1_sub_hg19_final.bam" \
#  "Tcell_Cpf1_ng_nD_L2_sub_hg19_final.bam"
#samtools index "Tcell_Cpf1_ng_nD_sub_hg19_merged.bam"


samtools merge \
  "HEK_WT_KU_r1_hg19_merged.bam" \
  "HEK_WT_KU_r1_L1_hg19_final.bam" \
  "HEK_WT_KU_r1_L1_hg19_final.bam"
samtools index "HEK_WT_KU_r1_hg19_merged.bam"
samtools merge \
  "HEK_WT_KU_r2_hg19_merged.bam" \
  "HEK_WT_KU_r2_L1_hg19_final.bam" \
  "HEK_WT_KU_r2_L1_hg19_final.bam"
samtools index "HEK_WT_KU_r2_hg19_merged.bam"
samtools merge \
  "HEK_WT_nD_r1_hg19_merged.bam" \
  "HEK_WT_nD_r1_L1_hg19_final.bam" \
  "HEK_WT_nD_r1_L1_hg19_final.bam"
samtools index "HEK_WT_nD_r1_hg19_merged.bam"
samtools merge \
  "HEK_WT_nD_r2_hg19_merged.bam" \
  "HEK_WT_nD_r2_L1_hg19_final.bam" \
  "HEK_WT_nD_r2_L1_hg19_final.bam"
samtools index "HEK_WT_nD_r2_hg19_merged.bam"

samtools merge \
  "K562_WT_KU_r1_hg19_merged.bam" \
  "K562_WT_KU_r1_L1_hg19_final.bam" \
  "K562_WT_KU_r1_L1_hg19_final.bam"
samtools index "K562_WT_KU_r1_hg19_merged.bam"
samtools merge \
  "K562_WT_KU_r2_hg19_merged.bam" \
  "K562_WT_KU_r2_L1_hg19_final.bam" \
  "K562_WT_KU_r2_L1_hg19_final.bam"
samtools index "K562_WT_KU_r2_hg19_merged.bam"
samtools merge \
  "K562_WT_nD_r1_hg19_merged.bam" \
  "K562_WT_nD_r1_L1_hg19_final.bam" \
  "K562_WT_nD_r1_L1_hg19_final.bam"
samtools index "K562_WT_nD_r1_hg19_merged.bam"
samtools merge \
  "K562_WT_nD_r2_hg19_merged.bam" \
  "K562_WT_nD_r2_L1_hg19_final.bam" \
  "K562_WT_nD_r2_L1_hg19_final.bam"
samtools index "K562_WT_nD_r2_hg19_merged.bam"

samtools merge \
  "mice_WT_KU_r1_hg19_merged.bam" \
  "mice_WT_KU_r1_L1_hg19_final.bam" \
  "mice_WT_KU_r1_L1_hg19_final.bam"
samtools index "mice_WT_KU_r1_hg19_merged.bam"
samtools merge \
  "mice_WT_KU_r2_hg19_merged.bam" \
  "mice_WT_KU_r2_L1_hg19_final.bam" \
  "mice_WT_KU_r2_L1_hg19_final.bam"
samtools index "mice_WT_KU_r2_hg19_merged.bam"
samtools merge \
  "mice_WT_KU_r3_hg19_merged.bam" \
  "mice_WT_KU_r3_L1_hg19_final.bam" \
  "mice_WT_KU_r3_L1_hg19_final.bam"
samtools index "mice_WT_KU_r3_hg19_merged.bam"
samtools merge \
  "mice_WT_nD_r1_hg19_merged.bam" \
  "mice_WT_nD_r1_L1_hg19_final.bam" \
  "mice_WT_nD_r1_L1_hg19_final.bam"
samtools index "mice_WT_nD_r1_hg19_merged.bam"
samtools merge \
  "mice_WT_nD_r2_hg19_merged.bam" \
  "mice_WT_nD_r2_L1_hg19_final.bam" \
  "mice_WT_nD_r2_L1_hg19_final.bam"
samtools index "mice_WT_nD_r2_hg19_merged.bam"
samtools merge \
  "mice_WT_nD_r3_hg19_merged.bam" \
  "mice_WT_nD_r3_L1_hg19_final.bam" \
  "mice_WT_nD_r3_L1_hg19_final.bam"
samtools index "mice_WT_nD_r3_hg19_merged.bam"
