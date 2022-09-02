#!/usr/bin/env bash

cd ../blender_rz
mm10_bowtie2="/mnt/c/Users/rzou4/bioinformatics/mm10_bowtie2/mm10.fa"

##############################################################
## Run BLENDER on samples with Cas9 targeting PCSK9 (mouse) ##
## in mouse livers using mm10, cutoff threshold 3.          ##
##############################################################

mice_ctrl="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq_public/SRR8553804_mm10_final.bam"
mice_P9_nD_r2="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A08_mPCSK9_nd24h_r2_sub_mm10_final.bam"
mice_P9_KU_r2="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A11_mPCSK9_KU24h_r2_sub_mm10_final.bam"
mice_P9_nD_r3="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A09_mPCSK9_nd24h_r3_sub_mm10_final.bam"
mice_P9_KU_r3="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A12_mPCSK9_KU24h_r3_sub_mm10_final.bam"
mice_P9_nD_r4="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A05_mPCSK9_nd24h_r4_sub_mm10_final.bam"
mice_P9_KU_r4="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A06_mPCSK9_KU24h_r4_sub_mm10_final.bam"
mice_P9_nD_r0="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A02_mPCSK9_nd24h_r0_sub_mm10_final.bam"
mice_P9_KU_r0="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A03_mPCSK9_KU24h_r0_sub_mm10_final.bam"
mice_P9_nD_r1="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A07_mPCSK9_nd24h_r1_sub_mm10_final.bam"
mice_P9_KU_r1="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A10_mPCSK9_KU24h_r1_sub_mm10_final.bam"
mice_P9_nD_m="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/mPCSK9_nD24h_merged.bam"
mice_P9_KU_m="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/mPCSK9_KU24h_merged.bam"
P9_gRNA="AGCAGCAGCGGCGGCAACAG"
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r0_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r2_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r3_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r4_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r0_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r2_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r3_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r4_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_merged_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_merged_c3

# Mouse PCSK9 replicate 0 (no drug)
sh run_blender.sh $mm10_bowtie2 $mice_P9_nD_r0 $mice_ctrl $P9_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r0_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 0 (w/ KU-60648)
sh run_blender.sh $mm10_bowtie2 $mice_P9_KU_r0 $mice_ctrl $P9_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r0_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 1 (no drug)
sh run_blender.sh $mm10_bowtie2 $mice_P9_nD_r1 $mice_ctrl $P9_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r1_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 1 (w/ KU-60648)
sh run_blender.sh $mm10_bowtie2 $mice_P9_KU_r1 $mice_ctrl $P9_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r1_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 2 (no drug)
sh run_blender.sh $mm10_bowtie2 $mice_P9_nD_r2 $mice_ctrl $P9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r2_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 2 (w/ KU-60648)
sh run_blender.sh $mm10_bowtie2 $mice_P9_KU_r2 $mice_ctrl $P9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r2_c3 "-c=3"


# Mouse PCSK9 replicate 3 (no drug)
sh run_blender.sh $mm10_bowtie2 $mice_P9_nD_r3 $mice_ctrl $P9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r3_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 3 (w/ KU-60648)
sh run_blender.sh $mm10_bowtie2 $mice_P9_KU_r3 $mice_ctrl $P9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r3_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 4 (no drug)
sh run_blender.sh $mm10_bowtie2 $mice_P9_nD_r4 $mice_ctrl $P9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r4_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 4 (w/ KU-60648)
sh run_blender.sh $mm10_bowtie2 $mice_P9_KU_r4 $mice_ctrl $P9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r4_c3 "-c=3"


# Mouse PCSK9 merged (w/ or without KU-60648)
sh run_blender.sh $mm10_bowtie2 $mice_P9_nD_m $mice_ctrl $P9_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_merged_c3 "-c=3" \
& \
sh run_blender.sh $mm10_bowtie2 $mice_P9_KU_m $mice_ctrl $P9_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_merged_c3 "-c=3"
