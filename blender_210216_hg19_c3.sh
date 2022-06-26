#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting VEGFAs2, VEGFAs3, and HEKs4
# in HEK293T cells using hg19, cutoff threshold 3.

cd /home/roger/bioinformatics/blender
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_NU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD04_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD12_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD24_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU04_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU12_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU24_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_KU_r1_c3
hg19_bowtie2="/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"
HEK_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam"
HEK_Vs3_nD="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A04_sub_hg19_final.bam"
HEK_Vs3_KU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A05_sub_hg19_final.bam"
HEK_Vs3_NU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A06_sub_hg19_final.bam"
HEK_Hs4_nD="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A11_sub_hg19_final.bam"
HEK_Hs4_KU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A14_sub_hg19_final.bam"
HEK_Hs4_nD_04h="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A10_sub_hg19_final.bam"
HEK_Hs4_nD_24h="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A12_sub_hg19_final.bam"
HEK_Hs4_KU_04h="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A13_sub_hg19_final.bam"
HEK_Hs4_KU_24h="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A15_sub_hg19_final.bam"
HEK_Vs2_nD="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A17_sub_hg19_final.bam"
HEK_Vs2_KU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A20_sub_hg19_final.bam"
Vs2_gRNA="GACCCCCTCCACCCCGCCTC"
Vs3_gRNA="GGTGAGTGAGTGTGTGCGTG"
Hs4_gRNA="GGCACTGCGGCTGGAGGTGG"


# HEK293T VEGFAs3 replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_nD $HEK_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_nD_r1_c3 "-c=3" \
& \
# HEK293T VEGFAs3 replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_KU $HEK_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_KU_r1_c3 "-c=3" \
& \
# HEK293T VEGFAs3 replicate 1 (w/ NU-7026)
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_NU $HEK_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_NU_r1_c3 "-c=3"


# HEK293T HEKs4 replicate 1 (no drug) 04h
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_nD_04h $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD04_r1_c3 "-c=3" \
& \
# HEK293T HEKs4 replicate 1 (no drug) 12h
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_nD $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD12_r1_c3 "-c=3" \
& \
# HEK293T HEKs4 replicate 1 (no drug) 24h
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_nD_24h $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD24_r1_c3 "-c=3"


# HEK293T HEKs4 replicate 1 (w/ KU-60648) 04h
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_KU_04h $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU04_r1_c3 "-c=3" \
& \
# HEK293T HEKs4 replicate 1 (w/ KU-60648) 12h
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_KU $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU12_r1_c3 "-c=3" \
& \
# HEK293T HEKs4 replicate 1 (w/ KU-60648) 24h
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_KU_24h $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU24_r1_c3 "-c=3"


# HEK293T VEGFAs2 replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $HEK_Vs2_nD $HEK_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_nD_r1_c3 "-c=3" \
& \
# HEK293T VEGFAs2 replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $HEK_Vs2_KU $HEK_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_KU_r1_c3 "-c=3
