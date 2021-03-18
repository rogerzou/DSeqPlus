#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting VEGFAs2, VEGFAs3, and HEKs4
# in HEK293T cells using hg19, cutoff threshold 2.

cd /mnt/c/users/Roger/bioinformatics/blender
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_nD_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_KU_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_NU_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD04_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD12_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD24_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU04_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU12_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU24_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_nD_r1_c2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_KU_r1_c2

# HEK293T VEGFAs3 replicate 1 (no drug)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A04_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGTGAGTGAGTGTGTGCGTG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_nD_r1_c2 "-c=2" \
& \
# HEK293T VEGFAs3 replicate 1 (w/ KU-60648)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A05_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGTGAGTGAGTGTGTGCGTG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_KU_r1_c2 "-c=2" \
& \
# HEK293T VEGFAs3 replicate 1 (w/ NU-7026)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A06_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGTGAGTGAGTGTGTGCGTG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs3_NU_r1_c2 "-c=2"


# HEK293T HEKs4 replicate 1 (no drug) 04h
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A10_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGCACTGCGGCTGGAGGTGG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD04_r1_c2 "-c=2" \
& \
# HEK293T HEKs4 replicate 1 (no drug) 12h
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A11_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGCACTGCGGCTGGAGGTGG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD12_r1_c2 "-c=2" \
& \
# HEK293T HEKs4 replicate 1 (no drug) 24h
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A12_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGCACTGCGGCTGGAGGTGG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_nD24_r1_c2 "-c=2"


# HEK293T HEKs4 replicate 1 (w/ KU-60648) 04h
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A13_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGCACTGCGGCTGGAGGTGG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU04_r1_c2 "-c=2" \
& \
# HEK293T HEKs4 replicate 1 (w/ KU-60648) 12h
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A14_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGCACTGCGGCTGGAGGTGG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU12_r1_c2 "-c=2" \
& \
# HEK293T HEKs4 replicate 1 (w/ KU-60648) 24h
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A15_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GGCACTGCGGCTGGAGGTGG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Hs4_KU24_r1_c2 "-c=2"


# HEK293T VEGFAs2 replicate 1 (no drug)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A17_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GACCCCCTCCACCCCGCCTC \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_nD_r1_c2 "-c=2" \
& \
# HEK293T VEGFAs2 replicate 1 (w/ KU-60648)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A20_sub_hg19_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam \
    GACCCCCTCCACCCCGCCTC \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/HEK_Vs2_KU_r1_c2 "-c=2
