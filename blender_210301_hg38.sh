#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting VEGFAs2 and FANCFs2,
# in K562 and/or iPSC cells using hg38

cd /mnt/c/users/Roger/bioinformatics/blender
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_nD_r1
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_KU_r1
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_nD_r1
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_KU_r1
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_nD_r1
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_KU_r1

# iPSC VEGFAs2 replicate 1 (no drug)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A13_hg38_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg38_final.bam \
    GACCCCCTCCACCCCGCCTC \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_nD_r1 \
& \
# iPSC VEGFAs2 replicate 1 (w/ KU-60648)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A14_hg38_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg38_final.bam \
    GACCCCCTCCACCCCGCCTC \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_KU_r1 \
& \
# K562 VEGFAs2 replicate 1 (no drug)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A15_hg38_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg38_final.bam \
    GACCCCCTCCACCCCGCCTC \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_nD_r1 \
& \
# K562 VEGFAs2 replicate 1 (w/ KU-60648)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A16_hg38_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg38_final.bam \
    GACCCCCTCCACCCCGCCTC \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_KU_r1 \
& \
# K562 FANCFs2 replicate 1 (no drug)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A17_hg38_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg38_final.bam \
    GCTGCAGAAGGGATTCCATG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_nD_r1 \
& \
# K562 FANCFs2  replicate 1 (w/ KU-60648)
sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A18_hg38_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg38_final.bam \
    GCTGCAGAAGGGATTCCATG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_KU_r1
