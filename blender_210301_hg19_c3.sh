#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting VEGFAs2 and FANCFs2,
# in K562 and/or iPSC cells using hg19, cutoff threshold 3.

cd /home/roger/bioinformatics/blender
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_KU_r1_c3
hg19_bowtie2="/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"
K562_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg19_final.bam"
K562_Vs2_nD="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A15_sub_hg19_final.bam"
K562_Vs2_KU="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A16_sub_hg19_final.bam"
K562_Fs2_nD="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A17_sub_hg19_final.bam"
K562_Fs2_KU="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A18_sub_hg19_final.bam"
HEK_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam"
iPSC_Vs2_nD="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A13_sub_hg19_final.bam"
iPSC_Vs2_KU="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A14_sub_hg19_final.bam"
Vs2_gRNA="GACCCCCTCCACCCCGCCTC"
Fs2_gRNA="GCTGCAGAAGGGATTCCATG"

# iPSC VEGFAs2 replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $iPSC_Vs2_nD $HEK_ctrl $Vs2_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_nD_r1_c3 "-c=3" \
& \
# iPSC VEGFAs2 replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 iPSC_Vs2_KU $HEK_ctrl $Vs2_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_KU_r1_c3 "-c=3" \
& \
# K562 VEGFAs2 replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_nD $K562_ctrl $Vs2_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_nD_r1_c3 "-c=3" \
& \
# K562 VEGFAs2 replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_KU $K562_ctrl $Vs2_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_KU_r1_c3 "-c=3" \
& \
# K562 FANCFs2 replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_nD $K562_ctrl $Fs2_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_nD_r1_c3 "-c=3" \
& \
# K562 FANCFs2 replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_KU $K562_ctrl $Fs2_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_KU_r1_c3 "-c=3"
