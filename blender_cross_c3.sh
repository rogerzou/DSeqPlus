#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting VEGFAs2 and FANCFs2,
# in K562 and/or iPSC cells using hg19, cutoff threshold 3.

cd /home/roger/bioinformatics/blender
hg19_bowtie2="/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"
K562_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg19_final.bam"
K562_Vs2_nD="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A15_sub_hg19_final.bam"
K562_Vs2_KU="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A16_sub_hg19_final.bam"
K562_Fs2_nD="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A17_sub_hg19_final.bam"
K562_Fs2_KU="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A18_sub_hg19_final.bam"
HEK_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam"
HEK_Vs3_nD="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A04_sub_hg19_final.bam"
HEK_Vs3_KU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A05_sub_hg19_final.bam"
HEK_Vs3_NU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A06_sub_hg19_final.bam"
HEK_Hs4_nD="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A11_sub_hg19_final.bam"
HEK_Hs4_KU="/mnt/z/rzou4/NGS_data/4_damage/210216_Dseq+/A14_sub_hg19_final.bam"
Vs2_gRNA="GACCCCCTCCACCCCGCCTC"
Fs2_gRNA="GCTGCAGAAGGGATTCCATG"
Vs3_gRNA="GGTGAGTGAGTGTGTGCGTG"
Hs4_gRNA="GGCACTGCGGCTGGAGGTGG"

# FANCF site 2 off-target search
  # Cas9 targeting VEGFA site 2
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_nD $K562_ctrl $Fs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Fs2g-Vs2s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_KU $K562_ctrl $Fs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Fs2g-Vs2s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting HEK site 4
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_nD $HEK_ctrl $Fs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Fs2g-Hs4s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_KU $HEK_ctrl $Fs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Fs2g-Hs4s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting VEGFA site 3
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_nD $HEK_ctrl $Fs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Fs2g-Vs3s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_KU $HEK_ctrl $Fs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Fs2g-Vs3s_KU_r1_c3 "-c=3"


# VEGFA site 2 off-target search
  # Cas9 targeting FANCF site 2
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_nD $K562_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs2g-Fs2s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_KU $K562_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs2g-Fs2s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting HEK site 4
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_nD $HEK_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs2g-Hs4s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_KU $HEK_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs2g-Hs4s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting VEGFA site 3
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_nD $HEK_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs2g-Vs3s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_KU $HEK_ctrl $Vs2_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs2g-Vs3s_KU_r1_c3 "-c=3"


# HEK site 4 off-target search
  # Cas9 targeting VEGFA site 2
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_nD $K562_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Hs4g-Vs2s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_KU $K562_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Hs4g-Vs2s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting FANCF site 2
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_nD $K562_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Hs4g-Fs2s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_KU $K562_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Hs4g-Fs2s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting VEGFA site 3
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_nD $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Hs4g-Vs3s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $HEK_Vs3_KU $HEK_ctrl $Hs4_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Hs4g-Vs3s_KU_r1_c3 "-c=3"


# VEGFA site 3 off-target search
  # Cas9 targeting FANCF site 2
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_nD $K562_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs3g-Fs2s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $K562_Fs2_KU $K562_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs3g-Fs2s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting HEK site 4
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_nD $HEK_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs3g-Hs4s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $HEK_Hs4_KU $HEK_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs3g-Hs4s_KU_r1_c3 "-c=3" \
& \
  # Cas9 targeting VEGFA site 2
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_nD $K562_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs3g-Vs2s_nD_r1_c3 "-c=3" \
& \
sh run_blender.sh $hg19_bowtie2 $K562_Vs2_KU $K562_ctrl $Vs3_gRNA \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/Vs3g-Vs2s_KU_r1_c3 "-c=3"
