#!/usr/bin/env bash

cd ../blender_rz
hg19_bowtie2="/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"
mm10_bowtie2="/mnt/c/Users/rzou4/bioinformatics/mm10_bowtie2/mm10.fa"

#########################################
## Run BLENDER on samples without Cas9 ##
## using hg19, cutoff threshold 3      ##
#########################################

HEK_WT_nD_r1="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/HEK_WT_nD_r1_hg19_merged.bam"
HEK_WT_KU_r1="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/HEK_WT_KU_r1_hg19_merged.bam"
HEK_WT_nD_r2="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/HEK_WT_nD_r2_hg19_merged.bam"
HEK_WT_KU_r2="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/HEK_WT_KU_r2_hg19_merged.bam"
K562_WT_nD_r1="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/K562_WT_nD_r1_hg19_merged.bam"
K562_WT_KU_r1="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/K562_WT_KU_r1_hg19_merged.bam"
K562_WT_nD_r2="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/K562_WT_nD_r2_hg19_merged.bam"
K562_WT_KU_r2="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/K562_WT_KU_r2_hg19_merged.bam"
mice_WT_nD_r1="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/mice_WT_nD_r1_hg19_merged.bam"
mice_WT_KU_r1="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/mice_WT_KU_r1_hg19_merged.bam"
mice_WT_nD_r2="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/mice_WT_nD_r2_hg19_merged.bam"
mice_WT_KU_r2="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/mice_WT_KU_r2_hg19_merged.bam"
mice_WT_nD_r3="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/mice_WT_nD_r3_hg19_merged.bam"
mice_WT_KU_r3="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/mice_WT_KU_r3_hg19_merged.bam"
Fs2_gRNA="GCTGCAGAAGGGATTCCATG"
Hs4_gRNA="GGCACTGCGGCTGGAGGTGG"
Vs2_gRNA="GACCCCCTCCACCCCGCCTC"
Vs3_gRNA="GGTGAGTGAGTGTGTGCGTG"
mP9_gRNA="AGCAGCAGCGGCGGCAACAG"
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Fs2_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Fs2_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Hs4_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Hs4_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs2_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs2_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs3_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs3_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Fs2_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Fs2_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Hs4_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Hs4_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs2_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs2_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs3_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs3_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Fs2_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Fs2_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Hs4_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Hs4_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Vs2_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Vs2_KU_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Vs3_nD_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-Vs3_KU_c3


#############
## HEK293T ##
#############

# HEK WT replicate 1 (w/ no drug) - Fs2
sh run_blender.sh $hg19_bowtie2 $HEK_WT_nD_r1 $HEK_WT_nD_r2 $Fs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Fs2_nD_c3 "-c=3" \
& \
# HEK WT replicate 1 (w/ KU-60648) - Fs2
sh run_blender.sh $hg19_bowtie2 $HEK_WT_KU_r1 $HEK_WT_KU_r2 $Fs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Fs2_KU_c3 "-c=3" \
& \
# HEK WT replicate 1 (w/ no drug) - Hs4
sh run_blender.sh $hg19_bowtie2 $HEK_WT_nD_r1 $HEK_WT_nD_r2 $Hs4_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Hs4_nD_c3 "-c=3" \
& \
# HEK WT replicate 1 (w/ KU-60648) - Hs4
sh run_blender.sh $hg19_bowtie2 $HEK_WT_KU_r1 $HEK_WT_KU_r2 $Hs4_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Hs4_KU_c3 "-c=3"

# HEK WT replicate 1 (w/ no drug) - Vs2
sh run_blender.sh $hg19_bowtie2 $HEK_WT_nD_r1 $HEK_WT_nD_r2 $Vs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs2_nD_c3 "-c=3" \
& \
# HEK WT replicate 1 (w/ KU-60648) - Vs2
sh run_blender.sh $hg19_bowtie2 $HEK_WT_KU_r1 $HEK_WT_KU_r2 $Vs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs2_KU_c3 "-c=3" \
& \
# HEK WT replicate 1 (w/ no drug) - Vs3
sh run_blender.sh $hg19_bowtie2 $HEK_WT_nD_r1 $HEK_WT_nD_r2 $Vs3_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs3_nD_c3 "-c=3" \
& \
# HEK WT replicate 1 (w/ KU-60648) - Vs3
sh run_blender.sh $hg19_bowtie2 $HEK_WT_KU_r1 $HEK_WT_KU_r2 $Vs3_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/HEK_WT-Vs3_KU_c3 "-c=3"


##########
## K562 ##
##########

# K562 WT replicate 1 (w/ no drug) - Fs2
sh run_blender.sh $hg19_bowtie2 $K562_WT_nD_r1 $K562_WT_nD_r2 $Fs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Fs2_nD_c3 "-c=3" \
& \
# K562 WT replicate 1 (w/ KU-60648) - Fs2
sh run_blender.sh $hg19_bowtie2 $K562_WT_KU_r1 $K562_WT_KU_r2 $Fs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Fs2_KU_c3 "-c=3" \
& \
# K562 WT replicate 1 (w/ no drug) - Hs4
sh run_blender.sh $hg19_bowtie2 $K562_WT_nD_r1 $K562_WT_nD_r2 $Hs4_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Hs4_nD_c3 "-c=3" \
& \
# K562 WT replicate 1 (w/ KU-60648) - Hs4
sh run_blender.sh $hg19_bowtie2 $K562_WT_KU_r1 $K562_WT_KU_r2 $Hs4_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Hs4_KU_c3 "-c=3"

# K562 WT replicate 1 (w/ no drug) - Vs2
sh run_blender.sh $hg19_bowtie2 $K562_WT_nD_r1 $K562_WT_nD_r2 $Vs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs2_nD_c3 "-c=3" \
& \
# K562 WT replicate 1 (w/ KU-60648) - Vs2
sh run_blender.sh $hg19_bowtie2 $K562_WT_KU_r1 $K562_WT_KU_r2 $Vs2_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs2_KU_c3 "-c=3" \
& \
# K562 WT replicate 1 (w/ no drug) - Vs3
sh run_blender.sh $hg19_bowtie2 $K562_WT_nD_r1 $K562_WT_nD_r2 $Vs3_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs3_nD_c3 "-c=3" \
& \
# K562 WT replicate 1 (w/ KU-60648) - Vs3
sh run_blender.sh $hg19_bowtie2 $K562_WT_KU_r1 $K562_WT_KU_r2 $Vs3_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/K562_WT-Vs3_KU_c3 "-c=3"


##########
## mice ##
##########

# mice WT replicate 1 (w/ no drug) - Fs2
sh run_blender.sh $mm10_bowtie2 $mice_WT_nD_r1 $mice_WT_nD_r2 $mP9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-mP9_nD_c3 "-c=3" \
& \
# mice WT replicate 1 (w/ KU-60648) - Fs2
sh run_blender.sh $mm10_bowtie2 $mice_WT_KU_r1 $mice_WT_KU_r2 $mP9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_cross_c3/mice_WT-mP9_KU_c3 "-c=3"
