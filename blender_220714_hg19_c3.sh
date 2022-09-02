#!/usr/bin/env bash

cd ../blender_rz
hg19_bowtie2="/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"

##########################################################################
## Run BLENDER on samples with Cas9 or Cpf1 targeting TRA,              ##
## in primary human T-cells using hg19, replicate 2, cutoff threshold 3 ##
##########################################################################

Tcell_Cas9_g_KU="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cas9_g_KU_r2_sub_hg19_final.bam"
Tcell_Cas9_g_nD="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cas9_g_nD_r2_sub_hg19_final.bam"
Tcell_Cas9_ng_KU="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cas9_ng_KU_r2_sub_hg19_final.bam"
Tcell_Cas9_ng_nD="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cas9_ng_nD_r2_sub_hg19_final.bam"
Tcell_Cpf1_g_KU="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cpf1_g_KU_r2_sub_hg19_final.bam"
Tcell_Cpf1_g_nD="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cpf1_g_nD_r2_sub_hg19_final.bam"
Tcell_Cpf1_ng_KU="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cpf1_ng_KU_r2_sub_hg19_final.bam"
Tcell_Cpf1_ng_nD="/mnt/z/rzou4/NGS_data/4_damage/220714_Dseq+/Tcell_Cpf1_ng_nD_r2_sub_hg19_final.bam"
Cas9_gRNA="AGAGTCTCTCAGCTGGTACA"
Cpf1_gRNA="GAGTCTCTCAGCTGGTACAC"
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_g_KU_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_g_nD_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_ng_KU_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_ng_nD_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_g_KU_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_g_nD_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_ng_KU_r2
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_ng_nD_r2

# Tcell Cas9 gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_g_KU $Tcell_Cas9_ng_nD $Cas9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_g_KU_r2 "-c=3" \
& \
# Tcell Cas9 gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_g_nD $Tcell_Cas9_ng_nD $Cas9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_g_nD_r2 "-c=3" \
& \
# Tcell Cas9 no gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_ng_KU $Tcell_Cas9_ng_nD $Cas9_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_ng_KU_r2 "-c=3" \
& \
# Tcell Cas9 no gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_ng_nD $Tcell_Cas9_ng_KU $Cas9_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cas9_ng_nD_r2 "-c=3"


# Tcell Cpf1 gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cpf1_g_KU $Tcell_Cpf1_ng_nD $Cpf1_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_g_KU_r2 "-c=3 -p=TTT -n=Cpf1" \
& \
# Tcell Cpf1 gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cpf1_g_nD $Tcell_Cpf1_ng_nD $Cpf1_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_g_nD_r2 "-c=3 -p=TTT -n=Cpf1" \
& \
# Tcell Cpf1 no gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cpf1_ng_KU $Tcell_Cpf1_ng_nD $Cpf1_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_ng_KU_r2 "-c=3 -p=TTT -n=Cpf1" \
& \
# Tcell Cpf1 no gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cpf1_ng_nD $Tcell_Cpf1_ng_KU $Cpf1_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220714_hg19_c3/Tcell_Cpf1_ng_nD_r2 "-c=3 -p=TTT -n=Cpf1"
