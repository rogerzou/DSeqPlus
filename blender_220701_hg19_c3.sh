#!/usr/bin/env bash

cd ../blender_rz
hg19_bowtie2="/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"

##########################################################################
## Run BLENDER on samples with Cas9 or Cpf1 targeting TRA,              ##
## in primary human T-cells using hg19, replicate 1, cutoff threshold 3 ##
##########################################################################

Tcell_Cas9_g_KU="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_g_KU_sub_hg19_merged.bam"
Tcell_Cas9_g_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_g_nD_sub_hg19_merged.bam"
Tcell_Cas9_ng_KU="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_ng_KU_sub_hg19_merged.bam"
Tcell_Cas9_ng_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_ng_nD_sub_hg19_merged.bam"
Tcell_Cpf1_g_KU="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_g_KU_sub_hg19_merged.bam"
Tcell_Cpf1_g_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_g_nD_sub_hg19_merged.bam"
Tcell_Cpf1_ng_KU="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_ng_KU_sub_hg19_merged.bam"
Tcell_Cpf1_ng_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_ng_nD_sub_hg19_merged.bam"
Cas9_gRNA="AGAGTCTCTCAGCTGGTACA"
Cpf1_gRNA="GAGTCTCTCAGCTGGTACAC"
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_g_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_g_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_ng_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_ng_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_g_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_g_nD_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_ng_KU_r1_c3
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_ng_nD_r1_c3

# Tcell Cas9 gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_g_KU $Tcell_Cas9_ng_nD $Cas9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_g_KU_r1_c3 "-c=3" \
& \
# Tcell Cas9 gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_g_nD $Tcell_Cas9_ng_nD $Cas9_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_g_nD_r1_c3 "-c=3" \
& \
# Tcell Cas9 no gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_ng_KU $Tcell_Cas9_ng_nD $Cas9_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_ng_KU_r1_c3 "-c=3" \
& \
# Tcell Cas9 no gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cas9_ng_nD $Tcell_Cas9_ng_KU $Cas9_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cas9_ng_nD_r1_c3 "-c=3"


# Tcell Cpf1 gRNA replicate 1 (w/ KU-60648)
sh run_blender_draw_only.sh $hg19_bowtie2 $Tcell_Cpf1_g_KU $Tcell_Cpf1_ng_nD $Cpf1_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_g_KU_r1_c3 "-c=3 -p=TTT -n=Cpf1" \
& \
# Tcell Cpf1 gRNA replicate 1 (no drug)
sh run_blender_draw_only.sh $hg19_bowtie2 $Tcell_Cpf1_g_nD $Tcell_Cpf1_ng_nD $Cpf1_gRNA \
  /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_g_nD_r1_c3 "-c=3 -p=TTT -n=Cpf1" \
& \
# Tcell Cpf1 no gRNA replicate 1 (w/ KU-60648)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cpf1_ng_KU $Tcell_Cpf1_ng_nD $Cpf1_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_ng_KU_r1_c3 "-c=3 -p=TTT -n=Cpf1" \
& \
# Tcell Cpf1 no gRNA replicate 1 (no drug)
sh run_blender.sh $hg19_bowtie2 $Tcell_Cpf1_ng_nD $Tcell_Cpf1_ng_KU $Cpf1_gRNA \
 /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/Tcell_Cpf1_ng_nD_r1_c3 "-c=3 -p=TTT -n=Cpf1"


##############################################################
## Publicly available DISCOVER-Seq files from Wienert et al ##
##############################################################

public_K562_Cpf1_bfp="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq_public/SRR8550683_hg19_final.bam"
public_K562_Cpf1_08h="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq_public/SRR8550685_hg19_final.bam"
public_K562_Cpf1_24h="/mnt/z/rzou4/NGS_data/4_damage/211213_Dseq_public/SRR8550686_hg19_final.bam"
DNMT1_gRNA_all="CTGATGGTCCATGTCTGTTACTC"
DNMT1_gRNA="CTGATGGTCCATGTCTGTTA"
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/public_K562_Cpf1_08h
mkdir -p /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/public_K562_Cpf1_24h

# public K562 Cpf1 DMNT3 24bp gRNA 8h
sh run_blender.sh $hg19_bowtie2 $public_K562_Cpf1_08h $public_K562_Cpf1_bfp $DNMT1_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/public_K562_Cpf1_08h "-c=3 -p=TTT -n=Cpf1" \
& \
# public K562 Cpf1 DMNT3 24bp gRNA 24h
sh run_blender.sh $hg19_bowtie2 $public_K562_Cpf1_24h $public_K562_Cpf1_bfp $DNMT1_gRNA \
   /mnt/z/rzou4/NGS_data/4_damage/Dseq+/220701_hg19_c3/public_K562_Cpf1_24h "-c=3 -p=TTT -n=Cpf1"
