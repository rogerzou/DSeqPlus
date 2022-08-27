#!/usr/bin/env bash

# Use macs2 to find all ChIP-seq peaks for wild type samples aligned to hg19.
# After finding the peaks, I go through one by one, making sure each region is
outf="/mnt/z/rzou4/NGS_data/4_damage/Dseq_macs/"

HEK_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A19_hg19_final.bam"
K562_ctrl="/mnt/z/rzou4/NGS_data/4_damage/210301_Dseq+/A20_hg19_final.bam"
Tcell_Cas9_ng_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_ng_nD_sub_hg19_merged.bam"
Tcell_Cpf1_ng_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_ng_nD_sub_hg19_merged.bam"

 macs2 callpeak -t $HEK_ctrl --outdir $outf --name "HEK_ctrl" -f BAMPE -g hs \
 & \
 macs2 callpeak -t $K562_ctrl --outdir $outf --name "K562_ctrl" -f BAMPE -g hs \
 & \
 macs2 callpeak -t $Tcell_Cas9_ng_nD --outdir $outf --name "Tcell_Cas9_ctrl" -f BAMPE -g hs \
 & \
 macs2 callpeak -t $Tcell_Cpf1_ng_nD --outdir $outf --name "Tcell_Cpf1_ctrl" -f BAMPE -g hs


Tcell_Cas9_g_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_g_nD_sub_hg19_merged.bam"
Tcell_Cpf1_g_nD="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_g_nD_sub_hg19_merged.bam"
Tcell_Cas9_g_KU="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cas9_g_KU_sub_hg19_merged.bam"
Tcell_Cpf1_g_KU="/mnt/z/rzou4/NGS_data/4_damage/220701_Dseq+/Tcell_Cpf1_g_KU_sub_hg19_merged.bam"

macs2 callpeak -t $Tcell_Cas9_g_nD -c $Tcell_Cas9_ng_nD --outdir $outf --name "Tcell_Cas9_g_nD" -f BAMPE -g hs \
& \
macs2 callpeak -t $Tcell_Cpf1_g_nD -c $Tcell_Cpf1_ng_nD --outdir $outf --name "Tcell_Cpf1_g_nD" -f BAMPE -g hs \
& \
macs2 callpeak -t $Tcell_Cas9_g_KU -c $Tcell_Cas9_ng_nD --outdir $outf --name "Tcell_Cas9_g_KU" -f BAMPE -g hs \
& \
macs2 callpeak -t $Tcell_Cpf1_g_KU -c $Tcell_Cpf1_ng_nD --outdir $outf --name "Tcell_Cpf1_g_KU" -f BAMPE -g hs
