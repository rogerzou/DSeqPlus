#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting PCSK9 (mouse)
# in mouse livers using mm10, cutoff threshold 3.

########### USER ENTRY SECTION ###########
cd /path/to/folder/blender
mkdir -p /path/to/output/mouse_PCSK9_nD_merged_c3
mkdir -p /path/to/output/mouse_PCSK9_KU_merged_c3
##########################################

 # Mouse PCSK9 replicate 3 (DISCOVER-Seq, i.e. no drug)
 sh run_blender.sh /home/roger/bioinformatics/mm10/mm10.fa \
    /mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A09_mPCSK9_nd24h_r3_sub_mm10_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/211213_Dseq_public/SRR8553804_mm10_final.bam \
    AGCAGCAGCGGCGGCAACAG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_nD_r3_c3 "-c=3" \
 & \
 # Mouse PCSK9 replicate 3 (DISCOVER-Seq+, i.e. w/ KU-60648)
 sh run_blender.sh /home/roger/bioinformatics/mm10/mm10.fa \
    /mnt/z/rzou4/NGS_data/4_damage/211213_Dseq+/A12_mPCSK9_KU24h_r3_sub_mm10_final.bam \
    /mnt/z/rzou4/NGS_data/4_damage/211213_Dseq_public/SRR8553804_mm10_final.bam \
    AGCAGCAGCGGCGGCAACAG \
    /mnt/z/rzou4/NGS_data/4_damage/Dseq+/mouse_PCSK9_KU_r3_c3 "-c=3"
