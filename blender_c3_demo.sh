#!/usr/bin/env bash

# Run BLENDER on samples with Cas9 targeting PCSK9 (mouse)
# in mouse livers using mm10, cutoff threshold 3.

########### USER ENTRY SECTION ###########
cd /mnt/c/users/Roger/bioinformatics/blender
mm10path="/home/roger/bioinformatics/mm10/mm10.fa"
mkdir -p /mnt/c/Users/rzou4/Downloads/demo/mouse_PCSK9_nD_c3
mkdir -p /mnt/c/Users/rzou4/Downloads/demo/mouse_PCSK9_KU_c3
##########################################

# Mouse PCSK9 replicate 3 (DISCOVER-Seq, i.e. no drug)
sh run_blender.sh $mm10path \
  /mnt/z/rzou4/Downloads/mm_mP9_nD24h_r3_sub_mm10_final.bam \
  /mnt/z/rzou4/Downloads/mm_negctrl_r1_mm10_final.bam \
  AGCAGCAGCGGCGGCAACAG \
  /mnt/c/Users/rzou4/Downloads/demo/mouse_PCSK9_nD_c3 "-c=3" \
& \
# Mouse PCSK9 replicate 3 (DISCOVER-Seq+, i.e. w/ KU-60648)
sh run_blender.sh /home/roger/bioinformatics/mm10/mm10.fa \
  /mnt/z/rzou4/Downloads/mm_mP9_KU24h_r3_sub_mm10_final.bam \
  /mnt/z/rzou4/Downloads/mm_negctrl_r1_mm10_final.bam \
  AGCAGCAGCGGCGGCAACAG \
  /mnt/c/Users/rzou4/Downloads/demo/mouse_PCSK9_KU_c3 "-c=3"
