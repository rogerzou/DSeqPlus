#!/usr/bin/env bash

cd ../blender_rz

################################################################
## Generate updated BLENDER target site list after removal of ##
## false positives and duplicates.                            ##
################################################################

bdir="/mnt/z/rzou4/NGS_data/4_damage/Dseq_out/3_minus/"

Vs2_gRNA="GACCCCCTCCACCCCGCCTC"
Vs3_gRNA="GGTGAGTGAGTGTGTGCGTG"
P9_gRNA="AGCAGCAGCGGCGGCAACAG"

sh draw_blender.sh $Vs3_gRNA "${bdir}/HEK_Vs3_nD_final.txt" "${bdir}/HEK_Vs3_nD_final"
sh draw_blender.sh $Vs3_gRNA "${bdir}/HEK_Vs3_KU_final.txt" "${bdir}/HEK_Vs3_KU_final"
sh draw_blender.sh $Vs2_gRNA "${bdir}/iPSC_Vs2_nD_final.txt" "${bdir}/iPSC_Vs2_nD_final"
sh draw_blender.sh $Vs2_gRNA "${bdir}/iPSC_Vs2_KU_final.txt" "${bdir}/iPSC_Vs2_KU_final"
sh draw_blender.sh $Vs2_gRNA "${bdir}/K562_Vs2_nD_final.txt" "${bdir}/K562_Vs2_nD_final"
sh draw_blender.sh $Vs2_gRNA "${bdir}/K562_Vs2_KU_final.txt" "${bdir}/K562_Vs2_KU_final"
sh draw_blender.sh $P9_gRNA "${bdir}/mice_mP9_KU_final_merge_rmdup.txt" "${bdir}/mice_mP9_KU_final_merge_rmdup"
sh draw_blender.sh $P9_gRNA "${bdir}/mice_mP9_KU_final_merge_rmdup.txt" "${bdir}/mice_mP9_KU_final_merge_rmdup"
