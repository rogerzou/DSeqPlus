"""
Script for:
(1)
"""
from Bio import SeqIO
import src.amplicon as amp
import src.getgenome as gg
import src.chipseq as c
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    genome_savepath = "/mnt/c/Users/rzou4/Desktop/"         # Directory to store indexed genome
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38/hg38.fa"]
    hg19 = ['hg19', "/mnt/c/Users/rzou4/bioinformatics/hg19/hg19.fa"]
    mm10 = ['mm10', "/mnt/c/Users/rzou4/bioinformatics/mm10/mm10.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    genome_savepath = "/Users/rogerzou/Desktop/"            # Directory to store indexed genome
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"]
    mm10 = ['mm10', "/Users/rogerzou/bioinformatics/mm10_bowtie2/mm10.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
HEK_Vs3_nD_r1_c3_txt = datadir + "Dseq+/HEK_Vs3_nD_r1_c3/filtered_blender_hits.txt"
HEK_Vs3_KU_r1_c3_txt = datadir + "Dseq+/HEK_Vs3_KU_r1_c3/filtered_blender_hits.txt"
HEK_Hs4_nD12_r1_c3_txt = datadir + "Dseq+/HEK_Hs4_nD12_r1_c3/filtered_blender_hits.txt"
HEK_Hs4_KU12_r1_c3_txt = datadir + "Dseq+/HEK_Hs4_KU12_r1_c3/filtered_blender_hits.txt"
K562_Vs2_nD_r1_c3_txt = datadir + "Dseq+/K562_Vs2_nD_r1_c3/filtered_blender_hits.txt"
K562_Vs2_KU_r1_c3_txt = datadir + "Dseq+/K562_Vs2_KU_r1_c3/filtered_blender_hits.txt"
K562_Fs2_nD_r1_c3_txt = datadir + "Dseq+/K562_Fs2_nD_r1_c3/filtered_blender_hits.txt"
K562_Fs2_KU_r1_c3_txt = datadir + "Dseq+/K562_Fs2_KU_r1_c3/filtered_blender_hits.txt"
Vs2g_Fs2s_nD_c3_txt = datadir + "Dseq+/Vs2g-Fs2s_nD_r1_c3/filtered_blender_hits.txt"
Vs2g_Fs2s_KU_c3_txt = datadir + "Dseq+/Vs2g-Fs2s_KU_r1_c3/filtered_blender_hits.txt"
Vs2g_Hs4s_nD_c3_txt = datadir + "Dseq+/Vs2g-Hs4s_nD_r1_c3/filtered_blender_hits.txt"
Vs2g_Hs4s_KU_c3_txt = datadir + "Dseq+/Vs2g-Hs4s_KU_r1_c3/filtered_blender_hits.txt"
Vs2g_Vs3s_nD_c3_txt = datadir + "Dseq+/Vs2g-Vs3s_nD_r1_c3/filtered_blender_hits.txt"
Vs2g_Vs3s_KU_c3_txt = datadir + "Dseq+/Vs2g-Vs3s_KU_r1_c3/filtered_blender_hits.txt"

""" gRNA sequences """
gRNA_Vs2 = "GACCCCCTCCACCCCGCCTC"
gRNA_Fs2 = "GCTGCAGAAGGGATTCCATG"
gRNA_Vs3 = "GGTGAGTGAGTGTGTGCGTG"
gRNA_Hs4 = "GGCACTGCGGCTGGAGGTGG"

""" Set analysis path """
ana = datadir + "Dseq_out/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "6_cross/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None


""" ############################################################################################ """
""" (Fig. 1l) """
# K562 FANCFs2 no drug - KU60648
print("K562 FANCFs2")
a1 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Fs2_nD_r1_c3_txt, 2000, hg19, gRNA_Fs2),
                                      c.blender_gen(K562_Fs2_KU_r1_c3_txt, 2000, hg19, gRNA_Fs2))]
a2 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Fs2_KU_r1_c3_txt, 2000, hg19, gRNA_Fs2),
                                      c.blender_gen(K562_Fs2_nD_r1_c3_txt, 2000, hg19, gRNA_Fs2))]

# VEGFAs2 no drug - KU60648
print("K562 VEGFAs2")
a3 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Vs2_nD_r1_c3_txt, 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(K562_Vs2_KU_r1_c3_txt, 2000, hg19, gRNA_Vs2))]
a4 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Vs2_KU_r1_c3_txt, 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(K562_Vs2_nD_r1_c3_txt, 2000, hg19, gRNA_Vs2))]

# HEK VEGFAs3 no drug - KU60648
print("HEK VEGFAs3")

a5 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Vs3_nD_r1_c3_txt, 2000, hg19, gRNA_Vs3),
                                      c.blender_gen(HEK_Vs3_KU_r1_c3_txt, 2000, hg19, gRNA_Vs3))]
a6 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Vs3_KU_r1_c3_txt, 2000, hg19, gRNA_Vs3),
                                      c.blender_gen(HEK_Vs3_nD_r1_c3_txt, 2000, hg19, gRNA_Vs3))]

# HEK HEKs4 no drug - KU60648
print("HEK HEKs4")
a7 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Hs4_nD12_r1_c3_txt, 2000, hg19, gRNA_Hs4),
                                      c.blender_gen(HEK_Hs4_KU12_r1_c3_txt, 2000, hg19, gRNA_Hs4))]
a8 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Hs4_KU12_r1_c3_txt, 2000, hg19, gRNA_Hs4),
                                      c.blender_gen(HEK_Hs4_nD12_r1_c3_txt, 2000, hg19, gRNA_Hs4))]
