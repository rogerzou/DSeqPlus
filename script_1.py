"""
Script for:
(1)
"""

import src.chipseq as c
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
iPSC_Vs2_nD_r1_bam = datadir + "210301_Dseq+/A13_hg38_final.bam"
iPSC_Vs2_KU_r1_bam = datadir + "210301_Dseq+/A14_hg38_final.bam"
K562_Vs2_nD_r1_bam = datadir + "210301_Dseq+/A15_hg38_final.bam"
K562_Vs2_KU_r1_bam = datadir + "210301_Dseq+/A16_hg38_final.bam"
K562_Fs2_nD_r1_bam = datadir + "210301_Dseq+/A17_hg38_final.bam"
K562_Fs2_KU_r1_bam = datadir + "210301_Dseq+/A18_hg38_final.bam"
iPSC_Vs2_nD_r1_c2_txt = datadir + "Dseq+/iPSC_Vs2_nD_r1_c2/filtered_blender_hits.txt"
iPSC_Vs2_KU_r1_c2_txt = datadir + "Dseq+/iPSC_Vs2_KU_r1_c2/filtered_blender_hits.txt"
K562_Vs2_nD_r1_c2_txt = datadir + "Dseq+/K562_Vs2_nD_r1_c2/filtered_blender_hits.txt"
K562_Vs2_KU_r1_c2_txt = datadir + "Dseq+/K562_Vs2_KU_r1_c2/filtered_blender_hits.txt"
K562_Fs2_nD_r1_c2_txt = datadir + "Dseq+/K562_Fs2_nD_r1_c2/filtered_blender_hits.txt"
K562_Fs2_KU_r1_c2_txt = datadir + "Dseq+/K562_Fs2_KU_r1_c2/filtered_blender_hits.txt"

""" gRNA sequences """
gRNA_Vs2 = "GACCCCCTCCACCCCGCCTC"
gRNA_Fs2 = "GCTGCAGAAGGGATTCCATG"

""" Set analysis path """
ana = datadir + "Dseq_out/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" For all on- and off-target sites from BLENDER, determine paired-end read subsets """
c.read_subsets(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
               iPSC_Vs2_nD_r1_bam, ana_1 + "iPSC_Vs2_nD_r1_c2")
c.read_subsets(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
               iPSC_Vs2_KU_r1_bam, ana_1 + "iPSC_Vs2_KU_r1_c2")
c.read_subsets(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
               K562_Vs2_nD_r1_bam, ana_1 + "K562_Vs2_nD_r1_c2")
c.read_subsets(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
               K562_Vs2_KU_r1_bam, ana_1 + "K562_Vs2_KU_r1_c2")
c.read_subsets(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg38, gRNA_Fs2),
               K562_Fs2_nD_r1_bam, ana_1 + "K562_Fs2_nD_r1_c2")
c.read_subsets(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg38, gRNA_Fs2),
               K562_Fs2_KU_r1_bam, ana_1 + "K562_Fs2_KU_r1_c2")


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from
    Cas9 and MRE11 ChIP-seq. (Fig. 1E-F, S2) """
c.peak_profile_bp_resolution(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
                             iPSC_Vs2_nD_r1_bam, ana_2 + "iPSC_Vs2_nD_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
                             iPSC_Vs2_KU_r1_bam, ana_2 + "iPSC_Vs2_KU_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
                             K562_Vs2_nD_r1_bam, ana_2 + "K562_Vs2_nD_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg38, gRNA_Vs2),
                             K562_Vs2_KU_r1_bam, ana_2 + "K562_Vs2_KU_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg38, gRNA_Fs2),
                             K562_Fs2_nD_r1_bam, ana_2 + "K562_Fs2_nD_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg38, gRNA_Fs2),
                             K562_Fs2_KU_r1_bam, ana_2 + "K562_Fs2_KU_r1_c2")
