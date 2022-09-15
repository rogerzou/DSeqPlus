"""
Script for:
(1) Calculating MRE11 enrichment at all putative off-target sites
"""
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
HEK_Vs3_nD_r1_bam = datadir + "210216_Dseq+/A04_sub_hg19_final.bam"
HEK_Vs3_KU_r1_bam = datadir + "210216_Dseq+/A05_sub_hg19_final.bam"
HEK_Hs4_nD12_r1_bam = datadir + "210216_Dseq+/A11_sub_hg19_final.bam"
HEK_Hs4_KU12_r1_bam = datadir + "210216_Dseq+/A14_sub_hg19_final.bam"
iPSC_Vs2_nD_r1_bam = datadir + "210301_Dseq+/A13_sub_hg19_final.bam"
iPSC_Vs2_KU_r1_bam = datadir + "210301_Dseq+/A14_sub_hg19_final.bam"
K562_Vs2_nD_r1_bam = datadir + "210301_Dseq+/A15_sub_hg19_final.bam"
K562_Vs2_KU_r1_bam = datadir + "210301_Dseq+/A16_sub_hg19_final.bam"
K562_Fs2_nD_r1_bam = datadir + "210301_Dseq+/A17_sub_hg19_final.bam"
K562_Fs2_KU_r1_bam = datadir + "210301_Dseq+/A18_sub_hg19_final.bam"
HEK_Vs3_nD_r1_c2_txt = datadir + "Dseq+/HEK_Vs3_nD_r1_c2/filtered_blender_hits.txt"
HEK_Vs3_KU_r1_c2_txt = datadir + "Dseq+/HEK_Vs3_KU_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_nD12_r1_c2_txt = datadir + "Dseq+/HEK_Hs4_nD12_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_KU12_r1_c2_txt = datadir + "Dseq+/HEK_Hs4_KU12_r1_c2/filtered_blender_hits.txt"
iPSC_Vs2_nD_r1_c2_txt = datadir + "Dseq+/iPSC_Vs2_nD_r1_c2/filtered_blender_hits.txt"
iPSC_Vs2_KU_r1_c2_txt = datadir + "Dseq+/iPSC_Vs2_KU_r1_c2/filtered_blender_hits.txt"
K562_Vs2_nD_r1_c2_txt = datadir + "Dseq+/K562_Vs2_nD_r1_c2/filtered_blender_hits.txt"
K562_Vs2_KU_r1_c2_txt = datadir + "Dseq+/K562_Vs2_KU_r1_c2/filtered_blender_hits.txt"
K562_Fs2_nD_r1_c2_txt = datadir + "Dseq+/K562_Fs2_nD_r1_c2/filtered_blender_hits.txt"
K562_Fs2_KU_r1_c2_txt = datadir + "Dseq+/K562_Fs2_KU_r1_c2/filtered_blender_hits.txt"
mm10_mP9_nD_r0_bam = datadir + "211213_Dseq+/A02_mPCSK9_nd24h_r0_sub_mm10_final.bam"
mm10_mP9_nD_r1_bam = datadir + "211213_Dseq+/A07_mPCSK9_nd24h_r1_sub_mm10_final.bam"
mm10_mP9_nD_r2_bam = datadir + "211213_Dseq+/A08_mPCSK9_nd24h_r2_sub_mm10_final.bam"
mm10_mP9_nD_r3_bam = datadir + "211213_Dseq+/A09_mPCSK9_nd24h_r3_sub_mm10_final.bam"
mm10_mP9_nD_r4_bam = datadir + "211213_Dseq+/A05_mPCSK9_nd24h_r4_sub_mm10_final.bam"
mm10_mP9_KU_r0_bam = datadir + "211213_Dseq+/A03_mPCSK9_KU24h_r0_sub_mm10_final.bam"
mm10_mP9_KU_r1_bam = datadir + "211213_Dseq+/A10_mPCSK9_KU24h_r1_sub_mm10_final.bam"
mm10_mP9_KU_r2_bam = datadir + "211213_Dseq+/A11_mPCSK9_KU24h_r2_sub_mm10_final.bam"
mm10_mP9_KU_r3_bam = datadir + "211213_Dseq+/A12_mPCSK9_KU24h_r3_sub_mm10_final.bam"
mm10_mP9_KU_r4_bam = datadir + "211213_Dseq+/A06_mPCSK9_KU24h_r4_sub_mm10_final.bam"
Tcel_Cas9_ng_nD_r1_bam = datadir + "220701_Dseq+/Tcell_Cas9_ng_nD_r1_sub_hg19_merged.bam"
Tcel_Cas9_ng_KU_r1_bam = datadir + "220701_Dseq+/Tcell_Cas9_ng_KU_r1_sub_hg19_merged.bam"
Tcel_Cas9_g_nD_r1_bam = datadir + "220701_Dseq+/Tcell_Cas9_g_nD_r1_sub_hg19_merged.bam"
Tcel_Cas9_g_KU_r1_bam = datadir + "220701_Dseq+/Tcell_Cas9_g_KU_r1_sub_hg19_merged.bam"
mm10_mP9_nD_r0_c3_txt = datadir + "Dseq+/mouse_PCSK9_nD_r0_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r1_c3_txt = datadir + "Dseq+/mouse_PCSK9_nD_r1_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r2_c3_txt = datadir + "Dseq+/mouse_PCSK9_nD_r2_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r3_c3_txt = datadir + "Dseq+/mouse_PCSK9_nD_r3_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r4_c3_txt = datadir + "Dseq+/mouse_PCSK9_nD_r4_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r0_c3_txt = datadir + "Dseq+/mouse_PCSK9_KU_r0_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r1_c3_txt = datadir + "Dseq+/mouse_PCSK9_KU_r1_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r2_c3_txt = datadir + "Dseq+/mouse_PCSK9_KU_r2_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r3_c3_txt = datadir + "Dseq+/mouse_PCSK9_KU_r3_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r4_c3_txt = datadir + "Dseq+/mouse_PCSK9_KU_r4_c3/filtered_blender_hits.txt"
Tcel_Cas9_ng_nD_r1_txt = datadir + "Dseq+/Tcell_Cas9_ng_nD_r1_c3/filtered_blender_hits.txt"
Tcel_Cas9_ng_KU_r1_txt = datadir + "Dseq+/Tcell_Cas9_ng_KU_r1_c3/filtered_blender_hits.txt"
Tcel_Cas9_g_nD_r1_txt = datadir + "Dseq+/Tcell_Cas9_g_nD_r1_c3/filtered_blender_hits.txt"
Tcel_Cas9_g_KU_r1_txt = datadir + "Dseq+/Tcell_Cas9_g_KU_r1_c3/filtered_blender_hits.txt"

""" gRNA sequences """
gRNA_Vs2 = "GACCCCCTCCACCCCGCCTC"
gRNA_Fs2 = "GCTGCAGAAGGGATTCCATG"
gRNA_Vs3 = "GGTGAGTGAGTGTGTGCGTG"
gRNA_Hs4 = "GGCACTGCGGCTGGAGGTGG"
gRNA_mP9 = "AGCAGCAGCGGCGGCAACAG"
gRNA_Cas9 = "AGAGTCTCTCAGCTGGTACA"
gRNA_Cpf1 = "GAGTCTCTCAGCTGGTACAC"

""" Set analysis path """
ana = datadir + "Dseq_out/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" For all on- and off-target sites from BLENDER, determine paired-end read subsets.
    (analysis1-2a) """
c.read_subsets(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
               iPSC_Vs2_nD_r1_bam, ana_1 + "iPSC_Vs2_nD_r1_c2")
c.read_subsets(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
               iPSC_Vs2_KU_r1_bam, ana_1 + "iPSC_Vs2_KU_r1_c2")
c.read_subsets(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
               K562_Vs2_nD_r1_bam, ana_1 + "K562_Vs2_nD_r1_c2")
c.read_subsets(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
               K562_Vs2_KU_r1_bam, ana_1 + "K562_Vs2_KU_r1_c2")
c.read_subsets(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
               K562_Fs2_nD_r1_bam, ana_1 + "K562_Fs2_nD_r1_c2")
c.read_subsets(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
               K562_Fs2_KU_r1_bam, ana_1 + "K562_Fs2_KU_r1_c2")

c.read_subsets(c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, mm10, gRNA_mP9),
               mm10_mP9_nD_r3_bam, ana_1 + "mm10_mP9_nD_r3_c3")
c.read_subsets(c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, mm10, gRNA_mP9),
               mm10_mP9_KU_r3_bam, ana_1 + "mm10_mP9_KU_r3_c3")

c.read_subsets(c.blender_gen(Tcel_Cas9_g_KU_r1_txt, 2000, hg19, gRNA_Cas9),
               Tcel_Cas9_g_nD_r1_bam, ana_1 + "Tcel_Cas9_g_nD_r1_c3")
c.read_subsets(c.blender_gen(Tcel_Cas9_g_KU_r1_txt, 2000, hg19, gRNA_Cas9),
               Tcel_Cas9_g_KU_r1_bam, ana_1 + "Tcel_Cas9_g_KU_r1_c3")
c.read_subsets(c.blender_gen(Tcel_Cas9_g_KU_r1_txt, 2000, hg19, gRNA_Cas9),
               Tcel_Cas9_ng_nD_r1_bam, ana_1 + "Tcel_Cas9_ng_nD_r1_c3")
c.read_subsets(c.blender_gen(Tcel_Cas9_g_KU_r1_txt, 2000, hg19, gRNA_Cas9),
               Tcel_Cas9_ng_KU_r1_bam, ana_1 + "Tcel_Cas9_ng_KU_r1_c3")


""" ############################################################################################ """
""" For all on- and off-target sites from BLENDER, determine paired-end read subsets at a location
    10kb away (to show that DNA-PKcs inhibition does not globally increase background).
    (analysis2b) """
c.read_subsets(c.shift_gen(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2)),
               iPSC_Vs2_nD_r1_bam, ana_1 + "Shift_iPSC_Vs2_nD_r1_c2")
c.read_subsets(c.shift_gen(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2)),
               iPSC_Vs2_KU_r1_bam, ana_1 + "Shift_iPSC_Vs2_KU_r1_c2")
c.read_subsets(c.shift_gen(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2)),
               K562_Vs2_nD_r1_bam, ana_1 + "Shift_K562_Vs2_nD_r1_c2")
c.read_subsets(c.shift_gen(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2)),
               K562_Vs2_KU_r1_bam, ana_1 + "Shift_K562_Vs2_KU_r1_c2")
c.read_subsets(c.shift_gen(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2)),
               K562_Fs2_nD_r1_bam, ana_1 + "Shift_K562_Fs2_nD_r1_c2")
c.read_subsets(c.shift_gen(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2)),
               K562_Fs2_KU_r1_bam, ana_1 + "Shift_K562_Fs2_KU_r1_c2")


""" ############################################################################################ """
""" Determine peak profiles at all discovered target sites. """
c.peak_profile_bp_resolution(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                             iPSC_Vs2_nD_r1_bam, ana_2 + "iPSC_Vs2_nD_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                             iPSC_Vs2_KU_r1_bam, ana_2 + "iPSC_Vs2_KU_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                             K562_Vs2_nD_r1_bam, ana_2 + "K562_Vs2_nD_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                             K562_Vs2_KU_r1_bam, ana_2 + "K562_Vs2_KU_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                             K562_Fs2_nD_r1_bam, ana_2 + "K562_Fs2_nD_r1_c2")
c.peak_profile_bp_resolution(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                             K562_Fs2_KU_r1_bam, ana_2 + "K562_Fs2_KU_r1_c2")
