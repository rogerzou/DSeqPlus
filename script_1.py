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
mm10_mP9_nD_merge_bam = datadir + "211213_Dseq+/mPCSK9_nD24h_merged.bam"
mm10_mP9_KU_merge_bam = datadir + "211213_Dseq+/mPCSK9_KU24h_merged.bam"
mm10_mP9_nD_merge_txt = datadir + "Dseq+/mouse_PCSK9_nD_merged_c3/filtered_blender_hits.txt"
mm10_mP9_KU_merge_txt = datadir + "Dseq+/mouse_PCSK9_KU_merged_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r2_c2_txt = datadir + "Dseq+/mouse_PCSK9_nD_r2_c2/filtered_blender_hits.txt"
mm10_mP9_KU_r2_c2_txt = datadir + "Dseq+/mouse_PCSK9_KU_r2_c2/filtered_blender_hits.txt"
mm10_mP9_nD_r3_c2_txt = datadir + "Dseq+/mouse_PCSK9_nD_r3_c2/filtered_blender_hits.txt"
mm10_mP9_KU_r3_c2_txt = datadir + "Dseq+/mouse_PCSK9_KU_r3_c2/filtered_blender_hits.txt"

ampngs1 = datadir + "210720_ampNGS/"
ampngs2 = datadir + "210912_ampNGS/"
ampngs3 = datadir + "210912_ampNGS/"
ampngs4 = datadir + "220118_ampNGS/"

""" gRNA sequences """
gRNA_Vs2 = "GACCCCCTCCACCCCGCCTC"
gRNA_Fs2 = "GCTGCAGAAGGGATTCCATG"
gRNA_Vs3 = "GGTGAGTGAGTGTGTGCGTG"
gRNA_Hs4 = "GGCACTGCGGCTGGAGGTGG"
gRNA_mP9 = "AGCAGCAGCGGCGGCAACAG"
gRNA_Fs2_OFF1 = "GCTGCAAAAAGGATTCCAGG"
gRNA_Fs2_OFF2 = "GACGCAGAAGGGACTCCATG"
gRNA_Fs2_OFF3 = "GCTGGAGGAAGGATTCCATG"
gRNA_Fs2_OFF4 = "GGTACAGAAGGGCTTCCATG"
gRNA_Fs2_OFF0 = "GCTGCAGAAGGGATTCCAAG"
gRNA_Fs2_OFF5 = "TGTGAAGAAGGGTTTCCATG"
gRNA_Fs2_OFF6 = "ACTGCAAAATGGATTCCATG"
gRNA_Fs2_OFF7 = "ACTGCTGAGAGGATTCCATG"
gRNA_Fs2_OFF8 = "GCTGCTGAACAGATTCCATG"
gRNA_Fs2_OFF9 = "GCTGAAGAAGGAATTTCATG"
gRNA_Fs2_OFF10 = "GTTGCAGCAAGAATTCCATG"
gRNA_Fs2_OFF11 = "CAAGCAGAAGGGATTCCACA"
gRNA_Fs2_OFF12 = "GCCACAGAAGGGATTCTATG"

""" Set analysis path """
ana = datadir + "Dseq_out/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_minus/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_amplicon/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None


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


""" ############################################################################################ """
""" (Fig. 1l) """
# K562 FANCFs2 no drug - KU60648
print("K562 FANCFs2")
a1 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                                      c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2))]
a2 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                                      c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2))]

# VEGFAs2 no drug - KU60648
print("K562 VEGFAs2")
a3 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2))]
a4 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(K562_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2))]

# HEK VEGFAs3 no drug - KU60648
print("HEK VEGFAs3")

a5 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Vs3_nD_r1_c2_txt, 2000, hg19, gRNA_Vs3),
                                      c.blender_gen(HEK_Vs3_KU_r1_c2_txt, 2000, hg19, gRNA_Vs3))]
a6 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Vs3_KU_r1_c2_txt, 2000, hg19, gRNA_Vs3),
                                      c.blender_gen(HEK_Vs3_nD_r1_c2_txt, 2000, hg19, gRNA_Vs3))]

# HEK HEKs4 no drug - KU60648
print("HEK HEKs4")
a7 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Hs4_nD12_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                                      c.blender_gen(HEK_Hs4_KU12_r1_c2_txt, 2000, hg19, gRNA_Hs4))]
a8 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(HEK_Hs4_KU12_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                                      c.blender_gen(HEK_Hs4_nD12_r1_c2_txt, 2000, hg19, gRNA_Hs4))]


""" ############################################################################################ """
""" nature16525 is Kleinstiver et al., 2016 | nbt3117 is Tsai et al., 2015 """
""" aav9023 is Wienert/Wyman et al., 2020 """
print("K562 FANCFs2 comparison to Kleinstiver et al.")
gen = c.gen_subtr_approx(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                         c.nature16525_gen(gRNA_Fs2))
c.save_gen(gen, ana_3 + "Vs2_KU_DseqP-Gseq")
gen = c.gen_subtr_approx(c.nature16525_gen(gRNA_Fs2),
                         c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2))
c.save_gen(gen, ana_3 + "Vs2_KU_Gseq-DseqP")
gen = c.gen_union_approx(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                         c.nature16525_gen(gRNA_Fs2))
c.save_gen(gen, ana_3 + "Vs2_KU_DseqP-U-Gseq")
gen = c.gen_union_approx(c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                         c.nature16525_gen(gRNA_Fs2))
c.save_gen(gen, ana_3 + "Vs2_nD_DseqP-U-Gseq")

print("K562 VEGFAs2 comparison to Tsai et al.")
b1 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                                       c.nbt3117_gen(gRNA_Vs2), approx=500)]
b2 = [g[0] for g in c.gen_subtr_approx(c.nbt3117_gen(gRNA_Vs2),
                                       c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                                       approx=500)]

print("Mouse PCSK9 pooled nD comparison to KU")
b3 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_merge_txt, 2000, mm10, gRNA_mP9),
                                      c.blender_gen(mm10_mP9_KU_merge_txt, 2000, mm10, gRNA_mP9))]
b4 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_merge_txt, 2000, mm10, gRNA_mP9),
                                      c.blender_gen(mm10_mP9_nD_merge_txt, 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 pooled nD comparison to pooled Wienert/Wyman et al.")
b5 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(mm10_mP9_nD_merge_txt, 2000, mm10, gRNA_mP9),
                                       c.aav9023_gen(gRNA_mP9, choice=0))]
b6 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=0),
                                       c.blender_gen(mm10_mP9_nD_merge_txt, 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 rep3 nD comparison to KU")
c1 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                      c.blender_gen(mm10_mP9_KU_r3_c2_txt, 2000, mm10, gRNA_mP9))]
c2 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                      c.blender_gen(mm10_mP9_nD_r3_c2_txt, 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 rep3 nD comparison to individual Wienert/Wyman et al.")
c3 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(mm10_mP9_nD_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                       c.aav9023_gen(gRNA_mP9, choice=1))]
c4 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                       c.blender_gen(mm10_mP9_nD_r3_c2_txt, 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 rep3 KU comparison to individual Wienert/Wyman et al.")
c5 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(mm10_mP9_KU_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                       c.aav9023_gen(gRNA_mP9, choice=1))]
c6 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                       c.blender_gen(mm10_mP9_KU_r3_c2_txt, 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 individual Wienert/Wyman SUBTRACT rep3 KU UNION rep3 nD")
d1 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                       c.gen_union_approx(c.blender_gen(mm10_mP9_KU_r3_c2_txt,
                                                                        2000, mm10, gRNA_mP9),
                                                          c.blender_gen(mm10_mP9_nD_r3_c2_txt,
                                                                        2000, mm10, gRNA_mP9)))]
print("Mouse PCSK9 rep3 KU SUBTRACT individual Wienert/Wyman UNION rep3 nD")
d2 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(mm10_mP9_KU_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                       c.gen_union_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                                          c.blender_gen(mm10_mP9_nD_r3_c2_txt,
                                                                        2000, mm10, gRNA_mP9)))]
print("Mouse PCSK9 rep3 nD SUBTRACT individual Wienert/Wyman UNION rep3 KU")
d3 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(mm10_mP9_nD_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                       c.gen_union_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                                          c.blender_gen(mm10_mP9_KU_r3_c2_txt,
                                                                        2000, mm10, gRNA_mP9)))]
print("Mouse PCSK9 rep3 nD SUBTRACT individual Wienert/Wyman INTERSECTION rep3 KU")
d4 = [g[0] for g in c.gen_inter_approx(c.blender_gen(mm10_mP9_nD_r3_c2_txt, 2000, mm10, gRNA_mP9),
                                       c.gen_inter_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                                          c.blender_gen(mm10_mP9_KU_r3_c2_txt,
                                                                        2000, mm10, gRNA_mP9)))]


""" ############################################################################################ """
""" Comparison between target sites from replicates pooled vs individual replicates """
gens = [c.blender_gen(mm10_mP9_nD_merge_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_nD_r0_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_nD_r1_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_nD_r2_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_nD_r3_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_nD_r4_c3_txt, 2000, mm10, gRNA_mP9)]
c.union_approx(gens, ana_3 + "union_nD", approx=1000)

gens = [c.blender_gen(mm10_mP9_KU_merge_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_KU_r0_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_KU_r1_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_KU_r2_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, mm10, gRNA_mP9),
        c.blender_gen(mm10_mP9_KU_r4_c3_txt, 2000, mm10, gRNA_mP9)]
c.union_approx(gens, ana_3 + "union_KU", approx=1000)


""" ############################################################################################ """
""" Amplicon NGS of to check indels in cells without DNA-PKi, with/without Cas9 (3 replicates). """
gg.generate_ref(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                ana_4 + "REF_K562_Fs2_KU_r1_c2", hg19[0], genome_savepath)
REF = SeqIO.to_dict(SeqIO.parse(ana_4 + "REF_K562_Fs2_KU_r1_c2.fa", "fasta"))

amp.getIndels(ampngs1 + "merged/Fs2_ON_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct1.txt")
amp.getIndels(ampngs1 + "merged/Fs2_ON_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep1.txt")
amp.getIndels(ampngs2 + "merged/Fs2_ON_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct2.txt")
amp.getIndels(ampngs2 + "merged/Fs2_ON_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep2.txt")
amp.getIndels(ampngs3 + "merged/Fs2_ON_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct3.txt")
amp.getIndels(ampngs3 + "merged/Fs2_ON_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep3.txt")

amp.getIndels(ampngs1 + "merged/Fs2_OFF1_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct1.txt")
amp.getIndels(ampngs1 + "merged/Fs2_OFF1_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep1.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF1_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct2.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF1_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep2.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF1_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct3.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF1_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep3.txt")

amp.getIndels(ampngs1 + "merged/Fs2_OFF2_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct1.txt")
amp.getIndels(ampngs1 + "merged/Fs2_OFF2_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep1.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF2_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct2.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF2_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep2.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF2_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct3.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF2_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep3.txt")

amp.getIndels(ampngs1 + "merged/Fs2_OFF3_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct1.txt")
amp.getIndels(ampngs1 + "merged/Fs2_OFF3_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep1.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF3_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct2.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF3_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep2.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF3_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct3.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF3_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep3.txt")

amp.getIndels(ampngs1 + "merged/Fs2_OFF4_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct1.txt")
amp.getIndels(ampngs1 + "merged/Fs2_OFF4_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep1.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF4_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct2.txt")
amp.getIndels(ampngs2 + "merged/Fs2_OFF4_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep2.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF4_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct3.txt")
amp.getIndels(ampngs3 + "merged/Fs2_OFF4_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ep3.txt")

amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ct1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ep1.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ct2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ep2.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ct3.txt")
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ep3.txt")
