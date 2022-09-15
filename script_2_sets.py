"""
Script for:
(1) Subtracting false positive sites from putative off-target sites
(2) Finding overlap between DISCOVER-Seq, DISCOVER-Seq+, GUIDE-seq, etc.
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
# HEK and iPSC with neg ctrls
HEK_Hs4_nD04_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Hs4_nD04_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_KU04_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Hs4_KU04_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_nD12_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Hs4_nD12_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_KU12_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Hs4_KU12_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_nD24_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Hs4_nD24_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_KU24_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Hs4_KU24_r1_c2/filtered_blender_hits.txt"
HEK_Vs3_nD_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Vs3_nD_r1_c2/filtered_blender_hits.txt"
HEK_Vs3_KU_r1_c2_txt = datadir + "Dseq+/210216_hg19_c2/HEK_Vs3_KU_r1_c2/filtered_blender_hits.txt"
iPSC_Vs2_nD_r1_c2_txt = datadir + "Dseq+/210301_hg19_c2/iPSC_Vs2_nD_r1_c2/filtered_blender_hits.txt"
iPSC_Vs2_KU_r1_c2_txt = datadir + "Dseq+/210301_hg19_c2/iPSC_Vs2_KU_r1_c2/filtered_blender_hits.txt"
HEK_Hs4_Kn_c3_txt = datadir + "Dseq+/220701_cross_c3/HEK_WT-Hs4_KU-nD_c3/filtered_blender_hits.txt"
HEK_Hs4_nn_c3_txt = datadir + "Dseq+/220701_cross_c3/HEK_WT-Hs4_nD-nD_c3/filtered_blender_hits.txt"
HEK_Vs3_Kn_c3_txt = datadir + "Dseq+/220701_cross_c3/HEK_WT-Vs3_KU-nD_c3/filtered_blender_hits.txt"
HEK_Vs3_nn_c3_txt = datadir + "Dseq+/220701_cross_c3/HEK_WT-Vs3_nD-nD_c3/filtered_blender_hits.txt"
HEK_Vs2_Kn_c3_txt = datadir + "Dseq+/220701_cross_c3/HEK_WT-Vs2_KU-nD_c3/filtered_blender_hits.txt"
HEK_Vs2_nn_c3_txt = datadir + "Dseq+/220701_cross_c3/HEK_WT-Vs2_nD-nD_c3/filtered_blender_hits.txt"
# K562 with neg ctrls
K562_Fs2_nD_r1_c2_txt = datadir + "Dseq+/210301_hg19_c2/K562_Fs2_nD_r1_c2/filtered_blender_hits.txt"
K562_Fs2_KU_r1_c2_txt = datadir + "Dseq+/210301_hg19_c2/K562_Fs2_KU_r1_c2/filtered_blender_hits.txt"
K562_Vs2_nD_r1_c2_txt = datadir + "Dseq+/210301_hg19_c2/K562_Vs2_nD_r1_c2/filtered_blender_hits.txt"
K562_Vs2_KU_r1_c2_txt = datadir + "Dseq+/210301_hg19_c2/K562_Vs2_KU_r1_c2/filtered_blender_hits.txt"
K562_Fs2_Kn_c3_txt = datadir + "Dseq+/220701_cross_c3/K562_WT-Fs2_KU-nD_c3/filtered_blender_hits.txt"
K562_Fs2_nn_c3_txt = datadir + "Dseq+/220701_cross_c3/K562_WT-Fs2_nD-nD_c3/filtered_blender_hits.txt"
K562_Vs2_Kn_c3_txt = datadir + "Dseq+/220701_cross_c3/K562_WT-Vs2_KU-nD_c3/filtered_blender_hits.txt"
K562_Vs2_nn_c3_txt = datadir + "Dseq+/220701_cross_c3/K562_WT-Vs2_nD-nD_c3/filtered_blender_hits.txt"
# mm10 with neg ctrls
mm10_mP9_nD_r0_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_nD_r0_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r1_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_nD_r1_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r2_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_nD_r2_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r3_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_nD_r3_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r4_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_nD_r4_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r0_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_KU_r0_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r1_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_KU_r1_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r2_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_KU_r2_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r3_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_KU_r3_c3/filtered_blender_hits.txt"
mm10_mP9_KU_r4_c3_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_KU_r4_c3/filtered_blender_hits.txt"
mm10_mP9_nD_merge_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_nD_merged_c3/filtered_blender_hits.txt"
mm10_mP9_KU_merge_txt = datadir + "Dseq+/211213_mm10_c3/mouse_PCSK9_KU_merged_c3/filtered_blender_hits.txt"
mm10_mP9_Kn_c3_txt = datadir + "Dseq+/220701_cross_c3/mice_WT-mP9_KU-nD_c3/filtered_blender_hits.txt"
mm10_mP9_nn_c3_txt = datadir + "Dseq+/220701_cross_c3/mice_WT-mP9_nD-nD_c3/filtered_blender_hits.txt"

""" gRNA sequences """
gRNA_Vs2 = "GACCCCCTCCACCCCGCCTC"
gRNA_Fs2 = "GCTGCAGAAGGGATTCCATG"
gRNA_Vs3 = "GGTGAGTGAGTGTGTGCGTG"
gRNA_Hs4 = "GGCACTGCGGCTGGAGGTGG"
gRNA_mP9 = "AGCAGCAGCGGCGGCAACAG"

""" Set analysis path """
ana = datadir + "Dseq_out/"
ana_3 = ana + "3_minus/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None


""" ############################################################################################ """
""" Remove false positive sites from initial putative DISCOVER-Seq+ off-target site detection """

print("HEK HEKs4 04h remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Hs4_nD04_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                        c.blender_gen(HEK_Hs4_nn_c3_txt, 2000, hg19, gRNA_Hs4))
c.save_gen(gen, ana_3 + "HEK_Hs4_nD04_final")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Hs4_KU04_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                        c.blender_gen(HEK_Hs4_Kn_c3_txt, 2000, hg19, gRNA_Hs4))
c.save_gen(gen, ana_3 + "HEK_Hs4_KU04_final")

print("HEK HEKs4 12h remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Hs4_nD12_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                        c.blender_gen(HEK_Hs4_nn_c3_txt, 2000, hg19, gRNA_Hs4))
c.save_gen(gen, ana_3 + "HEK_Hs4_nD12_final")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Hs4_KU12_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                        c.blender_gen(HEK_Hs4_Kn_c3_txt, 2000, hg19, gRNA_Hs4))
c.save_gen(gen, ana_3 + "HEK_Hs4_KU12_final")

print("HEK HEKs4 24h remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Hs4_nD24_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                        c.blender_gen(HEK_Hs4_nn_c3_txt, 2000, hg19, gRNA_Hs4))
c.save_gen(gen, ana_3 + "HEK_Hs4_nD24_final")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Hs4_KU24_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                        c.blender_gen(HEK_Hs4_Kn_c3_txt, 2000, hg19, gRNA_Hs4))
c.save_gen(gen, ana_3 + "HEK_Hs4_KU24_final")

print("HEK VEGFAs3 remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Vs3_nD_r1_c2_txt, 2000, hg19, gRNA_Vs3),
                        c.blender_gen(HEK_Vs3_nn_c3_txt, 2000, hg19, gRNA_Vs3))
c.save_gen(gen, ana_3 + "HEK_Vs3_nD_final")
gen = c.gen_subtr_exact(c.blender_gen(HEK_Vs3_KU_r1_c2_txt, 2000, hg19, gRNA_Vs3),
                        c.blender_gen(HEK_Vs3_Kn_c3_txt, 2000, hg19, gRNA_Vs3))
c.save_gen(gen, ana_3 + "HEK_Vs3_KU_final")

print("iPSC VEGFAs2 remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(iPSC_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(HEK_Vs2_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "iPSC_Vs2_nD_final")
gen = c.gen_subtr_exact(c.blender_gen(iPSC_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(HEK_Vs2_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "iPSC_Vs2_KU_final")

print("K562 FANCs2 remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                        c.blender_gen(K562_Fs2_nn_c3_txt, 2000, hg19, gRNA_Fs2))
c.save_gen(gen, ana_3 + "K562_Fs2_nD_final")
gen = c.gen_subtr_exact(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                        c.blender_gen(K562_Fs2_Kn_c3_txt, 2000, hg19, gRNA_Fs2))
c.save_gen(gen, ana_3 + "K562_Fs2_KU_final")

print("K562 VEGFAs2 remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(K562_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(K562_Vs2_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "K562_Vs2_nD_final")
gen = c.gen_subtr_exact(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(K562_Vs2_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "K562_Vs2_KU_final")

print("Mouse mPCSK9 remove false positives")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r0_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_r0")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r0_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_r0")

gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r1_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_r1")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r1_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_r1")

gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r2_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_r2")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r2_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_r2")

gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r3_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_r3")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_r3")

gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r4_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_r4")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r4_c3_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_r4")

gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_merge_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_nn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_merge")
gen = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_merge_txt, 2000, hg19, gRNA_Vs2),
                        c.blender_gen(mm10_mP9_Kn_c3_txt, 2000, hg19, gRNA_Vs2))
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_merge")
gen = c.gen_rmdup(c.blender_gen(ana_3 + "mice_mP9_nD_final_merge.txt", 2000, hg19, gRNA_Vs2), dist=1000)
c.save_gen(gen, ana_3 + "mice_mP9_nD_final_merge_rmdup")
gen = c.gen_rmdup(c.blender_gen(ana_3 + "mice_mP9_KU_final_merge.txt", 2000, hg19, gRNA_Vs2), dist=1000)
c.save_gen(gen, ana_3 + "mice_mP9_KU_final_merge_rmdup")


""" ############################################################################################ """
""" Determine overlap between no drug and Ku-60648 samples (Fig. 2j) """
print("HEK HEKs4 overlap")
a1 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "HEK_Hs4_nD12_final.txt", 2000, hg19, gRNA_Hs4),
                                      c.blender_gen(ana_3 + "HEK_Hs4_KU12_final.txt", 2000, hg19, gRNA_Hs4))]
a2 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "HEK_Hs4_KU12_final.txt", 2000, hg19, gRNA_Hs4),
                                      c.blender_gen(ana_3 + "HEK_Hs4_nD12_final.txt", 2000, hg19, gRNA_Hs4))]
print("HEK VEGFAs3 overlap")
a3 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "HEK_Vs3_nD_final.txt", 2000, hg19, gRNA_Vs3),
                                      c.blender_gen(ana_3 + "HEK_Vs3_KU_final.txt", 2000, hg19, gRNA_Vs3))]
a4 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "HEK_Vs3_KU_final.txt", 2000, hg19, gRNA_Vs3),
                                      c.blender_gen(ana_3 + "HEK_Vs3_nD_final.txt", 2000, hg19, gRNA_Vs3))]
print("iPSC VEGFAs2 overlap")
a5 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "iPSC_Vs2_nD_final.txt", 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(ana_3 + "iPSC_Vs2_KU_final.txt", 2000, hg19, gRNA_Vs2))]
a6 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "iPSC_Vs2_KU_final.txt", 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(ana_3 + "iPSC_Vs2_nD_final.txt", 2000, hg19, gRNA_Vs2))]
print("K562 FANCFs2 overlap")
a7 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "K562_Fs2_nD_final.txt", 2000, hg19, gRNA_Fs2),
                                      c.blender_gen(ana_3 + "K562_Fs2_KU_final.txt", 2000, hg19, gRNA_Fs2))]
a8 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "K562_Fs2_KU_final.txt", 2000, hg19, gRNA_Fs2),
                                      c.blender_gen(ana_3 + "K562_Fs2_nD_final.txt", 2000, hg19, gRNA_Fs2))]
print("K562 VEGFAs2 overlap")
a9 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "K562_Vs2_nD_final.txt", 2000, hg19, gRNA_Vs2),
                                      c.blender_gen(ana_3 + "K562_Vs2_KU_final.txt", 2000, hg19, gRNA_Vs2))]
a10 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "K562_Vs2_KU_final.txt", 2000, hg19, gRNA_Vs2),
                                       c.blender_gen(ana_3 + "K562_Vs2_nD_final.txt", 2000, hg19, gRNA_Vs2))]


""" ############################################################################################ """
""" nature16525 is Kleinstiver et al., 2016 | nbt3117 is Tsai et al., 2015
    Determine overlap between DISCOVER-Seq+ and GUIDE-seq """

print("K562 FANCFs2 comparison to Kleinstiver et al.")
b1 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(ana_3 + "K562_Fs2_KU_final.txt", 2000, hg19, gRNA_Fs2),
                                       c.nature16525_gen(gRNA_Fs2))]
b2 = [g[0] for g in c.gen_subtr_approx(c.nature16525_gen(gRNA_Fs2),
                                       c.blender_gen(ana_3 + "K562_Fs2_KU_final.txt", 2000, hg19, gRNA_Fs2))]

print("K562 VEGFAs2 comparison to Tsai et al.")
b4 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(ana_3 + "K562_Vs2_KU_final.txt", 2000, hg19, gRNA_Vs2),
                                       c.nbt3117_gen(gRNA_Vs2),
                                       approx=500)]
b5 = [g[0] for g in c.gen_subtr_approx(c.nbt3117_gen(gRNA_Vs2),
                                       c.blender_gen(ana_3 + "K562_Vs2_KU_final.txt", 2000, hg19, gRNA_Vs2),
                                       approx=500)]


""" ############################################################################################ """
""" aav9023 is Wienert/Wyman 2020
    Determine overlap between DISCOVER-Seq+ and DISCOVER-Seq in mice """

print("Mouse PCSK9 rep3 nD comparison to KU")
c1 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt", 2000, mm10, gRNA_mP9),
                                      c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt", 2000, mm10, gRNA_mP9))]
c2 = [g[0] for g in c.gen_subtr_exact(c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt", 2000, mm10, gRNA_mP9),
                                      c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt", 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 rep3 nD comparison to individual Wienert/Wyman et al.")
c3 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt", 2000, mm10, gRNA_mP9),
                                       c.aav9023_gen(gRNA_mP9, choice=1))]
c4 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                       c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt", 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 rep3 KU comparison to individual Wienert/Wyman et al.")
c5 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt", 2000, mm10, gRNA_mP9),
                                       c.aav9023_gen(gRNA_mP9, choice=1))]
c6 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                       c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt", 2000, mm10, gRNA_mP9))]

print("Mouse PCSK9 individual Wienert/Wyman SUBTRACT rep3 KU UNION rep3 nD")
d1 = [g[0] for g in c.gen_subtr_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                       c.gen_union_approx(c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt",
                                                                        2000, mm10, gRNA_mP9),
                                                          c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt",
                                                                        2000, mm10, gRNA_mP9)))]
print("Mouse PCSK9 rep3 KU SUBTRACT individual Wienert/Wyman UNION rep3 nD")
d2 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt", 2000, mm10, gRNA_mP9),
                                       c.gen_union_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                                          c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt",
                                                                        2000, mm10, gRNA_mP9)))]
print("Mouse PCSK9 rep3 nD SUBTRACT individual Wienert/Wyman UNION rep3 KU")
d3 = [g[0] for g in c.gen_subtr_approx(c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt", 2000, mm10, gRNA_mP9),
                                       c.gen_union_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                                          c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt",
                                                                        2000, mm10, gRNA_mP9)))]
print("Mouse PCSK9 rep3 nD SUBTRACT individual Wienert/Wyman INTERSECTION rep3 KU")
d4 = [g[0] for g in c.gen_inter_approx(c.blender_gen(ana_3 + "mice_mP9_nD_final_r3.txt", 2000, mm10, gRNA_mP9),
                                       c.gen_inter_approx(c.aav9023_gen(gRNA_mP9, choice=1),
                                                          c.blender_gen(ana_3 + "mice_mP9_KU_final_r3.txt",
                                                                        2000, mm10, gRNA_mP9)))]
