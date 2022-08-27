"""
Script for:
(2)
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
mm10_mP9_nD_merge_txt = datadir + "Dseq+/mouse_PCSK9_nD_merged_c3/filtered_blender_hits.txt"
mm10_mP9_KU_merge_txt = datadir + "Dseq+/mouse_PCSK9_KU_merged_c3/filtered_blender_hits.txt"
mm10_mP9_nD_r2_c2_txt = datadir + "Dseq+/mouse_PCSK9_nD_r2_c2/filtered_blender_hits.txt"
mm10_mP9_KU_r2_c2_txt = datadir + "Dseq+/mouse_PCSK9_KU_r2_c2/filtered_blender_hits.txt"
mm10_mP9_nD_r3_c2_txt = datadir + "Dseq+/mouse_PCSK9_nD_r3_c2/filtered_blender_hits.txt"
mm10_mP9_KU_r3_c2_txt = datadir + "Dseq+/mouse_PCSK9_KU_r3_c2/filtered_blender_hits.txt"

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




print("Mouse PCSK9 rep3 nD comparison to KU")
g = c.gen_subtr_exact(c.blender_gen(mm10_mP9_nD_r3_c3_txt, 2000, mm10, gRNA_mP9),
                       c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, mm10, gRNA_mP9))
c.save_gen(g, ana_3 + "mice_P9_subtr_nD-KU")
g = c.gen_subtr_exact(c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, mm10, gRNA_mP9),
                      c.blender_gen(mm10_mP9_nD_r3_c3_txt, 2000, mm10, gRNA_mP9))
c.save_gen(g, ana_3 + "mice_P9_subtr_KU-nD")

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
