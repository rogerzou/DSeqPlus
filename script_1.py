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
    genome_savepath = "/mnt/c/Users/Roger/Desktop/"         # Directory to store indexed genome
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"]
    mm10 = ['mm10', "/mnt/c/Users/Roger/bioinformatics/mm10_bowtie2/mm10.fa"]
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

ampngs = datadir + "210720_ampNGS/"

""" gRNA sequences """
gRNA_Vs2 = "GACCCCCTCCACCCCGCCTC"
gRNA_Fs2 = "GCTGCAGAAGGGATTCCATG"
gRNA_Vs3 = "GGTGAGTGAGTGTGTGCGTG"
gRNA_Hs4 = "GGCACTGCGGCTGGAGGTGG"
gRNA_Fs2_OFF1 = "GCTGCAAAAAGGATTCCAGG"
gRNA_Fs2_OFF2 = "GACGCAGAAGGGACTCCATG"
gRNA_Fs2_OFF3 = "GCTGGAGGAAGGATTCCATG"
gRNA_Fs2_OFF4 = "GGTACAGAAGGGCTTCCATG"

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
""" For all on- and off-target sites from BLENDER, determine paired-end read subsets """
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


""" ############################################################################################ """
""" For all on- and off-target sites from BLENDER, determine paired-end read subsets at a location
    10kb away (to show that DNA-PKcs inhibition does not globally increase background). """
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
"""  """
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
"""  """


def makegen_k562_fs2():
    return c.gen_subtract(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                          c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2))


def makegen_k562_vs2():
    return c.gen_subtract(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                          c.blender_gen(K562_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2))


# FANCFs2 KU60648 - no drug
c.save_subtract(makegen_k562_fs2(), ana_3 + "K562_Fs2_r1_c2")
c.peak_profile_bp_resolution(makegen_k562_fs2(), K562_Fs2_KU_r1_bam, ana_3 + "K562_Fs2_KU_r1_bam")
c.peak_profile_bp_resolution(makegen_k562_fs2(), K562_Fs2_nD_r1_bam, ana_3 + "K562_Fs2_nD_r1_bam")
# VEGFAs2 KU60648 - no drug
c.save_subtract(makegen_k562_vs2(), ana_3 + "K562_Vs2_r1_c2")
c.peak_profile_bp_resolution(makegen_k562_vs2(), K562_Vs2_KU_r1_bam, ana_3 + "K562_Vs2_KU_r1_bam")
c.peak_profile_bp_resolution(makegen_k562_vs2(), K562_Vs2_nD_r1_bam, ana_3 + "K562_Vs2_nD_r1_bam")


# FANCFs2 no drug - KU60648
gen = c.gen_subtract(c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                     c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2))
[g[0] for g in gen]
gen = c.gen_subtract(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                     c.blender_gen(K562_Fs2_nD_r1_c2_txt, 2000, hg19, gRNA_Fs2))
[g[0] for g in gen]

# VEGFAs2 no drug - KU60648
gen = c.gen_subtract(c.blender_gen(K562_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                     c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2))
[g[0] for g in gen]
gen = c.gen_subtract(c.blender_gen(K562_Vs2_KU_r1_c2_txt, 2000, hg19, gRNA_Vs2),
                     c.blender_gen(K562_Vs2_nD_r1_c2_txt, 2000, hg19, gRNA_Vs2))
[g[0] for g in gen]

# VEGFAs3 no drug - KU60648
gen = c.gen_subtract(c.blender_gen(HEK_Vs3_nD_r1_c2_txt, 2000, hg19, gRNA_Vs3),
                     c.blender_gen(HEK_Vs3_KU_r1_c2_txt, 2000, hg19, gRNA_Vs3))
[g[0] for g in gen]
gen = c.gen_subtract(c.blender_gen(HEK_Vs3_KU_r1_c2_txt, 2000, hg19, gRNA_Vs3),
                     c.blender_gen(HEK_Vs3_nD_r1_c2_txt, 2000, hg19, gRNA_Vs3))
[g[0] for g in gen]

# VEGFAs2 no drug - KU60648
gen = c.gen_subtract(c.blender_gen(HEK_Hs4_nD12_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                     c.blender_gen(HEK_Hs4_KU12_r1_c2_txt, 2000, hg19, gRNA_Hs4))
[g[0] for g in gen]
gen = c.gen_subtract(c.blender_gen(HEK_Hs4_KU12_r1_c2_txt, 2000, hg19, gRNA_Hs4),
                     c.blender_gen(HEK_Hs4_nD12_r1_c2_txt, 2000, hg19, gRNA_Hs4))
[g[0] for g in gen]



""" ############################################################################################ """
"""  """

gg.generate_ref(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                ana_4 + "REF_K562_Fs2_KU_r1_c2", hg19[0], genome_savepath)
REF = SeqIO.to_dict(SeqIO.parse(ana_4 + "REF_K562_Fs2_KU_r1_c2.fa", "fasta"))

amp.getIndels(ampngs + "merged/Fs2_ON_ct.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct.txt")
amp.getIndels(ampngs + "merged/Fs2_ON_ep.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep.txt")

amp.getIndels(ampngs + "merged/Fs2_OFF1_ct.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct.txt")
amp.getIndels(ampngs + "merged/Fs2_OFF1_ep.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep.txt")

amp.getIndels(ampngs + "merged/Fs2_OFF2_ct.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct.txt")
amp.getIndels(ampngs + "merged/Fs2_OFF2_ep.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep.txt")

amp.getIndels(ampngs + "merged/Fs2_OFF3_ct.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct.txt")
amp.getIndels(ampngs + "merged/Fs2_OFF3_ep.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep.txt")

amp.getIndels(ampngs + "merged/Fs2_OFF4_ct.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct.txt")
amp.getIndels(ampngs + "merged/Fs2_OFF4_ep.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep.txt")
