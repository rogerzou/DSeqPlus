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

ampngs1 = datadir + "210720_ampNGS/"
ampngs2 = datadir + "210912_ampNGS/"
ampngs3 = datadir + "210912_ampNGS/"
ampngs4 = datadir + "220118_ampNGS/"
ampngs5 = datadir + "220817_ampNGS/"

""" gRNA sequences """
gRNA_Fs2 = ["GCTGCAGAAGGGATTCCATG", None, None]
gRNA_Fs2_OFF1 = ["GCTGCAAAAAGGATTCCAGG", None, None]
gRNA_Fs2_OFF2 = ["GACGCAGAAGGGACTCCATG", None, None]
gRNA_Fs2_OFF3 = ["GCTGGAGGAAGGATTCCATG", None, None]
gRNA_Fs2_OFF4 = ["GGTACAGAAGGGCTTCCATG", None, None]
gRNA_Fs2_OFF0 = ["GCTGCAGAAGGGATTCCAAG", None, None]
gRNA_Fs2_OFF5 = ["TGTGAAGAAGGGTTTCCATG", None, None]
gRNA_Fs2_OFF6 = ["ACTGCAAAATGGATTCCATG", None, None]
gRNA_Fs2_OFF7 = ["ACTGCTGAGAGGATTCCATG", None, None]
gRNA_Fs2_OFF8 = ["GCTGCTGAACAGATTCCATG", None, None]
gRNA_Fs2_OFF9 = ["GCTGAAGAAGGAATTTCATG", None, None]
gRNA_Fs2_OFF10 = ["GTTGCAGCAAGAATTCCATG", None, None]
gRNA_Fs2_OFF11 = ["CAAGCAGAAGGGATTCCACA", None, None]
gRNA_Fs2_OFF12 = ["GCCACAGAAGGGATTCTATG", None, None]

gRNA_mP9 = ["AGCAGCAGCGGCGGCAACAG", "chr4", 106463862]
gRNA_mP9_OFF0 = ["AGCAGCAGCAGCAGCAACAG", "chr10", 37138630]
gRNA_mP9_OFF1 = ["AGCAGCAGCAGCAGCAACAG", "chr5", 103773150]
gRNA_mP9_OFF2 = ["AGCAGTAGCAGCGGCAACAA", "chr7", 19749248]
gRNA_mP9_OFF3 = ["AGCAGCAGCAGCAGCAACAG", "chr1", 182610256]
gRNA_mP9_OFF4 = ["AGCAGCAGCAGCAGCAACAG", "chr2", 35977759]
gRNA_mP9_OFF5 = ["AGCAACAGGGGCGGCAACAG", "chr11", 99649783]
gRNA_mP9_OFF6 = ["AGCAGCAGCAGCAGCAACAG", "chr16", 96068871]
gRNA_mP9_OFF7 = ["AGCAGCAGCAGCAGCAACAG", "chr19", 26631019]

""" Set analysis path """
ana = datadir + "Dseq_out/"
ana_4 = ana + "4_amplicon/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None


""" ############################################################################################ """
""" Amplicon NGS to check indels in cells without DNA-PKi, with/without Cas9 (3 replicates). """
gg.generate_ref(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                ana_4 + "REF_K562_Fs2_KU_r1_c2", hg19[0], genome_savepath)
REF = SeqIO.to_dict(SeqIO.parse(ana_4 + "REF_K562_Fs2_KU_r1_c2.fa", "fasta"))

amp.getIndels(ampngs1 + "merged/Fs2_ON_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct1.txt", stats=True)
amp.getIndels(ampngs1 + "merged/Fs2_ON_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep1.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_ON_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct2.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_ON_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep2.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_ON_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ct3.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_ON_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ep3.txt", stats=True)

amp.getIndels(ampngs1 + "merged/Fs2_OFF1_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct1.txt", stats=True)
amp.getIndels(ampngs1 + "merged/Fs2_OFF1_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep1.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF1_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct2.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF1_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep2.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF1_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ct3.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF1_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ep3.txt", stats=True)

amp.getIndels(ampngs1 + "merged/Fs2_OFF2_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct1.txt", stats=True)
amp.getIndels(ampngs1 + "merged/Fs2_OFF2_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep1.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF2_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct2.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF2_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep2.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF2_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ct3.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF2_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ep3.txt", stats=True)

amp.getIndels(ampngs1 + "merged/Fs2_OFF3_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct1.txt", stats=True)
amp.getIndels(ampngs1 + "merged/Fs2_OFF3_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep1.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF3_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct2.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF3_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep2.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF3_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ct3.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF3_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ep3.txt", stats=True)

amp.getIndels(ampngs1 + "merged/Fs2_OFF4_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct1.txt", stats=True)
amp.getIndels(ampngs1 + "merged/Fs2_OFF4_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep1.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF4_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct2.txt", stats=True)
amp.getIndels(ampngs2 + "merged/Fs2_OFF4_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep2.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF4_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ct3.txt", stats=True)
amp.getIndels(ampngs3 + "merged/Fs2_OFF4_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF0_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF5_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF7_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF8_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF9_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF10_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF11_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ep3.txt", stats=True)

amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ct1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ct1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ep1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ep1.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ct2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ct2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ep2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ep2.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ct3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ct3.txt", stats=True)
amp.getIndels(ampngs4 + "merged/Fs2_OFF12_ep3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ep3.txt", stats=True)


""" ############################################################################################ """
""" Amplicon NGS to check indels in cells with DNA-PKi, with Cas9 (3 replicates). """
gg.generate_ref(c.blender_gen(K562_Fs2_KU_r1_c2_txt, 2000, hg19, gRNA_Fs2),
                ana_4 + "REF_K562_Fs2_KU_r1_c2", hg19[0], genome_savepath)
REF = SeqIO.to_dict(SeqIO.parse(ana_4 + "REF_K562_Fs2_KU_r1_c2.fa", "fasta"))

amp.getIndels(ampngs5 + "merged/Fs2_ON_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_ON_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_ON_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2, ana_4 + "Fs2_ON_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF0_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF0_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF0_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF0, ana_4 + "Fs2_OFF0_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF1_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF1_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF1_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF1, ana_4 + "Fs2_OFF1_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF2_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF2_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF2_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF2, ana_4 + "Fs2_OFF2_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF3_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF3_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF3_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF3, ana_4 + "Fs2_OFF3_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF4_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF4_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF4_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF4, ana_4 + "Fs2_OFF4_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF5_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF5_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF5_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF5, ana_4 + "Fs2_OFF5_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF7_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF7_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF7_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF7, ana_4 + "Fs2_OFF7_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF8_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF8_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF8_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF8, ana_4 + "Fs2_OFF8_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF9_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF9_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF9_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF9, ana_4 + "Fs2_OFF9_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF10_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF10_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF10_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF10, ana_4 + "Fs2_OFF10_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF11_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF11_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF11_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF11, ana_4 + "Fs2_OFF11_ku3.txt", stats=True)

amp.getIndels(ampngs5 + "merged/Fs2_OFF12_ku1.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ku1.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF12_ku2.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ku2.txt", stats=True)
amp.getIndels(ampngs5 + "merged/Fs2_OFF12_ku3.extendedFrags.fastq", REF,
              gRNA_Fs2_OFF12, ana_4 + "Fs2_OFF12_ku3.txt", stats=True)


# """ ############################################################################################ """
# """ Amplicon NGS to check indels in mice, with/without Cas9 (3 replicates). """
# gg.generate_ref(c.blender_gen(mm10_mP9_KU_r3_c3_txt, 2000, mm10, gRNA_mP9[0]),
#                 ana_4 + "REF_mm10_mP9_KU_r3_c3", mm10[0], genome_savepath)
# REF = SeqIO.to_dict(SeqIO.parse(ana_4 + "REF_mm10_mP9_KU_r3_c3.fa", "fasta"))
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_ON_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9, ana_4 + "mP9_C9_ON_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_ON_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9, ana_4 + "mP9_C9_ON_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_ON_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9, ana_4 + "mP9_nC_ON_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_ON_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9, ana_4 + "mP9_nC_ON_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF0_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF0, ana_4 + "mP9_C9_OFF0_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF0_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF0, ana_4 + "mP9_C9_OFF0_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF0_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF0, ana_4 + "mP9_nC_OFF0_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF0_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF0, ana_4 + "mP9_nC_OFF0_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF1_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF1, ana_4 + "mP9_C9_OFF1_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF1_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF1, ana_4 + "mP9_C9_OFF1_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF1_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF1, ana_4 + "mP9_nC_OFF1_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF1_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF1, ana_4 + "mP9_nC_OFF1_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF2_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF2, ana_4 + "mP9_C9_OFF2_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF2_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF2, ana_4 + "mP9_C9_OFF2_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF2_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF2, ana_4 + "mP9_nC_OFF2_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF2_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF2, ana_4 + "mP9_nC_OFF2_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF3_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF3, ana_4 + "mP9_C9_OFF3_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF3_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF3, ana_4 + "mP9_C9_OFF3_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF3_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF3, ana_4 + "mP9_nC_OFF3_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF3_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF3, ana_4 + "mP9_nC_OFF3_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF4_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF4, ana_4 + "mP9_C9_OFF4_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF4_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF4, ana_4 + "mP9_C9_OFF4_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF4_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF4, ana_4 + "mP9_nC_OFF4_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF4_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF4, ana_4 + "mP9_nC_OFF4_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF5_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF5, ana_4 + "mP9_C9_OFF5_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF5_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF5, ana_4 + "mP9_C9_OFF5_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF5_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF5, ana_4 + "mP9_nC_OFF5_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF5_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF5, ana_4 + "mP9_nC_OFF5_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF6_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF6, ana_4 + "mP9_C9_OFF6_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF6_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF6, ana_4 + "mP9_C9_OFF6_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF6_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF6, ana_4 + "mP9_nC_OFF6_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF6_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF6, ana_4 + "mP9_nC_OFF6_ku1.txt", stats=True)
#
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF7_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF7, ana_4 + "mP9_C9_OFF7_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_C9_OFF7_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF7, ana_4 + "mP9_C9_OFF7_ku1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF7_ct1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF7, ana_4 + "mP9_nC_OFF7_ct1.txt", stats=True)
# amp.getIndels(ampngs5 + "merged/mP9_nC_OFF7_ku1.extendedFrags.fastq", REF,
#               gRNA_mP9_OFF7, ana_4 + "mP9_nC_OFF7_ku1.txt", stats=True)
