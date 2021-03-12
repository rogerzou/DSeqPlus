
main() {
  iPSC_Vs2_r1 & K562_Vs2_r1 & K562_Fs2_r1
}

iPSC_Vs2_r1() {
    # iPSC VEGFAs2 replicate 1 (no drug)
    mkdir mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_nD_r1
    sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38 \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A13_hg38_final.bam \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A19_hg38_final.bam \
        GACCCCCTCCACCCCGCCTC \
        mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_nD_r1
    # iPSC VEGFAs2 replicate 1 (w/ KU-60648)
    mkdir mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_KU_r1
    sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38 \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A14_hg38_final.bam \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A19_hg38_final.bam \
        GACCCCCTCCACCCCGCCTC \
        mnt/z/rzou4/NGS_data/4_damage/Dseq+/iPSC_Vs2_KU_r1
}

K562_Vs2_r1() {
    # K562 VEGFAs2 replicate 1 (no drug)
    mkdir mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_nD_r1
    sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38 \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A15_hg38_final.bam \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A20_hg38_final.bam \
        GACCCCCTCCACCCCGCCTC \
        mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_nD_r1
    # K562 VEGFAs2 replicate 1 (w/ KU-60648)
    mkdir mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_KU_r1
    sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38 \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A16_hg38_final.bam \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A20_hg38_final.bam \
        GACCCCCTCCACCCCGCCTC \
        mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Vs2_KU_r1
}

K562_Fs2_r1() {
    # K562 FANCFs2 replicate 1 (no drug)
    mkdir mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_nD_r1
    sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38 \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A17_hg38_final.bam \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A20_hg38_final.bam \
        GACCCCCTCCACCCCGCCTC \
        mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_nD_r1
    # K562 FANCFs2  replicate 1 (w/ KU-60648)
    mkdir mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_KU_r1
    sh run_blender.sh /mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38 \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A18_hg38_final.bam \
        mnt/z/rzou4/NGS_data/4_damage/210301_chipseq/A20_hg38_final.bam \
        GACCCCCTCCACCCCGCCTC \
        mnt/z/rzou4/NGS_data/4_damage/Dseq+/K562_Fs2_KU_r1
}

main "$@"; exit
