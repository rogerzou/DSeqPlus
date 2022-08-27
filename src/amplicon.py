
from Bio import SeqIO
from . import chipseq as c
import numpy as np
import subprocess as sp


def call_bowtie2(inf, hg38):
    """ Run bowtie2 to align 'inf' fasta file to hg38, result in 'outfile' sam file """
    sp.run(['bowtie2', '-f', '-x', hg38[:-3], '-U', inf + ".fa", '-S', inf + ".sam"])


def subset_reads(fq_file, out_file, i5, i7):
    """ Demultiplex reads by i5 and i7 indices. """
    infile = SeqIO.parse(fq_file, "fastq")
    out_reads = []
    for read in infile:
        if i5 == getI5(read) and i7 == getI7(read):
            out_reads.append(read)
    SeqIO.write(out_reads, out_file, "fastq")


def getI7(read):
    return read.description.split(" ")[1][6:14]


def getI5(read):
    return read.description.split(" ")[1][15:23]


def getIndels(f_fastq, ref, target, outpath, stats=False, maxcount=False):
    exp_seq = list(SeqIO.parse(f_fastq, "fastq"))
    # get index for sequence of interest in ref
    locus = None
    for key in ref.keys():
        chr_i, coo_i, seq_i = key.split('_')[0:3]
        if chr_i == target[1] and int(coo_i) == target[2]:
            locus = key
            break
    if not locus:
        for key in ref.keys():
            seq_i = key.split('_')[2]
            if seq_i == target[0] or seq_i == c.get_reverse_complement(target[0]):
                locus = key
                break
    ref_seq = ref[locus].seq.upper()
    ind = max(ref_seq.find(target[0]), ref_seq.find(c.get_reverse_complement(target[0])))
    f_ind = max(0, ind - 40)
    r_ind = ind + 40
    fwd = str(ref_seq[f_ind:f_ind + 20])
    rev = str(ref_seq[r_ind:r_ind + 20])
    ref_region = str(ref_seq[ref_seq.find(fwd) + 20: ref_seq.find(rev)])
    readcounter = 0
    indelcounter = 0
    indel_list = []
    if maxcount is not False:
        exp_seq = exp_seq[:min(len(exp_seq), maxcount)]
    for exp in exp_seq:
        find0 = max(exp.seq.find(fwd), exp.seq.find(c.get_reverse_complement(fwd)))
        find1 = max(exp.seq.find(rev), exp.seq.find(c.get_reverse_complement(rev)))
        fwd_find = min(find0, find1)
        rev_find = max(find0, find1)
        if fwd_find != -1 and rev_find != -1:
            exp_regraw = exp[fwd_find + 20: rev_find]
            exp_region = str(exp_regraw.seq)
            if np.mean(exp_regraw.letter_annotations["phred_quality"]) >= 20:
                readcounter += 1
                if len(ref_region) != len(exp_region):
                    indelcounter += 1
                    indel_list.append(exp_region)
    if readcounter > 0:
        printstring = "Detected %f indels (%i indels / %i total)\n" % \
                      (indelcounter * 1.0 / readcounter, indelcounter, readcounter)
        print(printstring)
        with open(outpath, 'w') as filehandle:
            filehandle.write(printstring)
            filehandle.writelines("%s\n" % indel for indel in indel_list)
        filehandle.close()
    if readcounter == 0:
        print("No reads passed filter for this sample :(\n")
    if stats:
        pluscounter = [0 for i in range(51)]
        minuscounter = [0 for i in range(51)]
        wtcounter = 0
        for indel in indel_list:
            diff = len(indel) - len(ref_region)
            if diff == 0:       wtcounter += 1
            if 50 > diff > 0:   pluscounter[diff-1] += 1
            if diff >= 50:      pluscounter[50] += 1
            if -50 < diff < 0:  minuscounter[abs(diff)-1] += 1
            if diff <= -50:     minuscounter[50] += 1
        totalcounter = minuscounter[::-1] + [readcounter - indelcounter] + pluscounter
        indexcounter = [i for i in range(-51, 52)]
        with open(outpath+".counter", 'w') as filehandle:
            filehandle.write("# of total reads: %i\n" % readcounter)
            filehandle.write("# of total indels: %i\n" % indelcounter)
            filehandle.write("# of +n insertions:\n")
            filehandle.writelines(["%i\t" % ind for ind in pluscounter])
            filehandle.write("\n# of +n deletions:\n")
            filehandle.writelines(["%i\t" % ind for ind in minuscounter])
            filehandle.write("\nIndel distribution:\n")
            filehandle.writelines(["%i\t" % ind for ind in indexcounter])
            filehandle.write("\n")
            filehandle.writelines(["%i\t" % ind for ind in totalcounter])
        filehandle.close()
