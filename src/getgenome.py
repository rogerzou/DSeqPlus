
from Bio import Entrez, SeqIO
from . import chipseq as c
import pickle
import os
import csv
import re

GENOME, GENOME_LIST = "", ['hg38', 'hg19', 'mm10']
genome_size, genome_id, genome_seq = None, None, None

CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
       'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
       'chr22', 'chrX', 'chrY']


def get_genome_seq(genome_str, savep=None):
    """ Return dict that holds entire sequence for each chromosome. Stores data with pickle to
        avoid multiple retrievals from NCBI.
    :param savep: save path of genome sequence using pickle. If not provided, will be saved to the
                  'lib' folder of this project
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :return: dict with keys as chromosomes, values as SeqIO sequence of each chromosome key
    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    mfile = savep + "/%s.pickle" % genome_str if savep else os.path.dirname(
        dirname) + "/lib/%s.pickle" % genome_str
    if os.path.isfile(mfile):
        print("get_genome_seq(): Loading existing %s sequence dict from %s." % (genome_str, mfile))
        return pickle.load(open(mfile, 'rb'))
    else:
        print("get_genome_seq(): Downloading %s sequence, saving to %s." % (genome_str, mfile))
        d = {}
        g_ids = get_genome_ids(genome_str)
        Entrez.email = "bob@example.org"
        for chr_i in g_ids.keys():
            handle = Entrez.efetch(db="nucleotide", id=g_ids[chr_i], rettype="fasta", strand=1)
            read_i = SeqIO.read(handle, "fasta")
            d[chr_i] = read_i
        pickle.dump(d, open(mfile, 'wb'))
        return d


def get_genome_ids(genome_str):
    """ Return dict that holds the NCBI GI number for each chromosome
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :return: dict with keys as chromosomes, values as NCBI GI number of each chromosome key
    """
    d = {}
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/%s_entrez.csv" % genome_str, 'r') as f:
        for row in csv.reader(f, delimiter=','):
            d[row[0]] = row[2].split(":")[1]
    return d


def genome_initialized(savepath, genome_str):
    """ Initialize downloading of genome GI indices, full genome, and sizes as global variables
    :param savepath: path to save the genome
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    """
    global GENOME, GENOME_LIST, genome_seq, genome_id, genome_size
    if genome_str in GENOME_LIST and genome_str != GENOME:
        GENOME = genome_str
        genome_seq = get_genome_seq(genome_str, savepath)
        genome_id = get_genome_ids(genome_str)
        genome_size = c.get_genome_dict(genome_str)


def generate_ref(generator, outfile, genome_str, savepath, rlen=300):
    """
    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param outfile: string path to output file (extension omitted)
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :param savepath: save path to genome sequence that will be downloaded from NCBI if not already
    :param rlen: number of bases to record left and right from the cut site.
    """
    genome_initialized(savepath, genome_str)
    with open(outfile + ".fa", 'w') as f1:
        for rs, cut, sen, pam, gui, mis, guide in generator:
            chr_i = re.split('[:-]', rs)[0]
            genome_i = genome_seq[chr_i]
            genome_i_len = len(genome_i)
            ind_lt, ind_rt = cut - rlen, cut + rlen
            if chr_i in CHR and ind_lt >= 0 and ind_rt < genome_i_len:
                f1.write(">%s_%i_%s_%s_%i\n%s\n" %
                         (chr_i, cut, gui, guide, mis, genome_i[ind_lt:ind_rt].seq))
