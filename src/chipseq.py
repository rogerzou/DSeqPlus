# -*- coding: utf-8 -*-
""" Useful functions for ChIP-seq analysis after Cas9 cleavage
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

from collections import defaultdict
import pysam
import os
import re
import numpy as np
import csv

OLD_CHAR = "ACGT"
NEW_CHAR = "TGCA"
WEIGHT = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 3]


def get_reverse_complement(seq):
    """ Return reverse complement of sequence string. """
    return seq.translate(str.maketrans(OLD_CHAR, NEW_CHAR))[::-1]


def status_statement(current, final, count, chromosome=None):
    """ Print progress statements for long processes

    :param current: current iteration number
    :param final: total number of iterations
    :param count: number of print statements total
    :param chromosome: print chromosome information of current progress

    """
    if current % int(final/count) == 0:
        if chromosome is None:
            print("Processed %i out of %i" % (current, final))
        else:
            print("Processed %i out of %i in %s" % (current, final, chromosome))


def region_string_split(rs):
    return re.split('[:-]', rs)


def get_genome_dict(genome_str):
    """ Return dict that holds the number of base pairs for each chromosome in 'hg38','hg19','mm10'.

    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :return: dict with keys as chromosomes, values as maximum coordinate of each chromosome key
    """
    d = {}
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/%s.sizes" % genome_str, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            d[row[0]] = int(row[1])
    return d


def get_genome_generator(genome_str):
    """ Generate the number of base pairs for each chromosome in 'hg38','hg19','mm10'.

    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :return generator that outputs the number of base pairs for each chromosome in hg38
            in the format: [chr7, 159345973]
    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/%s.sizes" % genome_str, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            yield row


def read_pair_generator(bam, region_string=None):
    """ Generate read pairs in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.

    :param bam: pysam AlignmentFile loaded with BAM file that contains paired-end reads
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :return generator of read pairs with the following format: [read1, read2]

    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_pair_align(read1, read2):
    """ Extract read pair locations as a fragment oriented in increasing chromosome coordinates

    :param read1: read #1 of pair in pysam AlignedSegment format
    :param read2: read #2 of pair in pysam AlignedSegment format
    :return 4-item array in the following format: [fragA-start, fragA-end, fragB-start, fragB-end]
            with monotonically increasing chromosome coordinates

    """
    r1pos = [x+1 for x in read1.positions]
    r2pos = [x+1 for x in read2.positions]
    if read1.mate_is_reverse and r1pos[0] < r2pos[0]:  # read1 is earlier
        read = [r1pos[0], r1pos[-1], r2pos[0], r2pos[-1]]
    elif read2.mate_is_reverse and r2pos[0] < r1pos[0]:  # read2 is earlier
        read = [r2pos[0], r2pos[-1], r1pos[0], r1pos[-1]]
    else:
        read = []
        # print("Skipping read pair from error in alignment.")
    # print("%s--%s>  <%s--%s" % tuple(read))
    return read


def load_nparray(array):
    """ Load dataset as numpy array. """
    return np.loadtxt(array, dtype=object, delimiter=',')


def load_npheader(array):
    """ Load the header (first line) of numpy array file. """
    with open(array) as f:
        head = f.readline()
    return head[2:-1]


def mismatch_filter_gen(generator, mismatch):
    """ Filter a cut site generator by the # of mismatches

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param mismatch: filter to only yield cut sites that have specified # of mismatches from
                     gRNA sequence
    """
    for rs, cut, sen, pam, gui, mis, tar, val in generator:
        if mis == mismatch:
            yield rs, cut, sen, pam, gui, mis, tar, val


def shift_gen(generator, dist=10000):
    """ Shift the coordinate position of input generator by bp distance 'dist'.
    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param dist: distance away from original position to shift
    :return: generator with each position shifted by dist
    """
    for rs, cut, sen, pam, gui, mis, guide, val_i in generator:
        rs_spt = re.split('[:-]', rs)
        chr_i, lt_i, rt_i = rs_spt[0], int(rs_spt[1]) + dist, int(rs_spt[2]) + dist
        cut += dist
        yield "%s:%i-%i" % (chr_i, lt_i, rt_i), cut, sen, pam, gui, mis, guide, val_i


def nbt3117_gen(expected_guide):
    """ Generator to yield all outputted sites from 41587_2015_BFnbt3117_MOESM22_ESM.xlsx,
        which is the Supplementary Table 2 of
        Tsai, Shengdar Q., et al. "GUIDE-seq enables genome-wide profiling of off-target cleavage
        by CRISPR-Cas nucleases." Nature biotechnology 33.2 (2015): 187-197.

    :param expected_guide: on-target protospacer sequence (no PAM)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i )
        ( region string, cut site, sense/antisense, PAM, discovered protospacer,
        # mismatches, non-mismatched protospacer, GUIDE-seq reads )
    """
    data = load_nparray('lib/nbt3117.csv')
    for d in data[1:]:
        try:
            span_rs = "%s:%s-%s" % (d[0], d[1], d[2])
            cut_i = int((int(d[1]) + int(d[2])) / 2)
            val_i = int(d[4])
            sen_i = d[5]
            pam_i = d[9][20:]
            gui_i = d[9][:20]
            mis_i = int(d[10])
            guide = d[8][:20]
            if guide == expected_guide:
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i
        except ValueError:
            print("nbt3117_gen(): Skipped entry due to incorrect formatting.")


def nature16525_gen(expected_guide, verbose=False):
    """ Generator to yield all outputted sites from 41586_2016_BFnature16526_MOESM91_ESM.xlsx,
        which is the Supplementary Table 4 of
        Kleinstiver, Benjamin P., et al. "High-fidelity CRISPR–Cas9 nucleases with no detectable
        genome-wide off-target effects." Nature 529.7587 (2016): 490-495.

    :param expected_guide: on-target protospacer sequence (no PAM)
    :param verbose: print ValueError statement
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i )
        ( region string, cut site, sense/antisense, PAM, discovered protospacer,
        # mismatches, non-mismatched protospacer, GUIDE-seq reads )
    """
    data = load_nparray('lib/nature16526.csv')
    for d in data[1:]:
        try:
            span_rs = "%s:%s-%s" % (d[0], d[1], d[2])
            cut_i = int((int(d[1]) + int(d[2])) / 2)
            val_i = int(d[4])
            sen_i = d[5]
            pam_i = d[10][20:]
            gui_i = d[10][:20]
            mis_i = int(d[11])
            guide = d[9][:20]
            if d[7] == 'wild-type' and guide == expected_guide:
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i
        except ValueError:
            if verbose:
                print("nature16525_gen(): Skipped entry due to incorrect formatting:.")
                print(d)


def aav9023_gen(expected_guide, choice=0, verbose=False):
    """ Generator to yield all outputted sites from aav9023_pooled.csv or aav9023_individual.csv,
        which is a Supplementary Table of
        Wienert, B., Wyman, S.K., et al. Unbiased detection of CRISPR off-targets in vivo using
        DISCOVER-Seq. Science 364, 286-289 (2019).

    :param expected_guide: on-target protospacer sequence (no PAM)
    :param choice: 0 for the pooled dataset, 1 for individual rep1, 2 for individual rep2
    :param verbose:
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i )
        ( region string, cut site, sense/antisense, PAM, discovered protospacer,
        # mismatches, non-mismatched protospacer, DISCO score )
    """
    if choice == 0:
        data = load_nparray('lib/aav9023_pooled.csv')
    elif choice == 1 or choice == 2:
        data = load_nparray('lib/aav9023_individual.csv')
    else:
        return ValueError("aav9023_gen(): Incorrect entry for choice variable.")
    for d in data[1:]:
        try:
            span_rs = d[5]
            cut_i = int(d[6])
            sen_i = '+' if d[7] == "sense" else '-'
            pam_i = d[8]
            gui_i = d[9]
            mis_i = int(d[10])
            guide = expected_guide
            val_i = int(d[12])
            if choice == 0:
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i
            elif choice == 1 and d[11] == '24h_1':
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i
            elif choice == 2 and d[11] == '24h_2':
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i
        except ValueError:
            if verbose:
                print("aav9023_gen(): Skipped entry due to incorrect formatting.")
                print(d)


def blender_gen(blender, span_r, genome, guide):
    """ Generator to yield all peaks from BLENDER output.

    :param blender: blender output file
    :param span_r: radius of window from peak center for downstream analysis
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param guide: on-target protospacer sequence (no PAM)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i )
        ( region string, cut site, sense/antisense, PAM, discovered protospacer,
        # mismatches, non-mismatched protospacer, DISCO score )
    """
    b_out = np.loadtxt(blender, dtype=object, skiprows=1, delimiter='\t', ndmin=2)
    numpeaks = b_out.shape[0]
    hgsize = get_genome_dict(genome[0])
    for i in range(numpeaks):
        chr_i = re.split('[:-]', b_out[i, 0])[0]
        if chr_i in hgsize:
            if b_out[i, 4] == 'sense':
                sen_i = '+'
                cut_i = int(b_out[i, 1])
            else:
                sen_i = '-'
                cut_i = int(b_out[i, 1]) + 1
            val_i = int(b_out[i, 2])
            pam_i = b_out[i, 5]
            gui_i = b_out[i, 6]
            mis_i = int(b_out[i, 7])
            span_sta = max(1, cut_i - span_r)
            span_end = min(hgsize[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide, val_i


def sub_findmis(s, matchstr, maxmismatch):
    """ Find the best protospacer (with mismatch tolerance) in a sequence range

    :param s: sequence string in which to find match
    :param matchstr: match sequence string
    :param maxmismatch: mismatch tolerance (int)
    """
    sublen = len(matchstr)
    len_full = sublen + 3
    candidates = []
    # check mismatches for original string
    for i in range(len(s)-sublen-6):
        i += 3
        # sense
        substr_1 = s[i:i + sublen]
        pam_1 = s[i + sublen:i + sublen + 3]
        cutsite_1 = i + sublen - 4
        str_1 = substr_1 + pam_1
        matchstr_1 = matchstr + 'NGG'
        wlist = [x + 1 if str_1[j] != matchstr_1[j] else x for j, x in enumerate([0] * len_full)]
        mismatches = sum(wlist) - 1
        if mismatches <= maxmismatch:
            score = sum([a*b for a, b in zip(wlist, WEIGHT)])
            candidates.append((score, mismatches, cutsite_1, substr_1, pam_1, 1))
        # anti-sense
        substr_0 = get_reverse_complement(s[i:i + sublen])
        pam_0 = get_reverse_complement(s[i - 3:i])
        cutsite_0 = i + 3
        str_0 = substr_0 + pam_0
        matchstr_0 = matchstr + 'NGG'
        wlist = [x + 1 if str_0[j] != matchstr_0[j] else x for j, x in enumerate([0] * len_full)]
        mismatches = sum(wlist) - 1
        if mismatches <= maxmismatch:
            score = sum([a*b for a, b in zip(wlist, WEIGHT)])
            candidates.append((score, mismatches, cutsite_0, substr_0, pam_0, 0))

    candidates.sort(key=lambda y: y[0])
    return candidates


def peak_profile_wide(generator, genome, bamfilein, fileout, norm_type=None,
                      span_rad=2000000, res=1000, wind_rad=10000):
    """ For each target location from generator, calculates enrichment at specified 'resolution'
        with sliding window of specified 'radius'. Outputs the enrichment from a BAM file as:
        (1) CSV file with each row one target location, column is enrichment values in a window
        centered at each target location, at a specific genomic resolution.
        (2) WIG file with the local (window-bounded) enrichment at each target location

    :param bamfilein: path to input BAM file
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param fileout: path to output file name (excludes extension)
    :param norm_type: If None (default), then normalize to RPM. If False, then no normalization.
                       Otherwise, if list, then assume list of region strings for normalization,
                       i.e. ['chr12:6532000-6536000', 'chr15:44709500-44713500']
                       Otherwise, assume it is an alternative BAM file to use for RPM normalization.
    :param span_rad: radius of analysis, centered at the cut site | default 2E6 bp
    :param res: resolution, i.e. bp to skip to calculate enrichment | default 1E3 bp
    :param wind_rad: radius of sliding window | default 1E4 bp

    Results in two files (WIG and CSV), described above.
    """
    hgsize = get_genome_dict(genome[0])
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    chr_old, csv_peaks = None, []
    wlist_all = []
    numrows = int(span_rad * 2 / res) + 1
    if norm_type is None:
        norm_num = bamin.mapped / 1E6
    elif not norm_type:
        norm_num = 1
    elif isinstance(norm_type, list):
        norm_num = sum(bamin.count(region=co) for co in norm_type) / len(norm_type) / 10
    else:
        bamalt = pysam.AlignmentFile(norm_type, 'rb')
        norm_num = bamalt.mapped / 1E6
        bamalt.close()
    for rs, cut, sen, pam, gui, mis, guide, val in generator:
        chr_i = re.split('[:-]', rs)[0]
        sta_i = cut - span_rad
        end_i = cut + span_rad
        if sta_i - wind_rad >= 0 and end_i + wind_rad < hgsize[chr_i]:
            wlist = [0] * numrows
            for row_i in range(numrows):
                center = sta_i + row_i * res
                rs_i = "%s:%i-%i" % (chr_i, center - wind_rad, center + wind_rad)
                wlist[row_i] = bamin.count(region=rs_i) / norm_num
            wlist_all.append([chr_i, sta_i] + wlist)
            csv_peaks.append([chr_i, cut, gui + pam, mis] + wlist)
    bamin.close()
    _peak_profile_helper(wlist_all, res, fileout)
    head = ",".join(["chr", "cut", "guide+PAM", "mismatches"] +
                    list(map(str, range(-span_rad, span_rad + 1, res))))
    np.savetxt(fileout + "_bpeaks.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',', header=head)


def peak_profile_bp_resolution(generator, bamfilein, fileout, norm_type=None):
    """ For each target location from generator, calculates enrichment at each base pair as the
        number of fragments that 'span' the base, i.e. the base is either (1) sequenced by either
        end of paired-end sequencing, or (2) not sequenced but spanned by the imputed DNA fragment.
        Output the enrichment from a BAM file as:
        (1) CSV file with each row one target location, column is enrichment values in a window
        centered at each target location
        (2) WIG file with the local (window-bounded) enrichment at each target location

    :param bamfilein: path to input BAM file
    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param fileout: path to output file name (excludes extension)
    :param norm_type: If None (default), then normalize to RPM. If False, then no normalization.
                       Otherwise, if list, then assume list of region strings for normalization,
                       i.e. ['chr12:6532000-6536000', 'chr15:44709500-44713500']
                       Otherwise, assume it is an alternative BAM file to use for RPM normalization.
    Results in two files (WIG and CSV), described above.
    """
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    chr_old, csv_peaks = None, []
    wlist_all = []
    if norm_type is None:
        norm_num = bamin.mapped / 1E6
    elif not norm_type:
        norm_num = 1
    elif isinstance(norm_type, list):
        norm_num = sum(bamin.count(region=co) for co in norm_type) / len(norm_type) / 10
    else:
        bamalt = pysam.AlignmentFile(norm_type, 'rb')
        norm_num = bamalt.mapped / 1E6
        bamalt.close()
    for rs, cut, sen, pam, gui, mis, guide, val in generator:
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        sta_i = int(sta_i)
        end_i = int(end_i)
        wlist = [0] * (end_i - sta_i + 1)
        for read1, read2 in read_pair_generator(bamin, rs):
            read = read_pair_align(read1, read2)
            if not read:
                continue
            wlist = [x + 1 if read[0] - sta_i <= j <= read[-1] - sta_i else x for j, x in
                     enumerate(wlist)]
        wlist = [x / norm_num for x in wlist]
        wlist_all.append([chr_i, sta_i] + wlist)
        wlist = wlist if sen == '+' else wlist[::-1]
        csv_peaks.append([chr_i, cut, gui + pam, mis] + wlist)
    bamin.close()
    _peak_profile_helper(wlist_all, 1, fileout)
    span_rad = int((len(wlist_all[0])-3) / 2)
    head = ",".join(
        ["chr", "cut", "guide+PAM", "mismatches"] + list(map(str, range(-span_rad, span_rad + 1))))
    np.savetxt(fileout + "_bpeaks.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',', header=head)


def _peak_profile_helper(wlist_all, resolution, fileout):
    """ Helper function for peak_profile_bp_resolution() or peak_profile_wide(). Writes all peak
        profiles to a wiggle file.

    :param wlist_all:
    :param resolution:
    :param fileout:
    """
    with open(fileout + "_bpeaks.wig", 'w') as wigout:
        wlist_all = np.asarray(wlist_all)
        wlist_all = wlist_all[np.lexsort((wlist_all[:, 1].astype(int), wlist_all[:, 0])), :]
        chr_prev = None
        for i in range(wlist_all.shape[0]):
            row = wlist_all[i, :]
            chr_i = row[0]
            if chr_prev != chr_i:
                chr_prev = chr_i
                wigout.write("variableStep\tchrom=%s\n" % chr_i)
            sta_i = int(row[1])
            wlist = row[2:].astype(float)
            for j, x in enumerate(wlist):
                wigout.write("%i\t%0.5f\n" % (sta_i + j * resolution, x))


def read_subsets(generator, filein, fileout):
    """ For each target location from generator, divides the reads from filein BAM file into ones
        spanning cut site, abutting cut site, or neither as separate BAM files. Also, calculates
        total enrichment in a window centered at each target location, as well as enrichment for
        different subsets of reads, including 'left' and 'right' of the cut site in the sense
        orientation.
        Determines if each target location resides in genes - if so, determines the direction, i.e.
        'left' or 'right' of the cut site in gene orientation.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param filein: path to input BAM file for analysis
    :param fileout: path to output file name (excluding extensions)

    OUTPUT: For all target sites, within window specified by generator, outputs
    (1) fileout_span.bam: BAM file with reads that span each target site
    (2) fileout_abut.bam: BAM file with reads that abut (i.e. within 5 bp) of each target site
    (3) fileout_else.bam: BAM file with reads that neither span nor abut
    (4) fileout.csv: CSV file with information on each target site, including target sequence,
    mismatch status, total enrichment, subset enrichment (span, abut, else, abut-left, abut-right),
    orientation relative to sense strand, presence on gene, orientation relative to transcription.
    (5) fileout_pambias.csv: CSV file with information about the number of reads that begin
    immediately PAM-distal or PAM-proximal from the cut site, with flexibility due to either
    staggered or blunt cut phenotype.
    """
    bamin = pysam.AlignmentFile(filein, 'rb')             # BAM file for analysis
    bamspan = pysam.AlignmentFile(fileout + "_span.bam", 'wb', template=bamin)
    bamabut = pysam.AlignmentFile(fileout + "_abut.bam", 'wb', template=bamin)
    bamelse = pysam.AlignmentFile(fileout + "_else.bam", 'wb', template=bamin)
    LS = []
    for rs, cut, sen, pam, gui, mis, tar, val in generator:
        ctM, ctN, ctL, ctR, ctT = 0, 0, 0, 0, 0
        for read1, read2 in read_pair_generator(bamin, rs):
            read = read_pair_align(read1, read2)
            if not read:
                continue
            ctT += 1
            if read[0] < cut < read[3]:         # fragments that span
                ctM += 1
                bamspan.write(read1)
                bamspan.write(read2)
            elif cut + 5 >= read[0] >= cut:     # fragments that begin 5bp of cleavage site
                ctR += 1
                ctN += 1
                bamabut.write(read1)
                bamabut.write(read2)
            elif cut - 5 <= read[-1] <= cut:    # fragments that end 5bp of cleavage site
                ctL += 1
                ctN += 1
                bamabut.write(read1)
                bamabut.write(read2)
            else:                               # other fragments
                bamelse.write(read1)
                bamelse.write(read2)
        fpm = bamin.mapped / 2E6      # fragments per millon
        ctM /= fpm                  # fragments that span
        ctN /= fpm                  # fragments that don't span
        ctL /= fpm                  # fragments that don't span, on left
        ctR /= fpm                  # fragments that don't span, on right
        ctT /= fpm                  # count total number of reads
        cPr, cDi = None, None
        if sen == '+':
            cDi, cPr = ctL, ctR
        if sen == '-':
            cDi, cPr = ctR, ctL
        LS.append((rs, cut, sen, gui, pam, tar, mis, ctT, ctM, ctN, cPr, cDi))
    bamin.close()
    bamspan.close()
    bamabut.close()
    bamelse.close()
    file_array = [fileout + "_span.bam", fileout + "_abut.bam", fileout + "_else.bam"]
    for file_i in file_array:
        pysam.sort("-o", file_i, file_i)
        os.system("samtools index " + file_i)
    LSheader = "region string, cut site, Cas9 sense, observed target sequence, PAM, " \
               "expected target sequence, mismatches, Ctotal, Cspan, Cend, Cproximal, Cdistal"
    LS = np.asarray(LS)
    np.savetxt(fileout + ".csv", LS, fmt='%s', delimiter=',', header=LSheader)


def read_counts(generator, filein, fileout=None):
    """ For each target location from generator, determine total enrichment within target-centered
        window.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param filein: path to input BAM file for analysis
    :param fileout: path to output file name (excluding extensions)

    For all target sites, within window specified by generator, outputs CSV file with information
    on each target site, including target sequence, mismatch status, and total enrichment in window.
    """
    bam = pysam.AlignmentFile(filein, 'rb')             # BAM file for analysis
    list_stat = []
    for rs, cut, sen, pam, gui, mis, tar, val in generator:
        ct_rpm = bam.count(region=rs) / bam.mapped * 1E6
        list_stat.append((rs, cut, sen, gui, mis, ct_rpm))
    bam.close()
    list_stat = np.asarray(list_stat)
    if fileout:
        header = 'region_string, cut site, sense, guide, mismatches, counts_RPM'
        np.savetxt(fileout, list_stat, fmt='%s', delimiter=',', header=header)
    return list_stat


def save_gen(gen, fileout):
    """ Save generator in the same format as BLENDER output
    :param gen: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param fileout: path to txt output, extension omitted.
    """
    outa = []
    [outa.append(g) for g in gen]
    head = 'Chr: Start - End\tCutsite\tDiscoscore\tCutsite Ends\tStrand\tPAM\tGuide sequence\tMismatches'
    out_array = np.asarray(outa)
    out_array = out_array[:, [0, 1, 7, 7, 2, 3, 4, 5]]
    antisense = out_array[:, 4] == '-'
    out_array[antisense, 4] = "antisense"
    out_array[~antisense, 4] = "sense"
    out_array[antisense, 2] = out_array[antisense, 2].astype(int) - 1
    np.savetxt(fileout + ".txt", out_array, fmt='%s', delimiter='\t', header=head)


def gen_rmdup(gen1, dist):
    """ Outputs generator for set of target sites from gen1, with duplicates removed.
    Duplicates defined as two target sites in the same chromosome, within 'dist' bp apart, same # mismatches
    and target sequence. Remove the duplicate site with lower enrichment value.
    :param gen1: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param dist: max distance to be considered a duplicate by distance criteria.
    :return: generator with duplicates removed.
    """
    gen_list = list(gen1)
    stop_while = False
    while not stop_while:
        stop_if = False
        for ind_i in range(len(gen_list) - 1):
            gen_1, gen_2 = gen_list[ind_i], gen_list[ind_i + 1]
            chr_1, chr_2 = re.split('[:-]', gen_1[0])[0], re.split('[:-]', gen_2[0])[0]
            cut_1, cut_2 = gen_1[1], gen_2[1]
            val_1, val_2 = gen_1[-1], gen_2[-1]
            mis_1, mis_2 = gen_1[-3], gen_2[-3]
            gui_1, gui_2 = gen_1[-2], gen_2[-2]
            if chr_1 == chr_2:  # if chr are the same and distance is close...remove site with lesser value
                if mis_1 == mis_2 and gui_1 == gui_2 and abs(cut_1 - cut_2) < dist:
                    gen_list.pop(ind_i + 1) if val_1 > val_2 else gen_list.pop(ind_i)
                    stop_if = True
                    print("gen_rmdup(): Removed %s" % gen_1[0])
            if stop_if:
                break
        if not stop_if:
            stop_while = True
    for gen_i in gen_list:
        yield gen_i


def gen_subtr_exact(gen1, gen2):
    """ Outputs generator for set of target sites from gen1, subtracted by gen2.
    :param gen1: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param gen2: same format as gen1
    :return: generator for set of target sites from gen1, subtracted by gen2
    """
    gen2list = []
    count2, count1, countdiff = 0, 0, 0
    for g in gen2:
        cutloc_i = (re.split('[:-]', g[0])[0], g[1])
        gen2list.append(cutloc_i)
        count2 += 1
    for g in gen1:
        cutloc_i = (re.split('[:-]', g[0])[0], g[1])
        count1 += 1
        if cutloc_i not in gen2list:
            countdiff += 1
            yield g
    print("# of target sites | gen1: %i | gen2: %i | gen1-gen2: %i" % (count1, count2, countdiff))


def gen_subtr_approx(gen1, gen2, approx=1):
    """ Outputs generator for set of target sites from gen1, subtracted by gen2.
    :param gen1: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param gen2: same format as gen1
    :param approx: max distance in coordinates between cut sites to be counted as the same one
    :return: generator for set of target sites from gen1, subtracted by gen2
    """
    gen2list = []
    count2, count1, countdiff = 0, 0, 0
    for g2 in gen2:
        cutloc_i = (re.split('[:-]', g2[0])[0], g2[1])
        gen2list.append(cutloc_i)
        count2 += 1
    for g1 in gen1:
        cutloc_1 = (re.split('[:-]', g1[0])[0], g1[1])
        count1 += 1
        bool_include = True
        for cutloc_2 in gen2list:
            if cutloc_1[0] == cutloc_2[0] and abs(cutloc_1[1] - cutloc_2[1]) <= approx:
                bool_include = False
        if bool_include:
            countdiff += 1
            yield g1
    print("# of target sites | gen1: %i | gen2: %i | gen1-gen2: %i" % (count1, count2, countdiff))


def gen_inter_approx(gen1, gen2, approx=1):
    """ Outputs generator for set of target sites that is intersection of gen1 and gen2
    :param gen1: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param gen2: same format as gen1
    :param approx: max distance in coordinates between cut sites to be counted as the same one
    :return: generator for set of target sites that is the intersection of gen1 and gen2
    """
    countuni = 0
    list1, list2 = list(gen1), list(gen2)
    for g1 in list1:
        cutloc_1 = (re.split('[:-]', g1[0])[0], g1[1])
        for g2 in list2:
            cutloc_2 = (re.split('[:-]', g2[0])[0], g2[1])
            if cutloc_1[0] == cutloc_2[0] and abs(cutloc_1[1] - cutloc_2[1]) <= approx:
                countuni += 1
                yield g2
                break
    count1, count2 = len(list1), len(list2)
    print("# of target sites | gen1: %i | gen2: %i | gen1 I gen2: %i" % (count1, count2, countuni))


def gen_union_approx(gen1, gen2, approx=1):
    """ Outputs generator for set of target sites that is union of gen1 and gen2
    :param gen1: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param gen2: same format as gen1
    :param approx: max distance in coordinates between cut sites to be counted as the same one
    :return: generator for set of target sites that is the union of gen1 and gen2
    """
    gen2list = []
    list1, list2 = list(gen1), list(gen2)
    countuni = 0
    for g2 in list2:
        cutloc_i = (re.split('[:-]', g2[0])[0], g2[1])
        gen2list.append(cutloc_i)
        countuni += 1
        yield g2
    for g1 in list1:
        cutloc_1 = (re.split('[:-]', g1[0])[0], g1[1])
        in_2 = False
        for cutloc_2 in gen2list:
            if cutloc_1[0] == cutloc_2[0] and abs(cutloc_1[1] - cutloc_2[1]) <= approx:
                in_2 = True
                break
        if not in_2:
            countuni += 1
            yield g1
    count1, count2 = len(list1), len(list2)
    print("# of target sites | gen1: %i | gen2: %i | gen1 U gen2: %i" % (count1, count2, countuni))


def union_approx(gens, fileout, approx=1000):
    """ For more than one generator of target sites, output the union wrt to the first generator.
    :param gens: list of generators of following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
                  val_i     =   enrichment value         (int)
    :param fileout: path to CSV output, extension omitted.
    :param approx: the same target site between different gen should be within 'approx' bp away.
    """
    out_all = []
    # Use the target sites from the first generator as point of comparison
    for g in gens[0]:
        out_all.append(g)
    out_all = np.asarray(out_all).astype(object)
    out_all[:, -1] = out_all[:, -1].astype(int)
    out_all = out_all[np.argsort(-out_all[:, -1])]
    # Iterate through target sites from remaining generators, match to those of first generator
    for gen in gens[1:]:        # iterate through each file (generator)
        out_i = np.zeros((len(out_all), 1)).astype(object)
        out_i[:] = ''
        for g in gen:           # iterate through each target site in specific generator
            for i, x in enumerate(out_all):     # compare to existing target sites
                # if { chr, cut, sen, pam, gui, mis } match
                if x[0].split(':')[0] == g[0].split(':')[0] and abs(int(x[1]) - g[1]) < approx and \
                        x[2] == g[2] and x[3] == g[3] and x[4] == g[4] and int(x[5]) == g[5]:
                    out_i[i] = g[-1]
                    break
        out_all = np.hstack((out_all, out_i))
    np.savetxt(fileout + ".csv", np.asarray(out_all), fmt='%s', delimiter=',')
