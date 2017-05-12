#!/usr/bin/python
# Lu Bingxin, 2017.3.1

'''
Given DNA sequence of a genome or a gene or a genomic region, analyze its GC content

Input:
1. fasta file

Output:
Measures of GC content:
GC content, GC skew, GC(k) content

'''

from __future__ import division
import itertools
import optparse
import random
from scipy.stats import chisquare


# Switch non-standard character to 'AGCT'
def switch_to_AGCT(Nuc):
    SWITCH = {'R': lambda: random.choice('AG'), \
           'Y': lambda: random.choice('CT'), \
           'M': lambda: random.choice('AC'), \
           'K': lambda: random.choice('GT'), \
           'S': lambda: random.choice('GC'), \
           'W': lambda: random.choice('AT'), \
           'H': lambda: random.choice('ATC'), \
           'B': lambda: random.choice('GTC'), \
           'V': lambda: random.choice('GAC'), \
           'D': lambda: random.choice('GAT'), \
           'N': lambda: random.choice('AGCT')}
    return SWITCH[Nuc]()


'''
Replace ambiguous bases with ATCG
'''
def standardize_DNASeq(genome):
    for i in ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']:
        if i in genome:
           genome = genome.replace(i, switch_to_AGCT(i))
    return genome


def isheader(line):
    return line[0] == '>'

'''
input: .ffn file, fasta format, DNA sequence for each gene
output: GC measures for each gene in a line
'''
def parse_genes(gfile, outfile):
    i = 0
    with open(gfile, 'rb') as fin, open(outfile, 'w') as fout:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
                #line = group.next()
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                gc = GC(sequence)
                gc1, gc2, gc3 = GC123(sequence)
                gc_skew = GC_skew(sequence)
                # line = '%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
                #     i, gc, gc1, gc2, gc3, gc_skew)
                line = '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
                    gc, gc1, gc2, gc3, gc_skew)
                fout.write(line)



'''
input: position for each segment (1-based), genome sequence (standardized)
output: GC measures for each segment in a line
'''
def parse_segs(pfile, gnome, outfile):
    i = 0
    with open(pfile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            sequence = gnome[start-1:end]
            # sequence = sequence.upper()
            # sequence = standardize_DNASeq(sequence)
            gc = GC(sequence)
            gc1, gc2, gc3 = GC123(sequence)
            gc_skew = GC_skew(sequence)
            # line = '%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
            #     i, gc, gc1, gc2, gc3, gc_skew)
            line = '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
                gc, gc1, gc2, gc3, gc_skew)
            fout.write(line)


'''
input: position for each segment (1-based), genome sequence (standardized)
output: GC measures for each segment in a line -- not meaningful to compute GC123 without knowing codons
'''
def parse_segs_dist(pfile, gnome):
    i = 0
    gc_genome = GC(gnome)
    gc1_genome, gc2_genome, gc3_genome = GC123(gnome)
    chi_dict = {}
    for i in range(4):
        chi_dict[i]=[]
    with open(pfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            sequence = gnome[start-1:end]
            # sequence = sequence.upper()
            # sequence = standardize_DNASeq(sequence)
            gc = GC(sequence)
            # gc1, gc2, gc3 = GC123(sequence)

            # multiple each number by 100 to avoid small value
            chi = chi_square(gc * 100, gc_genome * 100)
            chi_dict.setdefault(0,[]).append(chi)
            chi1 = chi_square(gc1 * 100, gc1_genome * 100)
            # chi_dict.setdefault(1,[]).append(chi1)
            # chi2 = chi_square(gc2 * 100, gc2_genome * 100)
            # chi_dict.setdefault(2,[]).append(chi2)
            # chi3 = chi_square(gc3 * 100, gc3_genome * 100)
            # chi_dict.setdefault(3,[]).append(chi3)
    # print chi_dict
    return chi_dict

def write_gc_chi(chi_dict, outfile):
    # new_dist_dict = {}
    # for key, values in chi_dict.items():
    #     maxv = max(values)
    #     minv = min(values)
    #     new_values = []
    #     for v in values:
    #         nv = (v-minv)/(maxv-minv)
    #         new_values.append(nv)
    #     new_dist_dict[key] = new_values
    # print new_dist_dict
    with open(outfile, 'w') as fout:
        for i in range(len(chi_dict[0])):
            dists=[]
            for k in range(1):
                dists.append(chi_dict[k][i])
            line = '%d' % i
            fout.write(line)
            for d in dists:
                line = "\t%0.3f" % d
                fout.write(line)
            fout.write('\n')



def get_genome(gfile):
    gnome = ''
    with open(gfile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                gnome += sequence
    return gnome


'''
input: .fna file, fasta format, DNA sequence for the whole genome
output: GC measures in a line
'''
def parse_genome(gfile, outfile):
    i = 0
    with open(gfile, 'rb') as fin, open(outfile, 'w') as fout:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
                #line = group.next()
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                gc = GC(sequence) * 100
                gc_skew = GC_skew(sequence) * 100
                # line = '%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
                #     i, gc, gc1, gc2, gc3, gc_skew)
                line = '%.2f\t%.2f\n' % (
                    gc, gc_skew)
                fout.write(line)


############################# compute GC measures ##############################

def GC(dna):
    gc = dna.count('G') + dna.count('C')
    return gc / (len(dna))


'''
Similar as in BioPython
http://biopython.org/DIST/docs/api/Bio.SeqUtils-pysrc.html
'''


def GC123(dna):
    d = {}
    for nt in ['A', 'T', 'G', 'C']:
        d[nt] = [0, 0, 0]

    remainder = len(dna)%3
    # exclude partial codons
    if remainder != 0:
        dna = dna[0:-remainder]
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        if len(codon) < 3:
            next
        for pos in range(0, 3):
            for nt in ['A', 'T', 'G', 'C']:
                if codon[pos] == nt:
                    d[nt][pos] += 1

    gc = {}
    for i in range(0, 3):
        try:
            n = d['G'][i] + d['C'][i] + d['T'][i] + d['A'][i]
            gc[i] = (d['G'][i] + d['C'][i]) / n
        except Exception:
            gc[i] = 0

    return gc[0], gc[1], gc[2]



def GC_skew(dna):
    c_count = dna.count('C')
    g_count = dna.count('G')
    gc_skew = (g_count - c_count) / (g_count + c_count)
    return gc_skew


################################# compute distances ################
def chi_square(val, avg):
    res = chisquare(val, f_exp=avg)
    return res[0]



if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-g", "--gfile", dest="gfile",
                      help="input gene file in fasta format")
    parser.add_option("-i", "--genome_file", dest="genome_file",
                      help="input genome file in fasta format")
    parser.add_option("-r", "--is_seg", dest="is_seg", action='store_true', default=False,
                      help="Analyze GC content in genomic segments rather than genes")
    parser.add_option("-m", "--is_metric", dest="is_metric", action='store_true', default=False,
                      help="Compute distance metrics")
    parser.add_option("-s", "--segfile", dest="segfile",
                      help="input file of segment position")
    parser.add_option("-n", "--is_genome", dest="is_genome", action='store_true', default=False,
                      help="input gene files in fasta format")
    parser.add_option("-o", "--output", dest="output",
                      help="output file of labeled intervals")

    options, args = parser.parse_args()
    if options.is_metric:
        gnome = get_genome(options.genome_file)
        chi_dict = parse_segs_dist(options.segfile, gnome)
        write_gc_chi(chi_dict, options.output)
    else:
        if options.is_genome:
            parse_genome(options.genome_file, options.output)
        else:
            if options.is_seg:
                gnome = get_genome(options.genome_file)
                parse_segs(options.segfile, gnome, options.output)
            else:
                parse_genes(options.gfile, options.output)
