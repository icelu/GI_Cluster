# Given DNA sequence of a genome or a gene or a genomic region, analyze its GC content
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# The fasta file of DNA sequence
#
# Output:
# Measures of GC content:
# GC content, GC skew, GC(k) content


from __future__ import division
import itertools
import optparse
import random
from scipy.stats import chisquare
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.parse_sequence import switch_to_AGCT, standardize_DNASeq, isheader


def parse_genes(gene_file, outfile):
    '''
    input: .ffn file, fasta format, DNA sequence for each gene
    output: GC measures for each gene in a line
    '''
    i = 0
    with open(gene_file, 'rb') as fin, open(outfile, 'w') as fout:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                gc = GC(sequence)
                gc1, gc2, gc3 = GC123(sequence)
                gc_skew = GC_skew(sequence)
                line = '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
                    gc, gc1, gc2, gc3, gc_skew)
                fout.write(line)


def parse_segs(pfile, gnome, outfile):
    '''
    input: position for each segment (1-based), genome sequence (standardized)
    output: GC measures for each segment in a line
    '''
    i = 0
    with open(pfile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            sequence = gnome[start - 1:end]
            # sequence = sequence.upper()
            # sequence = standardize_DNASeq(sequence)
            gc = GC(sequence)
            # gc1, gc2, gc3 = GC123(sequence)
            # gc_skew = GC_skew(sequence)
            # line = '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
            #     gc, gc1, gc2, gc3, gc_skew)
            line = '%.3f\n' % (gc)
            fout.write(line)            
            fout.write(line)


def parse_segs_dist(pfile, gnome):
    '''
    input: position for each segment (1-based), genome sequence (standardized)
    output: GC measures for each segment in a line -- not meaningful to compute GC123 without knowing codons
    '''
    i = 0
    gc_genome = GC(gnome)
    gc1_genome, gc2_genome, gc3_genome = GC123(gnome)
    chi_dict = {}
    for i in range(4):
        chi_dict[i] = []
    with open(pfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            sequence = gnome[start - 1:end]
            # sequence = sequence.upper()
            # sequence = standardize_DNASeq(sequence)
            gc = GC(sequence)
            # gc1, gc2, gc3 = GC123(sequence)

            # multiple each number by 100 to avoid small value
            chi = chi_square(gc * 100, gc_genome * 100)
            chi_dict.setdefault(0, []).append(chi)
    return chi_dict


def write_gc_chi(chi_dict, outfile):
    with open(outfile, 'w') as fout:
        for i in range(len(chi_dict[0])):
            dists = []
            for k in range(1):
                dists.append(chi_dict[k][i])
            line = '%d' % i
            fout.write(line)
            for d in dists:
                line = "\t%0.3f" % d
                fout.write(line)
            fout.write('\n')


# Read the genome sequence
def get_genome(gene_file):
    gnome = ''
    with open(gene_file, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                gnome += sequence
    return gnome


def parse_genome(gene_file, outfile):
    '''
    input: .fna file, fasta format, DNA sequence for the whole genome
    output: GC measures in a line
    '''
    i = 0
    with open(gene_file, 'rb') as fin, open(outfile, 'w') as fout:
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
                line = '%.2f\t%.2f\n' % (
                    gc, gc_skew)
                fout.write(line)


############################# compute GC measures ########################

def GC(dna):
    gc = dna.count('G') + dna.count('C')
    return gc / (len(dna))


def GC123(dna):
    '''
    Similar as in BioPython
    http://biopython.org/DIST/docs/api/Bio.SeqUtils-pysrc.html
    '''
    d = {}
    for nt in ['A', 'T', 'G', 'C']:
        d[nt] = [0, 0, 0]

    remainder = len(dna) % 3
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


######################### compute GC content for contigs #################
def get_contigs(infile):
    '''
    infile -- Input file containing the sequence of contigs
    '''
    contig_sequence = {}
    i = 0
    with open(infile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                contig_sequence[i] = sequence

    return contig_sequence


def parse_segs_contigs(pfile, contig_sequence, outfile):
    '''
    input: pfile -- position for each segment (1-based);  contig_sequence -- contig sequence (standardized)
    output: GC measures for each segment in a line
    '''
    i = 0
    with open(pfile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            p1 = fields[0]
            mark = p1.index('_')
            start = int(p1[mark + 1:])
            p2 = fields[1]
            end = int(p2[mark + 1:])
            contig_id = int(p1[0:mark])
            contig = contig_sequence[contig_id]
            sequence = contig[start - 1:end]
            gc = GC(sequence)
            line = '%.3f\n' % (gc)
            fout.write(line)


####################################################################
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-g", "--gene_file", dest="gene_file",
                      help="input gene file in fasta format")
    parser.add_option("-i", "--genome_file", dest="genome_file",
                      help="input genome file in fasta format")
    parser.add_option("-r", "--is_seg", dest="is_seg", action='store_true', default=False,
                      help="Analyze GC content in genomic segments rather than genes")
    parser.add_option("-m", "--is_metric", dest="is_metric", action='store_true', default=False,
                      help="Compute distance metrics")
    parser.add_option("-s", "--seg_file", dest="seg_file",
                      help="input file of segment position")
    parser.add_option("-n", "--is_genome", dest="is_genome", action='store_true', default=False,
                      help="input gene files in fasta format")
    parser.add_option("-o", "--output", dest="output",
                      help="output file of labeled intervals")
    parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true', default=False,
                      help="Analyze contigs from unassembled genomes")
    options, args = parser.parse_args()

    if options.is_metric:   # compute chisquare values
        gnome = get_genome(options.genome_file)
        chi_dict = parse_segs_dist(options.seg_file, gnome)
        write_gc_chi(chi_dict, options.output)
    else:
        if options.is_genome:
            parse_genome(options.genome_file, options.output)
        elif options.is_contig:
            contig_sequence = get_contigs(options.genome_file)
            # print contig_sequence
            parse_segs_contigs(options.seg_file, contig_sequence, options.output)
        else:
            if options.is_seg:
                gnome = get_genome(options.genome_file)
                parse_segs(options.seg_file, gnome, options.output)
            else:
                parse_genes(options.gene_file, options.output)
