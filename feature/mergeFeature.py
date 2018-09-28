#!/usr/bin/env python

# This script is used to find the features for each genomic segment
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# A list of segments
# A list of features for each gene/ORF
#
# Output:
# Group the genes by segments
#     format:
#     region position
#     genes inside this region along with their features
#
# Procedure:
# Store the features for each gene/ORF in a dicionary, and then write out the values in a file
#

import os
import optparse

import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.parse_sequence import get_contig_IDs, get_contig_gene_IDs


def get_gene_feature(infile):
    '''
    Note: the input genes are sorted and kept in this order
    Use a list of tuple to store a gene along with its features: (gid, gfeature)
    '''
    features = []
    with open(infile, 'r') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            id = fields[0]  # Sequential number from 1 to n
            feature = fields[1:]
            item = (id, feature)
            features.append(item)

    return features


# Find the positions of ncRNAs
def parse_cmscan(infile):
    hits = []
    with open(infile, 'r') as fin:
        for line in fin:
            if not line.startswith('#'):
                fields = line.strip().split()
                name = fields[1]
                accession = fields[2]
                start = int(fields[9])
                end = int(fields[10])
                strand = fields[11]
                if strand == '-':
                    tmp = start
                    start = end
                    end = tmp
                    strand = 'R'
                else:
                    strand = 'F'
                score = fields[16]
                evalue = fields[17]
                hit = (start, end, strand, name, accession, score, evalue)
                hits.append(hit)
    hits = sorted(hits)
    return hits



# Find the positions of ncRNAs when searching multiple sequences (contigs)
def parse_cmscan_multiple(infile):
    hits = []
    with open(infile, 'r') as fin:
        for line in fin:
            if not line.startswith('#'):
                fields = line.strip().split()
                name = fields[1]    # Name of RNA
                accession = fields[2]
                query = fields[3]   # Name of query sequence (different for different contigs)
                start = int(fields[9])
                end = int(fields[10])
                strand = fields[11]
                if strand == '-':
                    tmp = start
                    start = end
                    end = tmp
                    strand = 'R'
                else:
                    strand = 'F'
                score = fields[16]
                evalue = fields[17]
                hit = (start, end, strand, name, accession, score, evalue, query)
                hits.append(hit)
    hits = sorted(hits)
    return hits


# Output merged features for a single complete genome
def group_genes(infile, features, rnas, outfile):
    '''
    infile -- Input file containing segments of a genome. Each line contains the start and end positions of a segment
    features -- A list containing GI-related featues for each gene/ORF
    rnas -- A list of predicted ncRNAs
    outfile -- Output file containing ncRNAs and genes overlapping with each segment along with their featues
    '''
    fout = open(outfile, 'w')
    header = 'GeneID\tStart\tEnd\tStrand\tGC\tGC1\tGC2\tGC3\tGC_skew\tCUB\tAAB\tCHI\tCAI\tCBI\tFop\t2-mer\t3-mer\t4-mer\t5-mer\t6-mer\t7-mer\t8-mer\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\n'
    fout.write(header)
    # For each region, find the genes inside it
    i = 0
    with open(infile, 'r') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            i += 1
            oline = '>Region' + str(i) + '\t' + line
            fout.write(oline)
            # print 'start %d, end %d' % (start, end)
            ri = 0
            for hit in rnas:
                (rstart, rend, strand, name, accession, score, evalue) = hit
                # if ((rstart >= start and rstart <= end) or ( rend >= start
                # and rend <= end)):
                if not ((rstart > end) or (rend < start)):
                    # rmid = (rend + rstart) / 2
                    # if ((rstart >= start and rmid <= end) or (rmid >= start
                    # and rend <= end)):
                    ri += 1
                    id = 'r' + str(ri)
                    oline = id + '\t'
                    oline = oline + \
                        '\t'.join((str(rstart), str(rend), strand,
                                   name, accession, score, evalue)) + '\n'
                    fout.write(oline)
                    # a gene overlapping with two regions will be counted twice, since it is nontrivial to determine which region to include this gene
                    # rnas.remove(hit)

            for gid, gfeature in features:
                gstart = int(gfeature[0])
                gend = int(gfeature[1])
                gmid = (gend + gstart) / 2
                # print 'gstart %d, gend %d' % (gstart, gend)
                # genes strictly inside a region
                # if (gstart >= start and gend <= end):
                # include genes across the boundary (even with only 1 bp overlap), so one gene may be considerred in multiple segements
                # if ((gstart >= start and gstart <= end) or ( gend >= start and gend <= end)):
                if not ((gstart > end) or (gend < start)):
                    # since segments are non-overlapping, if a gene is across >1 segments, if half of the gene is inside a segment, this gene can be seen as inside this segment
                    # if ((gstart >= start and gmid <= end) or (gmid >= start and gend <= end)):
                    # write the output
                    items = [gid]
                    for v in gfeature:
                        items.append(v)
                    # print items
                    oline = '\t'.join(items) + '\n'
                    fout.write(oline)

                    # remove this item from the list since it will not be included in other regions
                    # a gene overlapping with two regions will be counted twice, since it is nontrivial to determine which region to include this gene
                    # features.remove((gid, gfeature))

    fout.close()


# Output merged features for multiple sequences
def group_genes_multiple(infile, features, rnas, outfile, contig_id_mapping, gene_id_mapping):
    '''
    infile -- Input file containing segments of a genome. Each line contains the start and end positions of a segment
    features -- A list containing GI-related featues for each gene/ORF
    rnas -- A list of predicted ncRNAs
    outfile -- Output file containing ncRNAs and genes overlapping with each segment along with their featues
    '''
    fout = open(outfile, 'w')
    header = 'GeneID\tStart\tEnd\tStrand\tGC\tGC1\tGC2\tGC3\tGC_skew\tCUB\tAAB\tCHI\tCAI\tCBI\tFop\t2-mer\t3-mer\t4-mer\t5-mer\t6-mer\t7-mer\t8-mer\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\n'
    fout.write(header)
    # For each region, find the genes inside it
    i = 0
    with open(infile, 'r') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            p1 = fields[0]
            mark = p1.index('_')
            start = int(p1[mark + 1:])
            p2 = fields[1]
            end = int(p2[mark + 1:])
            contig_id = int(p1[0:mark])
            coord = (start, end)

            i += 1
            oline = '>Region' + str(i) + '\t' + line
            fout.write(oline)
            # print 'start %d, end %d' % (start, end)
            ri = 0
            for hit in rnas:
                (rstart, rend, strand, name, accession, score, evalue, query) = hit
                # Get the ID of the query sequence via id_mapping
                query_ID = contig_id_mapping[query]
                if query_ID == contig_id and (not ((rstart > end) or (rend < start))):
                    ri += 1
                    id = 'r' + str(ri)
                    oline = id + '\t'
                    oline = oline + \
                        '\t'.join((str(rstart), str(rend), strand,
                                   name, accession, score, evalue)) + '\n'
                    fout.write(oline)

            for gid, gfeature in features:
                # Find the contig ID for the gene
                gene_name = gene_id_mapping[int(gid)]
                # Get contig name from gene name
                contig_name = gene_name[0:gene_name.rfind('_')]
                cid = contig_id_mapping[contig_name]
                if contig_id != cid:
                    continue
                gstart = int(gfeature[0])
                gend = int(gfeature[1])
                gmid = (gend + gstart) / 2
                if not ((gstart > end) or (gend < start)):
                    items = [gid]
                    for v in gfeature:
                        items.append(v)
                    # print items
                    oline = '\t'.join(items) + '\n'
                    fout.write(oline)

    fout.close()



if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-g", "--gene_file", dest="gene_file",
                      help="input file of genes and their features")
    parser.add_option("-s", "--seg_file", dest="seg_file",
                      help="input file of genomic segments")
    parser.add_option("-r", "--rna_file", dest="rna_file",
                      help="input file of predicted RNAs")
    parser.add_option("-m", "--genome_file", dest="genome_file",
                      help="input genome file in fasta format")
    parser.add_option("-d", "--gid_file", dest="gid_file",
                      help="input file of predicted gene IDs")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file")
    parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true', default=False,
                      help="Analyze contigs from unassembled genomes")
    options, args = parser.parse_args()

    features = get_gene_feature(options.gene_file)
    if options.is_contig:
        contig_id_mapping = get_contig_IDs(options.genome_file)
        gene_id_mapping = get_contig_gene_IDs(options.gid_file)
        rnas = parse_cmscan_multiple(options.rna_file)
        group_genes_multiple(options.seg_file, features, rnas, options.outfile, contig_id_mapping, gene_id_mapping)
    else:
        rnas = parse_cmscan(options.rna_file)
        group_genes(options.seg_file, features, rnas, options.outfile)
