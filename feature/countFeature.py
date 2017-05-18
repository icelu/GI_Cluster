# This script is used to convert feature values into values that can be used for clustering
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

from __future__ import division
import optparse
import numpy as np

import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.parse_sequence import get_rna_segment, get_repeat_segment



def get_gc_value(infile):
    gc_value = {}
    i = 0
    with open(infile, 'rb') as fin:
        for line in fin:
            i += 1
            fields = line.strip().split('\t')
            # Suppose the first column is ID
            gc_value[i] = [float(f) for f in fields[0:]]

    return gc_value


def get_kmer_dist(infile):
    kmer_dict = {}
    i = 0
    with open(infile, 'rb') as fin:
        for line in fin:
            i += 1
            fields = line.strip().split('\t')
            # Suppose the first column is ID
            kmer_dict[i] = [float(f) for f in fields[1:]]

    return kmer_dict


def get_feature_dict_segment(infile, rna_dict={}, repeat_dict={}, gc_value={}, kmer_dict={}):
    '''
    Summarize the features of the genes within a genomic segment
    Input:
    infile -- File containing the segments of a genome
    rna_dict -- A dicionary containing tRNAs around each segment
    repeat_dict -- A dicionary containing repeats around each segment
    gc_value -- A dicionary containing GC distance values (e.g. Chi-square) around each segment
    kmer_dict  -- A dicionary containing k-mer distance values (e.g. covariance) around each segment
    '''
    feature_dict = {}
    head_num = 0
    with open(infile, 'rb') as fin:
        for line in fin:
            head_num += 1
            fields = line.strip().split('\t')
            p1 = (fields[0])
            p2 = (fields[1])
            coord = (p1, p2)
            if '_' in p1:
                mark = p1.index('_')
                start = int(p1[mark + 1:])
                end = int(p2[mark + 1:])
            else:
                start = int(p1)
                end = int(p2)
            size = end - start + 1

            has_trna = 0
            has_repeat = 0
            if coord in rna_dict.keys():
                has_trna = 1
            if coord in repeat_dict.keys():
                has_repeat = 1

            gc = gc_value[head_num][0]
            i = 0
            kmer2 = kmer_dict[head_num][i]
            kmer3 = kmer_dict[head_num][i + 1]
            kmer4 = kmer_dict[head_num][i + 2]
            kmer5 = kmer_dict[head_num][i + 3]
            kmer6 = kmer_dict[head_num][i + 4]
            kmer7 = kmer_dict[head_num][i + 5]
            kmer8 = kmer_dict[head_num][i + 6]

            cluster = (head_num, p1, p2, size, has_trna, has_repeat)
            feature_list = (gc, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8)
            feature_dict[head_num] = (cluster, feature_list)

    return feature_dict



def get_feature_dict_gene(infile, rna_dict={}, repeat_dict={}, gc_value={}, kmer_dict={}):
    '''
    Summarize GI-related features of the genes within a genomic segment
    Input:
    infile -- File containing values of GI-related features for each gene
    rna_dict -- A dicionary containing tRNAs around each segment
    repeat_dict -- A dicionary containing repeats around each segment
    gc_value -- A dicionary containing GC distance values (e.g. Chi-square) around each segment
    kmer_dict  -- A dicionary containing k-mer distance values (e.g. covariance) around each segment
    '''
    head_num = 0    # The number of regions start with '>'
    feature_dict = {}
    cluster_set = []

    # The number of each feature in a genomic region
    gc = 0
    gc1 = 0
    gc2 = 0
    gc3 = 0
    gc_skew = 0

    cub = 0
    aab = 0
    chi = 0
    cai = 0
    cbi = 0
    fop = 0
    nc = 0

    kmer2 = 0
    kmer3 = 0
    kmer4 = 0
    kmer5 = 0
    kmer6 = 0
    kmer7 = 0
    kmer8 = 0

    mbt = 0
    phaget = 0
    mbpht = 0

    vft = 0
    art = 0
    vfart = 0

    ngt = 0
    vfarngt = 0

    # The percentage of each feature
    gcp = 0
    gc1p = 0
    gc2p = 0
    gc3p = 0
    gc_skewp = 0

    mbp = 0
    phagep = 0
    mbphp = 0

    vfp = 0
    arp = 0
    vfarp = 0

    ngp = 0
    vfarngp = 0

    cdst = 0
    cdsp = 0  # gene density, number of genes per kb
    rnat = 0
    rnap = 0

    with open(infile, 'rb') as fin:
        # Skip 1st line, which is header
        firstline = fin.readline()
        line = fin.readline()
        while line:
            if line.startswith('>'):
                head_num += 1
                fields = line.strip().split('\t')
                start = (fields[1])
                end = (fields[2])
                coord = (start, end)

                has_trna = 0
                has_repeat = 0
                if coord in rna_dict.keys():
                    has_trna = 1
                if coord in repeat_dict.keys():
                    has_repeat = 1

                if '_' in end:
                    size = int(end[end.index("_") + 1:]) - int(start[start.index("_") + 1:]) + 1
                else:
                    size = int(end) - int(start) + 1
                cluster = (head_num, start, end, size, has_trna, has_repeat)
                cluster_set.append(cluster)

                feature_list = [0] * 25
                feature_dict[head_num] = (cluster, tuple(feature_list))

                line = fin.readline()
            else:
                gc_list = []
                gc1_list = []
                gc2_list = []
                gc3_list = []
                gc_skew_list = []

                cub_list = []
                aab_list = []
                chi_list = []
                cai_list = []
                cbi_list = []
                fop_list = []
                nc_list = []

                gene_list = []
                while line and (not line.startswith('>')):
                    if (not line.startswith('r')):
                        cdst += 1
                        fields = line.strip().split('\t')

                        gstart = int(fields[1])
                        gend = int(fields[2])
                        gene_list.append((int(gstart), int(gend)))

                        i = 4
                        gc = float(fields[i])
                        gc1 = float(fields[i + 1])
                        gc2 = float(fields[i + 2])
                        gc3 = float(fields[i + 3])
                        gc_list.append(gc)
                        gc1_list.append(gc1)
                        gc2_list.append(gc2)
                        gc3_list.append(gc3)

                        i = 9
                        cub_list.append(float(fields[i]))
                        aab_list.append(float(fields[i + 1]))
                        chi_list.append(float(fields[i + 2]))
                        cai_list.append(float(fields[i + 3]))
                        cbi_list.append(float(fields[i + 4]))
                        fop_list.append(float(fields[i + 5]))

                        mb = fields[-5]
                        if 'No hits' not in mb:
                            mbt += 1
                        phage = fields[-4]
                        if 'No hits' not in phage:
                            phaget += 1
                        if 'No hits' not in mb or 'No hits' not in phage:
                            mbpht += 1
                        vf = fields[-3]
                        if 'No hits' not in vf:
                            vft += 1
                        ar = fields[-2]
                        if 'No hits' not in ar:
                            art += 1
                        if 'No hits' not in vf or 'No hits' not in ar:
                            vfart += 1
                        ng = fields[-1]
                        if 'No hits' in ng:
                            ngt += 1
                        if 'No hits' not in vf or 'No hits' not in ar or 'No hits' in ng:
                            vfarngt += 1

                    if (not line.startswith('>')) and (line.startswith('r')):
                        rfields = line.strip().split('\t')
                        rstart = rfields[1]
                        rend = rfields[2]
                        rfunc = rfields[4]
                        if 'rRNA' in rfunc or 'tRNA' in rfunc:
                            rnat += 1
                        gene_list.append((int(rstart), int(rend)))
                    line = fin.readline()

                # Summarize the features of a segment
                if cdst > 0:
                    gc = np.mean(gc_list)
                    gc1 = np.mean(gc1_list)
                    gc2 = np.mean(gc2_list)
                    gc3 = np.mean(gc3_list)
                    # gc_skew = np.mean(gc_skew_list)
                    # gc = gc_value[head_num][0]
                    # gc1 = gc_value[head_num][1]
                    # gc2 = gc_value[head_num][2]
                    # gc3 = gc_value[head_num][3]

                    cub = np.mean(cub_list)
                    aab = np.mean(aab_list)
                    chi = np.mean(chi_list)
                    cai = np.mean(cai_list)
                    cbi = np.mean(cbi_list)
                    fop = np.mean(fop_list)
                    #nc = np.mean(nc_list)

                    i = 0
                    kmer2 = kmer_dict[head_num][i]
                    kmer3 = kmer_dict[head_num][i + 1]
                    kmer4 = kmer_dict[head_num][i + 2]
                    kmer5 = kmer_dict[head_num][i + 3]
                    kmer6 = kmer_dict[head_num][i + 4]
                    kmer7 = kmer_dict[head_num][i + 5]
                    kmer8 = kmer_dict[head_num][i + 6]

                    mbp = mbt / cdst
                    phagep = phaget / cdst
                    mbphp = mbpht / cdst

                    vfp = vft / cdst
                    arp = art / cdst
                    vfarp = vfart / cdst

                    ngp = ngt / cdst
                    vfarngp = vfarngt / cdst

                    cdsp = (cdst + rnat) * 1000 / size
                    #cdsp = (cdst)*1000 / size
                if rnat > 0:
                    rnap = rnat / (cdst + rnat)
                if cdst > 1:
                    sorted_gene_list = sorted(
                        gene_list, key=lambda x: (int(x[0]), int(x[1])))
                    dist = 0
                    prev_gend = sorted_gene_list[0][1]
                    for gstart, gend in sorted_gene_list[1:]:
                        idist = gstart - prev_gend
                        if idist < 0:
                            idist = 0
                        dist += idist
                        prev_gend = gend
                    inter_gene_dist = dist / (cdst + rnat - 1)
                feature_dict[head_num] = (cluster, (cdst, rnat, gc, gc1, gc3, cub, aab, chi, cai, cbi, fop, kmer2,
                                                    kmer3, kmer4, kmer5, kmer6, kmer7, kmer8, mbp, phagep, vfp, arp, ngp, cdsp, inter_gene_dist))
                # reinitialize the variables after each loop
                # number of each feature
                mbt = 0
                phaget = 0
                mbpht = 0

                vft = 0
                art = 0
                vfart = 0

                ngt = 0
                vfarngt = 0

                # percentage of each feature
                mbp = 0
                phagep = 0
                mbphp = 0

                vfp = 0
                arp = 0
                vfarp = 0

                ngp = 0
                vfarngp = 0

                cdst = 0
                cdsp = 0
                rnat = 0
                rnap = 0

    return feature_dict


def write_feature_dict(output, feature_dict):
    with open(output, 'w') as fout:
        for cluster, feature in feature_dict.values():
            content = cluster + feature
            format = '%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d'
            for i in range(2, len(feature)):
                format += '\t%.3f'
            format += '\n'
            line = format % content
            fout.write(line)


def write_feature_dict_segment(output, feature_dict):
    with open(output, 'w') as fout:
        for cluster, feature in feature_dict.values():
            content = cluster + feature
            format = '%d\t%s\t%s\t%d\t%d\t%d'
            for i in range(0, len(feature)):
                format += '\t%.3f'
            format += '\n'
            line = format % content
            fout.write(line)


if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-i", "--input", dest="input",
                      help="input file of extracted features")
    parser.add_option("-o", "--output", dest="output",
                      help="output file of the number of each feature in a genomic interval")
    parser.add_option("-r", "--repeats", dest="repeats",
                      help="input file containing repeats overlapping with the input regions")
    parser.add_option("-g", "--gcfile", dest="gcfile",
                      help="input file containing gc content in a genomic interval")
    parser.add_option("-k", "--kmerfile", dest="kmerfile",
                      help="input file containing kmer distribution in a genomic interval")
    parser.add_option("-t", "--trnas", dest="trnas",
                      help="input file containing trnas overlapping with the input regions")
    parser.add_option("-a", "--has_gene", dest="has_gene", action='store_true', default=False,
                      help="The gene predictions are available")
    options, args = parser.parse_args()

    rna_dict = get_rna_segment(options.trnas)
    repeat_dict = get_repeat_segment(options.repeats)
    kmer_dict = get_kmer_dist(options.kmerfile)

    if options.has_gene:
        feature_dict = get_feature_dict_gene(
            options.input, rna_dict=rna_dict, repeat_dict=repeat_dict, gc_value={}, kmer_dict=kmer_dict)
        write_feature_dict(options.output, feature_dict)
    else:
        gc_value = get_gc_value(options.gcfile) # Only compute when gene predictions are not available
        feature_dict = get_feature_dict_segment(options.input, rna_dict=rna_dict, repeat_dict=repeat_dict, gc_value=gc_value, kmer_dict=kmer_dict)
        write_feature_dict_segment(options.output, feature_dict)
