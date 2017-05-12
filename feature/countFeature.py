#!/usr/bin/python

# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
# Version : 1.0

from __future__ import division
import optparse
import numpy as np


def get_rna(infile):
    rna_dict = {}
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[0]), int(fields[1]))
            num = int(fields[2])
            rnas = fields[3].split(';')
            rna_pos = []
            for rp in rnas:
                items = rp.split(',')
                #print items
                # remove brackets
                rpos = (int(items[0][1:]), int(items[1][:-1]))
                rna_pos.append(rpos)
            rna_dict[coord] = rna_pos
        #print rna_dict
    return rna_dict


def get_repeat(infile):
    repeat_dict = {}
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[0]), int(fields[1]))
            num = int(fields[2])
            repeats = fields[3].split(';')
            # keep only the position
            repeats_pos = []
            for rp in repeats:
                items = rp.split(',')
                #print items
                rpos = (int(items[1]), int(items[2]))
                repeats_pos.append(rpos)
            repeat_dict[coord] = repeats_pos
    #print repeat_dict
    return repeat_dict


def get_gc_dist(infile):
    gc_dict = {}
    i=0
    with open(infile, 'rb') as fin:
        for line in fin:
            i+=1
            fields = line.strip().split('\t')
            # 1st column is ID
            gc_dict[i] = [float(f) for f in fields[1:]]

    return gc_dict


def get_kmer_dist(infile):
    kmer_dict = {}
    i=0
    with open(infile, 'rb') as fin:
        for line in fin:
            i+=1
            fields = line.strip().split('\t')
            kmer_dict[i] = [float(f) for f in fields[1:]]

    return kmer_dict

'''
Summarize the features of the genes within a genomic segment
Use the indicator value for sequence compositional metrics (atypical outside 1.5 sd)
'''
def get_feature_dict_group(infile, rna_dict={}, repeat_dict={}, gc_dict={}, kmer_dict={}):
    head_num = 0
    feature_dict = {}
    cluster_set = []

    # number of each feature
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

    # percentage of each feature
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
        # skip 1st line, which is header
        firstline = fin.readline()
        line = fin.readline()
        while line:
            if line.startswith('>'):
                head_num += 1
                fields = line.strip().split('\t')
                # print 'fields in header: %s' % fields

                start = int(fields[1])
                end = int(fields[2])
                coord = (start, end)
                # has_bi = 0 # has boundary indicator: tRNA or short repeats
                has_trna = 0
                has_repeat = 0
                if coord in rna_dict.keys():
                    has_trna = 1
                if coord in repeat_dict.keys():
                    has_repeat = 1

                size = end - start + 1
                cluster = (head_num, start, end, size, has_trna, has_repeat)
                cluster_set.append(cluster)
                # number of features:  1 rna, 5 gc, 3 codon, 7 kmer, 5 function -> 2 function (combine mob and phage, combine VF and AR and novel gene), 2 gene density
                # 5  + 4 + 7 + 5 + 3 = 24 -> 5  + 3 + 7 + 2 + 3 = 20
                feature_list = [0]* 25
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

                # kmer2_list = []
                # kmer3_list = []
                # kmer4_list = []
                # kmer5_list = []
                # kmer6_list = []
                # kmer7_list = []
                # kmer8_list = []

                gene_list = []
                while line and (not line.startswith('>')):
                    if (not line.startswith('r')):
                        cdst += 1
                        fields = line.strip().split('\t')
                        # print 'fields',fields

                        gstart = int(fields[1])
                        gend = int(fields[2])
                        gene_list.append((int(gstart), int(gend)))

                        i=4
                        gc = float(fields[i])
                        gc1 = float(fields[i+1])
                        gc2 = float(fields[i+2])
                        gc3 = float(fields[i+3])
                        gc_list.append(gc)
                        gc1_list.append(gc1)
                        gc2_list.append(gc2)
                        gc3_list.append(gc3)
                        # if 'NORMAL' not in fields[i]: gc += 1
                        # if 'NORMAL' not in fields[i+1]: gc1 += 1
                        # if 'NORMAL' not in fields[i+2]: gc2 += 1
                        # if 'NORMAL' not in fields[i+3]: gc3 += 1
                        # # if 'NORMAL' not in fields[8]: gc_skew += 1
                        # gc_skew_list.append(float(fields[i+4]))
                        # if '*' not in fields[9]:
                        i=9
                        cub_list.append(float(fields[i]))
                        aab_list.append(float(fields[i+1]))
                        chi_list.append(float(fields[i+2]))
                        cai_list.append(float(fields[i+3]))
                        cbi_list.append(float(fields[i+4]))
                        fop_list.append(float(fields[i+5]))

                        # i=15
                        # kmer2_list.append(float(fields[i]))
                        # kmer3_list.append(float(fields[i+1]))
                        # kmer4_list.append(float(fields[i+2]))
                        # kmer5_list.append(float(fields[i+3]))
                        # kmer6_list.append(float(fields[i+4]))
                        # kmer7_list.append(float(fields[i+5]))
                        # kmer8_list.append(float(fields[i+6]))

                        mb = fields[-5]
                        if 'No hits' not in mb: mbt += 1
                        phage = fields[-4]
                        if 'No hits' not in phage: phaget += 1
                        if 'No hits' not in mb or 'No hits' not in phage: mbpht += 1
                        vf = fields[-3]
                        if 'No hits' not in vf: vft += 1
                        ar = fields[-2]
                        if 'No hits' not in ar: art += 1
                        if 'No hits' not in vf or 'No hits' not in ar: vfart += 1
                        ng = fields[-1]
                        if 'No hits' in ng: ngt += 1
                        if 'No hits' not in vf or 'No hits' not in ar or 'No hits' in ng: vfarngt += 1

                    if (not line.startswith('>')) and (line.startswith('r')):

                        rfields = line.strip().split('\t')
                        rstart = rfields[1]
                        rend = rfields[2]
                        rfunc = rfields[4]
                        if 'rRNA' in rfunc or 'tRNA' in rfunc:
                            rnat += 1
                        gene_list.append((int(rstart), int(rend)))
                    line = fin.readline()

                # summarize the features of a segment
                if cdst > 0:
                    gc = np.mean(gc_list)
                    gc1 = np.mean(gc1_list)
                    gc2 = np.mean(gc2_list)
                    gc3 = np.mean(gc3_list)
                    # gc_skew = np.mean(gc_skew_list)
                    # gc = gc_dict[head_num][0]
                    # gc1 = gc_dict[head_num][1]
                    # gc2 = gc_dict[head_num][2]
                    # gc3 = gc_dict[head_num][3]

                    cub = np.mean(cub_list)
                    aab = np.mean(aab_list)
                    chi = np.mean(chi_list)
                    cai = np.mean(cai_list)
                    cbi = np.mean(cbi_list)
                    fop = np.mean(fop_list)
                    #nc = np.mean(nc_list)

                    # kmer2 = np.mean(kmer2_list)
                    # kmer3 = np.mean(kmer3_list)
                    # kmer4 = np.mean(kmer4_list)
                    # kmer5 = np.mean(kmer5_list)
                    # kmer6 = np.mean(kmer6_list)
                    # kmer7 = np.mean(kmer7_list)
                    # kmer8 = np.mean(kmer8_list)

                    i=0
                    kmer2 = kmer_dict[head_num][i]
                    kmer3 = kmer_dict[head_num][i+1]
                    kmer4 = kmer_dict[head_num][i+2]
                    kmer5 = kmer_dict[head_num][i+3]
                    kmer6 = kmer_dict[head_num][i+4]
                    kmer7 = kmer_dict[head_num][i+5]
                    kmer8 = kmer_dict[head_num][i+6]

                    mbp = mbt / cdst
                    phagep = phaget / cdst
                    mbphp = mbpht / cdst

                    vfp = vft / cdst
                    arp = art / cdst
                    vfarp = vfart / cdst

                    ngp = ngt / cdst
                    vfarngp = vfarngt / cdst

                    cdsp = (cdst+rnat)*1000 / size
                    #cdsp = (cdst)*1000 / size
                if rnat > 0:
                    rnap = rnat / (cdst+rnat)
                if cdst > 1:
                    sorted_gene_list = sorted(gene_list, key=lambda x : (int(x[0]), int(x[1])))
                    dist = 0
                    prev_gend = sorted_gene_list[0][1]
                    for gstart, gend in sorted_gene_list[1:]:
                        idist = gstart-prev_gend
                        if idist<0:
                            idist = 0
                        dist += idist
                        prev_gend = gend
                    inter_gene_dist = dist/(cdst+rnat-1)
                # 'GeneID\tStart\tEnd\tStrand\tGC\tGC1\tGC2\tGC3\tGC_skew\tCAI\tCBI\tFop\tNc\t2-mer\t3-mer\t4-mer\t5-mer\t6-mer\t7-mer\t8-mer\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\tGene density\tRNA density\tIntergenetic distance\n'
                # feature_dict[head_num] = (cluster, (rnat, gc, gc1, gc2, gc3, gc_skew, cai, cbi, fop, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8, mbphp, vfarngp, cdsp, inter_gene_dist))
                # feature_dict[head_num] = (cluster, (rnat, gcp, gc1p, gc2p, gc3p, gc_skewp, cai, cbi, fop, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8, mbphp, vfarngp, cdsp, inter_gene_dist))
                feature_dict[head_num] = (cluster, (cdst, rnat, gc, gc1, gc3, cub, aab, chi, cai, cbi, fop, kmer2, kmer3, kmer4, kmer5, kmer6, kmer7, kmer8, mbp, phagep, vfp, arp, ngp, cdsp, inter_gene_dist))
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


if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-i", "--input", dest="input", help="input file of extracted gene features")
    parser.add_option("-o", "--output", dest="output", help="output file of the number of each feature in a genomic interval")
    parser.add_option("-r", "--repeats", dest="repeats",
                      help="input file containing repeats overlapping with the input regions")
    parser.add_option("-g", "--gcfile", dest="gcfile", help="input file containing gc content in a genomic interval")
    parser.add_option("-k", "--kmerfile", dest="kmerfile",
                      help="input file containing kmer distribution in a genomic interval")
    parser.add_option("-t", "--trnas", dest="trnas",
                      help="input file containing trnas overlapping with the input regions")
    options, args = parser.parse_args()

    # feature_dict = get_feature_dict_percentage(options.input, options.percentage, options.long)

    rna_dict = get_rna(options.trnas)
    repeat_dict = get_repeat(options.repeats)
    # gc_dict = get_gc_dist(options.gcfile)
    kmer_dict = get_kmer_dist(options.kmerfile)
    # feature_dict = get_feature_dict(options.input, rna_dict, repeat_dict)
    feature_dict = get_feature_dict_group(options.input, rna_dict=rna_dict, repeat_dict=repeat_dict, gc_dict={}, kmer_dict=kmer_dict)

    with open(options.output, 'w') as fout:
        for cluster, feature in feature_dict.values():
            content = cluster + feature
            format = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d'
            # more precision for gc
            # for i in range(2,5):
            #     format += '\t%.8f'
            for i in range(2,len(feature)):
                format += '\t%.3f'
            format += '\n'
            # print format
            # print content
            line = format % content
            fout.write(line)
