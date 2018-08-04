#!/usr/bin/env python

# Parse the output from COGcognitor to get genes assigned to "none" category
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# COG.csv --
#     gi_CDS_337-2799,820,1,820,102772,COG0460
# gene_locs --
#     >gi_CDS_337-2799
#
# Output:
# '*No hits*' or specific hit
#


import os
import optparse



def get_gene_locus(infile):
    gene_list = []
    with open(infile, 'rb') as fin:
        for line in fin:
            gene = line.strip().replace('>', '')
            # print gene
            gene_list.append(gene)
    return gene_list


def get_COG(infile):
    cog_dict = {}
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split(',')
            gid = fields[0]
            cog = fields[5]

            if gid not in cog_dict.keys():
                cog_dict[gid] = [cog]
            else:
                cog_dict.setdefault(gid, []).append(cog)
    # print cog_dict
    return cog_dict


def write_COG(gene_list, cog_dict, outfile):
    fout = open(outfile, 'w')
    for gene in gene_list:
        if gene not in cog_dict.keys() or '-1' in cog_dict[gene]:
            # print gene
            #line = '%s\n' % ("None category in COG")
            line = '%s\n' % ("*No hits*")
        else:
            #line = '%s\n' % ("Known category in COG")
            line = '%s\n' % (','.join(cog_dict[gene]))
        fout.write(line)
    fout.close()



if __name__ == '__main__':
        parser = optparse.OptionParser()

        parser.add_option("-g", "--gfile", dest="gfile", help="input file of genes")
        parser.add_option("-c", "--cfile", dest="cfile", help="input file of COGs assignment for each gene")
        parser.add_option("-o", "--outfile", dest="outfile", help="output file for the merged intervals")
        # parser.add_option("-c", "--cluster", dest="cluster", action="store_true", default=False, help="create cluster file from input interval file1")

        options, args = parser.parse_args()

        gene_list = get_gene_locus(options.gfile)
        cog_dict = get_COG(options.cfile)
        write_COG(gene_list, cog_dict, options.outfile)
