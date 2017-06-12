#!/usr/bin/env python

import os
import optparse


def read_loc_ncbi(locfile):
    '''
    find the location of all genes/ORFs in a genome from the annotations of NCBI
    '''
    loc_dict = {}
    i = 1
    with open(locfile, 'r') as fin:
        for line in fin:
            # multiple positions in a line
            # e.g. >gi|225869487|ref|NC_012471.1|:1177219-1177902,1177906-1178478 Streptococcus equi subsp. equi 4047, complete genome
            if ',' in line:
                if 'c' in line:
                   strand = 'R'
                   content = line.strip().replace('c', '')
                else:
                   strand = 'F'
                   content = line.strip()
                genes = content.split(',')
                poses = []
                for g in genes:
                    p = g.split('-')
                    intp = map(int, p)
                    poses.extend(intp)
                start = min(poses)
                end = max(poses)
            else:
                pos = line.strip().split('-')
                if 'c' in line:
                    start = pos[1]
                    end = pos[0][1:]
                    strand = 'R'
                else:
                    start = pos[0]
                    end = pos[1]
                    strand = 'F'
            loc_dict[i] = [start, end, strand]
            i += 1

    return loc_dict



def read_loc_pred(locfile, format=1):
    '''
    Find the location of all genes/ORFs in a genome from annotations
    format -- 0: ">gi_CDS_337-2796" (predictions from NCBI); 1: " 337 # 2796 # 1" or " 45473 # 46630 # -1" (predictions from Prodigal)
    '''
    loc_dict = {}
    i = 1

    if format == 0:     # for format ">gi_CDS_337-2796"
        with open(locfile, 'r') as fin:
            for line in fin:
                fields = line.strip().split('_')
                pos = fields[-1].split('-')
                p1 = int(pos[0])
                p2 = int(pos[1])
                if (p1 > p2):
                    start = p2
                    end = p1
                    strand = 'R'
                else:
                    start = p1
                    end = p2
                    strand = 'F'
                loc_dict[i] = [start, end, strand]
                i += 1

    elif format == 1:    # for format " 337 # 2796 # 1" or " 45473 # 46630 # -1"
        with open(locfile, 'r') as fin:
            for line in fin:
                fields = line.strip().split('#')
                start = int(fields[0].strip())
                end  = int(fields[1].strip())
                strand = int(fields[2].strip())
                if (strand == 1):
                    strand = 'F'
                else:
                    strand = 'R'
                loc_dict[i] = [start, end, strand]
                i += 1
    return loc_dict



def read_func(funcfile, loc_dict):
    i = 1
    with open(funcfile, 'r') as fin:
        for line in fin:
            func = line.strip()
            loc_dict.setdefault(i, []).append(func)
            i += 1

    return loc_dict


def write_gene(outfile, loc_dict):
    fout = open(outfile, 'w')
    for key, value in loc_dict.items():
        items = [str(key)]
        for v in value:
            items.append(str(v))
        # print items
        line = '\t'.join(items) + '\n'
        fout.write(line)
    fout.close()


if __name__ == '__main__':
        parser = optparse.OptionParser()

        parser.add_option("-l", "--locfile", dest="locfile", help="input file of gene/ORF locus")
        parser.add_option("-p", "--funcfile", dest="funcfile", help="input file of gene products/function")
        parser.add_option("-o", "--outfile", dest="outfile", help="output file")
        parser.add_option("-n", "--is_NCBI", dest="is_NCBI", action="store_true", default=False, help="the gene file is from ncbi")
        options, args = parser.parse_args()

        if options.is_NCBI:
            loc_dict = read_loc_ncbi(options.locfile)
        else:
            loc_dict = read_loc_pred(options.locfile, format=1)
        if options.funcfile:
            loc_dict = read_func(options.funcfile, loc_dict)
        write_gene(options.outfile, loc_dict)
