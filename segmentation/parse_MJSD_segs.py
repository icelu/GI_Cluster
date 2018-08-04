#!/usr/bin/env python

# Get the segments after the first step of MJSD
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# The output for the first step of MJSD
#
# Output:
# The positions for a set of intervals
#

import os
import optparse



# get segment boundaries
def parse_segs(fname, genome, minsize=5000):
    windows = []
    sids = []
    len_genome = len(genome)
    segs = []

    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    segs.append('0 1\n')
    try:
        for i, line in enumerate(lines):
            l = line.strip().split(':')
            if len(l) == 3 and line[:5] != 'Start':
                # print line
                # SEG POINT VALUE:50.7749SEG POINT LOCATION:4364008
                loc = line.strip().split(':')[2].strip()
                # val = line.strip().split(':')[1].split('SEG')[0]
                '''
                DMax Location 999889
                Rec Level 2
                end-start=1418913
                SHORTEST_LEN*2=100
                LOG(N)6.15196
                BETA: 0.855
                Neff: 30.8585
                POWER(202.257) = 1
                '''
                idx = lines.index('DMax Location ' + loc + '\n')
                val = lines[idx + 7].split(' = ')[1].strip()
                segs.append(loc + ' ' + val + '\n')
    except ValueError:
        print 'l%-ine', i
        print 'loc', loc
        print 'val', val
        raise Exception('Wrong')
    segs.append(str(len_genome) + ' 1\n')

    seg_sort = lambda x, y : int(x.split(' ')[0]) - int(y.split(' ')[0])
    segs.sort(seg_sort)
    start = -1
    size = 0
    for i in range(0, len(segs) - 1):
        loc1 = int(segs[i].split(' ')[0]) + 1
        loc2 = int(segs[i + 1].split(' ')[0])
        # val = segs[i + 1].split(' ')[1]
        size = loc2 - loc1 + 1
        # merge contiguous segments that are shorter than 5kb
        if size < minsize:
            if start < 0:
                start = loc1
            size += size
            continue
        if size >= minsize:
            if start > 0:
                print '%d\t%d\t%d' % (start, loc2, loc2 - start + 1)
                sids.append((start, loc2))
                # pos j is not included, use 0-based to extract substr
                slice = genome[start - 1:loc2].upper()
                # print 'slice: %s\n' % slice
                start = -1
                size = 0
            else:
                print '%d\t%d\t%d' % (loc1, loc2, loc2 - loc1 + 1)
                sids.append((loc1, loc2))
                # pos j is not included, use 0-based to extract substr
                slice = genome[loc1 - 1:loc2].upper()
            windows.append(slice)

    return segs



scale = 0.0
def mjsd_segs(segs, genome_sequence, order='2', use_genes=False, bgmult=1):
    ret = []
    yvalues = []
    guessed_positive_locs = []
    if not use_genes:
        seg_sort = lambda x, y : int(x.split(' ')[0]) - int(y.split(' ')[0])
        segs.sort(seg_sort)
    else:
        junk = 0
        # print 'USING GENES!'

    avg_seg_length = 0
    for i in range(1, len(segs)):
        loc1 = int(segs[i - 1].split(' ')[0])
        loc2 = int(segs[i].split(' ')[0])
        avg_seg_length += float(loc2 - loc1) / len(segs)


    pieces = []
    ymax = 0

    if use_genes:
        loop = range(0, len(segs) - 1, 2)
    else:
        loop = range(0, len(segs) - 1)
    for i in loop:
        loc1 = int(segs[i].split(' ')[0])
        loc2 = int(segs[i + 1].split(' ')[0])

        gene = genome_sequence[loc1:loc2]
        if len(gene) < 30:
            print 'gene is empty!'
            print 'gene is', gene
            print loc1
            print loc2
            print 'genome length is', len(genome_sequence)
            continue



    # global scale
    # scale = float(1500) / ymax
    for val in pieces:
        yvalue = val[0]
        loc1 = val[1]
        loc2 = val[2]
        yvalues.append(yvalue)
        guessed_positive_locs.append((loc1, loc2))
        line = str(loc1) + ' ' + '\n' + str(loc2) + ' ' + '\n'
        ret.append(line)

    return ret, guessed_positive_locs

def read_fasta(fname):
    f = open(fname)
    lines = f.readlines()
    f.close()
    i = 0
    for line in lines:
        if (line[0] == '>'):
            lines.pop(i)
            i -= 1
        i += 1
    seq = ''.join([line.strip() for line in lines])
    return seq


def get_final_segs(segs, genome_sequence, order):
    mjsd_regions_cp = []
    seg_sort = lambda x, y : int(x.split(' ')[0]) - int(y.split(' ')[0])
    segs.sort(seg_sort)
    mjsd_regions_cp, yvalues, guessed_positive_locs = mjsd_segs(segs, genome_sequence, order, use_genes, bgmult)


if __name__ == '__main__':
        parser = optparse.OptionParser()

        parser.add_option("-f", "--fastafile", dest="fastafile", help="input fasta file of genome sequence")
        parser.add_option("-s", "--segfile", dest="segfile", help="input file of genome segments")
        #parser.add_option("-o", "--output", dest="output", help="output file for the merged intervals")
        #parser.add_option("-c", "--cluster", dest="cluster", action="store_true", default=False, help="create cluster file from input interval file1")
        parser.add_option("-m", "--minl", dest="minl", type="int", default=5000, help="the minimum length of segments")

        options, args = parser.parse_args()
        genome_sequence = read_fasta(options.fastafile)
        len_genome = len(genome_sequence)
        parse_segs(options.segfile, genome_sequence, options.minl)
