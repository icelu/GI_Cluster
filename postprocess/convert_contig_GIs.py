#!/usr/bin/env python

# This script is used to convert the GIs predicted from contigs
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# id_start id_end
#
# Output:
# id start end size
#

import optparse

import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.parse_sequence import get_contig_IDs_rev
from util.interval_operations import get_intervals_contigs




if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--infile", dest="infile", help="input file of predicted GIs")
    parser.add_option("-g", "--genome_file", dest="genome_file", default="", help="input genome file in fasta format")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file")
    (options, args) = parser.parse_args()

    orig_intervals = get_intervals_contigs(options.infile)
    id_mapping = get_contig_IDs_rev(options.genome_file)

    with open(options.outfile, 'w') as fout:
        for id, intervals in orig_intervals.items():
            for start, end in intervals:
                size = end - start + 1
                line = '%s\t%d\t%d\t%d\t%d\n' % (id_mapping[id], id, start, end, size)
                fout.write(line)
