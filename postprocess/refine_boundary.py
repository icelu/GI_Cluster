# This script is used to refine the boundary of intervals based on tRNAs or repeats overlapping with an interval
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# A set of intervals
# tRNAs or repeats overlapping with an interval
#
# Output:
# Invervals with new start and end position


from __future__ import division
import optparse

import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.parse_sequence import get_rna_segment, get_repeat_segment
from util.interval_operations import get_orig_intervals


def refine_boundary(intervals, rna_dict, repeat_dict, threshold=500):
    new_intervals = []
    for coord in intervals:
        (p1, p2) = coord
        # Parse the postions to get integer coordinates
        contig_id = 0
        if '_' in p1:
            mark = p1.index('_')
            contig_id = int(p1[0:mark])
            start = int(p1[mark + 1:])
            end = int(p2[mark + 1:])
        else:
            start = int(p1)
            end = int(p2)

        if coord in rna_dict.keys():
            rnas = rna_dict[coord]
        else:
            rnas = []
        if coord in repeat_dict.keys():
            repeats = repeat_dict[coord]
        else:
            repeats = []
        total = rnas + repeats
        if len(total) > 0:
            # print total
            sorted_total = sorted(total, key=lambda x: (int(x[0]), int(x[1])))
            # Compare the distance between the interval
            rs = sorted_total[0][0]
            re = sorted_total[-1][1]
            # New boundary: the outmost position of trna or repeat
            if(abs(start - rs) < threshold):
                ns = rs
            else:
                ns = start
            if(abs(end - re) < threshold):
                ne = re
            else:
                ne = end
            if ns == rs or ne == re:
                print '%d\t%d\t--\t%d\t%d' % (start, end, ns, ne)

            if contig_id != 0:
                ns = str(contig_id) + '_' + str(ns)
                ne = str(contig_id) + '_' + str(ne)
            new_intervals.append((ns, ne))
        else:   # Not update intervals
            new_intervals.append((p1, p2))
    return new_intervals


def write_output(intervals, outfile):
    with open(outfile, 'w') as fout:
        for s, e in intervals:
            line = '%s\t%s\n' % (s, e)
            fout.write(line)


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input",
                      help="input interval file (predicted GIs)")
    parser.add_option("-r", "--repeats", dest="repeats",
                      help="input file containing repeats overlapping with the input regions")
    parser.add_option("-t", "--trnas", dest="trnas",
                      help="input file containing trnas overlapping with the input regions")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="int", default="500",
                      help="the cutoff number of base distances between tRNA/repeat and a region")
    parser.add_option("-o", "--output", dest="output",
                      help="output file of intervals with new boundaries")
    # parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true', default=False, help="Analyze contigs from unassembled genomes")
    options, args = parser.parse_args()

    rna_dict = get_rna_segment(options.trnas)
    repeat_dict = get_repeat_segment(options.repeats)

    intervals = get_orig_intervals(options.input)
    new_intervals = refine_boundary(intervals, rna_dict, repeat_dict, threshold=options.cutoff)

    write_output(new_intervals, options.output)
