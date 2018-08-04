#!/usr/bin/env python

# This script is used to find the intersection between two interval files
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

from __future__ import division
import optparse
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.interval_operations import find, get_window_tree, get_overlap


################################ parse input ############################################
def getIntervals(intervalfile):
    '''
    return a list of interval tuples with additional information (score)
    suppose 1-based coordinates
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                # coord = (int(fields[0]), int(fields[1]), float(fields[2]))
                coord = (int(fields[0]), int(fields[1]))
            else:  # len(fields)==2
                start = int(fields[0])
                end = int(fields[1])
                coord = (start, end)
            intervals.append(coord)

    return intervals


def getIntervals_0based(intervalfile):
    '''
    previous run use 0-based coords,
    this will affect the computation of boundary
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            # suppose 3rd column is a score
            if len(fields) == 3:
                coord = (int(fields[0]) + 1, int(fields[1]), float(fields[2]))
            else:  # len(fields)==2
                coord = (int(fields[0]) + 1, int(fields[1]))
            intervals.append(coord)

    return intervals


#################################### output formating ####################################################

def writeListTofile(filename, list):
    outfile = open(filename, 'w')

    for tuple in list:
        line = '%s\t%s\n' % (tuple[0], tuple[1])
        outfile.write(line)

    outfile.close()


def writeListOfTupleToFile(filename, list):
    outfile = open(filename, 'w')

    for value in list:
        if len(value) == 3:
            line = '%s\t%s\t%d\n' % (str(value[0]), str(value[1]), (value[2]))
        else:
            line = '%s\t%s\n' % (str(value[0]), str(value[1]))
        outfile.write(line)

    outfile.close()



def getCommonIntervals(interval1, interval2, fraction):
    '''
    find common intervals in both lists of intervals
    extract the common regions
    '''
    if len(interval1) == 0 or len(interval2) == 0:
        return []

    tree = get_window_tree(interval2)

    overlap_intervals = []
    overlap_interval1 = []
    overlap_interval2 = []
    # The number of reference intervals with no overlapping query intervals
    num_nooverlap = 0

    # for start, end, other in interval1:
    for start, end in interval1:
        # find all genes in a reference intervals
        # find all query intervals overlapping with the reference
        overlap = find(start, end, tree)
        if len(overlap) == 0:  # no query intervals overlapping with the reference interval
           num_nooverlap += 1

        for interval in overlap:
            o_start = max (interval[0], start)
            o_end = min(interval[1], end)
            o_size = o_end - o_start + 1
            overlap_intervals.append((o_start, o_end, o_size))
            ref_size = interval[1] - interval[0] + 1
            if o_size > fraction * ref_size:
                # not store size to facilitate set intersection
                overlap_interval1.append((start, end))
                overlap_interval2.append((interval[0], interval[1]))

    return set(overlap_intervals), set(overlap_interval1), set(overlap_interval2)



if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option("", "--interval1", dest="interval1", help="input interval file1")
        parser.add_option("", "--interval2", dest="interval2", help="input interval file2")
        parser.add_option("-o", "--output", dest="output", help="output file of common intervals")
        parser.add_option("-f", "", dest="fraction", type="float", default="0.5", help="the fraction of interval2 regions covered by interval1")
        parser.add_option("-z", dest="zero", default=False, action="store_true", help="the query interval is 0-based (default=false)")

        options, args = parser.parse_args()

        interval1 = getIntervals(options.interval1)
        interval2 = getIntervals(options.interval2)

        overlap_intervals, overlap_interval1, overlap_interval2 = getCommonIntervals(interval1, interval2, options.fraction)
        writeListOfTupleToFile(options.output, overlap_intervals)
        # the intervals in interval1 that overlap with interval2
        writeListOfTupleToFile(options.output + '_interval1', overlap_interval1)
        # the intervals in interval2 that overlap with interval1
        writeListOfTupleToFile(options.output + '_interval2', overlap_interval2)

        nonoverlap_interval1 = list(set(interval1) - set(overlap_interval1))
        nonoverlap_interval2 = list(set(interval2) - set(overlap_interval2))
        writeListOfTupleToFile(options.output + '_noninterval1', nonoverlap_interval1)
        writeListOfTupleToFile(options.output + '_noninterval2', nonoverlap_interval2)
