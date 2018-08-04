#!/usr/bin/env python

# For genomes with a large set of known GIs, assign labels to predicted regions overlapping with reference GIs for evaluation.
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# 1. predicted GIs with all features
# 2. predicted GIs that are overlapping with reference GIs
#
# Output:
# predicted GIs with additional labels


from __future__ import division
import os
import optparse


################################ parse input ############################################

'''
return a list of interval tuples
'''
def getIntervals(intervalfile):
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
                # size = end - start + 1
                # add size so that being consistent with overlap intervals
                # coord = (start, end, size)
                coord = (start, end)
            intervals.append(coord)

    # print len(intervals)
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



def assign_cluster(gifile, overlap_intervals, outfile):
    with open(gifile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[1])
            end = int(fields[2])
            coord = (start, end)
            if coord in overlap_intervals:
                oline = line.strip()+'\t'+'1\n'
                fout.write(oline)
            else:
                oline = line.strip()+'\t'+'0\n'
                fout.write(oline)

# assign clusters based on C-data set
def assign_cluster_cmp(gifile, overlap_pos, overlap_neg, outfile):
    with open(gifile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[1])
            end = int(fields[2])
            coord = (start, end)
            if coord in overlap_pos:
                oline = line.strip()+'\t'+'1\n'
                fout.write(oline)
            elif coord in overlap_neg:
                oline = line.strip()+'\t'+'0\n'
                fout.write(oline)
            else:
                oline = line.strip()+'\t'+'2\n'
                fout.write(oline)


if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option("", "--interval1", dest="interval1", help="input interval file1 (predicted GIs)")
        parser.add_option("", "--interval2", dest="interval2", help="input interval file2 (predicted GIs overlapping with reference GIs)")
        parser.add_option("", "--interval3", dest="interval3", help="input interval file3 (predicted GIs not overlapping with reference GIs)")
        parser.add_option("-o", "--output", dest="output", help="output file of labeled intervals")
        parser.add_option("-c", "--isCdata", dest="isCdata", action='store_true', default=False, help="assign cluster based on C-data set")

        options, args = parser.parse_args()

        if options.isCdata:
            overlap_pos = getIntervals(options.interval2)
            overlap_neg = getIntervals(options.interval3)
            assign_cluster_cmp(options.interval1, overlap_pos, overlap_neg, options.output)
        else:
            overlap_intervals = getIntervals(options.interval2)
            assign_cluster(options.interval1, overlap_intervals, options.output)
