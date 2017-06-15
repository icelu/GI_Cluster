# For comparison of tools for predicting genomic islands on reference datasets obtained from comparative genomics:
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
#
# Input:
# The benchmark dataset for reference positives,
# The benchmark dataset for reference negatives,
# The predicted intervals
#
# Output:
# Recall, Precision, Fmeasure:
#


from __future__ import division
import os
import optparse
import copy



def getOverlap(interval_list1, interval_list2):
    num_overlap = 0
    # Sort the list to facilitate searching
    # Suppose the intervals are not overlapping
    interval_list1 = sorted(interval_list1, key=lambda x : (int(x[0]), int(x[1])))
    interval_list2 = sorted(interval_list2, key=lambda x : (int(x[0]), int(x[1])))
    for i1 in interval_list1:
        for i2 in interval_list2:
            # when = holds, i2[1] = i1[1] = i2[0], not likely, as i2[1]>i2[0]
            overlap = max(0, min(i1[1], i2[1]) - max(i1[0], i2[0]) + 1)
            # print 'overlap: %s\t' % overlap
            num_overlap += overlap
    return num_overlap


def getIntervalSize(interval_list):
    size = 0
    for i in interval_list:
        size += i[1] - i[0] + 1
    return size


def getIntervals(intervalfile):
    '''
    return a list of interval tuples
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(float(fields[0])), int(float(fields[1])))
            intervals.append(coord)

    # print len(intervals)

    return intervals


def fmeasure(recall, precision):
    total = recall + precision
    # print '%d %d %d' % (recall, precision, total)
    fmeasure = 0
    if total != 0:
        fmeasure = (2 * recall * precision) / total
    # print 'after %d %d %d' % (recall, precision, total)
    return fmeasure


if __name__ == '__main__':
        parser = optparse.OptionParser()

        parser.add_option("-q", "--query", dest="query", help="query file")
        parser.add_option("-p", "--positive", dest="positive", help="positive file")
        parser.add_option("-n", "--negative", dest="negative", help="negative file")

        options, args = parser.parse_args()

        query_interval = getIntervals(options.query)
        positive_interval = getIntervals(options.positive)
        negative_interval = getIntervals(options.negative)

        tp = getOverlap(positive_interval, query_interval)
        fp = getOverlap(negative_interval, query_interval)
        real = getIntervalSize(positive_interval)

        recall = tp / real
        precision = tp / (tp + fp)
        fmeasure = fmeasure(recall, precision)

        print 'tp: %s\tfp: %s\treal: %s\n' % (tp, fp, real)
        print 'recall: %.3f\tprecision: %.3f\tfmeasure: %.3f\n' % (recall, precision, fmeasure)
