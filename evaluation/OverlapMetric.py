from __future__ import division
import os
import optparse
import copy


# For comparison of GI tools:
# Given the benchmark dataset R, the predicted dataset P,
# TP: nuleotides in both R and P
# FP: nuleotides in P not R
#
# Input:
# 2 sets of intervals
# R, P
#
# Output:
# overlap between R and P --> TP;
#
# Calculation of TP :
# do interval intersections of 2 sets of intervals, find the number of overlapped bases.
# a straightforward method is 2-level iteration as implemented here
#
# another way: (as IntervalOverlap.py)
# 1. find overlapped intervals
# 2. compute overlapping bases



def getOverlap(interval_list1, interval_list2):
    num_overlap = 0
    # sort the list to facilitate searching
    interval_list1 = sorted(interval_list1, key=lambda x : (int(x[0]), int(x[1])))
    interval_list2 = sorted(interval_list2, key=lambda x : (int(x[0]), int(x[1])))
    for i1 in interval_list1:
        for i2 in interval_list2:
            overlap = max(0, min(i1[1], i2[1]) - max(i1[0], i2[0]) + 1)
            # print 'overlap: %s\t' % overlap
            num_overlap += overlap
    return num_overlap


def getIntervalSize(interval_list):
    size = 0
    for i in interval_list:
        size += i[1] - i[0] + 1
    return size

'''
return a list of interval tuples
'''
def getIntervals(intervalfile):
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


        options, args = parser.parse_args()

        query_interval = getIntervals(options.query)
        positive_interval = getIntervals(options.positive)

        tp = getOverlap(positive_interval, query_interval)

        real = getIntervalSize(positive_interval)
        predicted = getIntervalSize(query_interval)

        # same result as IntervalOverlap.py, seems a bit slower
        recall = tp / real
        precision = tp / predicted
        fmeasure = fmeasure(recall, precision)

        print 'tp: %s\tpredicted: %s\treal: %s\n' % (tp, predicted, real)
        print 'recall: %s\tprecision: %s\tfmeasure: %s\n' % (recall, precision, fmeasure)
