#!/usr/bin/env python

# This script is used to evaluate the predicted genomic islands (GIs) on contigs
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# Predicted GIs: cid start end (size)
# Reference GIs: cid start end
# Output:
# Overlap between predictions and references
#

from __future__ import division
import optparse
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.interval_operations import find, get_window_tree, get_overlap, get_overlap_interval, get_interval_length

def get_intervals_dict(intervalfile):
    '''
    Input:  intervalfile -- file containing a list of intervals with different IDs
    Return a list of intervals for each ID
    '''
    intervals = {}
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            id = int(fields[0])
            start = int(fields[1])
            end = int(fields[2])
            size = end - start + 1
            coord = (start, end, size)
            if id not in intervals.keys():
                intervals[id] = [coord]
            else:
                intervals.setdefault(id, []).append(coord)
    return intervals

def getAvgBoundaryError(extentions):
    '''
    Input: a list of offset tuple (abs(left_ext), abs(right_ext))
    '''
    if len(extentions) <= 0:
        return (0, 0, 0)
    sum_left = 0
    sum_right = 0

    for (left, right) in extentions:
        sum_left += left
        sum_right += right

    total = len(extentions)
    avg_left = sum_left / total
    avg_right = sum_right / total
    avg_offset = (avg_left + avg_right) / 2

    # avg_offset2 = (sum_left + sum_right) / (total * 2)
    # # There may be some errors due to float precision
    # if not avg_offset == avg_offset2:
    #     print 'avg_offset:%s\tavg_offset2:%s' % (avg_offset, avg_offset2)

    return (avg_left, avg_right, avg_offset)

# Suppose intervals are not overlapping
def getOverlapIntervalSize(interval_list):
    size = 0
    for i in interval_list:
        size += i[1] - i[0] + 1
    return size



def get_overlap_statistics(ref_intervals, query_intervals, genelist, options):
    '''
    For each reference interval, check the query intervals that overlap with it
    '''
    overlap_intervals = {}
    unique_intervals = {}
    extentions_dict = {}
    # The number of reference intervals with no overlapping query intervals
    num_nooverlap = 0

    ref_totalSize = 0
    overlap_totalSize = 0
    tp_interval = 0

    for id, r_intervals in ref_intervals.items():
        # Find query intervals with the same ID
        if id not in query_intervals.keys():
            continue
        q_intervals = query_intervals[id]
        tree = get_window_tree(q_intervals)
        overlap_intervals[id] = []
        unique_intervals[id] = []
        extentions = []
        for start, end, size in r_intervals:
            # Find all query intervals overlapping with the reference
            overlap = find(start, end, tree)
            if len(overlap) == 0:  # There is no query intervals overlapping with a reference interval
               num_nooverlap += 1
            # Update overlap_intervals for all reference intervals
            overlap_intervals[id].extend(overlap)

            ################################## boundary offset ######################################################
            sorted_overlap = sorted(overlap, key=lambda x : (int(x[0]), int(x[1])))
            left_ext = 0
            right_ext = 0
            if(len(sorted_overlap) > 0):
                left_ext = sorted_overlap[0][0] - start
                right_ext = end - sorted_overlap[-1][1]
                extentions.append((abs(left_ext), abs(right_ext)))

            ##################################### base metrics ####################################################
            # Check the coverage of the overlap
            overlap_bases = 0
            # Record the intervals of reference and query to get the union of overlapping bases
            overlap_interval_list = []
            for interval in overlap:
                overlap_interval = get_overlap_interval(interval, (start, end))
                if overlap_interval is not None:
                    overlap_interval_list.append(overlap_interval)

            if len(overlap_interval_list) > 0:
                # Suppose intervals in overlap_interval_list do not overlap
                overlap_bases = get_interval_length(overlap_interval_list)

            refsize = int(end) - int(start) + 1
            ref_totalSize += refsize
            # only count overlap if the reference interval are counted as found
            if overlap_bases > options.cutoff_base * refsize:
                overlap_totalSize += overlap_bases
                tp_interval += 1
            # not enough overlap
            else:
                for interval in overlap:
                    overlap_intervals[id].remove(interval)

            # The fraction of a reference interval covered by all the query interval overlapping with it
            overlap_percentage = round(overlap_bases * 100 / refsize, 3)

            suffix_str = '\t%s'
            line = [id, start, end, (end - start + 1), left_ext, right_ext]
            overlap_len = len(sorted_overlap)
            line_str = '%d\t(%s, %s, %s)\t%d\t%d'

            # Compute the gap between intervals when there are more than two intervals
            for overlap in sorted_overlap:
                line_str += suffix_str
                line.append(overlap)

            gaps = []
            if overlap_len > 1:
                p1 = sorted_overlap[0][1]
                for overlap in sorted_overlap[1:]:
                    p2 = overlap[0]
                    gap = p2 - p1 - 1
                    gaps.append(gap)
                    p1 = overlap[1]
            # Find all genes in a reference interval
            if options.genefile:
                line_str += '\t%s\t%s\t%s\t%s\t%s\t%s'
                line.extend([overlap_percentage, num_refgenes, num_predictedgenes, num_overlapgenes, overlap_gene_percentage, gaps])
            else:
                line_str += '\t%s\t%s'
                line.extend([overlap_percentage, gaps])
            print line_str % tuple(line)

        unique_intervals[id] = list(set(query_intervals[id]) - set(overlap_intervals[id]))
        extentions_dict[id] = extentions

    # For all the contigs
    # The number of all the predicted intervals overlapping with the reference
    # Note: all the intervals are dictionaries
    num_overlap = sum(len(v) for v in overlap_intervals.itervalues())
    num_ref = sum(len(v) for v in ref_intervals.itervalues())
    num_pred = sum(len(v) for v in query_intervals.itervalues())
    for id in query_intervals.keys():
        if id not in unique_intervals.keys():
            unique_intervals[id] = query_intervals[id]
    num_uniq = sum(len(v) for v in unique_intervals.itervalues())
    print '\nThe number of predicted intervals: %s' % num_pred
    print 'The number of reference intervals: %s' % num_ref
    print 'The number of predicted reference intervals (TPs): %s' % tp_interval
    print 'The number of unpredicted reference intervals (FNs): %s' % (num_ref - tp_interval)
    print 'The number of reference intervals not overlapping with predictions: %s' % num_nooverlap
    # Some intervals may be overlapped with different reference intervals, so this number may be overestimated
    print 'The number of predicted intervals overlapping with the reference: %s' % num_overlap
    print 'The number of predicted intervals not overlapping with the reference (FPs): %s' % num_uniq
    #
    # ############################## PR in #overlap bases ###############################################

    tp = overlap_totalSize
    real = 0
    for id in ref_intervals.keys():
        real += get_interval_length(ref_intervals[id])

    # Two alternative ways to get the number of bases in reference intervals
    assert real == ref_totalSize
    # Merge before counting as query_intervals may be overlapping
    predicted = 0
    for id in query_intervals.keys():
        predicted += getOverlapIntervalSize(query_intervals[id])

    recall = tp / real
    precision = tp / predicted
    if recall != 0 and precision != 0:
        fmeasure = 2 * recall * precision / (recall + precision)
    else:
        fmeasure = 0
    # Append offset error at the end of last line
    lavg = 0
    ravg = 0
    oavg = 0
    count_ext = 0
    for id, extentions in extentions_dict.items():
        (avg_left, avg_right, avg_offset) = getAvgBoundaryError(extentions)
        lavg += avg_left
        ravg += avg_right
        oavg += avg_offset
        count_ext += 1
    avg_offset_all = oavg/count_ext
    avg_left_all = lavg/count_ext
    avg_right_all = ravg/count_ext
    avg_size = predicted/num_pred
    print 'The number of reference bases: %d' % real
    format_str = 'Bases Recall: %.3f\tPrecision: %.3f\tF-measure: %.3f\tLeft offset: %d\tRight: %d\tPredicted bases: %d\tOverlap bases: %d\tPredicted intervals: %d\tAverage interval size: %d\tAvg offset: %d'
    print format_str % (recall, precision, fmeasure, avg_left_all, avg_right_all, predicted, overlap_totalSize, num_pred,  avg_size, avg_offset_all)

    # # Output FP predictions
    # if options.overlap:
    #      # Sort intervals
    #      overlap_intervals = sorted(set(overlap_intervals), key=lambda x : (int(x[0]), int(x[1])))
    #      writeListOfTupleToFile(options.overlap, overlap_intervals)
    # if options.output:
    #      unique_intervals = sorted(unique_intervals, key=lambda x : (int(x[0]), int(x[1])))
    #      writeListOfTupleToFile(options.output, unique_intervals)
    #


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-q", "--query", dest="query", help="query file")
    parser.add_option("-r", "--reference", dest="reference", help="reference file")
    parser.add_option("-o", "--output", dest="output", help="output file of FP intervals")
    parser.add_option("-l", "--overlap", dest="overlap", help="output file of overlapping intervals")
    parser.add_option("-c", "", dest="cutoff_gene", type="float", default="0", help="the fraction of genes covered")
    parser.add_option("-b", "", dest="cutoff_base", type="float", default="0.4", help="the fraction of bases covered in a GI")
    parser.add_option("-g", "--genefile", dest="genefile", help="input file containing the position of genes")
    # parser.add_option("-z", dest="zero", default=False, action="store_true", \
    #       help="the query interval is 0-based (default=false)")
    options, args = parser.parse_args()

    # Suppose input intervals are 1-based
    query_intervals = get_intervals_dict(options.query)
    reference_intervals = get_intervals_dict(options.reference)
    # if options.genefile:
    #     genelist = (getGenelocList(options.genefile))
    # else:
    genelist = []

    get_overlap_statistics(reference_intervals, query_intervals, genelist, options)
