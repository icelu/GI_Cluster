# This script is used to evaluate the predicted genomic islands (GIs) on contigs
#
# Input:
# Predicted GIs: cid start end (size)
# Reference GIs: cid start end
# Output:
# Overlap between predictions and references
#


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
            coord = (start, end)
            if id not in intervals.keys():
                intervals[id] = [coord]
            else:
                intervals.setdefault(id, []).append(coord)
    return intervals



def get_overlap_statistics(ref_intervals, query_intervals, genelist, options):
    '''
    For each reference interval, check the query intervals that overlap with it
    '''
    overlap_intervals = {}
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

        for start, end in r_intervals:
            extentions = []
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

    # # For all the contigs
    # # The number of all the predicted intervals overlapping with the reference
    # num_overlap = len(set(overlap_intervals))
    # print '\nThe number of predicted intervals: %s' % len(query_intervals)
    # print 'The number of reference intervals: %s' % len(ref_intervals)
    # print 'The number of predicted reference intervals (TPs): %s' % tp_interval
    # print 'The number of unpredicted reference intervals (FNs): %s' % (len(ref_intervals) - tp_interval)
    # print 'The number of reference intervals not overlapping with predictions: %s' % num_nooverlap
    # # Some intervals may be overlapped with different reference intervals, so this number may be overestimated
    # print 'The number of predicted intervals overlapping with the reference: %s' % num_overlap
    # # The query intervals not overlapping with the reference
    # unique_intervals = set(query_intervals) - set(overlap_intervals)
    # print 'The number of predicted intervals not overlapping with the reference (FPs): %s' % len(unique_intervals)
    #
    # ############################## PR in #overlap bases ###############################################
    # '''
    # This is in term of overlapping bases.
    # '''
    # tp = overlap_totalSize
    # real = get_interval_length(ref_intervals)
    # # Two alternative ways to get the number of bases in reference intervals
    # assert real == ref_totalSize
    # # Merge before counting as query_intervals may be overlapping
    # predicted = getOverlapIntervalSize(query_intervals)
    #
    # recall = tp / real
    # precision = tp / predicted
    # if recall != 0 and precision != 0:
    #     fmeasure = 2 * recall * precision / (recall + precision)
    # else:
    #     fmeasure = 0
    # # Append offset error at the end of last line
    # (avg_left, avg_right, avg_offset) = getAvgBoundaryError(extentions)
    # print 'The number of reference bases: %d' % real
    # format_str = 'Bases Recall: %.3f\tPrecision: %.3f\tF-measure: %.3f\tLeft offset: %d\tRight: %d\tPredicted bases: %d\tOverlap bases: %d\tPredicted intervals: %d\tAverage interval size: %d\tAvg offset: %d'
    # print format_str % (recall, precision, fmeasure, avg_left, avg_right, predicted, overlap_totalSize, len(query_intervals), avg_query_len, avg_offset)
    #
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
