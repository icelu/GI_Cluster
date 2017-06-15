# The program for evaluation of GI prediction methods by comparing with reference GIs:
# Given a set of reference intervals (all reported GIs with high reliability), and a set of query intervals (predicted GIs from programs),
# find the overlap intervals (TPs) and unique regions in the query (FPs)
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Steps:
# 1. Find the overlap regions for each reference interval
# 2. The regions in query intervals with no overlap with reference intervals are FPs
# For accurate comparison, the intervals should be 1-based.
#
# Differences between 0-based and 1-based intervals:
# 0-based (12,15), (15,20) --> actually means 1-based (13,15), (16,20).
# For 0-based intervals, the condition for merging is st <= saved[1],
# but for 1-based intervals, the condition for merging is st-1 <= saved[1].
# As reference intervals are 1-based, here query intervals should also be 1-based.
#
# Use quicksect to find all the query intervals overlapping with a reference interval.
# Because quicksect use no equal comparison,
# given two intervals in 1-based query [(13, 15), (18, 25)] and 1-based reference interval (15,20),
# we cannot find overlap of 1bp. Namely, (13, 15) will not be reported.
# But if the query is 0-based, it will become [(12,15), (14,20)].
# Then if reference intervals are still 1-based, 1bp overlap can still not be found.
# If reference intervals are 0-based (14,20), the overlap of 1bp can be found, namely (12,15) will be reported.
# Since 1bp overlap is minimal compared with the size of a genomic island, we use 1-based query for consistency.
#
# Input:
#  Predicted genomic islands
#  Reference genomic islands
#  Gene positions
#
# Output:
# false postives
# multiple metrics (output to the screen, which can be redirected to a file)
# *_offset: the predicted regions overlapping with each reference region
# Predicted intervals that overlap with reference intervals (optional)
#
# Command sample:
# python IntervalOverlap.py -c 0.5 -q predicted_GI -i reference_GI -p reference.ptt -o  eval_FPs_GI > eval_std_GI


from __future__ import division
import os
import optparse
import math

import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.interval_operations import find, get_window_tree, get_overlap, get_overlap_interval, get_interval_length, merge, merge_ref, merge_score


################################ parse input ############################################

def getGeneCoords(intervalfile):
    '''
    Return an array of tuples
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(float(fields[0])), int(float(fields[1])))
            intervals.append(coord)

    return intervals



def getIntervals(intervalfile):
    '''
    Return a list of interval tuples with additional information (score)
    Assume 1-based coordinates
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                coord = (int(fields[0]), int(fields[1]), float(fields[2]))
            else:  # len(fields)==2
                start = int(fields[0])
                end = int(fields[1])
                size = end - start + 1
                # add size so that being consistent with overlap intervals
                coord = (start, end, size)
            intervals.append(coord)

    return intervals


def getIntervals_0based(intervalfile):
    '''
    Return a list of interval tuples with additional information (score)
    Assume 0-based coordinates
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                coord = (int(fields[0]) + 1, int(fields[1]), float(fields[2]))
            else:  # len(fields)==2
                coord = (int(fields[0]) + 1, int(fields[1]))
            intervals.append(coord)

    return intervals



############################## interval operation ########################################
def getOverlapIntervalSize(interval_list):
    '''
    The interval_list may contain overlapping regions
    Use merge() to get only unique regions
    '''
    intervals = []
    if len(interval_list[0]) == 3:
        intervals = merge_score(interval_list)
    else:
        intervals = merge(interval_list)
    size = 0
    for i in intervals:
        size += i[1] - i[0] + 1
    return size


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
            line = '%d\t%d\t%d\n' % (value[0], value[1], value[2])
        else:
            line = '%d\t%d\n' % (value[0], value[1])
        outfile.write(line)

    outfile.close()


#################################### GI prediction metrics #####################################################

def getGenelocList(infile):
    '''
    Protein-coding gene position can be read from NCBI ptt file
    Input format:
    Salmonella enterica subsp. enterica serovar Typhi str. CT18, complete genome - 1..4809037
    4111 proteins
    Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
    190..255        +       21      16758994        -       STY0001 -       -       thr operon leader peptide
    '''
    genelist = []
    with open(infile, 'rb') as fin:
        # skip the first 3 lines
        for _ in xrange(3):
            next(fin)
        for line in fin:
            fields = line.strip().split('\t')
            loc = fields[0].split('..')
            start = int(loc[0])
            end = int(loc[1])
            genelist.append((start, end))
    return genelist



def getGenesInInterval(interval, genelist, alpha=0):
    '''
    alpha: the degree of overlap
    '''
    start = interval[0]
    end = interval[1]
    genes = []
    for g_start, g_end in genelist:
        # check whether gene is included in the interval
        if end < g_start or start > g_end:
            continue

        # completely included genes
        if start <= g_start and end >= g_end:
            genes.append((g_start, g_end))
        # partly covered genes
        else:
            cov = get_overlap((g_start, g_end), interval)
            gene_len = g_end - g_start + 1
            if cov > alpha * gene_len:
                genes.append((g_start, g_end))

    return list(set(genes))


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




def getMetric_base(ref_intervals, query_intervals, genelist, options):
    '''
    For evaluation of the predictions of intervals based on nucleotides or genes.
    A given percentage is used to determine whether a reference interval (genomic island or gene) is detected.
    If only a base of an interval is detected, it is meaningless to classify this interval as predicted.
    '''
    print '(start, end, size)\tleft_offset\tright_offset\toverlap_regions\tpredicted_size\toverlap_percentage\toverlap_percentage_pred\tnum_reference_genes\tnum_predicted_genes\tnum_overlap_genes\toverlap_gene_percentage\toverlap_gene_percentage_pred\tgaps'
    tree = get_window_tree(query_intervals)
    avg_query_len = get_interval_length(query_intervals) / len(query_intervals)
    # The list of query intervals overlapping with the reference
    # Use list for convenience in inserting, convert to set later in case some intervals overlapping with multiple reference intervals
    overlap_intervals = []
    # The number of reference intervals with no overlapping query intervals
    num_nooverlap = 0

    # For evaluation based on bases
    # The total size of all the reference intervals
    ref_totalSize = 0
    overlap_totalSize = 0
    # To record the boundary offset for each reference interval
    extentions = []
    foffset = ''
    if options.output:
        offsetfile = options.output + '_offset'
        foffset = open(offsetfile, 'w')

    # For evaluation based on genes
    if options.pttfile:
        overlap_total_genes = set()
        ref_total_genes = set()
    tp_interval = 0

    for start, end in ref_intervals:
        # Find all query intervals overlapping with the reference (even 1 bp overlap is counted here)
        overlap = find(start, end, tree)
        if len(overlap) == 0:  # No query intervals overlapping with the reference interval
           num_nooverlap += 1
        # Update overlap_intervals for all reference intervals
        overlap_intervals.extend(overlap)

        ################################## boundary offset ######################################################
        offset_line_str = ''
        offset_line = []
        '''
        The predicted interval may be much larger than the reference interval.
        Compute the extended region for FPs to check boundary accuracy.
        There maybe several predicted intervals, only the boundary intervals need to be checked
        '''
        sorted_overlap = sorted(overlap, key=lambda x : (int(x[0]), int(x[1])))
        left_ext = 0
        right_ext = 0
        if(len(sorted_overlap) > 0):
            # negative number represents predicted interval crossing the reference region
            left_ext = sorted_overlap[0][0] - start
            right_ext = end - sorted_overlap[-1][1]
            # If we only count boundary error when the overlap is large enough, comment out this statement
            extentions.append((abs(left_ext), abs(right_ext)))
            offset_line_str += '%s\t%s\t%s\t%s'
            offset_line.extend([start, end, left_ext, right_ext])


        ##################################### gene metrics ####################################################
        if options.pttfile:
            ref_genes = getGenesInInterval((start, end), genelist, options.cutoff_gene)
            num_refgenes = len(ref_genes)
            ref_total_genes.update(ref_genes)
            # There may be overlap if some genes are counted in both intervals
            overlap_genes = set()
            # The number of genes in the predicted intervals, to get an idea of over/under estimation
            predicted_genes = set()

            for interval in overlap:
                genes = getGenesInInterval(interval, genelist, options.cutoff_gene)
                # Output genes across the boundary of all predicted intervals
                if len(genes) > 0:
                    sorted_genes = sorted(genes, key=lambda x : (int(x[0]), int(x[1])))
                    offset_line.extend([sorted_genes[0], sorted_genes[-1]])
                    offset_line_str += '\t%s\t%s'
                overlap_gene = set(ref_genes).intersection(genes)
                # For this ref interval
                predicted_genes.update(genes)
                overlap_genes.update(overlap_gene)
                # For all the ref intervals
                overlap_total_genes.update(overlap_genes)

            num_overlapgenes = len(overlap_genes)
            # num_predictedgenes > num_overlapgenes, since genes in the predicted interval may not overlap with reference intervals
            num_predictedgenes = len(predicted_genes)

            if num_refgenes > 0:
                overlap_gene_percentage = round(num_overlapgenes * 100 / num_refgenes, 3)
            else:
                overlap_gene_percentage = 0

            if num_predictedgenes > 0:
                overlap_gene_percentage_pred = round(num_overlapgenes * 100 / num_predictedgenes, 3)
            else:
                overlap_gene_percentage_pred = 0

         ##################################### base metrics ####################################################
        # Check the coverage of the overlap
        overlap_bases = 0
        # Record the intervals of reference and query to get the union of overlapping bases
        overlap_interval_list = []
        for interval in overlap:
            # overlap_interval is a tuple
            overlap_interval = get_overlap_interval(interval, (start, end))
            if overlap_interval is not None:
                overlap_interval_list.append(overlap_interval)

        if foffset != '':
            offset_line_str += '\n'
            foffset.write(offset_line_str % tuple(offset_line))

        if len(overlap_interval_list) > 0:
            overlap_bases = getOverlapIntervalSize(overlap_interval_list)

        refsize = int(end) - int(start) + 1
        ref_totalSize += refsize
        # only count overlap if the reference interval are counted as found
        if overlap_bases > options.cutoff_base * refsize:
            overlap_totalSize += overlap_bases
            tp_interval += 1
        # not enough overlap
        else:
            for interval in overlap:
                overlap_intervals.remove(interval)

        # The fraction of a reference interval covered by all the query interval overlapping with it
        overlap_percentage = round(overlap_bases * 100 / refsize, 3)

        line = [start, end, (end - start + 1), left_ext, right_ext]
        overlap_len = len(sorted_overlap)
        line_str = '(%s, %s, %s)\t%d\t%d'

        # Compute the gap between intervals when there are more than two intervals
        suffix_str = ';%s'
        query_size = 0
        for i, overlap in enumerate(sorted_overlap):
            if i==0:
                line_str += '\t%s'
            else:
                line_str += suffix_str
            line.append(overlap)
            query_size += overlap[1]-overlap[0]+1

        line_str += '\t%s'
        line.append(query_size)
        # The fraction of a reference interval covered by all the query interval overlapping with it
        if query_size > 0:
            overlap_percentage_pred = round(overlap_bases * 100 / query_size, 3)
        else:
            overlap_percentage_pred = 0

        gaps = []
        if overlap_len > 1:
            p1 = sorted_overlap[0][1]
            for overlap in sorted_overlap[1:]:
                p2 = overlap[0]
                # 2-1 = 1, but they are adjacent
                gap = p2 - p1 - 1
                gaps.append(gap)
                p1 = overlap[1]
        # Find all genes in a reference interval
        if options.pttfile:
            # 8 fields
            line_str += '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'
            line.extend([overlap_percentage, overlap_percentage_pred, num_refgenes, num_predictedgenes, num_overlapgenes, overlap_gene_percentage, overlap_gene_percentage_pred, gaps])
        else:
            line_str += '\t%s\t%s\t%s'
            line.extend([overlap_percentage, overlap_percentage_pred, gaps])
        print line_str % tuple(line)

    # The number of all the predicted intervals overlapping with the reference
    num_overlap = len(set(overlap_intervals))
    print '\nThe number of predicted intervals: %s' % len(query_intervals)
    print 'The number of reference intervals: %s' % len(ref_intervals)
    print 'The number of predicted reference intervals (TPs, at least certain fraction (0.4 by default) of the reference interval is predicted): %s' % tp_interval
    print 'The number of unpredicted reference intervals (FNs): %s' % (len(ref_intervals) - tp_interval)
    print 'The number of reference intervals not overlapping with predictions: %s' % num_nooverlap
    # Some intervals may be overlapped with different reference intervals, so this number may be overestimated
    print 'The number of predicted intervals overlapping with the reference: %s' % num_overlap
    # The query intervals not overlapping with the reference
    unique_intervals = set(query_intervals) - set(overlap_intervals)
    print 'The number of predicted intervals not overlapping with the reference (FPs): %s' % len(unique_intervals)

    ############################## PR in #intervals ###############################################
    '''
    This is in term of #intervals. Not meaningful. For reference.
    '''
    # recall: TP/(TP+FN), precision= TP/(TP+FP)
    tp = tp_interval
    # The number of predicted intervals not overlapping with the reference, hard to be mapped to reference intervals
    # fp = len(query_intervals) - tp
    fp = len(unique_intervals)
    # By definition, FN should be the number of unpredicted intervals overlapping with the reference interval
    # Here, for convenience, it is the number of unpredicted reference intervals
    fn = len(ref_intervals) - tp
    real = len(ref_intervals)

    predicted = len(set(query_intervals))
    '''
    When using real, recall may be larger than 1, as many predicted intervals may overlap with the same GI.
    So we use tp+fn instead.
    '''
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    if recall != 0 and precision != 0:
        fmeasure = 2 * recall * precision / (recall + precision)
    else:
        fmeasure = 0
    print 'Interval Recall: %.3f\tPrecision: %.3f\tF-measure: %.3f\tPredicted intervals: %s\tTotal intervals: %s\t' % (recall, precision, fmeasure, predicted, real)

    ############################## PR in #overlap bases ###############################################
    '''
    This is in term of overlapping bases.
    '''
    tp = overlap_totalSize
    real = get_interval_length(ref_intervals)
    # Two alternative ways to get the number of bases in reference intervals
    assert real == ref_totalSize
    # Merge before counting as query_intervals may be overlapping
    predicted = getOverlapIntervalSize(query_intervals)

    recall = tp / real
    precision = tp / predicted
    if recall != 0 and precision != 0:
        fmeasure = 2 * recall * precision / (recall + precision)
    else:
        fmeasure = 0
    # Append offset error at the end of last line
    (avg_left, avg_right, avg_offset) = getAvgBoundaryError(extentions)

    print 'The number of reference bases: %d' % real

    ############################## PR in #overlap genes ###############################################
    '''
    This is in term of overlapping genes.
    '''
    if options.pttfile:
        tp = len(set(overlap_total_genes))
        g_real = len(set(ref_total_genes))
        query_total_genes = set()
        for query_interval in query_intervals:
            query_genes = getGenesInInterval(query_interval, genelist, options.cutoff_gene)
            query_total_genes.update(query_genes)
        g_predicted = len(query_total_genes)

        g_recall = tp / g_real
        g_precision = tp / g_predicted
        if g_recall != 0 and g_precision != 0:
            g_fmeasure = 2 * g_recall * g_precision / (g_recall + g_precision)
        else:
            g_fmeasure = 0
        diff_fmeasure = g_fmeasure - fmeasure
        print 'The number of reference genes: %d' % g_real

        # Other measures
        fn = g_real - tp
        fp = g_predicted - tp
        neg = len(genelist) - g_real
        tn = neg - fp
        tnr = (tn)/(tn+fp)
        oacc = (tp+tn)/(tp+tn+fn+fp)
        acc = (g_recall+tnr)/2
        mcc = (tp*tn - fp*fn)/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

    if options.pttfile:
        format_str = 'Bases Recall: %.3f\tPrecision: %.3f\tF-measure: %.3f\tF-measure Difference: %.3f\tLeft offset: %d\tRight: %d\tPredicted bases: %d\tOverlap bases: %d\tAverage interval size: %d\tOverlap genes: %d\tPredicted genes: %d\tPredicted intervals: %d\tGene Recall: %.3f\tPrecision: %.3f\tF-measure: %.3f\tAvg offset: %d\tTNR: %.3f\tOACC: %.3f\tACC: %.3f\tMCC: %.3f'
        print format_str % (recall, precision, fmeasure, diff_fmeasure, avg_left, avg_right, predicted, overlap_totalSize, avg_query_len, len(overlap_total_genes), g_predicted, len(query_intervals), g_recall, g_precision, g_fmeasure, avg_offset, tnr, oacc, acc, mcc)
    else:
        format_str = 'Bases Recall: %.3f\tPrecision: %.3f\tF-measure: %.3f\tLeft offset: %d\tRight: %d\tPredicted bases: %d\tOverlap bases: %d\tPredicted intervals: %d\tAverage interval size: %d\tAvg offset: %d'
        print format_str % (recall, precision, fmeasure, avg_left, avg_right, predicted, overlap_totalSize, len(query_intervals), avg_query_len, avg_offset)

    # Output FP predictions
    if options.overlap:
        # Sort intervals
        overlap_intervals = sorted(set(overlap_intervals), key=lambda x : (int(x[0]), int(x[1])))
        writeListOfTupleToFile(options.overlap, overlap_intervals)
    if options.output:
        unique_intervals = sorted(unique_intervals, key=lambda x : (int(x[0]), int(x[1])))
        writeListOfTupleToFile(options.output, unique_intervals)

    if foffset != '':
        foffset.close()

    return (recall, precision, fmeasure, options.output)



if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-i", "--interval", dest="interval", help="input interval file of gene locus")
    parser.add_option("-q", "--qinterval", dest="qinterval", help="input query interval file of gene locus")
    parser.add_option("-o", "--output", dest="output", help="output file of FP intervals")
    parser.add_option("-l", "--overlap", dest="overlap", help="output file of overlapping intervals")
    parser.add_option("-c", "", dest="cutoff_gene", type="float", default="0", help="the fraction of genes covered")
    parser.add_option("-b", "", dest="cutoff_base", type="float", default="0.4", help="the fraction of bases covered in a GI")
    parser.add_option("-p", "--pttfile", dest="pttfile", help="input ptt file")
    parser.add_option("-z", dest="zero", default=False, action="store_true", \
          help="the query interval is 0-based (default=false)")

    options, args = parser.parse_args()

    ref_intervals = getGeneCoords(options.interval)
    nonoverlap_ref_intervals = merge_ref(ref_intervals)

    if options.zero:
        query_intervals = getIntervals_0based(options.qinterval)
    else:
        query_intervals = getIntervals(options.qinterval)

    if options.pttfile:
        genelist = (getGenelocList(options.pttfile))
    else:
        genelist = []

    getMetric_base(nonoverlap_ref_intervals, query_intervals, genelist, options)
