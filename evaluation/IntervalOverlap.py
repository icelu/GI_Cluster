# The program for evaluation of GI prediction methods by comparing with reference GIs:
# Given a set of reference intervals (all reported GIs with high reliability), and a set of query intervals (predicted GIs from programs),
# find the overlap intervals (TPs) and unique regions in the query (FPs)
#
# Steps:
# 1. Find the overlap regions for each reference interval
# 2. The regions in query intervals with no overlap with reference intervals are FPs
# for accurate comparison, the intervals should be 1-based.
#
# Differences between 0-based and 1-based intervals:
# 0-based (12,15), (15,20) --> actually means 1-based (13,15), (16,20).
# For 0-based intervals, the condition for merging is st <= saved[1],
# but for 1-based intervals, the condition for merging is st-1 <= saved[1].
# As reference intervals are 1-based, here query intervals should also be 1-based.
#
# Use quicksect to find all the query intervals overlapping with a reference interval
# quicksect use no equal comparison,
# so given two intervals in 1-based query [(13, 15), (18, 25)], and 1-based reference interval (15,20),
# we cannot find overlap of 1bp. Namely, (13, 15) will not be reported.
# but if the query is 0-based, it will become [(12,15), (14,20)].
# Then if reference intervals are still 1-based, 1bp overlap can still not be found.
# if reference intervals are 0-based (14,20), the overlap of 1bp can be found, namely (12,15) will be reported.
# since 1bp overlap is minimal compared with the size of GI, we use 1-based query for consistency.
#
# Output file:
# *_offset: the predicted regions overlapping with each reference region
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
from util.interval_operations import find, get_window_tree, get_overlap


################################ parse input ############################################
'''
Return an array of tuples
'''
def getGeneCoords(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(float(fields[0])), int(float(fields[1])))
            intervals.append(coord)

    return intervals


'''
Return a list of interval tuples with additional information (score)
Assume 1-based coordinates
'''
def getIntervals(intervalfile):
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

    # print len(intervals)
    return intervals

'''
Return a list of interval tuples with additional information (score)
Assume 0-based coordinates
'''
def getIntervals_0based(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) == 3:
                coord = (int(fields[0]) + 1, int(fields[1]), float(fields[2]))
            else:  # len(fields)==2
                coord = (int(fields[0]) + 1, int(fields[1]))
            intervals.append(coord)

    # print len(intervals)
    return intervals



############################## interval operation ########################################
'''
Get the union of (overlap) intervals
'''
def getOverlapInterval(a, b):
    start = max(a[0], b[0])
    end = min(a[1], b[1])
    if end >= start:
        return (start, end)

'''
The interval_list may contain overlapping regions
Use merge() to get only unique regions
'''
def getOverlapIntervalSize(interval_list):
    intervals = []
    if len(interval_list[0]) == 3:
        intervals = merge_score(interval_list)
    else:
        intervals = merge(interval_list)
    size = 0
    for i in intervals:
        size += i[1] - i[0] + 1
    return size

'''
Assume intervals are not overlapping,
Add up the bases of all the intervals
'''
def getIntervalLen(intervals):
    length = 0
    # can not use len as variable name, or else it will be confused with the function len
    # TypeError: object of type 'generator' has no len()
    # if len(intervals) > 0:
    for interval in intervals:
        length += interval[1] - interval[0] + 1
    return length


'''
Get the average size of all the intervals
'''
def getIntervalAvg(intervals):
    length = 0

    for interval in intervals:
        length += interval[1] - interval[0] + 1
    return length / len(intervals)


def merge(intervals):
    if len(intervals) < 0: return
    intervals = sorted(intervals, key=lambda x : (int(x[0]), int(x[1])))
    saved = list(intervals[0])
    for st, en in intervals:
        if st - 1 <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)

'''
Given a list of intervals, merge overlapping ones, return a list of non-overlapping intervals
'''
def merge_ref(intervals):
    if len(intervals) < 0: return
    intervals = sorted(intervals, key=lambda x : (int(x[0]), int(x[1])))
    saved = list(intervals[0])
    merged_intervals = []
    for st, en in intervals:
        if st - 1 <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            merged_intervals.append(tuple(saved))
            saved[0] = st
            saved[1] = en
    merged_intervals.append(tuple(saved))

    return merged_intervals


'''
Merge regions with score
'''
def merge_score(intervals):
    intervals = sorted(intervals, key=lambda x : (int(x[0]), int(x[1])))
    merged_intervals = []
    saved = list(intervals[0])
    for st, en, score in intervals:
        if st - 1 <= saved[1]:
            saved[1] = max(saved[1], en)
            # average the score, watch out data format
            saved[2] = (saved[2] + score) / 2
        else:
            # only add the interval to the result when not overlapping with adjacent regions
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
            saved[2] = float(score)
    yield tuple(saved)



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
'''
Protein-coding gene position can be read from ptt file
input format:
Salmonella enterica subsp. enterica serovar Typhi str. CT18, complete genome - 1..4809037
4111 proteins
Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
190..255        +       21      16758994        -       STY0001 -       -       thr operon leader peptide
'''
def getGenelocList(infile):
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


'''
alpha: the degree of overlap
'''
def getGenesInInterval(interval, genelist, alpha=0):
    start = interval[0]
    end = interval[1]
    # store all genes covered in an interval in a dict
    # dict = {}
    # dict[(start, end)] = []
    genes = []
    for g_start, g_end in genelist:
        # check whether gene is included in the interval
        # not overlap
        if end < g_start or start > g_end:
            continue

        # completely included genes
        if start <= g_start and end >= g_end:
        # dict[(start, end)].append((g_start, g_end))
            genes.append((g_start, g_end))
        # partly covered genes
        else:
            cov = get_overlap((g_start, g_end), interval)
            gene_len = g_end - g_start + 1
            if cov > alpha * gene_len:
                genes.append((g_start, g_end))

    return list(set(genes))



'''
Input: a list of offset tuple (abs(left_ext), abs(right_ext))
'''
def getAvgBoundaryError(extentions):
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

    '''
    avg_offset2 = (sum_left + sum_right) / (total * 2)
    # there may be some errors due to float precision
    if not avg_offset == avg_offset2:
        print 'avg_offset:%s\tavg_offset2:%s' % (avg_offset, avg_offset2)
    '''
    return (avg_left, avg_right, avg_offset)



'''
For evaluation based on nucleotides or genes,
A given percentage is used to determine whether a GI or gene is detected.
If only a base of a GI or gene is detected, it is meaningless to classify the GI or gene as predicted.
If less than a given fraction of sequence is detected as transferred, this predicted GI is discarded.
'''
def getMetric_base(ref_intervals, query_intervals, genelist, options):
        print '(start, end, size)\tleft_offset\tright_offset\toverlap_regions\toverlap_percentage\tnum_reference_genes\tnum_predicted_genes\tnum_overlap_genes\toverlap_gene_percentage\tgaps'
        tree = get_window_tree(query_intervals)
        avg_query_len = getIntervalLen(query_intervals) / len(query_intervals)
        # list of query intervals overlapping with the reference, use set in case some intervals overlapping with multiple reference intervals
        # use list for convenience in inserting, convert to set later, same effect as use set update()
        overlap_intervals = []
        # The number of reference intervals with no overlapping query intervals
        num_nooverlap = 0

        # for base statistics
        # The total size of all the reference intervals
        ref_totalSize = 0
        overlap_totalSize = 0

        # to record the boundary offset for each reference interval
        extentions = []
        foffset = ''
        if options.output:
            offsetfile = options.output + '_offset'
            foffset = open(offsetfile, 'w')

        # for gene statistics
        if options.pttfile:
            # num_total_refgenes = 0
            overlap_total_genes = set()
            ref_total_genes = set()
        tp_interval = 0

        for start, end in ref_intervals:
            # find all query intervals overlapping with the reference
            overlap = find(start, end, tree)
            if len(overlap) == 0:  # no query intervals overlapping with the reference interval
               num_nooverlap += 1
            # update overlap_intervals for all reference intervals
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
                # negative number represents crossing the reference region
                left_ext = sorted_overlap[0][0] - start
                right_ext = end - sorted_overlap[-1][1]
                # if only count boundary error when the overlap is large enough, comment out this statement
                extentions.append((abs(left_ext), abs(right_ext)))
                offset_line_str += '%s\t%s\t%s\t%s'
                offset_line.extend([start, end, left_ext, right_ext])


            ##################################### gene metrics ####################################################
            if options.pttfile:
                ref_genes = getGenesInInterval((start, end), genelist, options.cutoff_gene)
                num_refgenes = len(ref_genes)
                ref_total_genes.update(ref_genes)
                # there may be overlap if some genes are counted in both intervals
                # num_total_refgenes += num_refgenes
                overlap_genes = set()
                # the number of genes in the predicted, to get an idea of over/under estimation
                predicted_genes = set()

                for interval in overlap:
                    genes = getGenesInInterval(interval, genelist, options.cutoff_gene)
                    # output genes across the boundary of all predicted interval
                    if len(genes) > 0:
                        sorted_genes = sorted(genes, key=lambda x : (int(x[0]), int(x[1])))
                        offset_line.extend([sorted_genes[0], sorted_genes[-1]])
                        offset_line_str += '\t%s\t%s'

                    overlap_gene = set(ref_genes).intersection(genes)
                    # for this ref interval
                    predicted_genes.update(genes)
                    overlap_genes.update(overlap_gene)
                    # for all the ref intervals
                    overlap_total_genes.update(overlap_genes)

                num_overlapgenes = len(overlap_genes)
                # num_predictedgenes > num_overlapgenes, genes in the predicted region may not overlap with reference GI
                num_predictedgenes = len(predicted_genes)

                if num_refgenes > 0:
                    overlap_gene_percentage = round(num_overlapgenes * 100 / num_refgenes, 3)
                else:
                    overlap_gene_percentage = 0

             ##################################### base metrics ####################################################
            # check the coverage of the overlap, as intervals in overlap may be overlapping
            overlap_bases = 0

            '''
            Record the intervals of ref and query to get the union of overlapping bases
            '''
            overlap_interval_list = []
            for interval in overlap:
                # overlap_interval is a tuple
                overlap_interval = getOverlapInterval(interval, (start, end))
                if overlap_interval is not None:
                    overlap_interval_list.append(overlap_interval)

            if foffset != '':
                offset_line_str += '\n'
                foffset.write(offset_line_str % tuple(offset_line))

            # overlap_interval_list may be none
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

            # print out the results
            suffix_str = '\t%s'
            line = [start, end, (end - start + 1), left_ext, right_ext]
            overlap_len = len(sorted_overlap)
            line_str = '(%s, %s, %s)\t%d\t%d'

            # compute the gap between intervals when there are more than two intervals
            for overlap in sorted_overlap:
                line_str += suffix_str
                line.append(overlap)

            gaps = []
            if overlap_len > 1:
                p1 = sorted_overlap[0][1]
                for overlap in sorted_overlap[1:]:
                    p2 = overlap[0]
                    # 2-1 =1, but they are adjacent
                    gap = p2 - p1 - 1
                    gaps.append(gap)
                    p1 = overlap[1]
            # find all genes in a reference interval
            if options.pttfile:
                # 6 fields
                line_str += '\t%s\t%s\t%s\t%s\t%s\t%s'
                line.extend([overlap_percentage, num_refgenes, num_predictedgenes, num_overlapgenes, overlap_gene_percentage, gaps])
            else:
                line_str += '\t%s\t%s'
                line.extend([overlap_percentage, gaps])
            print line_str % tuple(line)

        # number of all the predicted intervals overlapping with the reference, use set to remove redundancy
        num_overlap = len(set(overlap_intervals))
        print '\nThe number of predicted intervals: %s' % len(query_intervals)
        print 'The number of reference intervals: %s' % len(ref_intervals)
        print 'The number of predicted reference intervals (TPs): %s' % tp_interval
        print 'The number of unpredicted reference intervals (FNs): %s' % (len(ref_intervals) - tp_interval)
        print 'The number of reference intervals not overlapping with predictions: %s' % num_nooverlap
        # some intervals may be overlapped with different reference intervals, so this number may be overestimated
        print 'The number of predicted intervals overlapping with the reference: %s' % num_overlap
        # The query intervals not overlap with the reference
        unique_intervals = set(query_intervals) - set(overlap_intervals)
        print 'The number of predicted intervals not overlapping with the reference (FPs): %s' % len(unique_intervals)

        ############################## PR in #intervals ###############################################
        '''
        this is in term of #intervals. not meaningful
        in principle, should be based on reference intervals
        '''
        # recall: TP/(TP+FN), precision= TP/(TP+FP)
        tp = tp_interval
        # The number of predicted intervals not overlapping with the reference, hard to be mapped to reference intervals
        # fp = len(query_intervals) - tp
        fp = len(unique_intervals)
        # By definition, should be the number of unpredicted intervals overlapping with the reference interval
        # Here, for convenience, it is the number of unpredicted reference intervals
        fn = len(ref_intervals) - tp
        real = len(ref_intervals)

        predicted = len(set(query_intervals))
        '''
        if using real, recall may be larger than 1, as many predicted intervals may overlap with the same GI
        so use tp+fn instead
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
        this is in term of overlap bases.
        '''
        tp = overlap_totalSize
        real = getIntervalLen(ref_intervals)
        # two alternative ways to get the number of bases in reference intervals
        assert real == ref_totalSize
        # merge before counting as query_intervals may be overlapping?
        predicted = getOverlapIntervalSize(query_intervals)

        recall = tp / real
        precision = tp / predicted
        if recall != 0 and precision != 0:
            fmeasure = 2 * recall * precision / (recall + precision)
        else:
            fmeasure = 0
        # append offset error at the end of last line
        (avg_left, avg_right, avg_offset) = getAvgBoundaryError(extentions)

        print 'The number of reference bases: %d' % real

        ############################## PR in #overlap genes ###############################################
        '''
        this is in term of overlap genes.
        '''
        if options.pttfile:
            tp = len(set(overlap_total_genes))
            g_real = len(set(ref_total_genes))
            # print 'num_total_refgenes: %d; g_real: %d'% (num_total_refgenes,g_real)
            # two alternative ways to get the number of genes in reference intervals
            # assert g_real == num_total_refgenes
            # merge before counting as query_intervals may be overlapping?
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

            # other measures
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

        # output FP predictions
        if options.overlap:
            # sort intervals
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
        parser.add_option("-l", "--overlap", dest="overlap", help="output file of overlap intervals")
        # to get the size of genome for calculation of fraction of predicted bases
        parser.add_option("-g", "--genomefile", dest="genomefile", help="input genome file")
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
