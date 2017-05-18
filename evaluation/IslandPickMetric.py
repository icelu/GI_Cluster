# For comparison of GI tools:
# Given the benchmark dataset B+, B-, the predicted dataset P,
# TP: nuleotides in both B+ and P
# FP: nuleotides in both B- and P
#
# following rules from IslandPick
#
# created by lbx, 2015/4/17
#
# input: 3 sets of intervals
# B+, B-, P
#
# output:
# compare B+ and P --> TP;
# compare B- and P --> FP;
# get the size of all intervals in B+ --> TP+FN
#
#
# calculation of TP and FP are similar:
# do interval intersections of 2 sets of intervals, find the number of overlapped bases.
# a straightforward method is 2-level iteration
# another way:
# 1. find overlapped intervals
# 2. compute overlapping bases



from __future__ import division
import os
import optparse
# from pyroc import *
import copy



def getOverlap(interval_list1, interval_list2):
    num_overlap = 0
    # sort the list to facilitate searching
    # suppose the intervals are not overlapping
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



# def usePyroc(label_list):
#     roc = ROCData(label_list)
#     # roc.auc()
#     # roc.plot(title='ROC Curve')  # Create a plot of the ROC curve
#     # roc.confusion_matrix(0.5, True)
#     matrix = roc.confusion_matrix(1, True)
#     print matrix
#     roc.evaluateMetrics(matrix, do_print=True)


'''
def useSklearnRoc(label_list):
    y_true = []
    y_score = []
    for label in label_list:
       y_true.append(label[0])
       y_score.append(label[1])

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    # print fpr
    # print tpr
    # print thresholds
    # Plot of a ROC curve for a specific class
    plt.figure()
    plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic curve')
    plt.legend(loc="lower right")
    plt.show()
'''


'''
input: genome coordinates of predictions and reference
There are two benchmark datasets, one positive, one negative.
output: label file for counting
                data: The data is a list l of tuples t (l = [t_0,t_1,...t_n]) where:
                      t[0] = 1 for positive class and 0 for negative class
                      t[1] = a score
                      t[2] = any label (optional)

pos_points = [x for x in self.data if x[1] >= threshold]
neg_points = [x for x in self.data if x[1] < threshold]

""" Calculates the number of false positives, true positives, false negatives and true negatives """
tp_count = len([x for x in pos_data if x[0] == 1])
fp_count = len([x for x in pos_data if x[0] == 0])
fn_count = len([x for x in neg_data if x[0] == 1])
tn_count = len([x for x in neg_data if x[0] == 0])

'''
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
