#!/usr/bin/env python

# Commonly used functions for analyzing genomic intervals
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

############################## interval operation ########################
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.quicksect import IntervalNode


def get_orig_intervals(intervalfile):
    '''
    Return a list of interval tuples without type conversion
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = (fields[0])
            end = (fields[1])
            coord = (start, end)
            intervals.append(coord)
    return intervals


def get_intervals(intervalfile):
    '''
    Return a list of interval tuples
    '''
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            coord = (start, end)
            intervals.append(coord)
    return intervals


def get_intervals_contigs(intervalfile):
    '''
    Input:  intervalfile -- file containing a list of intervals, with row has format "id_start\tid_end", such as "1_1     1_5000"
    Return a list of intervals for each contig
    '''
    intervals = {}
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            p1 = fields[0]
            mark = p1.index('_')
            start = int(p1[mark + 1:])
            p2 = fields[1]
            end = int(p2[mark + 1:])
            contig_id = int(p1[0:mark])
            coord = (start, end)
            if contig_id not in intervals.keys():
                intervals[contig_id] = [coord]
            else:
                intervals.setdefault(contig_id, []).append(coord)
    return intervals


def find(start, end, tree):
    '''
    Find all the intervals overlapping with the query interval
    Returns a list with the overlapping intervals
    '''
    out = []
    tree.intersect(start, end, lambda x: out.append(x))
    # x.other may be none if no score is assigned to the interval
    return [(x.start, x.end, x.other) for x in out]


def get_window_tree(selectRes):
    '''
    Build an interval tree for all the query intervals
    '''
    if len(selectRes[0]) == 3:
        start, end, score = selectRes[0]
        tree = IntervalNode(start, end, other=score)
        # build an interval tree from the rest of the data
        for start, end, score in selectRes[1:]:
            tree = tree.insert(start, end, other=score)
    else:
        start, end = selectRes[0]
        tree = IntervalNode(start, end, other=(end - start + 1))
        # build an interval tree from the rest of the data
        for start, end in selectRes[1:]:
            # use size as the 3rd column
            tree = tree.insert(start, end, other=(end - start + 1))
    return tree


def get_overlap(a, b):
    '''
    From http://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
    Find the number of overlapping bases between two intervals
    Add by 1 for 1-based interval
    '''
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def get_overlap_interval(a, b):
    '''
    Get the overlap interval of two invervals
    '''
    start = max(a[0], b[0])
    end = min(a[1], b[1])
    if end >= start:
        return (start, end)
    else:
        return None


def get_interval_length(intervals):
    '''
    Assume intervals are not overlapping.
    Add up the bases of all the intervals.
    '''
    length = 0
    for interval in intervals:
        length += interval[1] - interval[0] + 1
    return length


def merge(intervals):
    '''
    Given a list of intervals, merge overlapping ones and return results by generator
    '''
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


def merge_ref(intervals):
    '''
    Given a list of intervals, merge overlapping ones and return a list of non-overlapping intervals
    '''
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


def merge_score(intervals):
    '''
    Merge regions with score
    '''
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
