# The adjacent candidate segments were joined and the boundary is relocated by the outmost ORF if an ORF intersects with a segment.
# Note: the segment will only be extended
# For features: sum the numbers, recompute the percentages


import os
import optparse
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.quicksect import IntervalNode


def get_input_interval(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[0]), int(fields[1]))
            intervals.append(coord)

    # print len(intervals)
    return intervals


# gene positions in different fields
def get_genes(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[1]), int(fields[2]))
            intervals.append(coord)

    # print len(intervals)
    return intervals


def merge_intervals(intervals):
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
merge overlapped regions or regions with small gap
'''
def merge_intervals_offset(intervals, allow_gap, gap_len):
    intervals = sorted(intervals, key=lambda x : (int(x[0]), int(x[1])))
    merged_intervals = []
    saved = list(intervals[0])

    offset = 1
    if allow_gap:
		offset = gap_len

    for st, en in intervals:
        if st - offset <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            # only add the interval to the result when not overlapping with adjacent regions
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)



def writeListOfTupleToFile(filename, list):
    outfile = open(filename, 'w')

    for value in list:
        if len(value) == 3:
            line = '%s\t%s\t%s\n' % (str(value[0]), str(value[1]), str(value[2]))
        else:
            line = '%s\t%s\n' % (str(value[0]), str(value[1]))
        outfile.write(line)

    outfile.close()


'''
Build an interval tree for all the query intervals
'''
def getWindowTree(selectRes):
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


'''
Find all the intervals overlapping with the query interval
'''
def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect(start, end, lambda x: out.append(x))
    # x.other may be none if no score is assigned to the interval
    return [ (x.start, x.end, x.other) for x in out ]


def extend_boundary(intervals, genes):
    # build an interval to facilitate querying
    tree = getWindowTree(genes)
    new_intervals = []
    for start, end in intervals:
        overlap = find(start, end, tree)
        if len(overlap) > 0:
            # find the boundary coordinates of the intervals
            sorted_overlap = sorted(overlap, key=lambda x : (int(x[0]), int(x[1])))
            ostart = sorted_overlap[0][0]
            oend = sorted_overlap[-1][1]
            intv_size = sorted_overlap[0][1] - sorted_overlap[0][0] + 1
            hang_size = start - ostart + 1
            intv_size1 = sorted_overlap[-1][1] - sorted_overlap[-1][0] + 1
            hang_size1 = oend - end + 1
            if ostart < start and hang_size < intv_size/2:
                nstart = ostart
            else:
                nstart = start
            if oend > end and hang_size1 < intv_size1/2:
                nend = oend
            else:
                nend = end
            new_intervals.append((nstart, nend))
#             intersects = []
#             for ol in overlap:
#                intersects.append(ol[0])
#                intersects.append(ol[1])
#             minCoord = min(intersects)
#             maxCoord = max(intersects)
        else:
            new_intervals.append((start, end))
    return new_intervals


if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-g", "--genefile", dest="genefile", help="input file of genes and their locations")
    parser.add_option("-i", "--gifile", dest="gifile", help="input file of predicted GIs")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file")
    parser.add_option("-a", dest="allow_gap", default=False, action="store_true", help="allow to combine adjacent intervals that are very close")
    parser.add_option("-l", "--gap_len", dest="gap_len", type='int', default=2500, help="threshold to merge adjacent intervals")
    (options, args) = parser.parse_args()

    directory = os.path.dirname(os.path.realpath(options.gifile))
    suffix = os.path.basename(options.gifile)
    # i = base.find('_')
    # suffix = base[i + 1:]
    # print directory, suffix

    origRes = get_input_interval(options.gifile)
    print 'The number of intevals before merging adjacent ones: %d' % len(origRes)

    mergedRes = list(merge_intervals(origRes))
    print 'The number of intevals after merging adjacent ones: %d' % len(mergedRes)

    genes = get_genes(options.genefile)
    new_intervals = extend_boundary(mergedRes, genes)

    merged_new_intervals = merge_intervals_offset(new_intervals, options.allow_gap, options.gap_len)

    # merge intervals again to avoid overlapping regions and combine close intervals
    writeListOfTupleToFile(os.sep.join([directory, 'merged_']) + suffix, merged_new_intervals)
