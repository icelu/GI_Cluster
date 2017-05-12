#!/usr/bin/python
##=========================================================
# Script for finding repeats in a genome
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
# Version : 1.0

import optparse
from quicksect import IntervalNode

'''
For
Input:
1. predicted GIs (genomic intervals)
2. repseek output

Sample command to run repseek:
output_dir=/home/b/bingxin/genome/StypiCT18
organism=NC_010161
repseek -l 15 -O 0 -r $output_dir/boundary/$organism.repseek $output_dir/$organism.fna

Output:
predicted GIs annotated with nearby repeats

'''


############################## interval operation ########################################
'''
Find all the intervals overlapping with the query interval
'''
def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect(start, end, lambda x: out.append(x))
    # x.other may be none if no score is assigned to the interval
    return [ (x.start, x.end, x.other) for x in out ]

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
From http://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
Find the number of overlapping bases between two intervals
Add by 1 for 1-based interval
'''
def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


################################### find repeats in a segment #######################################
'''
return a list of interval tuples
'''
def get_intervals(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            coord = (start, end)
            intervals.append(coord)
    return intervals

'''

'''
def parse_repseek_output(infile, intervals):
    tree = getWindowTree(intervals)
    repeat_dict = {}

    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            type = fields[0]
            pos1 = int(fields[1])
            pos2 = int(fields[2])
            len1 = int(fields[3])
            len2 = int(fields[4])
            spacer = int(fields[5])
            seed = fields[6]
            identity = float(fields[7])
            score = float(fields[8])
            meanR = float(fields[9])
            modeR = float(fields[10])
            fraction = float(fields[11])
            # for a distant repeat, find its neighbor GIs
            if spacer >= 4500 and spacer <= 600000:
                end = pos2+len2
                overlap = find(pos1, end, tree)
                #line='%s\t%s\t%s\n' % ()
                if len(overlap) > 0:
                    for intv in overlap:
                        istart, iend, iscore = intv
                        if abs(istart-pos1)<3000 and abs(iend-end)<3000:
                            if intv not in repeat_dict.keys():
                                repeat_dict[intv]=[(type, pos1, end, len1, len2)]
                            else:
                                repeat_dict.setdefault(intv,[]).append((type, pos1, end, len1, len2))
    return repeat_dict


def write_repeats(outfile, repeat_dict):
    with open(outfile, 'w') as fout:
        for key, value in repeat_dict.items():
            #line='%s\t%s\t%d\t%s\n' % (key[0], key[1], len(value),value)
            line='%s\t%s\t%d\t' % (key[0], key[1], len(value))
            fout.write(line)
            line=''
            for pos in value:
                line+='%s;' % (str(pos))
            line = line[:-1]  + '\n'
            fout.write(line)



if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option("-i", "--interval1", dest="interval1", help="input interval file1 (predicted GIs)")
        parser.add_option("-r", "--repeats", dest="repeats", help="input file containing repeats in a genome")
        parser.add_option("-o", "--output", dest="output", help="output file of intervals with overlapping repeats")

        options, args = parser.parse_args()

        intervals = get_intervals(options.interval1)
        repeat_dict = parse_repseek_output(options.repeats, intervals)
        write_repeats(options.output, repeat_dict)
