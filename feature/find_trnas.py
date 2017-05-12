#!/usr/bin/python
# Lu Bingxin, 2017.2.27
import optparse
from quicksect import IntervalNode

'''
For
Input:
1. predicted GIs (genomic intervals)
2. predicted trnas

Sample command to predict tRNA by tRNAscan-SE
output_dir=/home/b/bingxin/genome/StypiCT18
organism=NC_010161
tRNAscan-SE -B --frag $output_dir/boundary/$organism.trna_frag -o $output_dir/boundary/$organism.trna_pred -m $output_dir/boundary/$organism.trna_stat --brief $output_dir/$organism.fna
less $output_dir/boundary/$organism.trna_pred | cut -f3-4 > $output_dir/boundary/$organism.pred_trna

Output:
predicted GIs annotated with nearby trnas

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


################################### find trnas in a segment #######################################
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
def parse_trnas(infile, intervals, offset):
    tree = getWindowTree(intervals)
    trna_dict = {}

    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            pos1 = int(fields[0])
            pos2 = int(fields[1])
            if pos1 > pos2:
                tstart = pos2
                tend = pos1
            else:
                tstart = pos1
                tend = pos2
            # extend trna in case it is nearby
            overlap = find(tstart-offset, tend+offset, tree)
            if len(overlap) > 0:
                for intv in overlap:
                    istart, iend, iscore = intv
                    if intv not in trna_dict.keys():
                            trna_dict[intv]=[(pos1, pos2)]
                    else:
                            trna_dict.setdefault(intv,[]).append((pos1, pos2))
    return trna_dict


def write_trnas(outfile, trna_dict):
    with open(outfile, 'w') as fout:
        for key, value in trna_dict.items():
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
        parser.add_option("-t", "--trnas", dest="trnas", help="input file containing trnas in a genome")
        parser.add_option("-s", "--offset", dest="offset", type=int, default=100, help="finding tRNAs certain bases away from a position")
        parser.add_option("-o", "--output", dest="output", help="output file of intervals with overlapping trnas")

        options, args = parser.parse_args()

        intervals = get_intervals(options.interval1)
        trna_dict = parse_trnas(options.trnas, intervals, options.offset)
        write_trnas(options.output, trna_dict)
