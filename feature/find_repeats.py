# For finding repeats around each genomic interval
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# 1. predicted GIs (genomic intervals)
# 2. repseek output
#
# Output:
# genomic interval annotated with nearby repeats
#


import optparse
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.interval_operations import get_intervals, find, get_window_tree, get_overlap


################################### find repeats in a segment #######################################
def parse_repseek_output(infile, intervals):
    tree = get_window_tree(intervals)
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
