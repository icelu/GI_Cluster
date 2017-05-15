
# Lu Bingxin, 2017.2.12

from __future__ import division
import os
import optparse

'''
input:
for each predicted region:
    DNA sequence + upstream + downstream
estimate the change point by Bayesian methods

Refine the boundary of a segment by trna and repeats
for each region, check the trna and repeat overlapping with it.


output:
new start and end position
'''


def get_input_interval(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[0]), int(fields[1]))
            intervals.append(coord)

    # print len(intervals)

    return intervals


def get_rna(infile):
    rna_dict = {}
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[0]), int(fields[1]))
            num = int(fields[2])
            rnas = fields[3].split(';')
            rna_pos = []
            for rp in rnas:
                items = rp.split(',')
                #print items
                # remove brackets
                rpos = (int(items[0][1:]), int(items[1][:-1]))
                rna_pos.append(rpos)
            rna_dict[coord] = rna_pos
        #print rna_dict
    return rna_dict


def get_repeat(infile):
    repeat_dict = {}
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[0]), int(fields[1]))
            num = int(fields[2])
            repeats = fields[3].split(';')
            # keep only the position
            repeats_pos = []
            for rp in repeats:
                items = rp.split(',')
                #print items
                rpos = (int(items[1]), int(items[2]))
                repeats_pos.append(rpos)
            repeat_dict[coord] = repeats_pos
    #print repeat_dict
    return repeat_dict


'''

'''
def refine_boundary(intervals, rna_dict, repeat_dict, threshold=500):
    new_intervals = []
    for coord in intervals:
        (start, end) = coord
        if coord in rna_dict.keys():
            rnas = rna_dict[coord]
        else:
            rnas = []
        if coord in repeat_dict.keys():
            repeats = repeat_dict[coord]
        else:
            repeats = []
        total = rnas + repeats
        if len(total) > 0:
            #print total
            sorted_total = sorted(total, key=lambda x: (int(x[0]), int(x[1])))
            # compare the distance between the interval
            rs = sorted_total[0][0]
            re = sorted_total[-1][1]
            # new boundary: the outmost position of trna or repeat
            if(abs(start - rs) < threshold):
                ns = rs
            else:
                ns = start
            if(abs(end - re) < threshold):
                ne = re
            else:
                ne = end
            if ns == rs or ne == re:
                print '%d\t%d\t--\t%d\t%d' % (start, end, ns, ne)
            new_intervals.append((ns, ne))
        else:
            new_intervals.append((start, end))
    return new_intervals


def write_output(intervals, outfile):
    with open(outfile, 'w') as fout:
        for s, e in intervals:
            line = '%d\t%d\n' % (s, e)
            fout.write(line)


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--interval", dest="interval",
                      help="input interval file1 (predicted GIs)")
    parser.add_option("-r", "--repeats", dest="repeats",
                      help="input file containing repeats overlapping with the input regions")
    parser.add_option("-t", "--trnas", dest="trnas",
                      help="input file containing trnas overlapping with the input regions")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="int", default="500",
                      help="the cutoff number of base distances between tRNA/repeat and a region")
    parser.add_option("-o", "--output", dest="output",
                      help="output file of intervals with new boundaries")

    options, args = parser.parse_args()

    intervals = get_input_interval(options.interval)
    rna_dict = get_rna(options.trnas)
    repeat_dict = get_repeat(options.repeats)
    new_intervals = refine_boundary(
        intervals, rna_dict, repeat_dict, threshold=options.cutoff)
    write_output(new_intervals, options.output)
