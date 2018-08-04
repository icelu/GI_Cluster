#!/usr/bin/env python

# Merge adjacent or overlapping intervals
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

import optparse
import os


def merge(intevals):
    intevals = sorted(intevals, key=lambda x : (int(x[0]), int(x[1])))
    saved = list(intevals[0])
    for st, en in intevals:
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)



def getIntervals(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(float(fields[0])), int(float(fields[1])))
            intervals.append(coord)
    return intervals


def writeListOfTupleToFile(filename, list):
    outfile = open(filename, 'w')

    for value in list:
        if len(value) == 3:
            line = '%s\t%s\t%s\n' % (str(value[0]), str(value[1]), str(value[2]))
        else:
            line = '%s\t%s\n' % (str(value[0]), str(value[1]))
        outfile.write(line)

    outfile.close()



if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--interval", dest="interval", help="input interval file of genomic regions")
    (options, args) = parser.parse_args()

    directory = os.path.dirname(os.path.realpath(options.interval))
    base = os.path.basename(options.interval)
    # i = base.find('_')
    # suffix = base[i + 1:]

    selectRes = getIntervals(options.interval)
    print 'intevals length before merge: %s' % str(len(selectRes))

    mergedRes = list(merge(selectRes))
    print 'intevals length after merge: %s' % str(len(mergedRes))
    writeListOfTupleToFile(os.sep.join([directory, 'merged_']) + base, mergedRes)
