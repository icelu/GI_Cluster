#!/usr/bin/env python

# Get the segments by calling GCProfile
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# The genome sequence in FASTA format
#
# Output:
# The positions for a set of intervals
#


import os
import optparse
import subprocess


def GetGenomeSize(genomefile):
    firstLine = open(genomefile).readline()
    assert ('>' in firstLine), "This is not a standard FASTA file!"
    genomesize = 0
    with open(genomefile, 'rb') as fin:
        # skip the first line
        fin.next()
        for line in fin:
            genomesize += len(line.strip())

    return genomesize



def Segment(tool, genomefile, halt=50, minsize=5000):
    Tag   = 1
    Start = 0
    Seg_Points = []
    Primary_Seg_Invals = []

    subprocess.call([tool + 'GCProfile', genomefile, '-t', str(halt), '-i', str(minsize)])

    Seg_File = genomefile.split('.')[0] + '_SegPoints.html'

    for Seg_Line in open(Seg_File):
        if '<TR align="center">' in Seg_Line:
           for Seg_Lines in Seg_Line.split('<TR align="center">')[1:]:
               Seg_Points.append(int(Seg_Lines.split('<div align="center">')[2].split('</div></td>')[0]))
    for Seg_Point in Seg_Points:
        if Tag != len(Seg_Points):
           Primary_Seg_Invals.append((Start + 1, Seg_Point))
           Start = Seg_Point
        else:
           Primary_Seg_Invals.append((Start + 1, Seg_Point))
           Start = Seg_Point
           Primary_Seg_Invals.append((Start + 1, GetGenomeSize(genomefile)))
        Tag += 1

    return  Primary_Seg_Invals


def WriteSegments(intervals, outfile):
    fout = open(outfile,'w')
    for start, end in intervals:
        line='%s\t%s\n' % (start, end)
        fout.write(line)

    fout.close()



if __name__ == '__main__':
        parser = optparse.OptionParser()

        parser.add_option("-f", "--fastafile", dest="fastafile", help="input fasta file of genome sequence")
        parser.add_option("-p", "--progdir", dest="progdir", help="the directory containing GCProfile")
        parser.add_option("-o", "--output", dest="output", help="output file for the merged intervals")
        parser.add_option("-t", "--halt", dest="halt", type="int", default=50, help="the halting parameter for segmentation")
        parser.add_option("-m", "--minl", dest="minl", type="int", default=5000, help="the minimum length of segments")

        options, args = parser.parse_args()

        intervals= Segment(options.progdir, options.fastafile, options.halt, options.minl)
        WriteSegments(intervals, options.output)
