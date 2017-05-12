#!/usr/bin/python

# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
# Version : 1.0

import optparse
import itertools
'''
Use a sliding window to generate segments along the genome
'''


def isheader(line):
    return line[0] == '>'

def get_genome_len(genome_file):
    with open(genome_file, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                sequence = ''.join(line.strip() for line in group)
                return len(sequence)


def sliding_window(outfile, seqlen, width=5000, step=5000):
    # for the last window if != winsize (mostly true),
	# take different start point to the end of the genome, but again with size=winsize
    with open(outfile,'w') as fout:
        for i in range(0, seqlen, step):
        	# last window
            if i + width > seqlen:
            	j = seqlen
            	i = seqlen - width
            else:
            	j = i + width

            line='%d\t%d\n' % ((i + 1, j))
            fout.write(line)

            if j == seqlen: break





if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--genome_file", dest="genome_file",
                      help="The input genome files in fasta format")
    parser.add_option("-w", "--width", dest="width", type=int, default=5000,
                      help="The width of sliding windows")
    parser.add_option("-s", "--step", dest="step", type=int, default=5000,
                      help="The step of sliding windows")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file of segment position")

    options, args = parser.parse_args()

    seqlen = get_genome_len(options.genome_file)
    sliding_window(options.output, seqlen, options.width, options.step)
