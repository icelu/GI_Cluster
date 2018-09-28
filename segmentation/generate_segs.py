#!/usr/bin/env python

# Use a sliding window to generate segments along the genome
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


import optparse
import itertools




def isheader(line):
    return line[0] == '>'


def get_genome_len(genome_file):
    with open(genome_file, 'r') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                sequence = ''.join(line.strip() for line in group)
                return len(sequence)


def sliding_window(outfile, seqlen, width=5000, step=5000):
    with open(outfile, 'w') as fout:
        for i in range(0, seqlen, step):
            # last window
            if i + width > seqlen:
                j = seqlen
                i = seqlen - width
            else:
                j = i + width

            line = '%d\t%d\n' % ((i + 1, j))
            fout.write(line)

            if j == seqlen:
                break


# Switch non-standard character to 'AGCT'
def switch_to_AGCT(Nuc):
    SWITCH = {'R': lambda: random.choice('AG'),
              'Y': lambda: random.choice('CT'),
              'M': lambda: random.choice('AC'),
              'K': lambda: random.choice('GT'),
              'S': lambda: random.choice('GC'),
              'W': lambda: random.choice('AT'),
              'H': lambda: random.choice('ATC'),
              'B': lambda: random.choice('GTC'),
              'V': lambda: random.choice('GAC'),
              'D': lambda: random.choice('GAT'),
              'N': lambda: random.choice('AGCT')}
    return SWITCH[Nuc]()


def standardize_DNASeq(genome):
    '''
    Replace ambiguous bases with ATCG
    '''
    for i in ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']:
        if i in genome:
            genome = genome.replace(i, switch_to_AGCT(i))
    return genome


def parse_contigs(infile, outfile, width=5000, step=5000):
    '''
    input: .fna file, fasta format, DNA sequence for each contigs
    output: The positions for a set of non-overlapping windows
    '''
    c = 0
    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                c += 1
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                seqlen = len(sequence)
                if (seqlen >= width):  # Use sliding window
                    for i in range(0, seqlen, step):
                        # last window
                        if i + width > seqlen:
                            j = seqlen
                            i = seqlen - width
                        else:
                            j = i + width

                        line = '%d_%d\t%d_%d\n' % ((c, i + 1, c, j))
                        fout.write(line)

                        if j == seqlen:
                            break
                else:   # Use original contig
                    line = '%d_%d\t%d_%d\n' % ((c, 1, c, seqlen))
                    fout.write(line)


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input",
                      help="The input DNA sequence file in fasta format")
    parser.add_option("-w", "--width", dest="width", type=int, default=5000,
                      help="The width of sliding windows")
    parser.add_option("-s", "--step", dest="step", type=int, default=5000,
                      help="The step of sliding windows")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file of segment position")
    parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true', default=False,
                      help="Analyze contigs from unassembled genomes")

    options, args = parser.parse_args()

    if options.is_contig:
        # read contigs and partition contigs of size > 'width' with sliding window
        parse_contigs(options.input, options.output, width=options.width, step=options.step)
    else:
        seqlen = get_genome_len(options.input)
        sliding_window(options.output, seqlen, options.width, options.step)
