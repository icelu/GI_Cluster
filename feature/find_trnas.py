# For finding tRNAs around a genomic interval
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# 1. predicted GIs (genomic intervals)
# 2. predicted trnas
#
# Output:
# predicted GIs annotated with nearby trnas
#

import optparse
import sys, os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.parse_sequence import switch_to_AGCT, standardize_DNASeq, isheader, get_contig_IDs
from util.interval_operations import get_intervals, get_intervals_contigs, find, get_window_tree, get_overlap
import itertools


################################### find trnas in a segment ##############
def parse_trnas(infile, intervals, offset):
    tree = get_window_tree(intervals)
    trna_dict = {}

    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            # Suppose the start and end positions are in 3rd and 4th columns
            pos1 = int(fields[2])
            pos2 = int(fields[3])
            if pos1 > pos2:
                tstart = pos2
                tend = pos1
            else:
                tstart = pos1
                tend = pos2
            # Extend search region in case a tRNA is nearby
            overlap = find(tstart - offset, tend + offset, tree)
            if len(overlap) > 0:
                for intv in overlap:
                    istart, iend, iscore = intv
                    if intv not in trna_dict.keys():
                        trna_dict[intv] = [(pos1, pos2)]
                    else:
                        trna_dict.setdefault(intv, []).append((pos1, pos2))
    return trna_dict


def write_trnas(outfile, trna_dict):
    with open(outfile, 'w') as fout:
        for key, value in trna_dict.items():
            line = '%s\t%s\t%d\t' % (key[0], key[1], len(value))
            fout.write(line)
            line = ''
            for pos in value:
                line += '%s;' % (str(pos))
            line = line[:-1] + '\n'
            fout.write(line)

# Find tRNAs around segments obtained from contigs
def parse_trnas_contigs(infile, id_mapping, intervals, offset):
    trna_dict = {}

    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            # Suppose the start and end positions are in 3rd and 4th columns
            name = fields[0].strip()
            id = id_mapping[name]
            pos1 = int(fields[2])
            pos2 = int(fields[3])
            if pos1 > pos2:
                tstart = pos2
                tend = pos1
            else:
                tstart = pos1
                tend = pos2

            tree = get_window_tree(intervals[id])
            # Extend search region in case a tRNA is nearby
            overlap = find(tstart - offset, tend + offset, tree)
            if len(overlap) > 0:
                for intv in overlap:
                    istart, iend, iscore = intv
                    key = (str(id) + '_' + str(istart), str(id) + '_' + str(iend))
                    if intv not in trna_dict.keys():
                        trna_dict[key] = [(pos1, pos2)]
                    else:
                        trna_dict.setdefault(key, []).append((pos1, pos2))

    return trna_dict


##########################################################################
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input",
                      help="Input file containing a list of intervals (predicted GIs)")
    parser.add_option("-t", "--trnas", dest="trnas",
                      help="Input file containing tRNAs in a genome")
    parser.add_option("-s", "--offset", dest="offset", type=int, default=100,
                      help="Finding tRNAs certain bases away from a position")
    parser.add_option("-o", "--output", dest="output",
                      help="Output file of intervals with overlapping trnas")
    parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true', default=False,
                      help="Analyze contigs from unassembled genomes")
    parser.add_option("-g", "--genome_file", dest="genome_file",
                      help="input genome file in fasta format")
    options, args = parser.parse_args()

    if options.is_contig:
        # Search on different contigs separately
        intervals = get_intervals_contigs(options.input)
        # print intervals
        id_mapping = get_contig_IDs(options.genome_file)
        # print id_mapping
        trna_dict = parse_trnas_contigs(
            options.trnas, id_mapping, intervals, options.offset)
        # print trna_dict
        write_trnas(options.output, trna_dict)
    else:
        intervals = get_intervals(options.input)
        trna_dict = parse_trnas(options.trnas, intervals, options.offset)
        write_trnas(options.output, trna_dict)
