# The adjacent candidate segments were joined and the boundary is relocated by the outmost ORF if an ORF intersects with a segment.
# Note: the segment will only be extended
# For features: sum the numbers, recompute the percentages
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

import os
import optparse

import sys
import os
parentdir = os.path.dirname(os.path.dirname(sys.argv[0]))
sys.path.insert(0, parentdir)
from util.quicksect import IntervalNode
from util.interval_operations import get_intervals, get_intervals_contigs, find, get_window_tree
from util.parse_sequence import get_contig_IDs, get_contig_gene_IDs

# Suppose the gene positions are at 2nd and 3rd columns
def get_genes(intervalfile):
    intervals = []
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            coord = (int(fields[1]), int(fields[2]))
            intervals.append(coord)

    return intervals


# Find genes in each sequence
def get_genes_contig(intervalfile, gene_id_mapping, contig_id_mapping):
    intervals = {}
    with open(intervalfile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            gid = int(fields[0])
            gene_name = gene_id_mapping[int(gid)]
            # Get contig name from gene name
            contig_name = gene_name[0:gene_name.rfind('_')]
            contig_id = contig_id_mapping[contig_name]
            coord = (int(fields[1]), int(fields[2]))
            if contig_id not in intervals.keys():
                intervals[contig_id] = [coord]
            else:
                intervals.setdefault(contig_id, []).append(coord)

    return intervals


def merge_intervals(intervals):
    if len(intervals) < 0:
        return
    intervals = sorted(intervals, key=lambda x: (int(x[0]), int(x[1])))
    saved = list(intervals[0])
    for st, en in intervals:
        if st - 1 <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)


# Merge overlapped regions or regions with small gap
def merge_intervals_offset(intervals, allow_gap, gap_len):
    intervals = sorted(intervals, key=lambda x: (int(x[0]), int(x[1])))
    merged_intervals = []
    saved = list(intervals[0])

    offset = 1
    if allow_gap:
        offset = gap_len

    for st, en in intervals:
        if st - offset <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)


# Merge overlapped regions or regions with small gap
def merge_intervals_contigs(intervals_dict, allow_gap, gap_len):
    merged_intervals = []
    for cid, intervals in intervals_dict.items():
        m_intervals = list(merge_intervals_offset(intervals, allow_gap, gap_len))
        for start, end in m_intervals:
            ns = '_'.join([str(cid), str(start)])
            ne = '_'.join([str(cid), str(end)])
            merged_intervals.append((ns, ne))
    return merged_intervals


def write2file(filename, list):
    outfile = open(filename, 'w')

    for value in list:
        if len(value) == 3:
            line = '%s\t%s\t%s\n' % (
                str(value[0]), str(value[1]), str(value[2]))
        else:
            line = '%s\t%s\n' % (str(value[0]), str(value[1]))
        outfile.write(line)

    outfile.close()


def extend_boundary(intervals, genes):
    # build an interval to facilitate querying
    tree = get_window_tree(genes)
    new_intervals = []
    for start, end in intervals:
        overlap = find(start, end, tree)
        if len(overlap) > 0:
            # find the boundary coordinates of the intervals
            sorted_overlap = sorted(
                overlap, key=lambda x: (int(x[0]), int(x[1])))
            ostart = sorted_overlap[0][0]
            oend = sorted_overlap[-1][1]
            intv_size = sorted_overlap[0][1] - sorted_overlap[0][0] + 1
            hang_size = start - ostart + 1
            intv_size1 = sorted_overlap[-1][1] - sorted_overlap[-1][0] + 1
            hang_size1 = oend - end + 1
            if ostart < start and hang_size < intv_size / 2:
                nstart = ostart
            else:
                nstart = start
            if oend > end and hang_size1 < intv_size1 / 2:
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


def extend_boundary_contig(intervals, genes_dict):
    # build an interval to facilitate querying
    tree_dict = {}
    for cid, genes in genes_dict.items():
        tree = get_window_tree(genes)
        tree_dict[cid] = tree
    new_intervals = {}  # Use dicionary to facilitate the merging of overlapping regions
    for p1, p2 in intervals:
        # Suppose the format for start/end is cid_coord
        mark = p1.index('_')
        start = int(p1[mark + 1:])
        end = int(p2[mark + 1:])
        contig_id = int(p1[0:mark])
        tree = tree_dict[contig_id]
        overlap = find(start, end, tree)
        if len(overlap) > 0:
            # find the boundary coordinates of the intervals
            print 'intervals with overlapping:'
            print p1
            print p2
            print overlap
            sorted_overlap = sorted(
                overlap, key=lambda x: (int(x[0]), int(x[1])))
            ostart = sorted_overlap[0][0]
            oend = sorted_overlap[-1][1]
            intv_size = sorted_overlap[0][1] - sorted_overlap[0][0] + 1
            hang_size = start - ostart + 1
            intv_size1 = sorted_overlap[-1][1] - sorted_overlap[-1][0] + 1
            hang_size1 = oend - end + 1
            if ostart < start and hang_size < intv_size / 2:    # More than half of the gene is outside this region
                nstart = ostart
            else:
                nstart = start
            if oend > end and hang_size1 < intv_size1 / 2:
                nend = oend
            else:
                nend = end
            coord = (nstart, nend)
            # ns = str(contig_id) + '_' + str(nstart)
            # ne = str(contig_id) + '_' + str(nend)
            # new_intervals.append((ns, ne))
        else:
            # new_intervals.append((p1, p2))
            coord = (start, end)
        if contig_id not in new_intervals.keys():
            new_intervals[contig_id] = [coord]
        else:
            new_intervals.setdefault(contig_id, []).append(coord)
    return new_intervals


if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-g", "--genefile", dest="genefile",
                      help="input file of genes and their locations")
    parser.add_option("-i", "--gifile", dest="gifile",
                      help="input file of predicted GIs")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file")
    parser.add_option("-p", dest="allow_gap", default=False, action="store_true",
                      help="allow to combine adjacent intervals that are very close")
    parser.add_option("-l", "--gap_len", dest="gap_len", type='int',
                      default=2500, help="threshold to merge adjacent intervals")
    parser.add_option("-a", "--has_gene", dest="has_gene", action='store_true', default=False,
                      help="The gene predictions are available")
    parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true',
                      default=False, help="Analyze contigs from unassembled genomes")
    parser.add_option("-m", "--genome_file", dest="genome_file",
                      help="input genome file in fasta format")
    parser.add_option("-d", "--gid_file", dest="gid_file",
                      help="input file of predicted gene IDs")
    (options, args) = parser.parse_args()

    directory = os.path.dirname(os.path.realpath(options.gifile))
    suffix = os.path.basename(options.gifile)

    if options.is_contig:
        orig_intervals = get_intervals_contigs(options.gifile)
        count_orig = 0
        count_merged = 0
        merged_intervals = []
        for id, intervals in orig_intervals.items():
            count_orig += len(intervals)
            m_intervals = list(merge_intervals(intervals))
            count_merged += len(m_intervals)
            for start, end in m_intervals:
                ns = '_'.join([str(id), str(start)])
                ne = '_'.join([str(id), str(end)])
                merged_intervals.append((ns, ne))
        print 'The number of intevals before merging adjacent ones: %d' % count_orig
        print 'The number of intevals after merging adjacent ones: %d' % count_merged
        if options.has_gene:
            contig_id_mapping = get_contig_IDs(options.genome_file)
            gene_id_mapping = get_contig_gene_IDs(options.gid_file)
            genes_dict = get_genes_contig(options.genefile, gene_id_mapping, contig_id_mapping)
            new_intervals_dict = extend_boundary_contig(merged_intervals, genes_dict)
            # Merge intervals again to avoid overlapping regions and combine close intervals
            merged_new_intervals = merge_intervals_contigs(
                new_intervals_dict, options.allow_gap, options.gap_len)
            print 'The number of intevals after merging close ones with gap %d bp: %d' % (options.gap_len, len(merged_new_intervals))
            write2file(os.sep.join([directory, 'merged_']
                                   ) + suffix, merged_new_intervals)
        else:
            write2file(os.sep.join(
                [directory, 'merged_']) + suffix, merged_intervals)
    else:
        orig_intervals = get_intervals(options.gifile)
        merged_intervals = list(merge_intervals(orig_intervals))
        print 'The number of intevals before merging adjacent ones: %d' % len(orig_intervals)
        print 'The number of intevals after merging adjacent ones: %d' % len(merged_intervals)

        if options.has_gene:
            genes = get_genes(options.genefile)
            new_intervals = extend_boundary(merged_intervals, genes)
            merged_new_intervals = list(merge_intervals_offset(
                new_intervals, options.allow_gap, options.gap_len))
            print 'The number of intevals after merging close ones with gap %d bp: %d' % (options.gap_len, len(merged_new_intervals))
            write2file(os.sep.join([directory, 'merged_']
                                   ) + suffix, merged_new_intervals)
        else:
            write2file(os.sep.join(
                [directory, 'merged_']) + suffix, merged_intervals)
