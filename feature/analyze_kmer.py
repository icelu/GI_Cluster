# Given DNA sequence of a genome or a gene or a genomic region, analyze its k-mer spectrum
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg

# Input:
# 1. gene sequence file, fasta format
# 2. genome sequence file, fasta format
#
# Output:
# Measures of k-mer:
# for each gene, the differences between its k-mer spectrum and the k-mer
# spectrum of whole genome


from __future__ import division
import itertools
import optparse
import re
import numpy as np
import math
import random


################################# k-mer operations #######################

# Conversion table form nucleotide to binary.
_nucleotide_to_binary = {
    'A': 0x00, 'a': 0x00,
    'C': 0x01, 'c': 0x01,
    'G': 0x02, 'g': 0x02,
    'T': 0x03, 't': 0x03
}

# Conversion table form binary to nucleotide.
_binary_to_nucleotide = {
    0x00: 'A',
    0x01: 'C',
    0x02: 'G',
    0x03: 'T'
}


############################# Most of  the following functions are from kM

###################### Extracting k-mers ####################
def dna_to_binary(sequence):
    """
    Convert a string of DNA to an integer.

    :arg str sequence: DNA sequence.

    :return: Binary representation of `sequence`.
    :rtype: int
    """
    result = 0x00

    for i in sequence:
        result <<= 2
        result |= _nucleotide_to_binary[i]

    return result


def binary_to_dna(number, length):
    """
    Convert an integer to a DNA string.

    :arg int number: Binary representation of a DNA sequence.

    :returns: DNA string corresponding to `number`.
    :rtype: str
    """
    sequence = ""

    for i in range(length):
        sequence += _binary_to_nucleotide[number & 0x03]
        number >>= 2

    return sequence[::-1]


def reverse_complement(number, length):
    """
    Calculate the reverse complement of a DNA sequence in a binary
    representation.

    :arg int number: Binary representation of a DNA sequence.

    :return: Binary representation of the reverse complement of the
      sequence corresponding to `number`.
    :rtype: int
    """
    number = ~number
    result = 0x00

    for i in range(length):
        result = (result << 2) | (number & 0x03)
        number >>= 2

    return result


def kmer_from_sequences(sequences, length):
    """
    Create a *k*-mer profile from `sequences` by counting all *k*-mers in
    each sequence.

    :arg sequences: An iterable of string sequences.
    :type sequences: iterator(str)
    :arg int length: Length of the *k*-mers.

    :return: A *k*-mer profile.
    """
    number = 4 ** length
    bitmask = number - 0x01
    counts = [0] * number
    alphabet = re.compile('[^' + ''.join(_nucleotide_to_binary) + ']')

    # for sequence in sequences:
    for part in alphabet.split(sequences):
        # print part
        if len(part) >= length:
            binary = 0x00

            # Calculate the binary representation of a k-mer.
            for i in part[:length]:
                binary = (binary << 2) | _nucleotide_to_binary[i]
            counts[binary] += 1

            # Calculate the binary representation of the next k-mer.
            for i in part[length:]:
                binary = ((binary << 2) |
                          _nucleotide_to_binary[i]) & bitmask
                counts[binary] += 1

    return np.array(counts, dtype='int64')


###################### distance metrics ####################
# Pairwise distance functions. Arguments should be of type `numpy.ndarray`.
pairwise = {
    "prod": lambda x, y: abs(x - y) / ((x + 1) * (y + 1)),
    "sum": lambda x, y: abs(x - y) / (x + y + 1)
}


def vector_length(vector):
    """
    Calculate the Euclidean length of a vector.

    :arg array_like vector: A vector.

    :return: The length of `vector`.
    :rtype: float
    """
    # Note: This seems faster than `numpy.linalg.norm`.
    return np.sqrt(np.dot(vector, vector))


def multiset(left, right, pairwise=pairwise['prod']):
    """
    Calculate the multiset distance between two vectors.

    :arg array_like left, right: Vector.
    :arg function pairwise: A pairwise distance function.

    :return: The multiset distance between `left` and `right`.
    :rtype: float

    Note that `function` must be vectorized, i.e., it is called directly on
    NumPy arrays, instead of on their pairwise elements. If your function only
    works on individual elements, convert it to a NumPy ufunc first. For
    example::

        >>> f = np.vectorize(f, otypes=['float'])
    """
    left = np.asanyarray(left)
    right = np.asanyarray(right)

    nonzero = np.where(np.logical_or(left, right))
    distances = pairwise(left[nonzero], right[nonzero])
    return distances.sum() / (len(distances) + 1)


def euclidean(left, right):
    """
    Calculate the Euclidean distance between two vectors.

    :arg array_like left, right: Vector.

    :return: The Euclidean distance between `left` and `right`.
    :rtype: float
    """
    return vector_length(np.subtract(left, right))


def cosine_similarity(left, right):
    """
    Calculate the Cosine similarity between two vectors.

    :arg array_like left, right: Vector.

    :return: The Cosine similarity between `left` and `right`.
    :rtype: float
    """
    return np.dot(left, right) / (vector_length(left) * vector_length(right))


def covariance(left, right):
    """
    Calculate the Cosine similarity between two vectors.

    :arg array_like left, right: Vector.

    :return: The Cosine similarity between `left` and `right`.
    :rtype: float
    """
    assert len(left) == len(right)
    dist = np.dot(left, right) * 100 / float(len(left))
    # print "size of kmer vector %d, distance %.3f" % (len(left), dist)
    return dist


def distance(left, right, pairwise=pairwise['prod'], distance_function=None):
    """
    Calculate the distance between two *k*-mer profiles.

    :arg left, right: Profiles to calculate distance
      between.

    :return: The distance between `left` and `right`.
    :rtype: float
    """
    if not distance_function:
        return multiset(left, right, pairwise)
    return distance_function(left, right)


############################# Most of the above functions are from kMer li

############################## parse files ##################################

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


def isheader(line):
    return line[0] == '>'

'''
input: .fna file, fasta format, genome sequence
output: k-mer spectrum
'''


def parse_genome_file(genome_file, klen):
    with open(genome_file, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                # merge multiple lines of sequence
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                glen = len(sequence)

                # To get normalized k-mer frequency
                # atcg_fraction= get_atcg_fraction(sequence)
                # base_frac = []
                # for nt in sequence:
                #     base_frac.append(atcg_fraction[nt])
                # # print 'base_frac %s' % base_frac
                # norm = np.prod(np.asarray(base_frac))
                # print 'norm for genome %.5f' % norm
                # if norm == 0:
                #     norm = 1

                kmer_count = kmer_from_sequences(sequence, klen)
                kmer_freq = kmer_count * 100 / float(glen - klen + 1)
                kmer_freq = kmer_freq / norm
                return kmer_freq


def get_genome_profiles(genome_file, klen):
    '''
    Compute k-mer frequency for a genome
    '''
    profiles = {}

    with open(genome_file, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                glen = len(sequence)
                # atcg_fraction = get_atcg_fraction(sequence)
                # base_frac = []
                # for nt in sequence:
                #     base_frac.append(atcg_fraction[nt])
                # # print 'base_frac %s' % base_frac
                # norm = np.prod(np.asarray(base_frac))
                # print 'norm for genome %.5f' % norm
                # if norm == 0:
                #     norm = 1
                for k in range(1, klen + 1):
                    p = kmer_from_sequences(sequence, k)
                    kmer_count = kmer_from_sequences(sequence, k)
                    kmer_freq = kmer_count * 100 / float(glen - k + 1)
                    # kmer_freq = kmer_freq / norm
                    profiles[k] = kmer_freq

                return sequence, glen, profiles



def get_atcg_fraction(dna):
    dict = {}
    total = float(len(dna))
    dict['A'] = dna.count('A') / total
    dict['T'] = dna.count('T') / total
    dict['C'] = dna.count('C') / total
    dict['G'] = dna.count('G') / total

    return dict


def parse_genes(gfile, genome_profiles, klen, morder=1, atcg_fraction={}, distance_function=None):
    '''
    input: .ffn file, fasta format, DNA sequence for each gene
    output: distance for each gene in a line
    '''
    i = 0
    dist_dict = {}
    for k in range(2, klen + 1):
        dist_dict[k] = []

    with open(gfile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                glen = len(sequence)
                dists = []
                gene_kmer_profiles = {}
                # base_frac = []
                # for nt in sequence:
                #     base_frac.append(atcg_fraction[nt])
                # norm = np.prod(np.asarray(base_frac))
                # print 'norm for gene %.5f' % norm
                # if norm == 0:
                #     norm = 1
                for k in range(2, klen + 1):
                    kmer_count = kmer_from_sequences(sequence, k)
                    kmer_freq = kmer_count * 100 / float(glen - k + 1)
                    # kmer_freq = kmer_freq / norm
                    gene_kmer_profiles[k] = kmer_freq

                for k in range(2, klen + 1):
                    gene_kmer = gene_kmer_profiles[k]
                    genome_kmer = genome_profiles[k]
                    # print 'gene_kmer %s' % gene_kmer
                    # print 'genome_kmer %s' % genome_kmer
                    dist = distance(gene_kmer, genome_kmer,
                                    distance_function=distance_function)
                    dist_dict.setdefault(k, []).append(dist)

    return dist_dict


def parse_segs(gfile, gnome, genome_profiles, klen, outfile, morder=1, atcg_fraction={}, distance_function=None):
    '''
    input: position for each segment (1-based)
    output: kmer measures for each segment in a line
    '''
    i = 0
    dist_dict = {}
    for k in range(2, klen + 1):
        dist_dict[k] = []
    with open(gfile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            sequence = gnome[start - 1:end]
            glen = len(sequence)
            dists = []
            seg_kmer_profiles = {}

            for k in range(2, klen + 1):
                kmer_count = kmer_from_sequences(sequence, k)
                kmer_freq = kmer_count * 100 / float(glen - k + 1)
                seg_kmer_profiles[k] = kmer_freq

            for k in range(2, klen + 1):
                seg_kmer = seg_kmer_profiles[k]
                genome_kmer = genome_profiles[k]
                dist = distance(seg_kmer, genome_kmer,
                                distance_function=distance_function)
                dist_dict.setdefault(k, []).append(dist)

    return dist_dict


def write_dists(dist_dict, klen, outfile):
    new_dist_dict = {}

    for key, values in dist_dict.items():
        maxv = max(values)
        minv = min(values)
        new_values = []
        for v in values:
            nv = (v - minv) / (maxv - minv)
            new_values.append(nv)
        new_dist_dict[key] = new_values

    with open(outfile, 'w') as fout:
        for i in range(len(new_dist_dict[2])):
            dists = []
            for k in range(2, klen + 1):
                dists.append(new_dist_dict[k][i])
            line = '%d' % i
            fout.write(line)
            for d in dists:
                line = "\t%0.3f" % d
                fout.write(line)
            fout.write('\n')


######################### compute k-mer frequency for contigs #################
def get_contigs_profile(genome_file, klen):
    profiles = {}

    for k in range(1, klen + 1):
        kmer_dict = {}
        number = 4 ** k
        bitmask = number - 0x01
        counts = [0] * number
        kmer_count = np.array(counts, dtype='int64')
        total_kmer = 0
        i = 0
        with open(genome_file, 'rb') as fin:
            for header, group in itertools.groupby(fin, isheader):
                if header:
                    i += 1
                else:
                    sequence = ''.join(line.strip() for line in group)
                    sequence = sequence.upper()
                    sequence = standardize_DNASeq(sequence)
                    glen = len(sequence)
                    p = kmer_from_sequences(sequence, k)
                    kmer_count += kmer_from_sequences(sequence, k)
                    total_kmer += glen - k + 1

            # print kmer_count
            # print total_kmer
            kmer_freq = kmer_count * 100 / float(total_kmer)
            profiles[k] = kmer_freq

    return profiles


def get_contigs(infile):
    '''
    infile -- Input file containing the sequence of contigs
    '''
    # id_mapping = {}
    contig_sequence = {}
    i = 0
    with open(infile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
            #     name = header.strip().split()
            #     id_mapping[i] = name
            else:
                sequence = ''.join(line.strip() for line in group)
                sequence = sequence.upper()
                sequence = standardize_DNASeq(sequence)
                contig_sequence[i] = sequence

    return contig_sequence


def parse_segs_contigs(pfile, contig_sequence, genome_profiles, klen, outfile, morder=1, atcg_fraction={}, distance_function=None):
    '''
    input: pfile -- position for each segment (1-based);  contig_sequence -- contig sequence (standardized)
    output: GC measures for each segment in a line
    '''
    dist_dict = {}
    for k in range(2, klen + 1):
        dist_dict[k] = []

    i = 0
    with open(pfile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            p1 = fields[0]
            mark = p1.index('_')
            start = int(p1[mark + 1:])
            p2 = fields[1]
            end = int(p2[mark + 1:])
            contig_id = int(p1[0:mark])
            contig = contig_sequence[contig_id]
            sequence = contig[start - 1:end]
            glen = len(sequence)
            dists = []
            seg_kmer_profiles = {}

            for k in range(2, klen + 1):
                kmer_count = kmer_from_sequences(sequence, k)
                kmer_freq = kmer_count * 100 / float(glen - k + 1)
                seg_kmer_profiles[k] = kmer_freq

            for k in range(2, klen + 1):
                seg_kmer = seg_kmer_profiles[k]
                genome_kmer = genome_profiles[k]
                dist = distance(seg_kmer, genome_kmer,
                                distance_function=distance_function)
                dist_dict.setdefault(k, []).append(dist)

    return dist_dict


####################################################################
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-g", "--gfile", dest="gfile",
                      help="The input gene files in fasta format")
    parser.add_option("-i", "--genome_file", dest="genome_file",
                      help="The input genome files in fasta format")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file of metrics for each gene")
    parser.add_option("-d", "--distance", dest="distance", default="covariance",
                      help="The distance function")
    parser.add_option("-k", "--klen", dest="klen", type=int, default=8,
                      help="The maximum k-mer size")
    parser.add_option("-m", "--morder", dest="morder", type=int, default=1,
                      help="The order of Markov model for background sequence")
    parser.add_option("-r", "--is_seg", dest="is_seg", action='store_true', default=False,
                      help="Analyze GC content in genomic segments rather than genes")
    parser.add_option("-s", "--seg_file", dest="seg_file",
                      help="input file of segment position")
    parser.add_option("-c", "--is_contig", dest="is_contig", action='store_true', default=False,
                      help="Analyze contigs from unassembled genomes")
    options, args = parser.parse_args()



    outfile = options.output + '.' + options.distance
    if options.is_seg:  # analyze k-mer for each interval
        gnome, genome_profiles = get_genome_profiles(options.genome_file, options.klen)
        dist_dict = parse_segs(options.seg_file, gnome, genome_profiles, options.klen,
                               outfile, morder=options.morder, atcg_fraction={}, distance_function=eval(options.distance))
        write_dists(dist_dict, options.klen, outfile)
    elif options.is_contig:
        genome_profiles = get_contigs_profile(options.genome_file, options.klen)
        contig_sequence = get_contigs(options.genome_file)
        dist_dict = parse_segs_contigs(options.seg_file, contig_sequence, genome_profiles, options.klen,
                               outfile, morder=options.morder, atcg_fraction={}, distance_function=eval(options.distance))
        write_dists(dist_dict, options.klen, outfile)
    else:
        gnome, genome_profiles = get_genome_profiles(options.genome_file, options.klen)
        dist_dict = parse_genes(options.gfile, genome_profiles, options.klen,
                                morder=options.morder, atcg_fraction={}, distance_function=eval(options.distance))
        write_dists(dist_dict, options.klen, outfile)
