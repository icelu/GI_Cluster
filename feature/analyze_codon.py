# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#


import optparse
import itertools
import math
import numpy as np
from scipy.stats import chisquare
import random

# 64 codons
CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


# this dictionary shows which codons encode the same AA
SynonymousCodons = {
    'CYS': ['TGT', 'TGC'],
    'ASP': ['GAT', 'GAC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'GLN': ['CAA', 'CAG'],
    'MET': ['ATG'],
    'ASN': ['AAC', 'AAT'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'LYS': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    'PHE': ['TTT', 'TTC'],
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'HIS': ['CAT', 'CAC'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'TRP': ['TGG'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
    'GLU': ['GAG', 'GAA'],
    'TYR': ['TAT', 'TAC']
}

Code = {}
#Standard codon table
Code[1] = {'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'C': ['TGT', 'TGC'], 'E': ['GAG', 'GAA'], 'D': ['GAC', 'GAT'], 'G': ['GGT', 'GGG', 'GGA', 'GGC'], 'F': ['TTT', 'TTC'], 'I': ['ATC', 'ATA', 'ATT'], 'H': ['CAT', 'CAC'], 'K': ['AAG', 'AAA'], 'M': ['ATG'], 'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'N': ['AAC', 'AAT'], 'Q': ['CAA', 'CAG'], 'P': ['CCT', 'CCG', 'CCA', 'CCC'], 'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'], 'R': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'], 'T': ['ACA', 'ACG', 'ACT', 'ACC'], 'W': ['TGG'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'Y': ['TAT', 'TAC']}
#Mycoplasma/Spiroplasma codon table
Code[2] = {'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'C': ['TGT', 'TGC'], 'E': ['GAG', 'GAA'], 'D': ['GAC', 'GAT'], 'G': ['GGT', 'GGG', 'GGA', 'GGC'], 'F': ['TTT', 'TTC'], 'I': ['ATC', 'ATA', 'ATT'], 'H': ['CAT', 'CAC'], 'K': ['AAG', 'AAA'], 'M': ['ATG'], 'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'N': ['AAC', 'AAT'], 'Q': ['CAA', 'CAG'], 'P': ['CCT', 'CCG', 'CCA', 'CCC'], 'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'], 'R': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'], 'T': ['ACA', 'ACG', 'ACT', 'ACC'], 'W': ['TGA', 'TGG'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'Y': ['TAT', 'TAC']}


def isheader(line):
    return line[0] == '>'


def complement_seq(GeneSeq):
    GeneSeq_complement_seq = ''
    for i in GeneSeq[::-1]:
        if   i == 'A': GeneSeq_complement_seq += 'T'
        elif i == 'T': GeneSeq_complement_seq += 'A'
        elif i == 'C': GeneSeq_complement_seq += 'G'
        elif i == 'G': GeneSeq_complement_seq += 'C'
    return GeneSeq_complement_seq


def vector_length(vector):
    """
    Calculate the Euclidean length of a vector.

    :arg array_like vector: A vector.

    :return: The length of `vector`.
    :rtype: float
    """
    # Note: This seems faster than `numpy.linalg.norm`.
    return np.sqrt(np.dot(vector, vector))


# Get codon/amino acid bias between individual genes and the overall genome based on cosine function
# <CODE_TYPE> either be 's' for standard codon table or 'm' for mycoplasma/spiroplasma codon table and default is 's'
# 's' -- code_type = 1, 'm' -- code_type = 2
# given genome sequence and gene position
def get_freq_mat(genes, genome, code_type=1):
    '''Calculate codon/amino acid frequencies of all genes'''
    codon_freq_mat = {}
    amino_freq_mat = {}
    for (start, end, strand) in genes:
        if start < end:
           if strand == '+':
              curr_gene = genome[start - 1: end]
           else:curr_gene = complement_seq(genome[start - 1: end])
        elif start >= end:
           if strand== '+':
              curr_gene = genome[start - 1: len(genome)] + genome[0: end]
           else:curr_gene = complement_seq(genome[start - 1: len(genome)] + genome[0: end])

        remainder = len(curr_gene)%3
        # exclude partial codons
        if remainder != 0:
            curr_gene = curr_gene[0:-remainder]

        Codon_List = [curr_gene[i: i + 3] for i in range(0, len(curr_gene), 3)]
        codon_freq = []
        amino_freq = []
        for i in Code[code_type]:
          TempFre = .0
          for j in Code[code_type][i]:
              codon_freq.append(Codon_List.count(j) *1.0 / len(Codon_List))
              TempFre += Codon_List.count(j) *1.0 / len(Codon_List)
          amino_freq.append(TempFre)
        Loci = (end + start) / 2.0
        codon_freq_mat[(start, end)] = codon_freq
        amino_freq_mat[(start, end)] = amino_freq

    return (codon_freq_mat, amino_freq_mat)

# Switch non-standard character to 'AGCT'
def switch_to_AGCT(Nuc):
    SWITCH = {'R': lambda: random.choice('AG'), \
           'Y': lambda: random.choice('CT'), \
           'M': lambda: random.choice('AC'), \
           'K': lambda: random.choice('GT'), \
           'S': lambda: random.choice('GC'), \
           'W': lambda: random.choice('AT'), \
           'H': lambda: random.choice('ATC'), \
           'B': lambda: random.choice('GTC'), \
           'V': lambda: random.choice('GAC'), \
           'D': lambda: random.choice('GAT'), \
           'N': lambda: random.choice('AGCT')}
    return SWITCH[Nuc]()


'''
Replace ambiguous bases with ATCG
'''
def standardize_DNASeq(genome):
    for i in ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']:
        if i in genome:
           genome = genome.replace(i, switch_to_AGCT(i))
    return genome


def get_freq_mat_file(gfile, code_type=1):
    '''Calculate codon/amino acid frequencies of all genes'''
    codon_freq_mat = {}
    amino_freq_mat = {}
    gene_id = 0
    with open(gfile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                gene_id += 1
                #line = group.next()
            else:
                sequence = ''.join(line.strip() for line in group)
                curr_gene = sequence.upper()
                curr_gene = standardize_DNASeq(curr_gene)
                remainder = len(curr_gene)%3
                # exclude partial codons
                if remainder != 0:
                    curr_gene = curr_gene[0:-remainder]
                # if len(curr_gene) % 3 == 0:
                Codon_List = [curr_gene[i: i + 3] for i in range(0, len(curr_gene), 3)]
                codon_freq = []
                amino_freq = []
                for i in Code[code_type]:
                  TempFre = .0
                  for j in Code[code_type][i]:
                      codon_freq.append(Codon_List.count(j) *1.0 / len(Codon_List))
                      TempFre += Codon_List.count(j) *1.0 / len(Codon_List)
                  amino_freq.append(TempFre)

                codon_freq_mat[gene_id] = codon_freq
                amino_freq_mat[gene_id] = amino_freq

    # for k,v in  codon_freq_mat.items():
    #     print 'key:%s, val:%s \n' % (k,v)
    # for k,v in  amino_freq_mat.items():
    #     print 'key:%s, val:%s \n' % (k,v)
    return (codon_freq_mat, amino_freq_mat)


def get_amino_bias(amino_freq_mat):
    '''Calculate average codon/amino acid frequencies of genome'''
    amino_bias={}
    Ave_Amino_Mat=[.0 for i in range(20)]

    for i in amino_freq_mat:
        for j in range(20):
            Ave_Amino_Mat[j] += amino_freq_mat[i][j] / len(amino_freq_mat)

    for i in amino_freq_mat:
        Cosine1 = np.dot(amino_freq_mat[i], Ave_Amino_Mat)
        Cosine2 = vector_length(amino_freq_mat[i])
        Cosine3 = vector_length(Ave_Amino_Mat)
        amino_bias[i] = 1 - Cosine1 / (Cosine2 * Cosine3)
        # Cosine1 = .0
        # Cosine2 = .0
        # Cosine3 = .0
        # for j in range(20):
        #     Cosine1 += amino_freq_mat[i][j] * Ave_Amino_Mat[j]
        #     Cosine2 += amino_freq_mat[i][j] * amino_freq_mat[i][j]
        #     Cosine3 += Ave_Amino_Mat[j] * Ave_Amino_Mat[j]
        # amino_bias[i] = 1 - Cosine1 / math.sqrt(Cosine2 * Cosine3)

    return amino_bias

def get_codon_bias(codon_freq_mat, codon_num):
    '''Get codon/amino acid bias of each gene'''
    codon_bias = {}
    Ave_Codon_Mat=[.0 for i in range(codon_num)]

    for i in codon_freq_mat:
        for j in range(codon_num):
            Ave_Codon_Mat[j] += codon_freq_mat[i][j] / len(codon_freq_mat)

    for i in codon_freq_mat:
        Cosine1 = np.dot(codon_freq_mat[i], Ave_Codon_Mat)
        Cosine2 = vector_length(codon_freq_mat[i])
        Cosine3 = vector_length(Ave_Codon_Mat)
        codon_bias[i] = 1 - Cosine1 / (Cosine2 * Cosine3)
        # Cosine1 = .0
        # Cosine2 = .0
        # Cosine3 = .0
        # for j in range(61):
        #     Cosine1 += codon_freq_mat[i][j] * Ave_Codon_Mat[j]
        #     Cosine2 += codon_freq_mat[i][j] * codon_freq_mat[i][j]
        #     Cosine3 += Ave_Codon_Mat[j] * Ave_Codon_Mat[j]
        # codon_bias[i] = 1 - Cosine1 / math.sqrt(Cosine2 * Cosine3)

    return codon_bias


def write_output(codon_bias, amino_bias, rcsu_list, rcsu_avg, outfile):
    with open(outfile, 'w') as fout:
        for i in range(1,len(codon_bias)+1):
            cb = codon_bias[i]
            ab = amino_bias[i]
            chi = chi_square(rcsu_list[i-1], rcsu_avg)
            line = '%.3f\t%.3f\t%.3f\n' % (cb, ab, chi)
            fout.write(line)


def count_codons(dna_sequence):
    # make the codon dictionary local
    codon_count = CodonsDict.copy()

    remainder = len(dna_sequence)%3
    # exclude partial codons
    if remainder != 0:
        dna_sequence = dna_sequence[0:-remainder]

    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            raise TypeError("illegal codon %s in gene: %s" % (codon, dna_sequence))
    return codon_count


def get_rcsu(dna_sequence):
    # count codon occurrences in the file.
    codon_count = count_codons(dna_sequence)

    # now to calculate the index we first need to sum the number of times
    # synonymous codons were used all together.
    rcsu = []
    for aa in SynonymousCodons:
        total = 0.0
          # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)
        codons = SynonymousCodons[aa]

        for codon in codons:
            total += codon_count[codon]

        # calculate the RSCU value for each of the codons
        for codon in codons:
            if total == 0.0:
                rcsu.append(0.0)
            else:
                denominator = float(total) / len(codons)
                rcsu.append(codon_count[codon] / denominator)
    # print len(rcsu) # 64 codons

    return rcsu


def chi_square(rcsu_gene, rcsu_avg):
    res = chisquare(rcsu_gene, f_exp=rcsu_avg)
    return res[0]

def get_codon_chi(gfile):
    rcsu_list=[]
    with open(gfile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if not header:
                #gene_id += 1
                #line = group.next()
            # else:
                sequence = ''.join(line.strip() for line in group)
                curr_gene = sequence.upper()
                curr_gene = standardize_DNASeq(curr_gene)

                rcsu = get_rcsu(curr_gene)
                rcsu_list.append(rcsu)

    rcsu_avg= np.mean(np.array(rcsu_list), axis=0)
    # print rcsu_avg

    return (rcsu_list, rcsu_avg)


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-g", "--gfile", dest="gfile",
                      help="input gene files in fasta format")
    parser.add_option("-o", "--output", dest="output",
                      help="output file of labeled intervals")
    parser.add_option("-c", "--code_type", dest="code_type", default=1, type=int,
                      help="either be '1' for standard codon table or '2' for mycoplasma/spiroplasma codon table and default is '1'")

    options, args = parser.parse_args()

    (codon_freq_mat, amino_freq_mat) = get_freq_mat_file(options.gfile, code_type=options.code_type)
    codons = []
    for a, c in Code[options.code_type].items():
        codons.extend(c)
    # print len(codons)
    codon_bias = get_codon_bias(codon_freq_mat, len(codons))
    amino_bias = get_amino_bias(amino_freq_mat)
    (rcsu_list, rcsu_avg) = get_codon_chi(options. gfile)

    write_output(codon_bias, amino_bias, rcsu_list, rcsu_avg, options.output)
