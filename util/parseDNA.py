# Functions to read and standardize fasta files

import itertools

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


def get_contig_IDs(infile):
    '''
    infile -- Input file containing the sequence of contigs or multiple sequences
    '''
    id_mapping = {}
    i = 0
    with open(infile, 'rb') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                i += 1
                name = ''.join(line.strip() for line in group)
                name = name.split()[0][1:]  # Ignore the begining '>'
                id_mapping[name] = i

    return id_mapping



def get_contig_gene_IDs(infile):
    '''
    infile -- Input file containing the sequence of contigs or multiple sequences
    '''
    id_mapping = {}
    i = 0
    with open(infile, 'rb') as fin:
        for line in fin:
            i += 1
            name = line.strip()[1:] # Remove leading '>'
            id_mapping[i] = name

    return id_mapping
