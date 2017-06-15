# Find the infernal hits for a genome
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#
# Input:
# infernal output, format --tblout
#
# Output:
# id start end name rfam_accession score E-value
# 



import optparse


def parse_cmscan(infile):
    hits = []
    with open(infile, 'rb') as fin:
        for line in fin:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                name = fields[1]
                accession = fields[2]
                start = fields[9]
                end = fields[10]
                strand = fields[11]
                score = fields[16]
                evalue = fields[17]
                hit = (name, accession, start, end, strand, score, evalue)
                hits.append(hit)
    return hits


def write_hit(outfile, hits):
    fout = open(outfile, 'w')
    # print gene_list
    for hit in hits:
        line = '\t'.join(hit)
        fout.write(line)
    fout.close()



if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-b", "--infernaloutput", dest="infernaloutput", help="input file of infernal hits")
    parser.add_option("-g", "--genefile", dest="genefile", help="input file of gene lists")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file")
    options, args = parser.parse_args()

    hit_dict = parse_infernal(options.infernaloutput)
    write_hit(options.outfile, gene_list, hit_dict)
