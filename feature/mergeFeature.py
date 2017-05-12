import os
import optparse


'''
Find the features for each genomic segment
Input:
A list of clusters (segments)
A list of features for each ORF
Output:
Group the ORFs by clusters
    format:
    region position
    genes inside this region along with their features

Procedure:
store the features for each ORF in a dicionary

'''

'''
Note: the input genes are sorted and kept in this order
Use a list of tuple to store a gene along with its features: (gid, gfeature)
'''
def get_gene_feature(infile):
    features = []
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            id = fields[0]
            feature = fields[1:]
            # features[id] = fields[1:]
            item = (id, feature)
            features.append(item)
    # sort the genes by position
    return features

def parse_cmscan(infile):
    hits = []
    with open(infile, 'rb') as fin:
        for line in fin:
            if not line.startswith('#'):
                fields = line.strip().split()
                name = fields[1]
                accession = fields[2]
                start = int(fields[9])
                end = int(fields[10])
                strand = fields[11]
                if strand == '-':
                    tmp = start
                    start = end
                    end = tmp
                    strand = 'R'
                else:
                    strand = 'F'
                score = fields[16]
                evalue = fields[17]
                hit = (start, end, strand, name, accession, score, evalue)
                hits.append(hit)
    hits = sorted(hits)
    return hits


'''
Output
'''
def group_genes(infile, features, rnas, outfile):
    fout = open(outfile, 'w')
    # header
    #header = 'GeneID\tStart\tEnd\tStrand\tGC content\tCodon usage\tk-mer frequency\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\n'
    #header = 'GeneID\tStart\tEnd\tStrand\tGC\tGC1\tGC2\tGC3\tGC_skew\tCAI\tCBI\tFop\tNc\t2-mer\t3-mer\t4-mer\t5-mer\t6-mer\t7-mer\t8-mer\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\n'
    # header = 'GeneID\tStart\tEnd\tStrand\tGC\tGC1\tGC2\tGC3\tGC_skew\tCAI\tCBI\tFop\t2-mer\t3-mer\t4-mer\t5-mer\t6-mer\t7-mer\t8-mer\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\n'
    header = 'GeneID\tStart\tEnd\tStrand\tGC\tGC1\tGC2\tGC3\tGC_skew\tCUB\tAAB\tCHI\tCAI\tCBI\tFop\t2-mer\t3-mer\t4-mer\t5-mer\t6-mer\t7-mer\t8-mer\tMobility gene\tPhage\tVirulence factor\tAntibiotic resistance gene\tNovel gene\n'
    fout.write(header)
    # for each region, find the genes inside it
    i = 0
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            i+=1
            oline = '>Region' + str(i)+'\t'+line
            fout.write(oline)
            # print 'start %d, end %d' % (start, end)
            ri = 0
            for hit in rnas:
                (rstart, rend, strand, name, accession, score, evalue) = hit
                # rstart = int(rstart)
                # rend = int(rend)
                # if ((rstart >= start and rstart <= end) or ( rend >= start and rend <= end)):
                if not ((rstart > end) or (rend < start)):
                #rmid = (rend + rstart) / 2
                #if ((rstart >= start and rmid <= end) or (rmid >= start and rend <= end)):
                    ri += 1
                    id = 'r' + str(ri)
                    oline = id + '\t'
                    oline = oline + '\t'.join((str(rstart), str(rend), strand, name, accession, score, evalue)) + '\n'
                    fout.write(oline)
                    # a gene overlapping with two regions will be counted twice, since it is nontrivial to determine which region to include this gene
                    # rnas.remove(hit)

            for gid, gfeature in features:
                gstart = int(gfeature[0])
                gend = int(gfeature[1])
                gmid = (gend + gstart) / 2
                # print 'gstart %d, gend %d' % (gstart, gend)
                # genes strictly inside a region
                # if (gstart >= start and gend <= end):
                # include genes across the boundary (even with only 1 bp overlap), so one gene may be considerred in multiple segements
                # if ((gstart >= start and gstart <= end) or ( gend >= start and gend <= end)):
                if not ((gstart > end) or (gend < start)):
                # since segments are non-overlapping, if a gene is across >1 segments, if half of the gene is inside a segment, this gene is seen as inside this segment
                # if ((gstart >= start and gmid <= end) or (gmid >= start and gend <= end)):
                    # write the output
                    items = [gid]
                    for v in gfeature:
                        items.append(v)
                    # print items
                    oline = '\t'.join(items) + '\n'
                    fout.write(oline)

                    # remove this item from the list since it will not be included in other regions
                    # a gene overlapping with two regions will be counted twice, since it is nontrivial to determine which region to include this gene
                    # features.remove((gid, gfeature))

    fout.close()

if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-g", "--genefile", dest="genefile", help="input file of genes and their features")
    parser.add_option("-s", "--segfile", dest="segfile", help="input file of genomic segments")
    parser.add_option("-r", "--rnafile", dest="rnafile", help="input file of predicted RNAs")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file")
    options, args = parser.parse_args()

    features = get_gene_feature(options.genefile)
    rnas = parse_cmscan(options.rnafile)
    group_genes(options.segfile, features, rnas, options.outfile)
