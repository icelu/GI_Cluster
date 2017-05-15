

import optparse


'''
find the blast hits for a predicted gene
Input:
blast output, format 6
gene locus for each gene
'''

def get_gene_id(infile):
    gene_list = []
    with open(infile, 'rb') as fin:
        for line in fin:
            gene = line.strip().replace('>', '')
            # print gene
            gene_list.append(gene)
    return gene_list


def parse_blast(infile):
    hit_dict = {}
    with open(infile, 'rb') as fin:
        for line in fin:
            fields = line.strip().split('\t')
            gene = fields[0]
            hit = fields[1]
            hit_dict[gene] = hit
    return hit_dict


def write_hit(outfile, gene_list, hit_dict):
    fout = open(outfile, 'w')
    #print gene_list
    for gene in gene_list:
        if gene in hit_dict.keys():
            line = '%s\n' % hit_dict[gene]
        else:
            line = '*No hits*\n'
        fout.write(line)    
    fout.close()
    
    

if __name__ == '__main__':    
    parser = optparse.OptionParser()
                 
    parser.add_option("-b", "--blastoutput", dest="blastoutput", help="input file of blast hits") 
    parser.add_option("-g", "--genefile", dest="genefile", help="input file of gene lists") 
    parser.add_option("-o", "--outfile", dest="outfile", help="output file") 
    options, args = parser.parse_args() 
        
    gene_list = get_gene_id(options.genefile) 
    hit_dict = parse_blast(options.blastoutput)
    write_hit(options.outfile, gene_list, hit_dict)
        
