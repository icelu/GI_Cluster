#!/usr/bin/python

import optparse

def get_gene_id(infile):
    gene_list = []
    with open(infile, 'rb') as fin:
        for line in fin:
            gene = line.strip().replace('>', '')
            # print gene
            gene_list.append(gene)
    return gene_list


def parse_hmmer(infile):
    hit_dict = {}
    hit_product = ''
    with open(infile, 'rb') as fin:
        for line in fin:
            if line.startswith('Description:'):
                # print line
                fields = line.strip().split(':')
                hit_product = fields[-1].strip()
                #print hit_product
            elif line.strip():                
                if line.strip().startswith('E-value'):
                    #print 'evalue line'
                    for line in fin:
                        if line.strip() and line.strip()[0].isdigit():
                            #print line.strip()
                            values = line.strip().split() 
                            #print values
                            gid = values[8] 
                            #evalue = values[3]
                            evalue = float(values[3])
                            #print 'gid', gid
                            #print 'evalue', evalue
                            if  evalue <= 0.0001:
                                if gid not in hit_dict.keys():
                                    hit_dict[gid] = hit_product
                                else:
                                    hit_dict[gid] = hit_dict[gid] + ', ' + hit_product
                        elif line.strip() and line.strip()[0] == '-':
                            # print line.strip()
                            continue
                        else:
                            hit_product = ''
                            break
    #print hit_dict
    return hit_dict                       


def write_hit(outfile, gene_list, hit_dict):
    fout = open(outfile, 'w')
    for gene in gene_list:
        if gene in hit_dict.keys():
            line = '%s\n' % hit_dict[gene]
        else:
            line = '*No hits*\n'
        fout.write(line)    
    fout.close()
    
    

if __name__ == '__main__':    
    parser = optparse.OptionParser()
                 
    parser.add_option("-d", "--dboutput", dest="dboutput", help="input file of database search hits") 
    parser.add_option("-g", "--genefile", dest="genefile", help="input file of gene lists") 
    parser.add_option("-o", "--outfile", dest="outfile", help="output file") 
    options, args = parser.parse_args() 
        
    gene_list = get_gene_id(options.genefile) 
    hit_dict = parse_hmmer(options.dboutput)
    write_hit(options.outfile, gene_list, hit_dict)    
        
        
