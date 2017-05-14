import os
import optparse

'''
find the location of all genes/ORFs in a genome
'''
def read_loc_ncbi(locfile):
    loc_dict = {}
    i = 1
    with open(locfile, 'r') as fin:
        for line in fin:
            # multiple positions in a line
            # e.g. >gi|225869487|ref|NC_012471.1|:1177219-1177902,1177906-1178478 Streptococcus equi subsp. equi 4047, complete genome
            if ',' in line:
                if 'c' in line:
                   strand = 'R'
                   content = line.strip().replace('c', '')
                else:
                   strand = 'F'
                   content = line.strip()
                genes = content.split(',')
                poses = []
                for g in genes:
                    p = g.split('-')
                    intp = map(int, p)
                    poses.extend(intp)
                start = min(poses)
                end = max(poses)
            else:
                pos = line.strip().split('-')
                if 'c' in line:
                    start = pos[1]
                    end = pos[0][1:]
                    strand = 'R'
                else:
                    start = pos[0]
                    end = pos[1]
                    strand = 'F'
            loc_dict[i] = [start, end, strand]
            i += 1

    return loc_dict


'''
find the location of all genes/ORFs in a genome
'''
def read_loc_pred(locfile):
    loc_dict = {}
    i = 1
    with open(locfile, 'r') as fin:
        for line in fin:
            fields = line.strip().split('_')
            pos = fields[-1].split('-')
            p1 = int(pos[0])
            p2 = int(pos[1])
            if (p1 > p2):
                start = p2
                end = p1
                strand = 'R'
            else:
                start = p1
                end = p2
                strand = 'F'
            loc_dict[i] = [start, end, strand]
            i += 1

    return loc_dict



def read_func(funcfile, loc_dict):
    i = 1
    with open(funcfile, 'r') as fin:
        for line in fin:
            func = line.strip()
            loc_dict.setdefault(i, []).append(func)
            i += 1

    return loc_dict


def write_gene(outfile, loc_dict):
    fout = open(outfile, 'w')
    for key, value in loc_dict.items():
        items = [str(key)]
        for v in value:
            items.append(str(v))
        # print items
        line = '\t'.join(items) + '\n'
        fout.write(line)
    fout.close()


if __name__ == '__main__':
        parser = optparse.OptionParser()

        parser.add_option("-l", "--locfile", dest="locfile", help="input file of gene/ORF locus")
        parser.add_option("-p", "--funcfile", dest="funcfile", help="input file of gene products/function")
        parser.add_option("-o", "--outfile", dest="outfile", help="output file")
        parser.add_option("-n", "--ncbi", dest="ncbi", action="store_true", default=False, help="the gene file is from ncbi")
        options, args = parser.parse_args()

        if options.ncbi:
            loc_dict = read_loc_ncbi(options.locfile)
        else:
            loc_dict = read_loc_pred(options.locfile)
        if options.funcfile:
            loc_dict = read_func(options.funcfile, loc_dict)
        write_gene(options.outfile, loc_dict)