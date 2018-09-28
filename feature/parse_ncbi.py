#!/usr/bin/env python

# This script is used to parse the files obtained from NCBI ftp
#
# Input:
# GCF_XXX_XXX_protein.faa
#
# Output:
# The location and ID of each protein
#

import re
import optparse
import itertools

def get_protIDs(infile):
    protIDs = []
    with open(infile,'r') as fin:
        for line in fin:
            field = line.strip()
            protIDs.append(field)
    return protIDs


def isheader(line):
    return line[0] == '>'




# Retrieve DNA-related information for a set of protein IDs
def parse_cds_IDs(infile, protIDs, outfile):
    id_pattern='\[protein_id=([\S]+)\]'+''
    loc_pattern='\[location=([\S]+)\]'+''
    ptn_id = re.compile(id_pattern)
    ptn_loc = re.compile(loc_pattern)
    with open(infile,'r') as fin, open(outfile,'w') as fout:
        for line in fin:
            fields = line.strip().split()
            gid = fields[0]
            match_id = ptn_id.search(line.strip())
            match_loc = ptn_loc.search(line.strip())
            if match_id:   # Pseudo genes do not have protein_id
                id = match_id.group(1)
                loc = match_loc.group(1)
                if 'complement' in loc:
                    strand = 'R'
                    loc.replace('complement(','').replace(')','')
                else:
                    strand = 'F'
                pos = loc.split('\.\.')
                line = '%s\n' % (id)
                fout.write(line)


# Retrieve DNA information for a set of protein IDs
# Output files:
# .ffn file -- DNA sequence for the set of protein IDs
# gene_id -- IDs for the retrieved DNA sequence
# glist -- locations for the retrieved DNA sequence
def retrieve_cds(infile, protIDs, dna_file, id_file, loc_file):
    id_pattern='\[protein_id=([\S]+)\]'+''
    loc_pattern='\[location=([\S]+)\]'+''
    coord_pattern='([0-9]+)\.\.([0-9]+)'
    ptn_id = re.compile(id_pattern)
    ptn_loc = re.compile(loc_pattern)
    ptn_coord = re.compile(coord_pattern)

    fout_dna = open(dna_file,'w')
    fout_id = open(id_file,'w')
    fout_loc = open(loc_file,'w')

    i = 0
    write_dna = False
    with open(infile,'r') as fin:
        for header, group in itertools.groupby(fin, isheader):
            if header:
                # print header
                # Check prot ID to see whether it is in protIDs
                line = group.next()
                # print line
                match_id = ptn_id.search(line.strip())
                if match_id:   # Pseudo genes do not have protein_id
                    id = match_id.group(1)
                    # Find other information
                    if id in protIDs:
                        i += 1
                        write_dna = True
                        match_loc = ptn_loc.search(line.strip())
                        loc = match_loc.group(1)
                        if 'complement' in loc:
                            strand = 'R'
                            loc.replace('complement(','').replace(')','')
                        else:
                            strand = 'F'
                        # print loc
                        match_pos = ptn_coord.search(loc)
                        # print match_pos.group()
                        start = int(match_pos.group(1))
                        end = int(match_pos.group(2))
                        loc_line = '%d\t%d\t%d\t%s\n' % (i, start, end, strand)
                        fout_loc.write(loc_line)

                        fields = line.strip().split()
                        # gid = fields[0]
                        gid = '>' + str(id)
                        fout_id.write(gid+'\n')

                        fout_dna.write(line)

                        # Remove this id from the list -- only consider one protein if there are two proteins with the same ID
                        protIDs.remove(id)

            else:
                if write_dna:
                    for line in group:
                        fout_dna.write(line)
                write_dna=False

    fout_loc.close()
    fout_id.close()
    fout_dna.close()




if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-i", "--infile", dest="infile", help="The input file")
    parser.add_option("-p", "--protID_file", dest="protID_file", help="The input file containing a list of protein IDs")
    parser.add_option("-o", "--outfile", dest="outfile", help="The output file")
    parser.add_option("-d", "--id_file", dest="id_file", help="The output file")
    parser.add_option("-l", "--loc_file", dest="loc_file", help="The output file")
    parser.add_option("-g", "--dna_file", dest="dna_file", help="The output file")
    options, args = parser.parse_args()

    # parse_cds_IDs(options.infile, options.outfile)
    protIDs = get_protIDs(options.protID_file)
    retrieve_cds(options.infile, protIDs, options.dna_file, options.id_file, options.loc_file)
