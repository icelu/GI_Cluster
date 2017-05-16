import os
import optparse
import re

# input file: Pfam-A.hmm, downloaded from ftp site of Pfam
# ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
#
# Extract Pfam HMMs representing mobility genes or transposable elements.
#
# PATTERN:
# less Pfam-A.hmm.dat | grep -E '[Tt]ransposase|[Ii]ntegrase|[Tt]ransposon|[Rr]esolvase|[Rr]elaxase|[Ii]nsertion element|[Rr]ecombinase|Bacterial mobilisation|[Tt]ransposable|Mu\b|[Tt]ransposition'
#
# Sample command:
# python extract_mobgene.py -i Pfam-A.hmm.dat -o mobgene.list
# hmmfetch -f Pfam-A.hmm mobgene.list > Pfam_mobgene.hmm


def extract_mob_domain(infile, outfile):
    with open(infile, 'rb') as fin, open(outfile, 'w') as fout:
        for line in fin:
            if line.startswith('#=GF AC'):
                fields = line.strip().split()
                acc = fields[-1]

                for line in fin:
                    content = line.strip()
                    if content.startswith('#=GF DE'):
                        pattern = re.compile('[Tt]ransposase|[Ii]ntegrase|[Tt]ransposon|[Rr]esolvase|[Rr]elaxase|[Ii]nsertion element|[Rr]ecombinase|Bacterial mobilisation|[Tt]ransposable|Mu[-\s]|[Tt]ransposition')
                        res = pattern.search(content)
                        if res != None:
                            line = "%s\n" % (acc)
                            fout.write(line)
                        # break to go into next section
                        break



if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option("-i", "--infile", dest="infile", help="input file of pfam HMM profile")
        parser.add_option("-o", "--outfile", dest="outfile", help="output file for the HMM profile related to mobility genes")
        options, args = parser.parse_args()

        extract_mob_domain(options.infile, options.outfile)
