import os
import optparse


def getGenomeSize(genomefile):
    firstLine = open(genomefile).readline()
    assert ('>' in firstLine), "This is not a standard FASTA file!"
    genomesize = 0
    with open(genomefile, 'rb') as fin:
        # skip the first line
        fin.next()
        for line in fin:
            genomesize += len(line.strip())

    return genomesize




def createGIFile(input_file, output_file, genomesize):
    '''
    given an input file with a list of start and end positions
    create a file: where the intervals not shown in the input has an additional column 0,
    while the intervals  shown in the input has an additional column 1.
    eg.
    Input:
    270001 275000
    560001 565000
    565001 570000

    Output:
    hs1 0 270000 0
    hs1 270001 275000 1
    hs1 275001 560000 0
    hs1 560001 565000 1
    '''
    gene_dict = {}
    product_dict = {}

    with open(input_file, 'rb') as infile, open(output_file, 'w') as outfile:
        # next(infile)
        last_end = 0
        for line in infile:
            # should always strip to trim the trailing special characters
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            # if start begins at 1, no line,
            # TODO: if end smaller than the total length, one 0 line to make genome connected
            if start != '1':
                nomark_line = "hs1\t%s\t%s\t0\n" % (last_end + 1, start - 1)
                outfile.write(nomark_line)
            marked_line = "hs1\t%s\t%s\t1\n" % (start, end)
            outfile.write(marked_line)
             # remember the last position of previous line
            last_end = end
        # if last_end < genome_len, should output additional line to complete the circle
        # print 'last_end', last_end
        if last_end < genomesize:
            nomark_line = "hs1\t%s\t%s\t0\n" % (last_end + 1, genomesize)
            outfile.write(nomark_line)


def createHightFile(input_file, output_file):
    gene_dict = {}
    product_dict = {}

    with open(input_file, 'rb') as infile, open(output_file, 'w') as outfile:
        # next(infile)
        last_end = 0
        for line in infile:
            # should always strip to trim the trailing special characters
            fields = line.strip().split('\t')
            start = int(fields[0])
            end = int(fields[1])
            marked_line = "hs1\t%s\t%s\n" % (start, end)
            outfile.write(marked_line)




def createGenomeFile(outfile, gname, genomesize):
    '''
    Output format:
    type parent name start end color
    chr - hs1 CT18 0 48090376 black
    '''
    with open(outfile, 'w') as fout:
        line = 'chr - hs1 %s 0 %d black\n' % (gname, genomesize)
        fout.write(line)

if __name__ == '__main__':
        parser = optparse.OptionParser()
        # parser.add_option("-l", "--length", dest="length", type="int", default="-1",  help="the size of the microbial genome")
        parser.add_option("-g", "--gfile", dest="gfile", help="input file containing the genome sequence")
        parser.add_option("-i", "--input", dest="input", help="input file containing GIs")
        parser.add_option("-o", "--output", dest="output", help="output file for visualizing genomic islands")
        parser.add_option("-c", "--cfile", dest="cfile", help="output file for visualing the whole genome")
        parser.add_option("-f", "--hfile", dest="hfile", help="output file for highlighting the genome islands")
        options, args = parser.parse_args()

        genomesize = getGenomeSize(options.gfile)
        createGIFile(options.input, options.output, genomesize)
        if options.cfile:
            gname = os.path.basename(options.gfile)
            createGenomeFile(options.cfile, gname, genomesize)
        if options.hfile:
            createHightFile(options.input, options.hfile)
