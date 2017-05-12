#!/usr/bin/python
# Lu Bingxin

import os
import optparse

'''
given an input file with a list of start and end positions
create a file: where the intervals not shown in the input has an additional column 0,
while the intervals  shown in the input has an additional column 1.
eg.
Input:
1       1       5000    5000    0       0       6       0       0.403   0.101   4.968   0.068   0.034   21.513  0.225   -0.075  0.376   0.697   0.701   0.712   0.707   0.709   0.705   0.696   0.167   0.333   0.000   0.000   0.000   1.200   68.400
for each column starting from 5 in the input
Output:
For each feature 4 columns
hs1 1 14708 0.0
'''


def createFeatureFile(input_file, outdir, index, id):
    gene_dict = {}
    product_dict = {}

    output_file = os.path.join(outdir, os.path.basename(input_file) + '_F' + str(id))

    with open(input_file, 'rb') as infile, open(output_file, 'w') as outfile:
        # next(infile)
        last_end = 0
        for line in infile:
            # should always strip to trim the trailing  special characters
            fields = line.strip().split('\t')
            start = int(fields[1])
            end = int(fields[2])

            feature = float(fields[index])
            line = "hs1\t%d\t%d\t%.3f\n" % (start, end, feature)
            outfile.write(line)


if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option("-d", "--index", dest="index", type="int", default="4", help="the index of the feature to extract in the input file")
        parser.add_option("-i", "--input", dest="input", help="input file containing GIs")
        parser.add_option("-o", "--outdir", dest="outdir", help="the output directory for files used in visualization")
        parser.add_option("-f", "--features", dest="features", default="9,15,20,25,26,27,28,29", help="the features used in visualization: a list of numbers separted by comma")
        parser.add_option("-a", "--all", dest="all", default=False, action="store_true", help="extract all the features in the input file[default=False]")

        options, args = parser.parse_args()
        if options.all:
            # GC, CAI, 4mer, MOB, PHAGE, VF, AR, NG
            # ft_range=[9,15,20,25,26,27,28,29]
            ft_range=options.features.strip().split(',')
            # for i in range(1, 9, 1):
            for i, ft in enumerate(ft_range):
                createFeatureFile(options.input, options.outdir, int(ft)-1, i+1)
        else:
            createFeatureFile(options.input, options.outdir, options.index, options.index)
