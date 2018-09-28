#!/usr/bin/env python

# Fill template file for visualizing GI-related features in Circos
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

import os
import optparse


def createConfig(input_file, name, output_file, prefix, oname):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if '$organism' in line:
                line=line.replace('$organism', name)
            if '$prefix' in line:
                line=line.replace('$prefix', prefix)
            if '$outfile' in line:
                line=line.replace('$outfile', oname)
            outfile.write(line)



if __name__ == '__main__':
        parser = optparse.OptionParser()
        parser.add_option("-n", "--name", dest="name", help="The name of the organism.")
        parser.add_option("-p", "--prefix", dest="prefix", default="feature.multi.percentage", help="The prefix of the feature file.")
        parser.add_option("-i", "--infile", dest="infile", help="The input template of the configuration file for circos.")
        parser.add_option("-o", "--outfile", dest="outfile", help="The output file for visualization.")
        parser.add_option("-f", "--oname", dest="oname", default="gifeature", help="The name of the output figure.")

        options, args = parser.parse_args()

        createConfig(options.infile, options.name, options.outfile, options.prefix, options.oname)
