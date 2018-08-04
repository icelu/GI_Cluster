#!/usr/bin/env python

# Fill template file for comparing different predictions of genomic islands in Circos
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

import os
import optparse


def createConfig(input_file, output_file, programs, oname, name, colors):
    num_prog = len(programs)
    prog_colors=colors[num_prog]
    radius = 0.9-(num_prog-1)*0.1
    i=0
    with open(input_file, 'rb') as fin, open(output_file, 'w') as fout:
        line = fin.readline()
        while line:
            if '# start prog' in line:    # Start of the settings for one program
                # Find the number of program
                begin = line.strip()
                print begin
                nums = [int(s) for s in begin.split() if s.isdigit()]
                print nums
                id = nums[0]
                if id + 1 > num_prog:
                    # delete the following lines until <highlights>
                    line = fin.readline()
                    while line:
                        if '# end prog' in line:
                            line = fin.readline()
                            break
                        line = fin.readline()
                else:
                    prog_id = '$prog'+str(id)
                    color_id = '$color'+str(id)
                    # Read in the modules for one program
                    line = fin.readline()
                    while line:
                        if  prog_id in line:
                            line=line.replace(prog_id, programs[id])
                        if color_id in line:
                            line=line.replace(color_id, prog_colors[id])
                        if '# end prog' in line:
                            line = fin.readline()
                            break
                        fout.write(line)
                        line = fin.readline()
            else:
                if '$outfile' in line:
                    line=line.replace('$outfile', oname)
                if '$radius' in line:
                    line=line.replace('$radius', str(radius)+"r")
                if '$organism' in line:
                    line=line.replace('$organism', name)
                fout.write(line)
                line = fin.readline()



if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-n", "--name", dest="name", help="The name of the organism.")
    parser.add_option("-p", "--programs", dest="programs", default="gicluster", help="The programs for comparison.")
    parser.add_option("-i", "--infile", dest="infile", help="The input template of the configuration file for circos.")
    parser.add_option("-o", "--outfile", dest="outfile", help="The output file for visualization.")
    parser.add_option("-f", "--oname", dest="oname", default="gifeature", help="The name of the output figure.")

    options, args = parser.parse_args()

    colors={}
    # The following color codes are from http://colorbrewer2.org
    colors[2]=['27,158,119','217,95,2']
    colors[3]=['27,158,119','217,95,2','117,112,179']
    colors[4]=['27,158,119','217,95,2','117,112,179', '231,41,138']
    colors[5]=['27,158,119','217,95,2','117,112,179', '231,41,138', '102,166,30']
    colors[6]=['27,158,119','217,95,2','117,112,179', '231,41,138', '102,166,30', '230,171,2']
    colors[7]=['27,158,119','217,95,2','117,112,179', '231,41,138', '102,166,30', '230,171,2', '166,118,29']
    colors[8]=['27,158,119','217,95,2','117,112,179', '231,41,138', '102,166,30', '230,171,2', '166,118,29', '102,102,102']
    colors[9]=['166,206,227','31,120,180','178,223,138', '51,160,44', '251,154,153', '227,26,28', '253,191,111', '255,127,0', '202,178,214']

    programs=options.programs.strip().split(",")
    print programs
    createConfig(options.infile, options.outfile, programs, options.oname, options.name, colors)
