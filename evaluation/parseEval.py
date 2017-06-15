# Extract metric values obtained by IntervalOverlap.py
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#

 
from os import listdir
from os.path import isfile, join
from optparse import OptionParser



def getMetricsFromEvalFile(folder, linenum, outfile, tag='eval_std_'):
    '''
    Parse eval_std_ to get metrics
    Only read the last line --> change to any line specified by linenum.
    '''
    # Find all files with name as eval_std_*
    file_list = [ join(folder, f) for f in listdir(folder) if isfile(join(folder, f)) and f.startswith(tag) ]
    base_metric_list = []
    for file in file_list:
        # start = file.index(tag[-1]) + 1
        # suffix = file[start:]
        base_res = getMetrics(file, linenum)
        if base_res:
            # Use abosolute path to faciliate merging multiple files
            base_metric_list.append((file, base_res))

    final_metric = join(folder, outfile)
    writeListOfTupleToFile(final_metric, base_metric_list)



def getMetrics(infile, linenum):
    last = ''
    with open(infile, 'r') as fh:
        lines = fh.readlines()
        # in case some output files are not complete
        if len(lines) > 0:
            last = lines[linenum]
    if last is not '':
        metric = []
        fields = last.strip().split('\t')
        for f in fields:
            value = f.split(':')[-1].strip()
            metric.append(value)

        return tuple(metric)


# tuple format: (v, (v,v,...))
def writeListOfTupleToFile(filename, list):
    outfile = open(filename, 'w')

    for key, values in list:
        suffix_str = '\t%s'
        line = [key]
        line_str = '%s'
        for value in values:
            line.append(value)
            line_str += suffix_str
        line_str += '\n'
        outfile.write(line_str % tuple(line))

    outfile.close()


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", "", dest="directory", help="home directory")
    parser.add_option("-t", "--tag", dest="tag", help="tag in the name output file. Default: 'eval_std_'")
    parser.add_option("-l", "--linenum", dest="linenum", type='int', default='-1', help="line number of the output to extract")
    parser.add_option("-o", "", dest="outfile", default='metric_list', help="The final list of metrics")
    (options, args) = parser.parse_args()

    if options.tag:
        getMetricsFromEvalFile(options.directory, options.linenum, options.outfile, tag=options.tag)
    else:
        getMetricsFromEvalFile(options.directory, options.linenum, options.outfile)
