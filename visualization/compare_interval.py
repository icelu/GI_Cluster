# This script is used to visualize several sets of intervals in a range
#
# Author: Bingxin Lu
# Affiliation : National University of Singapore
# E-mail : bingxin@comp.nus.edu.sg
#


from optparse import OptionParser
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams['font.size'] = 20
plt.rcParams["figure.figsize"] = (33,3)

def get_data(indir):
    data = []
    # linewidths = []
    for file in os.listdir(indir):
        # print file
        fname = os.path.join(indir, file)
        print fname
        positions = []
        width = []
        with open(fname, 'rb') as fin:
            for line in fin:
                location = line.split()
                start = int(location[0])
                end = int(location[1])
                positions.append((start, end))
                # width.append(int(location[1]) - int(location[0]) + 1)
        # linewidths.append(width)
        data.append(positions)

    return data


def plot_intervals(data, fig_file):
    # set different colors for each set of positions
    col = {'ggred':'#F8766d', 'ggblue':'#619CFF', 'ggpurple':'#C77CFF'}
    colors = [col['ggred'], col['ggblue'], col['ggpurple']]
    # colors = np.array([[1, 0, 0],
    #                 [0, 1, 0],
    #                 [0, 0, 1]])

    # set different line properties for each set of positions
    # note that some overlap
    lineoffsets = np.array([0.3, 0.2, 0.1])
    linewidths = [25, 25, 25]
    y=[0.04, 0.1, 0.2, 0.3, 0.36]
    ylabels = ['','Reference','IslandViewer','GI-Cluster','']
    plt.ylim(0.04, 0.36)
    for i, positions in enumerate(data):
        print positions
        for xmin, xmax in positions:
            print 'start:%d, end: %d' % (xmin, xmax)
            plt.hlines(lineoffsets[i], xmin, xmax, colors[i], linewidth=linewidths[i])
    # plt.axes().set_aspect(0.6)
    # plt.show()
    plt.yticks(y, ylabels)
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.savefig(fig_file, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", "", dest="directory", help="The directory containing all the interval files")
    parser.add_option("-o", "", dest="outfile", default='cmp_fig', help="The final plot of intervals")
    (options, args) = parser.parse_args()

    data = get_data(options.directory)
    plot_intervals(data, options.outfile)
