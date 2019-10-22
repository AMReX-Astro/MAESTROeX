"""
This script takes two or three plotfiles and a single variable 
as arguments and plots contours of the datasets on the same 
set of axes. 
"""

import yt 
import numpy as np 
import sys
import argparse
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams["mathtext.fontset"] = "stix"

def contour_compare(plotfiles, outputfile_name, var, norm_axis):

    # dictionary of coordinates so can map our norm axis
    geometries = {'cartesian':['x', 'y', 'z'], 
                  'cylindrical': ['r', 'z', 'theta'], 
                  'spherical':['r', 'theta', 'phi']}

    linestyles = ['-', ':', '--']

    fig, ax = plt.subplots()

    labels = []
    contours = []

    # for each plot, extract the data then plot as a contour 
    for i, pf in enumerate(plotfiles):

        # load data 
        ds = yt.load(pf)

        # first take a slice through the center
        coords = geometries[ds.geometry]
        norm_idx = coords.index(norm_axis)
        dims = ds.domain_dimensions

        dat = ds.all_data()[var]
        dat = dat.reshape(dims)

        if dims[2] == 1 or norm_idx == 2:
            dat = dat[:,:,dims[2]//2]

            X = ds.all_data()[coords[0]].reshape(dims)[:,:,dims[2]//2]
            Y = ds.all_data()[coords[1]].reshape(dims)[:,:,dims[2]//2]

            ax_labels = coords[:2]
        elif norm_idx == 1:
            dat = dat[:,dims[1]//2,:]

            X = ds.all_data()[coords[0]].reshape(dims)[:,dims[1]//2,:]
            Y = ds.all_data()[coords[2]].reshape(dims)[:,dims[1]//2,:]

            ax_labels = coords[::2]
        else:
            dat = dat[ds.dims[0]//2,:,:]

            X = ds.all_data()[coords[1]].reshape(dims)[ds.dims[0]//2,:,:]
            Y = ds.all_data()[coords[2]].reshape(dims)[ds.dims[0]//2,:,:]

            ax_labels = coords[1:]

        # draw the contours
        cntr = ax.contour(X.v, Y.v, dat.v, linestyles=linestyles[i], 
                          colors=f"C{i}", levels=3)
        c,_ = cntr.legend_elements()
        contours.append(c[0])
        labels.append(pf.split('/')[-1])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    ax.set_aspect('equal')
    plt.legend(contours, labels)

    fig.savefig(outputfile_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', type=str, default="plot.png", help="name of output file")
    parser.add_argument('variable', type=str, help="variable to plot")
    parser.add_argument('plotfiles', type=str, nargs='+', help="name of plotfile")
    parser.add_argument('-n', '--norm', type=str, default='z', help="Axis normal to the plot")

    args = parser.parse_args()

    contour_compare(args.plotfiles, args.outfile, args.variable, args.norm)