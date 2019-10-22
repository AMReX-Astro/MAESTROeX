"""
This script produces a slice plot for a given plotfile and variable. 
"""

import yt 
import numpy as np 
import sys
import argparse

def plot_single_var(plotfile_name, outputfile_name, var_name, use_log, norm_axis):
    # load data 
    ds = yt.load(plotfile_name)

    # make the slice plot 
    fig = yt.SlicePlot(ds, norm_axis, var_name)
    fig.set_log(var_name, use_log)

    fig.save(outputfile_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--log', action='store_const', const=True, help="plot with a log scale")
    parser.add_argument('-o', '--outfile', type=str, default="plot.png", help="name of output file")
    parser.add_argument('plotfile', type=str, help="name of plotfile")
    parser.add_argument('variable', type=str, help="variable to plot")
    parser.add_argument('-n', '--norm', type=str, default='z', help="Axis normal to the plot")

    args = parser.parse_args()

    use_log = args.log 
    plotfile_name = args.plotfile
    outputfile_name = args.outfile
    var_name = args.variable
    norm_axis = args.norm

    if use_log is None:
        use_log = False

    plot_single_var(plotfile_name, outputfile_name, var_name, use_log, norm_axis)