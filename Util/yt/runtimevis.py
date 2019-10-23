"""
This script is designed to be run from a submission script to produce plots from 
plotfiles as they are produced. It reads in variables and parameters from 
an inputs file.

If no output file name is provided, then it will append '.png' to the 
name of the plotfile and save it there.
"""

import yt 
import numpy as np 
import sys
import re
import argparse
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams["mathtext.fontset"] = "stix"

def read_inputs_file(inputs_file):
    varname_re = re.compile(r"\[([-\w()]+)\]")
    values_re = re.compile(r"(\w+)\s*=\s*(\w+)")

    # we're going to store all of the parameters in a dictionary of dictionaries
    variables = {}

    with open(inputs_file, 'r') as f:

        txt = f.read()

        # find all the variables
        matches = [m for m in re.finditer(varname_re, txt)]

        for i, m in enumerate(matches):
            var = m.group(1)
            variables[var] = {}

            # we only want to search up to either the end of file or to the 
            # position of the start of the next variable name
            if i == len(matches)-1:
                endpoint = -1
            else:
                endpoint = matches[i+1].start()

            # find all the parameters associated with that variable 
            for p in re.finditer(values_re, txt[m.end():endpoint]):

                param = p.group(1)
                value = p.group(2)

                variables[var][param] = value

    return variables
    

def runtime_vis(plotfile_name, outputfile_name, inputs_file):

    variables = read_inputs_file(inputs_file)

    if outputfile_name is None:
        if plotfile_name[-1] == '/':
            outputfile_name = plotfile_name[:-1] + '.png'
        else:
            outputfile_name = plotfile_name + '.png'
        
    # load data 
    ds = yt.load(plotfile_name)

    fig = plt.figure(figsize=(16,9))

    nvars = len(variables.keys())

    # calculate grid dimensions
    if nvars == 1:
        grid_dims = (1,1)
    elif nvars % 2 == 0:
        grid_dims = (nvars // 2, 2)
    elif nvars % 3 == 0:
        grid_dims = (nvars // 3, 3)
    else:
        # here we give up and just put them in a column
        grid_dims = (nvars, 1)

    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = grid_dims,
                axes_pad = 1.7,
                label_mode = "all",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

    plots = yt.SlicePlot(ds, 'z', list(variables.keys()))

    # iterate over the variables, assign them to the grid and update their 
    # parameters
    for i, (var, params) in enumerate(variables.items()):

        p = plots.plots[var]
        p.figure = fig 
        p.axes = grid[i].axes 
        p.cax = grid.cbar_axes[i]

        plots.set_zlim(var, params.get('min', 'min'), params.get('max', 'max'))

        use_log = params.get('log', '0')
        if use_log == '1':
            use_log = True
        else:
            use_log = False 

        plots.set_log(var, use_log)

    plots._setup_plots()
            
    fig.savefig(outputfile_name, bbox_inches='tight', dpi=80)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', type=str, help="name of output file")
    parser.add_argument('plotfile', type=str, help="name of plotfile")
    parser.add_argument('-i', '--inputs', type=str, default="vis.in", help="name of inputs file")
    args = parser.parse_args()

    runtime_vis(args.plotfile, args.outfile, args.inputs)