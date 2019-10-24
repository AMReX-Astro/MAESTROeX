"""
This script produces a slice plot for a given plotfile and variable. 
"""

import yt 
import sys
import argparse
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams["mathtext.fontset"] = "stix"

def plot_single_var(plotfile_name, outputfile_name, var_names, 
                    use_log, norm_axis, minimum, maximum):

    if outputfile_name is None:
        suffix = '_' + '_'.join(var_names) + '.png'
        if plotfile_name[-1] == '/':
            outputfile_name = plotfile_name[:-1] + suffix
        else:
            outputfile_name = plotfile_name + suffix
        
    # load data 
    ds = yt.load(plotfile_name)

    fig = plt.figure()

    nvars = len(var_names)

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

    # make the slice plot 
    plots = yt.SlicePlot(ds, norm_axis, var_names)
        
    if minimum is None:
        minimum = 'min'
    if maximum is None:
        maximum = 'max'

    # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
    # axes.
    for i, var in enumerate(var_names):

        p = plots.plots[var]
        p.figure = fig
        p.axes = grid[i].axes 
        p.cax = grid.cbar_axes[i]

        plots.set_zlim(var, minimum, maximum)
        plots.set_log(var, use_log)

    plots._setup_plots()
            
    fig.savefig(outputfile_name, bbox_inches='tight')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--log', action='store_const', const=True, help="plot with a log scale")
    parser.add_argument('-o', '--outfile', type=str, help="name of output file")
    parser.add_argument('plotfile', type=str, help="name of plotfile")
    parser.add_argument('variables', type=str, nargs='+', help="variable(s) to plot")
    parser.add_argument('-n', '--norm', type=str, default='z', help="Axis normal to the plot")
    parser.add_argument('-min', '--minimum', type=str, help="Set minimum data value")
    parser.add_argument('-max', '--maximum', type=str, help="Set maximum data value")

    args = parser.parse_args()

    use_log = args.log 

    if use_log is None:
        use_log = False

    plot_single_var(args.plotfile, args.outfile, args.variables, use_log, args.norm, args.minimum, args.maximum)