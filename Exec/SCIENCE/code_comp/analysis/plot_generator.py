#!/usr/bin/env python3

import yt
from yt.units import amu, cm

import os
import sys
import glob
import argparse
import numpy as np
import string

from collections import namedtuple
from functools import reduce

def parse_args():
    # Argument information
    description = """Generates plots of datasets using a specified yt plot function. Works with any slice or projection
            plot, as well as ParticlePlot."""
    datasets_help = "A list of datasets to be loaded by yt. Will be sorted by plot number by default."
    func_help = "The plotting function to use. SlicePlot by default."
    out_help = "The desired output directory for the image files."
    var_help = "The variable to plot. Set to 'Temp' by default."
    bounds_help = "The bounds for the colorbar."
    cmap_help = "The colormap for the variable to plot."
    log_help = "If provided, sets the plot to a logarithmic scale."
    linthresh_help = "If provided, sets the linear threshold for a symlog plot"
    time_help = "If provided, adds a timestamp to each plot with the given precision."
    ext_help = "The extension of the file format to save to. PNG by default."
    sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding the decimal point
            give the starting index, digits following the decimal point give the number of characters. Make negative for
            descending order."""
    xlim_help = "The x-axis limits."
    ylim_help = "The y-axis limits."
    zlim_help = "The z-axis limits."
    normal_help = "The normal direction"

    # Construct parser and parse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('datasets', nargs='*', help=datasets_help)
    parser.add_argument('-f', '--func', default='SlicePlot', help=func_help)
    parser.add_argument('-o', '--out', default='', help=out_help)
    parser.add_argument('-v', '--var', default='Temp', help=var_help)
    parser.add_argument('-b', '--bounds', nargs=2, type=float, metavar=('LOWER', 'UPPER'), help=bounds_help)
    parser.add_argument('-c', '--cmap', metavar=('NAME',), help=cmap_help)
    parser.add_argument('--log', action='store_true', help=log_help)
    parser.add_argument('--linthresh', type=float, help=linthresh_help)
    parser.add_argument('-t', '--time', type=int, metavar=('PRECISION',), help=time_help)
    parser.add_argument('-e', '--ext', type=lambda s: s.lower(), default='png', help=ext_help)
    parser.add_argument('-s', '--sort', type=float, default=0.0, help=sort_help)
    parser.add_argument('-x', '--xlim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=xlim_help)
    parser.add_argument('-y', '--ylim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=ylim_help)
    parser.add_argument('-z', '--zlim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=zlim_help)
    parser.add_argument('-n', '--normal', default='z', help=normal_help)

    return parser.parse_args(sys.argv[1:])


def plot_generator(args):
    coloropts = ['field_color', 'cmap', 'display_threshold', 'cbar']
    ColorOpt = namedtuple('ColorOpt', field_names=coloropts)
    optdict = dict(field_color=None, display_threshold=None, cmap=None, cbar=False)

    color_opt = ColorOpt(**optdict)

    # Make output directory
    if not args.out:
        args.out = os.getcwd()
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # Grab files from working directory if none were specified
    ts = args.datasets
    if not ts:
        ts = glob.glob('plt*')

    # Exit if nothing could be loaded
    if len(ts) < 1:
        sys.exit("No files were available to be loaded.")

    # Sort and load files
    desc = args.sort < 0
    start = abs(int(args.sort))
    nchars = int(str(args.sort).split('.')[1])

    if nchars == 0:
        key = lambda fname: fname[start:]
    else:
        key = lambda fname: fname[start:start + nchars]
    ts.sort(key=key, reverse=desc)

    tf = lambda file: yt.load(file.rstrip('/'))
    ts = list(map(tf, ts))
    print("Successfully loaded the following files: {}\n".format(ts))

    # Generate plots
    func = getattr(yt, args.func)
    field = args.var

    def get_width(ds, xlim=None, ylim=None, zlim=None):
        """ Get the width of the plot. """

        if xlim is None: 
            xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
        else: 
            xlim = xlim[0] * cm, xlim[1] * cm

        if ylim is None: 
            ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
        else: 
            ylim = ylim[0] * cm, ylim[1] * cm

        xwidth = (xlim[1] - xlim[0]).in_cgs()
        ywidth = (ylim[1] - ylim[0]).in_cgs()

        if ds.domain_dimensions[2] == 1:
            zwidth = 0.0
        else:
            if zlim is None: 
                zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
            else: 
                zlim = zlim[0] * cm, zlim[1] * cm

            zwidth = (zlim[1] - zlim[0]).in_cgs()

        return xwidth, ywidth, zwidth

    def get_center(ds, xlim=None, ylim=None, zlim=None):
        """ Get the coordinates of the center of the plot. """

        if xlim is None: 
            xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
        else: 
            xlim = xlim[0] * cm, xlim[1] * cm

        if ylim is None: 
            ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
        else: 
            ylim = ylim[0] * cm, ylim[1] * cm

        xctr = 0.5 * (xlim[0] + xlim[1])
        yctr = 0.5 * (ylim[0] + ylim[1])

        if ds.domain_dimensions[2] == 1:
            zctr = 0.0
        else:
            if zlim is None: 
                zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
            else: 
                zlim = zlim[0] * cm, zlim[1] * cm

            zctr = 0.5 * (zlim[0] + zlim[1])

        return xctr, yctr, zctr

    print("Generating...")

    # Loop and generate
    for ds in ts:
        
        settings = {}
        settings['center'] = get_center(ds, args.xlim, args.ylim, args.zlim)
        settings['width'] = get_width(ds, args.xlim, args.ylim, args.zlim)
        settings['normal'] = args.normal

        plot = func(ds, fields=field, **settings)
        if args.cmap: 
            plot.set_cmap(field=field, cmap=args.cmap)

        if args.linthresh:
            plot.set_log(field, args.log, linthresh=args.linthresh)
        else:
            plot.set_log(field, args.log)

        # print(args.bounds)

        # sys.exit()

        if args.bounds is not None:
            plot.set_zlim(field, *args.bounds)

        if args.time:
            
            time_format = f't = {{time:.{args.time}f}}{{units}}'
            
            plot.annotate_timestamp(corner='upper_left', time_format=time_format,
                    time_unit='s', draw_inset_box=True, inset_box_args={'alpha': 0.0})

        suffix = args.func.replace('Plot', '').lower()
        plot.save(os.path.join(args.out, f'{ds}_{field.translate(str.maketrans("","", string.punctuation))}_{suffix}.{args.ext}'))
        print()

    print("Task completed.")


if __name__ == "__main__":
    args = parse_args()
    plot_generator(args)