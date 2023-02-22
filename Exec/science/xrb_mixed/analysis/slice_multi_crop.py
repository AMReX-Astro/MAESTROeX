#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm, amu
from yt.frontends.boxlib.api import MaestroDataset

plotfile = sys.argv[1]
ds = MaestroDataset(plotfile)

xmin = ds.domain_left_edge[0]
xmax = ds.domain_right_edge[0]
xctr = 0.5*(xmin + xmax)
L_x = xmax - xmin


ymin = ds.domain_left_edge[1]
ymax = ds.domain_right_edge[1]
yctr = 0.5*(ymin + ymax)
L_y = ymax - ymin


fig = plt.figure()
fig.set_size_inches(12.0, 9.0)


fields = ["tfromp", "MachNumber", "Hnuc", "vort"] #, "density"]

grid = ImageGrid(fig, 111, nrows_ncols=(1, len(fields)),
                 axes_pad=1.0, cbar_pad=0.2, label_mode="L", cbar_mode="each")


for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "z", f, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="12")
    sp.set_buff_size((2400,2400))

    if f == "tfromp":
        sp.set_zlim(f, 5.e7, 1.e9)
        sp.set_cmap(f, "viridis")
    elif f == "Hnuc":
        sp.set_zlim(f, 1.e16, 1.e19)
        sp.set_cmap(f, "plasma_r")
    elif f == "vort":
        sp.set_zlim(f, -2.e6, 2.e6)
        sp.set_log(f, False)
        sp.set_cmap(f, "bwr")
    elif f == "MachNumber":
        sp.set_zlim(f, 1.e-4, 1.e-1)
        sp.set_log(f, True)
        sp.set_cmap(f, "cividis")

    #if f != "density":
    #    # now do a contour of density
    #    sp.annotate_contour("density", ncont=2, clim=(1.e2, 2.e6),
    #                        plot_args={"colors": "0.5", "linewidths": 1, "linestyle": ":"})

    sp.set_axes_unit("cm")

    #sp.annotate_text((0.05, 0.05), "{:8.5f} s".format(float(ds.current_time.in_cgs())),
    #                 coord_system="figure", text_args={"color": "black"})

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)

    sp._setup_plots()

fig.text(0.05, 0.05, "time = {:8.5f} s".format(float(ds.current_time)), transform=fig.transFigure)

fig.set_size_inches(19.2, 10.8)
plt.tight_layout()
plt.savefig("{}_slice.png".format(os.path.basename(plotfile)))
