#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams.update({
    "font.family": "serif",
    "xtick.direction" : "in",
    "ytick.direction" : "in",
    "xtick.top" : True,
    "ytick.right" : True,
})

# assume that our data is in CGS
from yt.units import cm, amu
from yt.frontends.boxlib.api import MaestroDataset

# Disable excessive printing
yt.set_log_level(50)

plotfile = sys.argv[1]
if plotfile[-1] == '/': plotfile=plotfile[:-1]
ds = yt.load(plotfile, hint="maestro")

w,h,_ = ds.domain_right_edge
nx,ny,_ = ds.domain_dimensions

#fig,ax = plt.subplots(1,1,figsize=(3,7))
fig,ax = plt.subplots(1,1)

slc = yt.SlicePlot(ds,'z',"vort")
frb = slc.data_source.to_frb(width=w,height=h,resolution=(nx,ny))

z = np.array(frb["vort"])
img = ax.imshow(z, origin="lower", extent=(-w/2,w/2,-h/2,h/2), cmap="PiYG", vmin=-1e5, vmax=1e5)

# cbar
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
divider = make_axes_locatable(ax)
height = axes_size.AxesX(ax, aspect=1/20)
pad = axes_size.Fraction(0.5, height)
cax = divider.append_axes("top", size=height, pad=pad)
cbar = plt.colorbar(img, cax=cax, orientation="horizontal", ticks=[-1e5,1e5])
cax.xaxis.set_ticks_position("top")
#cax.set_xlabel(r"$\nabla\times U$ (1/s)")
#cax.xaxis.set_label_coords(0.5, 6)

ax.set_xlabel(r"x (cm)")
ax.set_ylabel(r"y (cm)")

fig.suptitle(fr"time = {float(ds.current_time):8.5f} s")

fig.tight_layout()
basename = os.path.basename(plotfile)
plt.savefig("plt_vort/" + basename + ".png", bbox_inches="tight", dpi=500)
