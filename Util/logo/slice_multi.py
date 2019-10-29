#!/usr/bin/env python3

import os
import sys
import yt
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm

plotfile = sys.argv[1]

field = "vort"

# slice plot of temperature
ds = yt.load(plotfile)

xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

sp = yt.SlicePlot(ds, "z", field)
sp.set_buff_size((2000,2000))

sp.set_cmap(field, "bwr")
sp.set_log(field, False)
sp.set_axes_unit("cm")

sp.save("vort.png")
