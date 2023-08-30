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


W,H,_ = ds.domain_right_edge
nx,ny,_ = ds.domain_dimensions
x = np.linspace(-W.v/2,W.v/2,nx)
y = np.linspace(-H.v/2,H.v/2,ny)

# Create fig
fig,axes = plt.subplots(1,2,figsize=(5,4))
fig.subplots_adjust(wspace=0)

# del minus del ad
slc = yt.SlicePlot(ds,'z','ad_excess')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
z = frb['ad_excess']
img = axes[0].imshow(z, origin="lower", extent=(-W/2,0,-H/2,H/2), cmap="bwr", vmin=-0.5, vmax=0.5)

r = 1.0 # width ratio of cbar to plot
pos = axes[0].get_position()
x0,y0,h,w = pos.x0,pos.y0,pos.height,pos.width
cax1 = fig.add_axes([x0+(1-r)*w/2, y0+0.87*h, w*r, h/20])

cbar = plt.colorbar(img, cax=cax1, orientation="horizontal", ticks=[-1,1])
cax1.xaxis.set_ticks_position("top")
cax1.set_xlabel(r"$\nabla-\nabla_{\rm ad}$")
cax1.xaxis.set_label_coords(0.5, 4)
cax1.set_xticks([-0.5,0,0.5])

# heat flux
slc = yt.SlicePlot(ds,'z','rho')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
rho = frb['rho']

slc = yt.SlicePlot(ds,'z','vely')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
vy = frb['vely']

slc = yt.SlicePlot(ds,'z','tpert')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
dT = frb['tpert']

z = rho*vy*dT
# img = axes[1].imshow(z, origin="lower", extent=(0,W/2,-H/2,H/2), cmap="inferno", vmin=16, vmax=20)
img = axes[1].imshow(z, origin="lower", extent=(0,W/2,-H/2,H/2), cmap="coolwarm", vmin=-1e20, vmax=1e20)


pos = axes[1].get_position()
x0,y0 = pos.x0,pos.y0
cax2 = fig.add_axes([x0+(1-r)*w/2, y0+0.87*h, w*r, h/20])

cbar = plt.colorbar(img, cax=cax2, orientation="horizontal")
cax2.xaxis.set_ticks_position("top")
cax2.set_xlabel(r"$\rho v_y\;\delta T$")
cax2.xaxis.set_label_coords(0.5, 4)
cax2.set_xticks([-1e20,1e20],[r"$-10^{20}$","$10^{20}$"])

axes[0].set_ylim([-800,800])
axes[1].set_ylim([-800,800])
axes[0].set_yticks([-500,0,500],["-5","0","5"])
axes[1].set_yticks([])
axes[0].set_xlim([-1500,0])
axes[1].set_xlim([0,1500])
axes[0].set_xticks([-1500,-1000,-500],["-15","-10","-5"])
axes[1].set_xticks([0,500,1000],["0","5","10"])
axes[0].set_ylabel("y (meters)")
axes[0].set_xlabel("x (meters)")
axes[0].xaxis.set_label_coords(1., -0.15)


# Streamplot
slc = yt.SlicePlot(ds,'z','velx')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
vx = np.array(frb['velx'])
slc = yt.SlicePlot(ds,'z','vely')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
vy = np.array(frb['vely'])

V = np.sqrt(vx**2+vy**2)
lw = V/V.max()

axes[0].streamplot(x,y,vx,vy, color='k', density=2, linewidth=lw, arrowsize=0.7)
axes[1].streamplot(x,y,vx,vy, color='k', density=2, linewidth=lw, arrowsize=0.7)

fig.suptitle(fr"time = {float(ds.current_time):8.5f} s", y=1.05)

basename = os.path.basename(plotfile)
plt.savefig("plt2/" + basename + ".png", bbox_inches="tight", dpi=500)
