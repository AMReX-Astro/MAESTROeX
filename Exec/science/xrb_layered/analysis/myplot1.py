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

# h1
slc = yt.SlicePlot(ds,'z','X(h1)')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
z = np.log10(np.maximum(frb['X(h1)'],1e-99))
img = axes[0].imshow(z, origin="lower", extent=(-W/2,0,-H/2,H/2), cmap="Blues", vmin=-3, vmax=0)

r = 1.0 # width ratio of cbar to plot
pos = axes[0].get_position()
x0,y0,h,w = pos.x0,pos.y0,pos.height,pos.width
cax1 = fig.add_axes([x0+(1-r)*w/2, y0+0.87*h, w*r, h/20])

cbar = plt.colorbar(img, cax=cax1, orientation="horizontal", ticks=[-1,1])
cax1.xaxis.set_ticks_position("top")
cax1.set_xlabel(r"$\log X(^1{H})$")
cax1.xaxis.set_label_coords(0.5, 4)
cax1.set_xticks([-3,-2,-1,0])

# Hnuc
slc = yt.SlicePlot(ds,'z','Hnuc')
frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
z = np.log10(np.maximum(frb['Hnuc'],1e-99))
img = axes[1].imshow(z, origin="lower", extent=(0,W/2,-H/2,H/2), cmap="inferno", vmin=16, vmax=20)

pos = axes[1].get_position()
x0,y0 = pos.x0,pos.y0
cax2 = fig.add_axes([x0+(1-r)*w/2, y0+0.87*h, w*r, h/20])

cbar = plt.colorbar(img, cax=cax2, orientation="horizontal")
cax2.xaxis.set_ticks_position("top")
cax2.set_xlabel(r"log $H_{\rm nuc}$ (erg g$^{-1}$ s$^{-1}$)")
cax2.xaxis.set_label_coords(0.5, 4)
cax2.set_xticks([16,18,20])

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
# slc = yt.SlicePlot(ds,'z','velx')
# frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
# vx = np.array(frb['velx'])
# slc = yt.SlicePlot(ds,'z','vely')
# frb = slc.data_source.to_frb(width=W,height=H,resolution=(nx,ny))
# vy = np.array(frb['vely'])

# V = np.sqrt(vx**2+vy**2)
# lw = V/V.max()

# axes[0].streamplot(x,y,vx,vy, color='k', density=2, linewidth=lw, arrowsize=0.7)
# axes[1].streamplot(x,y,vx,vy, color='k', density=2, linewidth=lw, arrowsize=0.7)


# Swz Ledoux boundaries from conv grad averages file
avg_fname = ds.filename + "/convgrad/averages"
header = np.loadtxt(avg_fname, dtype='unicode', max_rows=1)
# print(header)
data = np.loadtxt(avg_fname, skiprows=1).T
z = data[0]
z = (z-z[-1]/2) # center and convert to meters
del_ = data[list(header).index('del')]
del_ad = data[list(header).index('del_ad')]
del_l = data[list(header).index('del_ledoux')]

iswz = np.argwhere(np.logical_and(z > -400 , del_-del_ad < -0.2))[0]
plt.axvline(z[iswz], color='k', ls=':')
iled = np.argwhere(np.logical_and(z > -400 , del_-del_l < -0.2))[0]


# # Plot Ledoux line?
# # Yes, if it's far enough from Schwarzschild. But also once we plot it, we always don't want it to disappear, so create a flag file in the directory
# # If flag is here, plot Ledoux
# plot_Ledoux = False
# if os.path.exists("flag_ledoux"):
#     plot_Ledoux = True
# elif abs(z[iled]-z[iswz])>50:
#     print("Start plotting Ledoux because %d cm away from Schwarzschild"%abs(z[iled]-z[iswz]))
#     open("flag_ledoux", "w")
#     plot_Ledoux = True

# This doesn't quite work because we're running these plot scripts in parallel
plot_Ledoux = True


for ax in axes:
    ax.axhline(z[iswz], color='r', ls='--', lw=0.7)
    if plot_Ledoux:
        ax.axhline(z[iled], color='g', ls='--', lw=0.7)

axes[1].text(+12e2,z[iswz]+90,r"$\bar\nabla=\bar\nabla_{\rm ad}$", color='r', ha='center', va='center')
if plot_Ledoux:
    axes[0].text(-12e2,z[iled]-90,r"$\bar\nabla=\bar\nabla_{\rm L}$", color='g', ha='center', va='center')


# Save
fig.suptitle(fr"time = {float(ds.current_time):8.5f} s", y=1.05)
basename = os.path.basename(plotfile)
plt.savefig("plt1/" + basename + ".png", bbox_inches="tight", dpi=500)
