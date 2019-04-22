#!/usr/bin/env python
import yt
import numpy as np
import matplotlib.pylab as pl

from yt.visualization.api import Streamlines
from yt.units import cm
from mpl_toolkits.mplot3d import Axes3D

# Load the dataset
ds = yt.load('wd_512_rhoc4-5_plt64636')

# Define c: the center of the box, N: the number of streamlines,
# scale: the spatial scale of the streamlines relative to the boxsize,
# and then pos: the random positions of the streamlines.
c = ds.domain_center
N = 100
scale = ds.domain_width[0]
pos_dx = np.random.random((N,3))*scale
pos = pos_dx

# Create streamlines of the 3D vector velocity and integrate them through
# the box defined above
streamlines = Streamlines(ds, pos, ('boxlib', 'x_vel'), ('boxlib', 'y_vel'), ('boxlib', 'z_vel'), get_magnitude=True, volume=ds.all_data())
streamlines.integrate_through_volume()

# Create a 3D plot, trace the streamlines throught the 3D volume of the plot
fig=pl.figure()
ax = Axes3D(fig)
for stream in streamlines.streamlines:
    stream = stream[np.all(stream != 0.0, axis=1)]
    ax.plot3D(stream[:,0], stream[:,1], stream[:,2], alpha=0.1)

# Save the plot to disk.
pl.savefig('vel_streamlines.png')
