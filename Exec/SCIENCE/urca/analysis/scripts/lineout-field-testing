#!/usr/bin/env python

import matplotlib.pyplot as plt
import yt
from yt import derived_field
import yt.visualization.volume_rendering.api as vr
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.api import Scene, VolumeSource, Camera, ColorTransferFunction
import numpy as np

# Open Dataset
ds = yt.load('wd_512_rhoc4-5_plt64636')

# Hack: because rendering likes log fields ...
## create positive_radial_velocity and negative_radial_velocity fields.
@derived_field(name="pos_radial_velocity", units="cm/s")
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')], 1.0e-99)

@derived_field(name="neg_radial_velocity", units="cm/s")
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')], 1.0e-99)
    
# Create lineout along x axis through the center
c = ds.domain_center
ax = 0 # take a line cut along the x axis

# cutting through the y0,z0 such that we hit the center
ray = ds.ortho_ray(ax, (c[1], c[2]))

# Sort the ray values by 'x' so there are no discontinuities
# in the line plot
srt = np.argsort(ray['x'])

plt.subplot(211)
plt.semilogy(np.array(ray['x'][srt]), np.array(ray['pos_radial_velocity'][srt]))
plt.ylabel('pos rad vel')
plt.ylim([1.0e-1, 1.0e6])
plt.subplot(212)
plt.semilogy(np.array(ray['x'][srt]), np.array(ray['neg_radial_velocity'][srt]))
plt.xlabel('x')
plt.ylim([1.0e-1, 1.0e6])
plt.ylabel('neg rad vel')

plt.savefig("rad_vel_xsweep.png")
