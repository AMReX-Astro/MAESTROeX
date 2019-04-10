#!/usr/bin/env python
import yt
from yt.units import dimensions
from yt import derived_field
import yt.visualization.volume_rendering.api as vr
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.api import Scene, VolumeSource, Camera, ColorTransferFunction
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-rup', '--rup', type=float, default=1.0e8, help='Maximum radius (cm). Default 1.0e8.')
parser.add_argument('-zoom', '--zoom', type=float, default=1.0, help='Camera zoom factor. Default 1.0.')
parser.add_argument('-dd', '--drawdomain', action='store_true', help='If supplied, draw the boundaries of the domain.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-da', '--drawaxes', action='store_true', help='If supplied, draw an axes triad.')
parser.add_argument('-alpha_ones', '--alpha_ones', action='store_true', help='If supplied, set the transfer function values to ones.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
args = parser.parse_args()

# Open Dataset
ds = yt.load(args.infile)
core = ds.sphere(ds.domain_center, (args.rup, 'cm'))

# Create Scene
sc = Scene()

# Create Sources
so_circum_vel = VolumeSource(core, ('boxlib', 'circum_velocity'))

mag_vel_bounds = np.array([1.0e1, 1.0e6])
mag_vel_sigma  = 0.08

nlayers = 6
if args.alpha_ones:
    alphavec = np.ones(nlayers)
else:
    alphavec = np.logspace(-5,0,nlayers)

tfh = TransferFunctionHelper(ds)
tfh.set_field(('boxlib', 'circum_velocity'))
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(nlayers, colormap='viridis', w=mag_vel_sigma**2, mi=1, ma=6, alpha=alphavec)
tfh.plot('{}_tfun_circum_vel.png'.format(args.infile))
so_circum_vel.transfer_function = tfh.tf

# Add sources to scene
sc.add_source(so_circum_vel)

# Add camera to scene
sc.add_camera()

# Set camera properties
sc.camera.focus = ds.domain_center
sc.camera.resolution = args.resolution
sc.camera.north_vector = [0, 0, 1]
sc.camera.position = ds.domain_center + [1.0, 1.0, 1.0] * ds.domain_width * args.rup/5.12e8
sc.camera.zoom(2.5*args.zoom)

# Annotate domain - draw boundaries
if args.drawdomain:
    sc.annotate_domain(ds, color=[1, 1, 1, 0.01])

# Annotate by drawing grids
if args.drawgrids:
    sc.annotate_grids(ds, alpha=0.01)

# Annotate by drawing axes triad
if args.drawaxes:
    sc.annotate_axes(alpha=0.01) 

# Render
sc.render()
sc.save('{}_rendering_circum-vel.png'.format(args.infile), sigma_clip=6)
