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
parser.add_argument('-cpos', '--camera_position', type=float, nargs=3, help='3-D Camera position in fractions of maximum radius (--rup).')
parser.add_argument('-cnorth', '--camera_north', type=float, nargs=3, help='Camera north vector (direction of up).')
parser.add_argument('-vmin', '--velocity_minimum', type=float, default=1.0e3, help='Minimum velocity for transfer function. (Default is 1.0e3 cm/s).')
parser.add_argument('-vmax', '--velocity_maximum', type=float, default=1.0e7, help='Maximum velocity for transfer function. (Default is 1.0e7 cm/s).')
parser.add_argument('-vsig', '--velocity_sigma', type=float, default=0.08, help='Velocity transfer function width parameter. (Default is 0.08).')
parser.add_argument('-n', '--num_layers', type=int, default=5, help='Number of layers for each of +/- velocity. (Default is 5).')
parser.add_argument('-dd', '--drawdomain', action='store_true', help='If supplied, draw the boundaries of the domain.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-da', '--drawaxes', action='store_true', help='If supplied, draw an axes triad.')
parser.add_argument('-alpha_ones', '--alpha_ones', action='store_true', help='If supplied, set the transfer function values to ones.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
parser.add_argument('-dry', '--dry_run', action='store_true', help='Plot only the transfer functions and quit.')
args = parser.parse_args()

# Hack: because rendering likes log fields ...
## create positive_radial_velocity and negative_radial_velocity fields.
## must do this before opening dataset
@derived_field(name='pos_radial_velocity', units='cm/s')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')], yt.YTQuantity(1.0e-99, 'cm/s'))
@derived_field(name='neg_radial_velocity', units='cm/s')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')], yt.YTQuantity(1.0e-99, 'cm/s'))

# Open Dataset
ds = yt.load(args.infile)
core = ds.sphere(ds.domain_center, (args.rup, 'cm'))

# Create Scene
sc = Scene()

# Create Sources
#so_enuc = VolumeSource(core, ('boxlib','enucdot'))
so_pos_vrad = VolumeSource(core, 'pos_radial_velocity')
so_neg_vrad = VolumeSource(core, 'neg_radial_velocity')

# Assign Transfer Functions to Sources
# tfh_en = TransferFunctionHelper(ds)
# tfh_en.set_field(('boxlib','enucdot'))
# tfh_en.set_log(True)
# tfh_en.set_bounds()
# tfh_en.build_transfer_function()
# tfh_en.tf.add_layers(10, colormap='black_green', w=0.01)
# tfh_en.grey_opacity = False
# tfh_en.plot('{}_tfun_enuc.png'.format(args.infile), profile_field=('boxlib','enucdot'))
# so_enuc.transfer_function = tfh_en.tf

mag_vel_bounds = np.array([args.velocity_minimum, args.velocity_maximum])
mag_vel_sigma  = args.velocity_sigma
log_min = np.log10(args.velocity_minimum)
log_max = np.log10(args.velocity_maximum)

nlayers = args.num_layers
if args.alpha_ones:
    alphavec = np.ones(nlayers)
else:
    alphavec = np.logspace(-3, 0, num=nlayers, endpoint=True)

tfh = TransferFunctionHelper(ds)
tfh.set_field('pos_radial_velocity')
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(nlayers, colormap='Blues', w=mag_vel_sigma**2, mi=log_min, ma=log_max, alpha=alphavec)
tfh.plot('{}_tfun_pos_vrad.png'.format(args.infile))
so_pos_vrad.transfer_function = tfh.tf

tfh = TransferFunctionHelper(ds)
tfh.set_field('neg_radial_velocity')
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(nlayers, colormap='Reds', w=mag_vel_sigma**2, mi=log_min, ma=log_max, alpha=alphavec)
tfh.plot('{}_tfun_neg_vrad.png'.format(args.infile))
so_neg_vrad.transfer_function = tfh.tf

if args.dry_run:
    exit()

# Add sources to scene
#sc.add_source(so_enuc)
sc.add_source(so_pos_vrad)
sc.add_source(so_neg_vrad)

# Add camera to scene
sc.add_camera()

# Set camera properties
sc.camera.focus = ds.domain_center
sc.camera.resolution = args.resolution
sc.camera.north_vector = yt.YTArray(args.camera_north, 'cm')
sc.camera.position = ds.domain_center + yt.YTArray(args.camera_position, 'cm') * args.rup
sc.camera.zoom(args.zoom)

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
sc.save('{}_rendering_rad-vel.png'.format(args.infile), sigma_clip=3)
