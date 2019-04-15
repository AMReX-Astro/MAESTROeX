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
parser.add_argument('-cpos', '--camera_position', type=float, nargs=3,
                    help='3-D Camera position in fractions of maximum radius (--rup).')
parser.add_argument('-cnorth', '--camera_north', type=float, nargs=3, help='Camera north vector (direction of up).')
parser.add_argument('-lepmax', '--log_epos_max', type=float, default=12, help='Log of maximum magnitude for positive enucdot.')
parser.add_argument('-lepmin', '--log_epos_min', type=float, default=10, help='Log of minimum magnitude for positive enucdot.')
parser.add_argument('-lemmax', '--log_eneg_max', type=float, default=12, help='Log of maximum magnitude for negative enucdot.')
parser.add_argument('-lemmin', '--log_eneg_min', type=float, default=10, help='Log of minimum magnitude for negative enucdot.')
parser.add_argument('-n', '--num_layers', type=int, default=5, help='Number of layers for each of +/- enucdot. (Default is 5).')
parser.add_argument('-dd', '--drawdomain', action='store_true', help='If supplied, draw the boundaries of the domain.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-da', '--drawaxes', action='store_true', help='If supplied, draw an axes triad.')
parser.add_argument('-alpha_ones', '--alpha_ones', action='store_true',
                    help='If supplied, set the transfer function values to ones.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for volume-rendering.')
args = parser.parse_args()

# Hack: because rendering likes log fields ...
## create positive_radial_velocity and negative_radial_velocity fields.
## must do this before opening dataset
@derived_field(name='pos_enucdot', units='erg/(g*s)')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','enucdot')], yt.YTQuantity(1.0e-99, 'erg/(g*s)'))
@derived_field(name='neg_enucdot', units='erg/(g*s)')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','enucdot')], yt.YTQuantity(1.0e-99, 'erg/(g*s)'))

# Open Dataset
ds = yt.load(args.infile)
core = ds.sphere(ds.domain_center, (args.rup, 'cm'))

# Create Scene
sc = Scene()

# Create Sources
so_pos_enuc = VolumeSource(core, 'pos_enucdot')
so_neg_enuc = VolumeSource(core, 'neg_enucdot')

# Get maximum, minimum values for color settings
if args.log_epos_max:
    pos_maxv = args.log_epos_max
else:
    pos_maxv = np.ceil(np.log10(core.max('pos_enucdot')))
if args.log_eneg_max:
    neg_maxv = args.log_eneg_max
else:
    neg_maxv = np.ceil(np.log10(core.max('neg_enucdot')))
rat_neg_pos = 10.0**(neg_maxv-pos_maxv)

if args.log_epos_min:
    pos_minv = args.log_epos_min
else:
    pos_minv = pos_maxv-args.num_layers+1
if args.log_eneg_min:
    neg_minv = args.log_eneg_min
else:
    neg_minv = neg_maxv-args.num_layers+1

# Assign Transfer Functions to Sources
mag_enuc_sigma  = 0.1

if args.alpha_ones:
    alphavec = np.ones(args.num_layers)
else:
    alphavec = np.logspace(-2,0,args.num_layers)

tfh = TransferFunctionHelper(ds)
tfh.set_field('pos_enucdot')
mag_enuc_bounds = np.array([10.0**(pos_maxv-2*(pos_maxv-pos_minv)-3), 10.0**pos_maxv])
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_enuc_bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(args.num_layers, colormap='PiYG', w=mag_enuc_sigma**2, mi=pos_minv, ma=pos_maxv, alpha=alphavec)
try:
    tfh.plot('{}_tfun_pos_enucdot.png'.format(args.infile))
except:
    print('could not plot pos_enucdot transfer function.')
so_pos_enuc.transfer_function = tfh.tf

tfh = TransferFunctionHelper(ds)
tfh.set_field('neg_enucdot')
mag_enuc_bounds = np.array([10.0**(neg_maxv-2*(neg_maxv-neg_minv)-3), 10.0**neg_maxv])
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_enuc_bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(args.num_layers, colormap='PiYG_r', w=mag_enuc_sigma**2, mi=neg_minv, ma=neg_maxv, alpha=rat_neg_pos*alphavec)
try:
    tfh.plot('{}_tfun_neg_enucdot.png'.format(args.infile))
except:
    print('could not plot neg_enucdot transfer function.')
so_neg_enuc.transfer_function = tfh.tf

# Add sources to scene
sc.add_source(so_pos_enuc)
sc.add_source(so_neg_enuc)

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
sc.save('{}_rendering_enucdot.png'.format(args.infile), sigma_clip=3)
