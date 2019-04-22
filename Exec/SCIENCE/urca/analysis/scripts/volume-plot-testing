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
args = parser.parse_args()

# Hack: because rendering likes log fields ...
## create positive_radial_velocity and negative_radial_velocity fields.
## must do this before opening dataset
@derived_field(name='pos_radial_velocity', units='cm/s')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')], 1.0e-99)
@derived_field(name='neg_radial_velocity', units='cm/s')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')], 1.0e-99)

# Open Dataset
ds = yt.load(args.infile)
core = ds.sphere(ds.domain_center, (0.25e8, 'cm'))

# Create Scene
sc = Scene()

# Create Sources
so_enuc = VolumeSource(core, ('boxlib','enucdot'))
so_pos_vrad = VolumeSource(core, 'pos_radial_velocity')
so_neg_vrad = VolumeSource(core, 'neg_radial_velocity')

# Assign Transfer Functions to Sources
tfh_en = TransferFunctionHelper(ds)
tfh_en.set_field(('boxlib','enucdot'))
tfh_en.set_log(True)
tfh_en.set_bounds()
tfh_en.build_transfer_function()
tfh_en.tf.add_layers(10, colormap='black_green', w=0.01)
tfh_en.grey_opacity = False
tfh_en.plot('{}_tfun_enuc.png'.format(args.infile), profile_field=('boxlib','enucdot'))
so_enuc.transfer_function = tfh_en.tf

# tfh = TransferFunctionHelper(ds)
# tfh.set_log(True)
# tfh.set_bounds()
# tfh.build_transfer_function()
# tfh.tf.add_layers(5, colormap='Blues')
# tfh.grey_opacity = False
# tfh.plot('tfun_pos_vrad.png', profile_field="pos_radial_velocity")
# so_pos_vrad.set_fields("pos_radial_velocity", no_ghost=False)
# so_pos_vrad.set_transfer_function(tfh.tf)

# tfh = TransferFunctionHelper(ds)
# tfh.set_field("neg_radial_velocity", no_ghost=False)
# tfh.set_log(True)
# tfh.set_bounds()
# tfh.build_transfer_function()
# tfh.tf.add_layers(5, colormap='Reds')
# tfh.grey_opacity = False
# tfh.plot('tfun_neg_vrad.png', profile_field="neg_radial_velocity")
# so_neg_vrad.transfer_function = tfh.tf

# tfh_vr = TransferFunctionHelper(ds)
# tfh_vr.set_field(('boxlib','radial_velocity'))
# tfh_vr.set_log(True)
# tfh_vr.set_bounds()
# tfh_vr.build_transfer_function()
# tfh_vr.tf.add_layers(5, colormap='coolwarm_r')
# tfh_vr.grey_opacity = False
# tfh_vr.plot('tfun_vrad.png', profile_field=('boxlib','radial_velocity'))
# so_vrad.transfer_function = tfh_vr.tf

mag_vel_bounds = np.array([1.0e4, 1.0e6])
mag_vel_sigma  = 0.08

tfh = TransferFunctionHelper(ds)
tfh.set_field('pos_radial_velocity')
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
nlayers = 3
tfh.tf.add_layers(nlayers, colormap='Blues', w=mag_vel_sigma**2, mi=4.5, ma=5.5, alpha=np.logspace(-2,0,nlayers))
tfh.plot('{}_tfun_pos_vrad.png'.format(args.infile))
so_pos_vrad.transfer_function = tfh.tf

tfh = TransferFunctionHelper(ds)
tfh.set_field('neg_radial_velocity')
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(mag_vel_bounds)
tfh.build_transfer_function()
nlayers = 3
tfh.tf.add_layers(nlayers, colormap='Reds', w=mag_vel_sigma**2, mi=4.5, ma=5.5, alpha=np.logspace(-2,0,nlayers))
tfh.plot('{}_tfun_neg_vrad.png'.format(args.infile))
so_neg_vrad.transfer_function = tfh.tf


# field = "neg_radial_velocity"
# vals = [1.e3, 5.e3, 1.e4, 5.e4, 1.e5]
# bounds = (min(vals), max(vals))
# sigma = 2.e2
# cm = "Reds"
# tf = ColorTransferFunction(np.log10(bounds))
# for v in vals:
#     tf.sample_colormap(np.log10(v), np.log10(sigma**2), colormap=cm)
# so_neg_vrad.transfer_function = tf
# so_neg_vrad.set_field(field)
# so_neg_vrad.bounds = bounds
# so_neg_vrad.set_log(True)
# so_neg_vrad.tfh.grey_opacity = False

# Add sources to scene
sc.add_source(so_enuc)
sc.add_source(so_pos_vrad)
sc.add_source(so_neg_vrad)

# Add camera to scene
sc.add_camera()

# Set camera properties
sc.camera.focus = ds.domain_center
sc.camera.resolution = 2048
sc.camera.north_vector = [0, 0, 1]
sc.camera.position = ds.domain_center + [1.0, 1.0, 1.0] * ds.domain_width * 1.0/5.12
sc.camera.zoom(2.5)

# Annotate domain - draw boundaries
#sc.annotate_domain(ds, color=[1, 1, 1, 0.01])
# Annotate by drawing grids
#sc.annotate_grids(ds, alpha=0.01)
# Annotate by drawing axes triad
#sc.annotate_axes(alpha=0.01) 

# Render
sc.render()
sc.save('{}_rendering.png'.format(args.infile), sigma_clip=6)
