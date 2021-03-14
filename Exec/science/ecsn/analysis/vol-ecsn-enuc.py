#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import sys

import yt
import numpy as np
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource


# this is for the wdconvect problem

def doit(plotfile):

    ds = yt.load(plotfile)
    ds.periodicity = (True, True, True)

    field = ('boxlib', 'Hnuc')
    ds._get_field_info(field).take_log = True
        
    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = VolumeSource(ds, field=field)
    sc.add_source(vol)


    # transfer function
    vals = [14, 14.5, 15, 15.5, 16]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cm = "viridis"

    for v in vals:
        if v < 15.5:
            alpha = 0.1
        else:
            alpha = 0.75

        tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cm)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")        
    cam.resolution = (1080, 1080)
    cam.position = 1.0*ds.domain_right_edge
    
    # look toward the center -- we are dealing with an octant
    center = ds.domain_left_edge
    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])
    cam.set_width(0.5*ds.domain_width)
    cam.zoom(1.5)
    sc.camera = cam
    #sc.annotate_axes(alpha=0.05)
    #sc.annotate_domain(ds, color=np.array([0.05, 0.05, 0.05, 0.05]))
    #sc.annotate_grids(ds, alpha=0.05)

    sc.render()
    sc.save("{}_Hnuc".format(plotfile), sigma_clip=4.0)
    sc.save_annotated("{}_Hnuc_annotated.png".format(plotfile), 
                      text_annotate=[[(0.05, 0.05), 
                                      "t = {}".format(ds.current_time.d),
                                      dict(horizontalalignment="left")],
                                     [(0.5,0.95), 
                                      "MAESTROeX simulation of ECSN convection",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)


        
