#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
import os
import argparse

"""
Make a turbulent KE power spectrum.  Since we are stratified, we use
a rho**(1/3) scaling to the velocity to get something that would
look Kolmogorov (if the turbulence were fully developed).

Ultimately, we aim to compute:

                      1  ^      ^*                                           
     E(k) = integral  -  V(k) . V(k) dS                                      
                      2                                                      
 
             n                                                              
where V = rho  U is the density-weighted velocity field.
 
(Note: sometimes we normalize by 1/volume to get a spectral
energy density spectrum).
"""

parser = argparse.ArgumentParser()
parser.add_argument('plotfile', type=str, help='Path to plotfile for which to compute the power spectrum.')
parser.add_argument('-l', '--level', type=int, default=-1, help='Level of refinement at which to compute the spectrum (starts at 0). Defaults to maximum level.')
parser.add_argument('-r', '--radius', type=float, help='Radius of the sphere inscribing the power spectrum domain in cm. Width of the power spectrum domain is 2*radius.')
args = parser.parse_args()

def doit(ds):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    if args.level >= 0 and args.level <= ds.index.max_level:
        max_level = args.level
    else:
        max_level = ds.index.max_level
    print('using covering grid at level {}'.format(max_level))

    if max_level == 0:
        ref = 1
    else:
        ref = int(np.product(ds.ref_factors[0:max_level]))

    if args.radius:
        L = yt.YTArray([2.0*args.radius, 2.0*args.radius, 2.0*args.radius], 'cm')
        low = ds.domain_center - 0.5*L
        L_domain = ds.domain_right_edge - ds.domain_left_edge
        dims_domain = ds.domain_dimensions*ref
        dims = np.ceil((dims_domain * L/L_domain).d)
    else:
        low = ds.domain_left_edge
        dims = ds.domain_dimensions*ref
        L = (ds.domain_right_edge - ds.domain_left_edge).d

    nx, ny, nz = dims

    nindex_rho = 1./3.

    Kk = np.zeros((int(np.floor(nx/2))+1, int(np.floor(ny/2))+1, int(np.floor(nz/2))+1))

    for vel in [("gas", "velocity_x"), ("gas", "velocity_y"), 
                ("gas", "velocity_z")]:

        Kk += 0.5*fft_comp(ds, ("gas", "density"), vel,
                           nindex_rho, max_level, low, dims)

    # wavenumbers
    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]
    
    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)
    
    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)
    
    E_spectrum = np.zeros(len(ncount)-1)

    for n in range(1,len(ncount)):
        E_spectrum[n-1] = np.sum(Kk.flat[whichbin==n])

    k = 0.5*(kbins[0:N-1] + kbins[1:N])
    E_spectrum = E_spectrum[1:N]

    # Sometimes there will be a spike at the largest wavenumber
    # Omit such a spike from this index calculation
    index = np.argmax(E_spectrum[:-5])
    kmax = k[index]
    Emax = E_spectrum[index]

    fig, ax = plt.subplots()
    ax.loglog(k, E_spectrum)
    ax.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$E(k)dk$")
    plt.savefig("{}_spectrum.png".format(os.path.basename(args.plotfile)), dpi=600)


def fft_comp(ds, irho, iu, nindex_rho, level, low, dimensions):

    print('calculating covering grid')
    cube = ds.covering_grid(level, left_edge=low,
                            dims=dimensions,
                            fields=[irho, iu])

    print('got covering grid')

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
    ru = np.fft.fftn(rho**nindex_rho * u)[0:int(np.floor(nx/2))+1,0:int(np.floor(ny/2))+1,0:int(np.floor(nz/2))+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2


if __name__ == "__main__":

    ds = yt.load(args.plotfile)
    doit(ds)
