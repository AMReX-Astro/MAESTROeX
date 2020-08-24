# Rotating fully convective star diagnostics

Process a fully convective star problem to produce a set of diagnostic measures
including:

- **Convection speed**: The root-mean-squared radial velocity
$$v_c(r) = \sqrt(\int_0^{2\pi}d\psi \int_0^{\pi}d\theta \sin\theta v_r^2)$$

- **Ratio of mean rotation rates to surface rotation**: 
$$\text{ratio}(r) = \dfrac{\int dV \Omega(r,\theta,\psi)}{\int dV \Omega(\text{surface},\theta,\psi)}$$
Note that the denominator is currently set at the surface, but it does not strictly need to be at the surface; it just needs to be at a fixed radius.

- **Meridional circulation**: There are two components, radial and polar circulation 
$$u_r(r,\theta,t) = \dfrac{1}{t} \int_0^t dt \int_0^{2\pi}d\psi v_r$$, 
$$u_{\theta}(r,\theta,t) = \dfrac{1}{t} \int_0^t dt \int_0^{2\pi}d\psi v_{\theta}$$
where $v_r$ and $v_\theta$ are the radial and polar velocities.
The instantaneous circulation components are computed by neglecting the time-averaging terms.

- **Latitudinal shear**: Expressed as the projection of the rotation profile onto the $(l=2,m=0)$ spherical harmonic, averaged over the surface. This can be computed in Cartesian coordinates 
$$\Omega_2(r) = \int dV \Omega(x,y,z) Y_{2,0}(x,y,z) K(r'-r,\delta r)$$
using a normalized Gaussian kernel 
$$K(r'-r,\delta r) = \dfrac{1}{sqrt{\pi\delta r}} \exp[-(r'-r)^2/\delta r]$$
where $r'=\sqrt{x^2+y^2+z^2}$ and kernel width is
$$\delta r = f(r) dx^2 \quad , \quad f(r) = \min(r/L, 0.5)$$
and $L$ is the half length of the problem domain.
Here,
$$Y_{2,0}(x,y,z) = \dfrac{1}{4}\sqrt{\dfrac{5}{\pi}} \dfrac{2z^2-x^2-y^2}{r'^2}.

- **Brünt-Väisälä frequency**:
$$N^2 = -\dfrac{\gamma-1}{\gamma}\mathbf{g}\cdot\nabla s$$
where $\mathbf{g}$ is the acceleration due to gravity and $s$ is a dimensionless entropy.


## Building & running

The n-dimensional diagnostic can be built by executing `make DIM=n`. This will
produce the executable `radial_nd.exe`. Note that `DIM=3` is default setting and
the only one that is currently supported. To run, the executable must be provided
with the name of the plotfile to be analyzed:
```
./radial_3d.exe infile=plotfile_name
```

Additional arguments are as follows:

- **model file diagnostics**: extra argument `modelfile=modelfile_name` needs to be provided
- **exact solution test problem**: if `plotfile_name` is not provided, the code will write and evaluate a test problem with an exact solution.
- **time-averaged values**: two extra arguments need to be provided
  1. `dt=timesteps_between_plotfiles`, which is an integer value
  2. `nfiles=num_plotfiles` indicates the number of *additional* plotfiles needed, excluding the `infile` plotfile.
