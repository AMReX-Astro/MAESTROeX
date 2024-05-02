**********
Unit Tests
**********

In addition to the MAESTROeX science problems, which use the full
capabilities of MAESTROeX, there are a number of unit tests that
exercise only specific components of the MAESTROeX solvers. These
tests have their own drivers (a custom varden.f90) that
initialize only the data needed for the specific test and call
specific MAESTROeX routines directly.

test_advect
===========

This test initializes a Gaussian density field (no other scalar
quantities are used) and a uniform velocity field in any one of the
coordinate directions. The Gaussian profile is advected through the
period domain exactly once and the error in the density profile (L2
norm) is computed. The driver for this problem does this for every
dimension twice (once with a positive velocity and once with a
negative velocity), and loops over all advection methods
(ppm_type = 0,1,2 and bds_type = 1). After
all coordinate directions are tested, the norms are compared to
ensure that the error does not show any directional bias.

Note: the BDS advection method does not guarantee that the error be
independent of the advection direction—small differences can
arise. What’s happening here is that within each cell BDS is trying
to define a tri-linear profile for rho subject to the constraints in
the BDS paper (:cite:`bds3d`) (constraints 1 and 2 on p. 2044 after eq. 3.4). We
do not solve the L2 minimization problem exactly so we iterate up to
3 times using a simple heuristic algorithm designed to work toward
the constraint. The iteration will give different results depending
on orientation since we work through the corners in arbitrary order.

.. figure:: dens_2d_orig_density.png
   :align: center
   :width: 80%

   Initial Gaussian density profile

.. figure:: dens_2d_ppm1_xp_final_density.png
   :align: center
   :width: 80%

   density profile after advecting to the right for one period

.. figure:: dens_2d_ppm1_xp_final_abserror.png
   :align: center
   :width: 80%

   absolute error between the final and initial density fields,
   showing the error in the advection scheme

test_average
============

This test initializes a 1D radial base state quantity with a
Gaussian distribution, maps it into the 3D domain (assuming a
spherical geometry) using the routines provided by
the fill_3d_module module, and then calls average to
put it back onto a 1D radial array. This way we test the accuracy
of our procedure to map between the 1D radial and 3D Cartesian
states. The output from this test was described in detail
in :cite:`multilevel`.

test_basestate
==============

This test initializes the base state to contain a hydrostatic
model and then evolves the state with heating to watch the
hydrostatic adjustment of the atmosphere. In particular,
the base state velocity, :math:`w_0`, is computed in response to
the heating and this is used to advect the base state density
and compute the new pressure, :math:`p_0`. An early version of
this routine was used for the plane-parallel expansion test
in :cite:`lowMach2`. This version of the test was also shown
for a spherical, self-gravitating star in :cite:`multilevel`.

test_diffusion
==============

This test initializes a Gaussian temperature profile and calls
the thermal diffusion routines in MAESTROeX to evolve the state
considering only diffusion. The driver estimates a timestep
based on the explicit thermal diffusion timescale and loops
over calls to the thermal diffusion solver. A Gaussian remains
Gaussian when diffusing, so an explicit error can be computed
by comparing to the analytic solution. This test is
described in :cite:`xrb`.

test_eos
========

This test sets up a 3-d cube with :math:`\rho` varying on one axis, :math:`T` on
another, and the composition on the third. The EOS is then called
in every zone, doing :math:`(\rho, T) \rightarrow  p, h, s, e` and stores those
quantities. Then it does each of the different EOS types to recover
either :math:`T` or :math:`\rho` (depending on the type), and stores the new :math:`T` (or
:math:`\rho`) and the relative error with the original value. A plotfile is
stored holding the results and errors. This allows us to determine
whether the EOS inversion routines are working right.

.. test_particles
.. ==============
..
.. This test exercises the particle advection routine. A simple
.. circular velocity field, with the magnitude increasing with radius
.. from the center is initialized. A number of particles are then
.. initialized at various radii from the center and they are advected
.. for one period. The particle paths should be perfect circles, and
.. the final particle position should overlap with the initial
.. position.
..
.. Particle data is stored separately from the fluid data. Instead
.. of being part of the plotfiles, the particle data is outputted
.. each timestep into files named ``timestamp_NN``, where
.. the number indicates which processor did the writing. These
.. particle files can be processed and the particle data plotted
.. using the python routines in ``data_processing/python/``.
..
.. The output from this test can be visualized with the script
.. plot.py in the test directory. The output shows the particle
.. paths (see below):
..
.. .. figure:: particle_paths.png
..    :align: center
..    :width: 80%
..
..    Particle paths for the test_particles problem. The initial
..    position of the particles is marked with an :math:`\times`.

test_projection
===============

This tests the projection routines in 2- and 3-d—either the hgprojection
(project_type = 1) or the MAC projection (project_type =
2). A divergence-free velocity field is initialized and then
“polluted” by adding the gradient of a scalar. The form of the
scalar differs depending on the boundary conditions (wall and
periodic are supported currently). Finally, the hgproject routine
is called to recover the initial divergence-free field.
The figures below show the initial field, polluted
field, and result of the projection for the hgproject case.

.. figure:: wall_u_init_x-velocity.png
   :align: center
   :width: 60%

   Initial divergence free velocity field (x-component)

.. figure:: wall_u_plus_grad_phi_x-velocity.png
   :align: center
   :width: 60%

   Velocity field plus gradient of a scalar (x-component)

.. figure:: wall_u_new_x-velocity.png
   :align: center
   :width: 60%

   Resulting velocity after projecting out the non-divergence free
   portion (x-component).

This is with slipwall boundary conditions on all sides, a 2-level grid
with the left half refined and right half coarse, and the hgprojection
tested.

.. figure:: test_project_3d.png
   :align: center
   :width: 80%

   Projection test in 3-d showing the x-velocity (left), y-velocity
   (middle), and z-velocity (right) initially (top row), after the
   gradient of a scalar is added (center row), and the resulting
   velocity after the projection. This is with slipwall boundary conditions
   on all sides, a 2-level grid with an octant refined, and the hgprojection.

test_react
==========

This simply tests the reaction network by calling
the MAESTROeX react_state routine directly. The network is
selected in the GNUmakefile by setting the ``NETWORK_DIR``
variable. A 3d cube is setup with density varying on one axis,
temperature varying on another, and the composition varying on the
third. The density and temperature ranges are set in the inputs
file. The composition is read in via an input file.

A good use of this test is to test whether a burner is threadsafe.
This is accomplished by compiling with OpenMP (setting OMP=t)
and the running with 1 thread and multiple threads (this can be done
by setting the environment variable ``OMP_NUM_THREADS`` to the
desired number of threads). Since each zone is independent of the
others, the results should be identical regardless of the number
of threads. This can be confirmed using the fcompare tool
in ``BoxLib/Tools/Postprocessing/F_Src/``.
