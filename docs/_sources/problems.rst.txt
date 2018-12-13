**************
Problem Setups
**************

‘Standard’ Test Problems
========================

Different problems in MAESTROeX are contained in one of three
sub-directories under MAESTROeX/Exec: (SCIENCE/,
TEST_PROBLEMS/, or UNIT_TESTS/). The GNUmakefile in each
problem directory lists the components of MAESTRO that are used
to build the executable. TEST_PROBLEMS/ contains simple
problems that were used in the development of MAESTROeX. Many
of these were featured in the papers describing the MAESTROeX algorithm.

Some of the test problems available are:

-  | double_bubble
   | A rising bubble problem where the bubble(s) can have a different gamma
     than the surrounding atmosphere. This uses the multigamma EOS.

-  | incomp_shear_jet
   | A simple pure-incompressible shear layer problem. This is the example
     problem used in :raw-latex:`\cite{bellcolellaglaz}`. This is useful to see how to
     use MAESTROeX as an incompressible solver.

-  | reacting_bubble
   | reacting_bubble places 3 hots spots in a plane-parallel atmosphere.
     Burning makes these bubbles buoyant, and then roll up. This problem was
     used in :raw-latex:`\cite{lowMach3}` to compare with compressible solvers.

   This problem can also be run adaptively. The tag_boxes.f90
   file in the problem directory tags cells for refinement if the
   perturbational temperature, :math:`T^\prime`, exceeds some threshold.

-  | rt
   | [-3mm]

   A Rayleigh-Taylor instability problem. There are two methods that the
   code is run here, in the standard (using inputs_2d), the base state
   has the stratified atmosphere and we introduce a velocity perturbation
   to start the instability. The alternate method, inputs_2d_SNe, uses
   the do_smallscale runtime parameter to eliminate the base state
   and instead use the incompressible constraint to evolve the system.

-  | test_convect
   | test_convect drives convection through a plane-parallel
     atmosphere using an externally-specified heat source. This problem
     was used to compare with compressible solvers in :raw-latex:`\cite{lowMach3}`
     and to test the multilevel algorithm in :raw-latex:`\cite{multilevel}`.

-  | test_spherical
   | This problem sets up an isentropically stratified star and stirs it up
     with a random velocity field. The low Mach number constraint is
     replaced with the anelastic constraint (through
     the beta_type runtime parameter). Analytically, under
     these conditions, the density of the star should not change. This
     test problem was discussed in Maestro paper IV :raw-latex:`\cite{lowMach4}`.

Distributed Science Problems
============================

The following problems were used for science studies. It is
anticipated that more will be made available with time.

-  | flame
   | A combustion-mode problem where we model a thermonuclear flame in a
     small domain. This enforces the low Mach combustion constraint
     divU = S. Hot ash and cool fuel are put into contact and a flame
     will ignite and propagate across the grid. Inflow boundary
     conditions are used to allow for an inflow velocity to be set to
     keep the laminar flame stationary.

   In this mode, MAESTROeX behaves like the code described
   in :raw-latex:`\cite{SNe}`, which was used for models of Rayleigh-Taylor
   unstable flames :raw-latex:`\cite{SNld,SNrt,SNrt3d}`.

-  | flame_1d
   | A 1-d version of the flame problem above. This uses a special
     elliptic solver in AMReX that only works for a single grid, so
     no parallel runs are allowed for this problem.

-  | toy_convect
   | A nova-like problem for studying convection. This problem has seen
     extensive use in understanding which prediction types are the best
     when we have sharp species gradients. See Mike Z or Ryan for details.

-  | wdconvect
   | Model convection leading up to ignition in the Chandraseskhar-mass SNe
     Ia progenitor model. This setup was the basis for the simulations
     presented in :raw-latex:`\cite{lowMach4,wdconvect,wdturb}`.

.. _sec:gettingstarted:runtime:
