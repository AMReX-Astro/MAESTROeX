**************
Problem Setups
**************

‘Standard’ Test Problems
========================

Different problems in MAESTROeX are contained in one of three
sub-directories under ``MAESTROeX/Exec/``: ``SCIENCE/``,
``TEST_PROBLEMS/``, or ``UNIT_TESTS/``).  The ``GNUmakefile`` in each
problem directory lists the components of MAESTROeX that are used
to build the executable. ``TEST_PROBLEMS/`` contains simple
problems that were used in the development of MAESTROeX. Many
of these were featured in the papers describing the original MAESTRO algorithm.

Some of the test problems available are:

* ``double_bubble`` : A rising bubble problem where the bubble(s) can
  have a different gamma than the surrounding atmosphere. This uses
  the multigamma EOS.

* ``incomp_shear_jet`` : A simple pure-incompressible shear layer
  problem. This is the example problem used in
  :cite:`bellcolellaglaz`. This is useful to see how to use MAESTROeX
  as an incompressible solver.

* ``reacting_bubble`` : The reacting bubble test places 3 hots spots
  in a plane-parallel atmosphere.  Burning makes these bubbles
  buoyant, and then roll up. This problem was used in :cite:`lowMach3`
  to compare with compressible solvers.

  This problem can also be run adaptively. The tagging.f90 file in the
  problem directory tags cells for refinement if the temperature
  exceeds some threshold.

* ``rt`` : A Rayleigh-Taylor instability problem.  Here, the base
  state has the stratified atmosphere and we introduce a velocity
  perturbation to start the instability.

* ``test_convect`` : This is a convection test that drives convection
  through a plane-parallel atmosphere using an externally-specified
  heat source. This problem was used to compare with compressible
  solvers in :cite:`lowMach3` and to test the multilevel
  algorithm in :cite:`multilevel`.


Distributed Science Problems
============================

* ``wdconvect`` : Model convection leading up to ignition in the
  Chandraseskhar-mass SNe Ia progenitor model. This setup was the
  basis for the simulations presented in
  :cite:`lowMach4,wdconvect,wdturb`.
