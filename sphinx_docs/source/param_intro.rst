.. _sec:runtime_parameters:

******************
Runtime Parameters
******************


Introduction to Runtime Parameters
==================================

MAESTROeX runtime parameters are set
in the inputs file and managed by the AMReX ``ParmParse``
class. For MAESTROeX-specific parameters, we list the runtime
parameters in a file ``Source/param/_cpp_parameters`` and generate the
C++ code and headers automatically at compilation time.

The behavior of the network, EOS, and other microphysics routines are
controlled by a different set of runtime parameters. These parameters are defined
in plain-text files ``_parameters`` located in the different
directories that hold the microphysics code. At compile time, a
script in the AMReX build system, ``findparams.py``, locates all
of the ``_parameters`` files that are needed for the given choice
of network, integrator, and EOS, and assembles all of the runtime
parameters into a module named ``extern_probin_module`` (using the
``write_probin.py`` script). The parameters are set in the ``&probin`` namelist
in the problem inputs file.

C++ parameter format
--------------------

The C parameters take the form of::

    # comment describing the parameter
    name   type   default   need in Fortran?   ifdef    fortran name    fortran type

Here,

  * ``name`` is the name of the parameter that will be looked for
    in the inputs file.

    The name can actually take the form of ``(a, b)``, where ``a`` is
    the name to be used in the inputs file where the parameter is set
    and ``b`` is the name used within the MAESTROeX C++ class.  It is not
    recommended to name new parameters with this functionality—this
    was implemented for backwards compatibility.


  * ``type`` is one of int, Real, or string

  * ``default`` is the default value of the parameter.

The next columns are outdated.

Finally, any comment (starting with ``#``) immediately before the
parameter definition will be used to generate the documentation
describing the parameters.

Microphysics/extern parameter format
------------------------------------

The microphysics/extern parameter definitions take the form of::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the ``priority`` is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.


Common Runtime Parameters
=========================

Controlling Timestepping and Output
-----------------------------------

Parameters that set the maximum time for the simulation to run
include:

-  ``stop_time`` is the maximum simulation time, in seconds,
   to evolve the system for.

-  ``max_step`` is the maximum number of steps to take.

Parameters affecting the size of the timestep include:

-  ``cfl`` is a multiplicative factor (:math:`\le 1`)
   applied to the advective CFL timestep

-  ``init_shrink`` is the factor (:math:`\le 1`) by which to reduce
   the initial timestep from the estimated first timestep.

Parameters affecting output and restart include:

-  ``restart_file`` tells MAESTROeX to restart from a checkpoint. The
   value of this parameter should be the name of the checkpoint file to
   restart from. For example, to restart from the checkpoint file ``chk00010``,
   you would set ``maestro.restart_file = chk00010``.

-  ``plot_int`` is the number of steps to take between
   outputting a plotfile

-  ``plot_deltat`` is the simulation time to evolve between
   outputting a plotfile. Note: to output only based on simulation
   time, set ``maestro.plot_int = -1``.

-  ``chk_int`` is the number of steps to take between
   outputting a checkpoint.

-  ``plot_base_name`` is the basename to use for the
   plotfile filename. The step number will be appended to
   this name.

Note that in addition to the normal plotfiles, there are *small* plotfiles
that store a small subset of the fields, intended to be output more frequently.
These are described in § \ `[vis:sec:miniplotfile] <#vis:sec:miniplotfile>`__.

Defining the Grid and Boundary Conditions
-----------------------------------------

Parameters that determine the spatial extent of the grid,
the types of boundaries, and the number of computational cells include:

-  ``amr.max_level`` is the maximum number of grid levels in the AMR
   hierarchy to use. ``amr.max_level = 0`` indicates running with only a
   single level spanning the whole domain.

-  ``amr.n_cell`` is the size of the
   base level in terms of number of cells. It takes a list of integer numbers,
   defining the size in the :math:`x`, :math:`y`, and :math:`z`
   coordinate directions, e.g. ``amr.n_cell = 64 64 64`` would define a domain
   of size :math:`64^3`

-  ``amr.max_grid_size`` is the maximum extent of a grid, in any
   coordinate direction, as measured in terms of number of cells.

   For multilevel problems, this parameter can take a list of values,
   each value determining the maximum extent at each level.

-  ``geometry.prob_lo`` is
   the physical coordinate of the lower extent of the domain boundary.
   It takes a list of values defining the coordinates
   in the :math:`x`, :math:`y`, and :math:`z` coordinate directions.

-  ``geometry.prob_hi`` is
   the physical coordinate of the upper extent of the domain boundary.
   It takes a list of values defining the coordinates
   in the :math:`x`, :math:`y`, and :math:`z` coordinate directions.

-  Boundary conditions are specified via integer flags.

   -  ``maestro.lo_bc`` and  ``maestro.hi_bc`` are the
      boundary condition types at the lower (‘lo’) and upper
      (‘hi’) domain boundaries in the :math:`x`, :math:`y`, and :math:`z`
      coordinate directions. The different types are set via integer
      flags listed in table :ref:`table-bc-flags`.

      .. _table-bc-flags:

      .. table:: Boundary condition flags

         +----------------------+--------------+
         | BC type              | integer flag |
         +======================+==============+
         | Interior             | :math:`0`    |
         +----------------------+--------------+
         | Inflow               | :math:`1`    |
         +----------------------+--------------+
         | Outflow              | :math:`2`    |
         +----------------------+--------------+
         | Symmetry             | :math:`3`    |
         +----------------------+--------------+
         | Slip wall            | :math:`4`    |
         +----------------------+--------------+
         | No-slip wall         | :math:`5`    |
         +----------------------+--------------+

   -  Periodic boundary conditions are set using the parameter
      ``geometry.is_periodic``. This takes a list of values, with ``1`` indicating
      periodic boundary conditions (and ``0`` non-periodic). So to set periodic
      boundary conditions in the :math:`x`-direction only, you would set
      ``geometry.is_periodic = 1 0 0``.

Note that grid cells must be square, i.e. :math:`\Delta x = \Delta y = \Delta z`
where :math:`\Delta x` on the base grid is computed as :math:`({\tt prob\_hi\_x}
- {\tt prob\_lo\_x})/{\tt n\_cellx}`. For multilevel problems, the effective
number of zones on the finest grid in the :math:`x` direction will be
:math:`{\tt n\_cellx} \cdot 2^{({\tt max\_levels} -1)}`.



Parameters by Namespace
=======================

The following tables list all of the general
MAESTROeX runtime parameters.
These tables are generated automatically using the
``rp.py`` script in ``MAESTROeX/sphinx_docs`` by parsing
the ``MAESTROeX/Source/param/_cpp_parameters`` file. The problem-specific parameters
are not shown here.

.. toctree::

   runtime_parameters
