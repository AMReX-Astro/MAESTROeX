*************************
Common Runtime Parameters
*************************

Controlling Timestepping and Output
===================================

Parameters that set the maximum time for the simulation to run
include:

-  stop_time is the maximum simulation time, in seconds,
   to evolve the system for.

-  max_step is the maximum number of steps to take.

Parameters affecting the size of the timestep include:

-  cflfac is a multiplicative factor (:math:`\le 1`)
   applied to the advective CFL timestep

-  init_shrink is the factor (:math:`\le 1`) by which to reduce
   the initial timestep from the estimated first timestep.

Parameters affecting output and restart include:

-  restart tells MAESTROeX to restart from a checkpoint. The
   value of this parameter should be the file number to restart from.
   For example, to restart from the checkpoint file chk00010,
   you would set restart = 10.

-  plot_int is the number of steps to take between
   outputting a plotfile

-  plot_deltat is the simulation time to evolve between
   outputting a plotfile. Note: to output only based on simulation
   time, set plot_int = -1.

-  check_int is the number of steps to take between
   outputting a checkpoint.

-  plot_base_name is the basename to use for the
   plotfile filename. The step number will be appended to
   this name.

Note that in addition to the normal plotfiles, there are *mini* plotfiles
that store a small subset of the fields, intended to be output more frequently.
These are described in § \ `[vis:sec:miniplotfile] <#vis:sec:miniplotfile>`__.

Defining the Grid and Boundary Conditions
=========================================

Parameters that determine the spatial extent of the grid,
the types of boundaries, and the number of computational cells include:

-  max_levs is the maximum number of grid levels in the AMR
   hierarchy to use. max_levs = 1 indicates running with only a
   single level spanning the whole domain.

-  n_cellx, n_celly, n_cellz the size of
   base level in terms of number of cells, in the :math:`x`, :math:`y`, and :math:`z`
   coordinate directions.

-  max_grid_size the maximum extend of a grid, in any
   coordinate direction, as measured in terms of number of cells.

   For multilevel problems, the parameter max_grid_size_1
   controls the maximum extent on level 1 (the base
   grid), max_grid_size_2 controls the maximum extent on
   level 2, and max_grid_size_3 controls the maximum extent on
   levels 3 and higher.

-  prob_lo_x, prob_lo_y, prob_lo_z is
   the physical coordinate of the lower extent of the domain boundary
   in the :math:`x`, :math:`y`, and :math:`z` coordinate directions.

-  prob_hi_x, prob_hi_y, prob_hi_z is
   the physical coordinate of the upper extent of the domain boundary
   in the :math:`x`, :math:`y`, and :math:`z` coordinate directions.

-  There are two ways to specify boundary conditions—via integer flags
   or descriptive string names. If the string names are present,
   then they will override the integer quantities in determining
   the boundary conditions.

   -  bcx_lo, bcy_lo, bcz_lo
      , bcx_hi, bcy_hi, bcz_hi are the
      boundary condition types at the lower (‘lo’) and upper
      (‘hi’) domain boundaries in the :math:`x`, :math:`y`, and :math:`z`
      coordinate directions. The different types are set via integer
      flags listed in table \ `[gs:table:bcflags] <#gs:table:bcflags>`__.

      .. table:: [gs:table:bcflags] Boundary condition flags

         +----------------------+--------------+
         | BC type              | integer flag |
         +======================+==============+
         | periodic             | :math:`-1`   |
         +----------------------+--------------+
         | inlet (user-defined) | :math:`11`   |
         +----------------------+--------------+
         | outlet               | :math:`12`   |
         +----------------------+--------------+
         | symmetry             | :math:`13`   |
         +----------------------+--------------+
         | slip wall            | :math:`14`   |
         +----------------------+--------------+
         | no-slip wall         | :math:`15`   |
         +----------------------+--------------+

   -  xlo_boundary_type, ylo_boundary_type, zlo_boundary_type, xhi_boundary_type, yhi_boundary_type, zhi_boundary_type
      are the boundary condition types at the lower and upper domain
      boundaries in the :math:`x`, :math:`y`, and :math:`z` coordinate directions. The
      boundary type is set by providing a string name—valid values are
      listed in table \ `[gs:table:bcstrings] <#gs:table:bcstrings>`__

      .. table:: [gs:table:bcstrings] Boundary condition string names

         +----------------------+----------------+
         | BC type              | integer flag   |
         +======================+================+
         | periodic             | “periodic”     |
         +----------------------+----------------+
         | inlet (user-defined) | “inlet”        |
         +----------------------+----------------+
         | outlet               | “outlet”       |
         +----------------------+----------------+
         | symmetry             | “symmetry”     |
         +----------------------+----------------+
         | slip wall            | “slip wall”    |
         +----------------------+----------------+
         | no-slip wall         | “no slip wall” |
         +----------------------+----------------+

   The string-based parameters are a newer option for specifying
   boundary conditions, and are preferred due to clarity. The
   conversion between the string names and integer flags is done
   by AMReX in the bc_module at the time of initializing
   the runtime parameters.

Note that grid cells must be square, i.e. :math:`\Delta x = \Delta y = \Delta z`
where :math:`\Delta x` on the base grid is computed as :math:`({\tt prob\_hi\_x}
- {\tt prob\_lo\_x})/{\tt n\_cellx}`. For multilevel problems, the effective
number of zones on the finest grid in the :math:`x` direction will be
:math:`{\tt n\_cellx} \cdot 2^{({\tt max\_levels} -1)}`.
