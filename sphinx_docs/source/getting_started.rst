***************
Getting Started
***************

In this chapter we give an overview of MAESTROeX, including some of the
standard problems, how to run the code, some basic runtime parameters,
and how to look at the output.

Quick Start
===========

Here we will run the standard reacting_bubble problem (three
reacting bubbles in a plane-parallel stratified atmosphere) on a
single processor [1]_. This test problem was shown in
paper 3.

#. *Get a copy of MAESTROeX*.

   If you don’t already have a copy of MAESTROeX, you can obtain one
   from the project’s github page:
   https://github.com/AMReX-Astro/MAESTROeX  . There are several
   options: you can fork it directly on github (recommended if
   you intend to develop the code) or clone it using git from the
   project page.

   MAESTROeX is under active development, so you will want to keep in
   sync with the changes by periodically pulling from the repository.
   Simply type

   ::

       git pull
         

   in the MAESTROeX/ directory.

#. *Get a copy of Microphysics*.

   MAESTROeX and its compressible counterpart CASTRO share a
   common-set of microphysics solvers (nuclear reaction networks and
   equations of state). These are kept in a separate repo.
   Microphysics is also available on github and can be obtained
   via:

   ::

       git clone https://github.com/starkiller-astro/Microphysics.git
         

   You will periodic want to update Microphysics by doing

   ::

       git pull
         

   in the Microphysics/ directory
   .

#. *Get a copy of AMReX*.

   MAESTROeX requires the AMReX library to manage the grids and
   parallelization. We also rely on the build system in AMReX to
   build a MAESTROeX executable. AMReX is also available on github
   and can be obtained via:

   ::

       git clone https://github.com/AMReX-Codes/amrex.git
         

   You will periodic want to update AMReX by doing

   ::

       git pull
         

   in the amrex/ directory.

#. *Setup your shell environment*.

   MAESTROeX needs to know where to find AMReX, by specifying the
   AMREX_HOME environment variable, and where
   to find Microphysics, bt specifying the
   MICROPHYSICS_HOME environment variable.

   If your shell is Bash, add

   ::

       export AMREX_HOME="/path/to/amrex/"
       export MICROPHYSICS_HOME="/path/to/Microphysics/"
         

   to your .bashrc.

   If your shell is Csh/Tcsh, add

   ::

       setenv AMREX_HOME /path/to/amrex/
       setenv MICROPHYSICS_HOME /path/to/Microphysics/
         

   to your .cshrc.

   Note: you must specify the full path to the AMReX/ and
   Microphysics/ directory. Do not use “:math:`\sim`” to refer to your
   home directory—the scripts used by the build system will not be
   able to process this.

#. *Setup the problem’s GNUmakefile*.

   In MAESTROeX, each problem lives under one of three sub-directories
   of MAESTROeX/Exec: SCIENCE/, TEST_PROBLEMS/, or
   UNIT_TESTS/. This problem sub-directory will contain any
   problem-specific files as well as the GNUmakefile that
   specifies how to build the executable. Note: we rely on features of
   GNU make. Full details of the GNUmakefile can be found
   in § \ `[sec:adding_problems] <#sec:adding_problems>`__. Here we will configure for a
   serial build.

   Change directory to
   MAESTROeX/Exec/TEST_PROBLEMS/reacting_bubble/.
   We only need to worry about the options at the very top of the
   GNUmakefile for now. These should be set as follows:

   -  DEBUG := TRUE

      This option determines whether we compile with support for
      less-optimized code with debugging runtime checks. Setting
      DEBUG := FALSE turns off debugging.

   -  DIM := 2

      The dimensionality of the problem must be specified at compile-time.

   -  COMP := gnu

      This option specifies the gnu compiler suite (g++/gfortran).
      We will use gnu, which is the preferred compiler suite for MAESTROeX.
      Specifying this compiler will automatically pull in the compiler
      settings as specified in AMREX_HOME/Tools/GNUMake/Make.defs.
      (Alternate compiler choices include
      intel, cray, and pgi.

   -  USE_MPI := TRUE

      This determines whether we are doing a parallel build, using the
      Message Passing Interface (MPI) library. If you set this option
      to FALSE, you will disable MPI
      and will build MAESTROeX in serial
      mode, so no MPI library needs to be present on the host system.

   -  USE_OMP := FALSE

      This determines whether we are using OpenMP to do parallelism
      within a shared memory node. OpenMP is used together with MPI,
      with MPI distributing the grids across the processors and within a
      shared-memory node, OpenMP allows many cores to operate on the
      same grid. For now, we leave this option as FALSE, disabling OpenMP.

   -  USE_REACT := TRUE

      Some test problems in MAESTROeX do not use reactions, so there is an
      option to disable the compilation of reaction-specific source code.

   -  TINY_PROFILE := FALSE

      Profiling tool that generates a text file with profiling information.
      Refer to the AMReX User’s Guide at
      https://amrex-codes.github.io/amrex/

   -  PROFILE := FALSE

      More advanced profiling tool that generates a text file with profiling
      information, or data files that can be interpreted with a special build of
      amrvis. Selecting TRUE overrides the TINY_PROFILE setting.
      Refer to the AMReX User’s Guide at https://amrex-codes.github.io/amrex/

#. *Build the executable*.

   Type make. The build system will first find the dependencies
   amongst all the source files and then build the executable. When
   finished, the executable will have a name like
   Maestro2d.gnu.DEBUG.MPI.ex, where the specific parts of the name
   depend on the options used in GNUmakefile.

   Note, at the end of the build process, a link will be made in the
   current directory to the data table needed for the equation of state
   (Microphysics/EOS/helmholtz/helm_table.dat).

#. *Run!*

   Each problem requires an input file. The inputs file
   consists of lines of the form *parameter = value*,
   where *parameter* is one of the many runtime parameters
   MAESTROeX knows, and *value* overrides the default value for
   that parameter. For the reacting_bubble problem, we will use
   the inputs file inputs_2d_C. An overview of some of the more
   common runtime parameters is given in
   § \ `5 <#sec:gettingstarted:runtime>`__, and a full list of all
   MAESTROeX runtime parameters and their default values is given in
   Chapter \ `[ch:runtimeparameters] <#ch:runtimeparameters>`__.

   MAESTROeX is run simply as:

   ::

         ./Maestro2d.gnu.DEBUG.MPI.ex inputs_2d_C
         

   or to run in parallel on a local workstation:

   ::

         mpiexec -n 4 ./Maestro2d.gnu.DEBUG.MPI.ex inputs_2d_C
         

   We can also override the default value of any runtime parameter by
   specifying it on the commandline as, e.g.,

   ::

         ./Maestro2d.gnu.DEBUG.MPI.ex inputs_2d_C maestro.max_step=0 amr.n_cell=192 320
         

   As the code runs, a lot of information will pass through the screen.
   For each timestep, each of the steps 1 through 12 shown in the
   MAESTROeX flowchart (Chapter `[ch:flowchart] <#ch:flowchart>`__) will be shown along
   with diagnostic information about the solution. Upon completion
   some memory usage information is printed.

#. *Examine the output*.

   As the code runs, it will output both plotfiles and checkpoints as
   well as one or more text diagnostic files (maestro_diag.out
   by default) with integral or extrema information (like maximum Mach
   number) from each timestep.

   By default, the plotfiles will be named plt\ *nnnnnnn*, where
   the number *nnnnnnn* is the timestep number when the file was
   outputted. Similarly, the checkpoints are named
   chk\ *nnnnnnn*. AMReX plotfiles and checkpoints are actually
   directories, with the data stored in sub-directories grouped by
   refinement level. Details of the simulation (build information,
   number of processors used, output date, output directory, runtime
   parameter values, ...) are stored in the plaintext job_info
   file in each plotfile and checkpoint directory.

   **Note: unless otherwise specified all quantities in
   MAESTROeX are assumed to be in CGS units.**

   Visualization of results is described in the next section.

Working with the Output
=======================

Visualization and analysis are done on the plotfiles. A number of
in-house and externally developed tools can work with AMReX-formatted
plotfiles [2]_.
An example plot of the reacting_bubble problem run above is
shown in Figure \ `[fig:gettingstarted:test2] <#fig:gettingstarted:test2>`__.

.. raw:: latex

   \centering

.. figure:: \gsfigpath/plt00133_tfromp
   :alt: [fig:gettingstarted:test2] Visualization of the
   final output of the reacting_bubble problem showing the temperature
   field (as derived from the pressure). This plot was done with
   the AmrPostprocessing tools.
   :width: 3in

   [fig:gettingstarted:test2] Visualization of the
   final output of the reacting_bubble problem showing the temperature
   field (as derived from the pressure). This plot was done with
   the AmrPostprocessing tools.

Amrvis
------

Amrvis is an easy-to-use visualization tool developed at LBL for
2- and 3D datasets which can plot slices through 3D datasets as well
as volume-renderings. It can also very easily extract 1D lines
through the dataset along any coordinate direction. It is distributed
separately from the MAESTROeX distribution.

Amrvis can be obtained via git from github as:

::

    git clone https://github.com/AMReX-Codes/Amrvis.git

Also, to build a 3D version of Amrvis you need to obtain volpack using:

::

    git clone https://ccse.lbl.gov/pub/Downloads/volpack.git

Amrvis is built in the C++ AMReX framework (instead of the Fortran
AMReX framework that MAESTROeX uses). The build systems are similar,
but differ in a few ways.

Amrvis uses the Motif library for defining the GUI. On a Linux
system, you may need to install the lesstif package and any
related development packages (e.g. lesstif-devel). Depending
on your Linux system, you may also need to install libXpm and
related development packages (e.g. libXpm-devel).

Further details on the C++ AMReX build system used by Amrvis can be found in the AMReX documentation.

AmrPostprocessing scripts
-------------------------

Several useful analysis scripts (written in Fortran 90) can be found
in amrex/Tools/Postprocessing/F_Src/.
The GNUmakefile there needs to be edited to
indicate which of the tools to build. For example, to extract the
density along a line from the center of a plotfile, plt00200, in
the :math:`y`-direction:

::

    fextract.Linux.gfortran.exe -d 2 -v "density" -p plt00200

These routines are described in § \ `[sec:analysis] <#sec:analysis>`__.

There is also a python visualization method in
AmrPostprocessing/python. This is described
in § \ `[sec:vis:python] <#sec:vis:python>`__.

VisIt
-----

VisIt is a powerful, DOE-supported visualization tool for 2- and 3D
datasets. It can do contouring, volume rendering, streamlines, ...  ,
directly from AMReX plotfiles. Details on
VisIt can be found at:
https://wci.llnl.gov/codes/visit/home.html .
The easiest way to get started with VisIt is to download a precompiled
binary from the VisIt webpage.

Once VisIt is installed, you can open a AMReX plotfile by pointing
VisIt to the Header file in the plotfile directory.

yt
--

yt (version 3.0 and later) can natively read the MAESTROeX plotfiles. See
the yt documentation or § \ `[sec:vis:yt] <#sec:vis:yt>`__.

Diagnostic Files
----------------

By default, MAESTROeX outputs global diagnostics each timestep into a
file called maestro_diag.out. This includes the maximum Mach
number, peak temperature, and peak nuclear energy generation rate.
Individual problems can provide their own diag.f90 file to
produce custom diagnostic output. This information can be plotted
directly with GNUplot, for example.

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

Common Runtime Parameters
=========================

Controlling Timestepping and Output
-----------------------------------

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
-----------------------------------------

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

Development Model
=================

When you clone MAESTROeX from github, you will be on the master
branch of the repo. New changes to MAESTROeX are first introduced
into the development branch in the MAESTROeX git repository.
Nightly regression tests are run on development to ensure that
our answers don’t change. Around the first work day of each month, we
merge from development :math:`\rightarrow` master (assuming
tests pass) and tag the state of the code with a date-based tag
YY-MM. We do this on all the other repos in the AMReX-ecosystem,
including amrex/, Microphysics/, and Castro/.

If you want to contribute to MAESTROeX’s development, issue a pull-request
through github onto the development branch.

Parallel Jobs
=============

To run in parallel with MPI, you would set MPI := t in your
GNUmakefile. For a machine with working MPI compiler wrappers
(mpif90 and mpicc), the build system should find these and
compile with MPI support automatically. This is the easiest way to do
a parallel build, and should work on most Linux systems.

More generally, the build system needs to know about your MPI
installation. For popular national computing facilities, this is
already setup, and the build system looks at the machine hostname to
set the proper libraries. For other machines, you may need to edit
the GMake.MPI file in the AMReX build files. See
§ \ `[ch:make] <#ch:make>`__ for more details.

OpenMP can be used to parallelize on shared-memory machines (i.e. within a node). OpenMP support is accomplished through the compiler.
Setting OMP := t in the GNUmakefile will enable the proper
compiler flags to build with OpenMP. Note: not all MAESTROeX modules
have OpenMP support. Microphysics routines need to be written in a
threadsafe manner. This can be tested via the test_react unit
test (see  § `[chapter:unit_tests] <#chapter:unit_tests>`__).

.. [1]
   In earlier versions of MAESTROeX this
   problem was called test2

.. [2]
   The plotfiles are in the same format as those made
   by the BoxLib library upon which MAESTROeX was previously based.
