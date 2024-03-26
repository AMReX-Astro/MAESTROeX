.. _sec:adding_problems:

**********************
Adding A New Problem
**********************

Introduction to Adding a New Problem
====================================

Different MAESTROeX problems are defined in sub-directories under
``Exec/`` in ``SCIENCE/``, ``TEST_PROBLEMS/``, or ``UNIT_TESTS/``.
To add a problem, start by creating a new sub-directory—this is
where you will compile your problem and store all the problem-specific
files.

The minimum requirement to define a new problem would be a
GNUmakefile which describes how to build the application and an
input file which lists the runtime parameters. The problem-specific
executable is built in the problem directory by typing ``make``.
Source files are found automatically by searching the directories
listed in the GNUmakefile. Customized versions of any source
files placed in the problem-directory override those with the same
name found elsewhere. Any unique source files (and not simply a
custom version of a file found elsewhere) needs to be listed in a file
called GPackage.mak in the problem-directory (and this needs to
be told to the build system—see below).

.. _sec:makefile:

The GNUmakefile
===============

A basic GNUmakefile begins with:

::

      DEBUG   := FALSE
      USE_MPI :=
      USE_OMP :=

Here, ``DEBUG`` is false if we are building an optimized executable.
Otherwise, the debug version is built—this typically uses less
optimization and adds various runtime checks through compiler flags.
``USE_MPI`` and ``USE_OMP`` are set to true if we want to use either MPI
or OpenMP for parallelization. If ``USE_MPI := TRUE``, you will need to
have the MPI libraries installed.

The next line sets the compiler to be used for compilation:

::

      COMP := gnu

The MAESTROeX build system knows what options to use for various
compiler families. The ``COMP`` flag specifies which compiler to
use. Allowed values include ``intel``, ``gnu``, ``pgi``, and ``cray``.

``USE_REACT`` set to true will turn on the reactions solvers.

::

      USE_REACT := TRUE

The next line defines where the top of the MAESTROeX source tree is located.

::

      MAESTROEX_HOME := ../../..

An MAESTROeX application is built from several packages (the
multigrid solver, an EOS, a reaction network, etc.). The core
MAESTROeX packages are always included, so a problem only needs
to define the EOS, reaction network, and conductivities to
use, as well as any extra, problem-specific files.

::

    EOS_DIR          := helmholtz
    CONDUCTIVITY_DIR := constant
    NETWORK_DIR      := ignition_simple

    Bpack := ./Make.package
    Blocs := .

Note that the microphysics packages are listed simply by the name of
the directory containing the specific implementation (e.g. ``helmholtz``).
By default, the build system will look in ``Microphysics/EOS/`` for
the EOS, ``Microphysics/conductivity/`` for the conductivity routine,
and ``Microphysics/networks/`` for the reaction network. To
override this default search path, you can set ``EOS_TOP_DIR``,
``CONDUCTIVITY_TOP_DIR``, and ``NETWORK_TOP_DIR`` respectively.

Generally, one does not need to include the problem directory itself
in ``Bpack`` and ``Blocs``, unless there are unique source files found there,
described in a ``Make.package`` file. These variables are
interpreted by the ``Make.Maestro`` file and used to build a main
list of packages called Bdirs. The build system will attempt
to build all of the files listed in the various ``Make.package``
files found in the Bdirs directories. Furthermore,
Bdirs will be will be added to the make VPATH, which
is the list of directories to search for source files. The problem
directory will always be put first in the VPATH, so any source
files placed there override those with the same name found elsewhere
in the source tree.

By default, some of the runtime parameters are listed in
``_parameters`` file in each problem directory.

::

    PROBIN_PARAMETER_DIRS := .

They are parsed at compile time and the file ``extern.F90``
is written and compiled. This is a Fortran module that holds the values of
the runtime parameters and makes them available to any routine via
``probin_module``.

The final line in the GNUmakefile includes the rules to actually
build the executable.

::

      include $(MAESTROEX_HOME)/Exec/Make.Maestro

Handling Problem-Specific Source Files
--------------------------------------

As mentioned above, any source files placed in the problem directory
override a files with the same name found elsewhere in the source
tree. This allows you to create a problem-specific version of any
routine. Source files that are unique to this problem (i.e. there is
no file with the same name elsewhere in the source tree) need to be
listed in a file ``Make.package`` in the problem directory, and
the problem-directory needs to be explicitly listed in the ``Bpack``
and ``Blocs`` lists in the GNUmakefile.

.. _sec:def_runtime_param:

Defining Runtime Parameters
===========================

The runtime parameters for the core MAESTROeX algorithm are listed in
``_cpp_parameters`` file in ``MAESTROeX/Source/param/`` folder.
That file is parsed at compile-time by the ``parse_maestro_params.py``
script (along with any problem-specific parameters) in the same folder.
The script outputs the ``extern.F90``source file.
Each line in the ``_cpp_parameters`` file has the form::

  *parameter*    *data-type*    *value*    *need in Fortran?*

where *parameter* is the name of the runtime parameter,
*data-type* is one of {string, Real, int, bool},
the *value* specifies the default value for the runtime parameter,
and *need in Fortran?* is marked *y* only if that parameter is
used in Fortran. Comments are indicated by a ‘#’ character and are
used to produce documentation about the available runtime parameters.
For the documentation, runtime parameters are grouped together
in the ``_cpp_parameters`` file into categories. The category headings
are defined by comments in the ``_cpp_parameters`` file and any comments
following that heading are placed into that category.

At runtime, the default values for the parameters can be overridden
either through the inputs file (by adding a line of the form:
``parameter = value``) or through a command-line argument (taking the
form: ``-parameter value``). The probin_module makes the
values of the runtime parameters available to the various functions
in the code (see § :ref:`sec:runtime_parameters`).

Problem-specific runtime parameters should be defined in the
problem-directory in a file called ``_parameters``. This file will
be automatically found at compile time.

.. _sec:initial_models:

Preparing the Initial Model
===========================

MAESTROeX models subsonic, non-hydrostatic flows as deviations from
a background state in hydrostatic equilibrium.
The solution in MAESTROeX is broken up into a 1D base state and the 2-
or 3D full state. The job of the 1D base state in the algorithm is
to represent the hydrostatic structure. The full, Cartesian state
carries the departures from hydrostatic equilibrium. The underlying
formulation of the low Mach number equations assumes that the base
state is in hydrostatic equilibrium. At the start of a simulation,
the initial model is read in and taken as the base state. Therefore,
any initial model needs to already be in hydrostatic equilibrium.

An initial model is generated by ``initial_models`` found in
the same github repo as MAESTROeX. You can obtain it via::

    git clone https://github.com/AMReX-Astro/initial_models

In general, there are two different proceduces that are
needed. The first type modify an existing 1D initial model produced
somewhere else (e.g. a 1D stellar evolution code), and map it onto a
uniform grid, at the desired resolution, using the desired
equation of state and discretization of hydrostatic
equilibrium. The second type generate the initial model internally,
by integrating the condition of hydrostatic equilibrium together with
a simplifying assumption on the energy (e.g. isothermal or
isentropic). In both cases hydrostatic equilibrium is enforced as:

.. math::

   \frac{p_{i+1} - p_i}{\Delta r} = \frac{1}{2} (\rho_i + \rho_{i+1})
   g_{i+1/2}

Here, :math:`g_{i+1/2}` is the edge-centered gravitational acceleration.

Full details on which initial model routine matches each problem and
how the initial models are used to initialize the full state data can
be found in § :ref:`sec:initial_models_main`.

Customizing the Initialization
==============================

The best way to customize the initialization (e.g. add perturbations)
is to copy from one of the existing problems.
The file ``initdata.f90`` controls initialization of scalars
(:math:`\rho`, :math:`\rho X_k`, :math:`\rho h`) and velocity field.
The ``reacting_bubble`` problem is a good
starting point for plane-parallel and ``wdconvect`` is a good
starting point for full spherical stars.
