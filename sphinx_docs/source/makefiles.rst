.. _ch:make:

**********************
MAESTROeX Build System
**********************

Build Process Overview
======================

The MAESTROeX build system uses features of GNU make (version 3.82 or
later), which is typically the default on systems. The MAESTROeX
executable is built in the problem’s directory (one of the directories
under ``SCIENCE/``, ``TEST_PROBLEMS/``, or ``UNIT_TESTS/``). This
directory will contain a makefile, ``GNUmakefile``, that includes all the
necessary information to build the executable.

We use the build system in AMReX in ``amrex/Tools/GNUmake``.  A set of
MAESTROeX-specific macros are defined in ``MAESTROeX/Exec/Make.Maestro``.

MAESTROeX gets the location of the AMReX library through the
``AMREX_HOME`` variable. This should be set as an environment
variable in your shell start-up files (e.g. ``.bashrc`` or
``.cshrc``).

The AMReX build system separates the compiler-specific information
from the machine-specific information—this allows for reuse of the
compiler information. The only machine-specific parts of the build
system are for the MPI library locations, contained in
``amrex/Tools/GNUmake/sites/``.  The compiler flags for the various
compilers are listed in the files in
``amrex/Tools/GNUmake/comps/``. The compiler is set via the ``COMP``
variable in the problem’s GNUmakefile.

There are several options in addition to the compiler that affect the
build process: ``USE_MPI``, ``USE_OMP``, ``USE_CUDA`` for the parallelism,
and ``DEBUG`` and ``TEST`` for diagnostic information.. Together,
these choices along with the compiler name are reflected in the name
of the executable.

When the ``make`` command is run, the object and module files are
written into sub-directories under ``tmp_build_dir/``. Separating each build into
a separate sub-directory under the problem directory allows for
multiple builds of MAESTROeX to exist side-by-side in the problem
directory.

Finding Source Files
--------------------

The MAESTROeX executable is built from source distributed across a
number of directories. In each of these directories containing source
files, there is a ``Make.package`` file. This file has a number of
lines of the form:

::

    CEXE_sources += file.cpp
    F90EXE_sources += file.F90

where ``file.cpp`` is a C++ source file and ``file.F90`` is a Fortran
source file that should be built when this directory is added to the
list of build directories.

The AMReX build system relies on the ``vpath`` functionality of
make. In a makefile, the ``vpath`` variable holds search path used to
locate source files. When a particular file is needed, the directories
listed in ``vpath`` are searched until the file is found. The first
instance of a file is used. We exploit this feature by always putting
the build directory first in ``vpath``.  This means that if a source
file is placed in the build directory, that copy will override any
other version in the source path.


Dependencies
------------

There is no need to explicitly define the dependencies between the
source files for Fortran modules. Scripts in
AMReX are run at the start of the build
process and parse all the source files and make an explicit list of
the dependency pairs.

Files Created at Compile-time
-----------------------------

Several files are created at build-time:

-  ``extern.F90``:

   This is the module that controls runtime parameters. This is
   created by the script ``write_probin.py`` in
   ``amrex/Tools/F_scripts/``. The makefile logic for creating it is
   in ``Make.Maestro``. At compile time, the problem, main MAESTROeX/,
   and any microphysics directories (set from the ``EOS_DIR``,
   ``CONDUCTIVITY_DIR``, and ``NETWORK_DIR`` parameters in the ``GNUmakefile`` are
   searched for ``_parameter`` files. These files are parsed and the
   ``extern.F90`` file is output containing the runtime parameters, their
   default values, and the logic for reading and interpreting the
   inputs file.

-  AMReX_buildInfo.cpp:

   This file contains basic information about the build
   environment (compiler used, build machine, build directory, compiler
   flags, etc.). This is created by the script ``makebuildinfo_C.py``
   in ``amrex/Tools/C_scripts/`` from
   ``Make.Maestro`` by passing in a number of makefile variables. This is
   rewritten every time the code is compiled. The primary use of this
   module is writing out the ``job_info`` file in the plotfiles.

-  ``actual_network.F90``:

   This is generated at compile time *only* for the
   ``general_null`` network. The ``general_null`` network allows the
   user to define a list of non-reacting species builds the
   ``actual_network.F90`` based on this list. The makefile logic for building
   this is in the ``Make.package`` in
   ``Microphysics/networks/general_null``. The script ``write_network.py``
   in that directory does the actual parsing of the species file and
   outputs the ``actual_network.f90``.

MAESTROex Problem Options
=========================

.. _sec:make:otherfiles:

Problem-specific Files
----------------------

If a problem has a unique file that is needed as part of the build,
then that file should be added to a ``Make.package`` file in the
problem’s directory.

Note that this is not necessary if you place a custom version of
a source file in the problem’s directory. Since that file is already
listed in the ``Make.package`` in its original location, the build
system will know that it needs to be built. Since the ``vpath``
variable puts the problem’s directory at the start of the search
path, the version of the file in the problem’s directory will be
found first.

Defining EOS, Network, and Conductivity Routines
------------------------------------------------

Each MAESTROeX problem needs to define an equation of state, a
reaction network, and a routine to compute the conductivities (for
thermal diffusion). This is true even if the problem is not doing
reactions of thermal diffusion. These definitions are specified
in the problem’s ``GNUmakefile``.

-  ``EOS_DIR``:

   This variable points to the directory (by default, relative to
   ``Microphysics/EOS/``) of the equation of state used for the build.
   Choices that work with MAESTROeX are:

   -  ``helmholtz``

   -  ``gamma_law``

   -  ``multigamma``


-  ``CONDUCTIVITY_DIR``:

   This variable points to the conductivity routine used for the build
   (by default, relative to ``Microphysics/conductivity/``). Choices
   that work with MAESTROeX are:

   -  ``constant``

   -  ``stellar``

   If diffusion is not being used for the problem, this should be set
   to constant.

-  ``NETWORK_DIR``:

   This variable points to the reaction network used for the build (by
   default, relative to ``Microphysics/networks/``). Several options
   are present in ``Microphysics/networks/``. A network is required even
   if you are not doing reactions, since the network defines the
   species that are advected and interpreted by the equation of state.

   A special choice, ``Microphysics/networks/general_null`` is a general
   network that simply defines the properties of one or more species.
   This requires an inputs file, specified by
   ``GENERAL_NET_INPUTS``. This inputs file is read at compile-time and
   used to build the ``actual_network.F90`` file that is compiled into the
   source.


Core MAESTROeX modules
----------------------

Several modules are included in all MAESTROeX builds by default.
In addition to the AMReX sources, we also include

-  ``MAESTROeX/Source``

-  ``Util/model_parser``

The microphysics used may bring in its own dependencies.

For each of these included directories, ``Make.Maestro`` adds the
list of source files defined in their ``Make.package`` to the list
of files to be compiled. It also adds each of these directories to
the ``vpath`` as a directory for the build process to search in for
source files.

Special Targets
===============


``print-*``
-----------

To see the contents of any variable in the build system, you can build
the special target ``print-varname``, where ``varname`` is the name of
the variable. For example, to see what the network directory is, you would
do::

    make print-NETWORK_DIR

This functionality is useful for debugging the makefiles.

``file_locations``
------------------

Source files are found by searching through the make
``vpath``. The first instance of the file found in the ``vpath``
is used in the build. To see which files are used and their locations,
do:

::

    make file_locations

This will also show any files that aren’t found. Some are expected
(e.g., ``extern.F90`` is created at compile time), but other files
that are not found could indicate an incomplete ``vpath``.

``clean`` and ``realclean``
---------------------------

Typing ``make clean`` deletes the object and module files for the
current build (i.e., the current choice of ``USE_MPI``, ``DEBUG``,
``COMP``, and ``USE_OMP``). This also removes any of the compile-time
generated source files. Any other builds are left unchanged.

Typing ``make realclean`` deletes the object and module files for
all builds—i.e., the entire ``tmp_build_dir/`` is removed.

.. _ch:makefiles:special:

Special Debugging Modes
=======================

AMReX has several options that produce executables that can help
track down memory issues, uninitialized variables, NaNs, etc.

-  ``DEBUG = TRUE`` :

   Setting ``DEBUG=TRUE`` on the ``make`` commandline or in the
   ``GNUmakefile`` generates an executable with debugging information
   included in the executable (e.g., to be interpreted by the
   debugger, gdb). This will usually add -g to the compile line and
   also lower the optimization. For gfortran it will add several
   options to catch uninitialize variables, bounds errors, etc.
   The resulting executable will have ``DEBUG`` in its name.

-  ``TEST = TRUE``

   Setting ``TEST=TRUE`` on the ``make`` commandline or in the
   ``GNUmakefile`` will enable routines in AMReX to initialize MultiFabs
   and arrays with signalliing NaNs. This
   behavior is the same as ``DEBUG=TRUE``, but ``TEST`` uses the same
   compiler optimizations as a normal build.

   This can be useful with compiler flags that trap floating point
   exceptions (FPEs), but checks on floating point exceptions can also
   be enabled through runtime parameters passed to AMReX’s
   backtrace functionlity:

   -  ``amrex.fpe_trap_invalid``: enabling FPE trapping for
      invalid operations (e.g. ``0 * inf``, ``sqrt(-1)``)

   -  ``amrex.fpe_trap_zero``: enable FPE trapping for
      divide-by-zero

   -  ``amrex.fpe_trap_overflow``: enable FPE trapping for
      overflow

-  backtracing

   When exception trapping is enabled (either via AMReX or the
   compiler), the code will abort, and the backtrace information will
   be output to a file ``Backtrace.N``, where N is the
   processor number. AMReX will also initialize multifabs with
   signaliing NaNs to help uncover any floating point issues.

   This is also useful to diagnose deadlocks in parallel regions.
   If the code is hanging, doing “control-C” will be intercepted
   and the code will generate a backtrace which will identify
   where in the code there was a deadlock.

   Behind the scenes, AMReX implements this capability via the
   Linux/Unix ``feenableexcept`` function (this is in
   ``backtrace_c.cpp`` in AMReX).

-  ``FSANITIZER``

   For gfortran, gcc, g++, setting ``FSANITIZER=TRUE``
   in ``GNUmakefile`` will enable the
   address sanitizer support built into GCC. This is enabled through
   integration with https://github.com/google/sanitizers in GCC.

   Note: you will need to have the libraries libasan and
   libubsan installed on your machine to use this functionality.

Extending the Build System
==========================

Adding a Compiler
-----------------

Properties for different compilers are already defined in
``${AMREX_HOME}/Tools/GNUmake``. Each compiler is given its
own file in the ``comps/`` sub-directory.  These
compiler files define the compiler flags for both optimized and debug
compiling.

Parallel (MPI) Builds
---------------------

When building with MPI, the build system needs to know about the
location of the MPI libraries. If your local MPI has the ``mpif90``
and ``mpicxx`` wrappers installed and working, then MAESTROeX will
attempt to use these. Otherwise, you will need to add a site to the
``sites/`` sub-directory in the AMReX build system specifying the
details of your environment.
