.. _ch:make:

**********************
MAESTROeX Build System
**********************

Build Process Overview
======================

The MAESTRO build system uses features of GNU make (version 3.81 or
later), which is typically the default on systems. The MAESTRO executable is built in the problem’s directory (one of the directories
under SCIENCE/, TEST_PROBLEMS, or UNIT_TESTS). This
directory will contain a makefile, GNUmakefile, that includes
all the necessary information to build the executable.

The main macros that define the build process are split across several
files. The 4 main files are:

-  ${AMREX_HOME}/Tools/F_mk/GMakedefs.mak:

   This setups the basic macros, includes the options for the selected
   compiler, builds the list of object and source files, and defines
   the compile and link command lines.

-  ${AMREX_HOME}/Tools/F_mk/GMakeMPI.mak:

   This implements any changes to the compiler names and library
   locations necessary to build a parallel executable with MPI.

-  ${AMREX_HOME}/Tools/F_mk/GMakerules.mak:

   This creates the various build targets and specifies the rules for
   building the object files, the list of dependencies, and some other
   lesser-used targets (tags for editors, documentation, etc.)

-  MAESTRO/GMaestro.mak:

   This is a MAESTRO-specific file that gathers all of the various
   modules that are used to build a typical MAESTRO application
   and integrates with the AMReX build system. Every MAESTRO problem’s GNUmakefile will include this file.

MAESTRO gets the location of the AMReX library through the
AMREX_HOME variable. This should be set as an environment
variable in your shell start-up files (e.g. .bashrc or
.cshrc).

The AMReX build system separates the compiler-specific information
from the machine-specific information—this allows for reuse of the
compiler information. The only machine-specific parts of the build system
are for the MPI library locations, contained in GMakeMPI.mak.
The compiler flags for the various compilers are listed in the
files in ${AMREX_HOME}/Tools/F_mk/comps/. The compiler
is set via the COMP variable in the problem’s GNUmakefile.

There are several options in addition to the compiler that affect the
build process: MPI, OMP, and NDEBUG—these turn on/off
MPI parallelization, OpenMP parallelization, and debugging. Together,
these choices along with the compiler name are reflected in the name
of the executable.

When the ‘make’ command is run, the object and module files are
written into a directory t/\ *OS*.\ *COMP*.\ *other*/,
where *OS* is the operating system detected by the build
system, *COMP* is the compiler used, and *other*
reflects any other build options (MPI, OpenMP, debugging, etc.) used.
Separating each build into a separate subdirectory under the problem
directory allows for multiple builds of MAESTRO to exist
side-by-side in the problem directory.

Finding Source Files
--------------------

The MAESTRO executable is built from source distributed across a
number of directories. In each of these directories containing source
files, there is a GPackage.mak file. This file has a number of
lines of the form:

::

    f90sources += file.f90

where file.f90 is a source file that should be built when this
directory is added to the list of build directories. For old
fixed-form Fortran files, the files should be added to the
fsources variable instead of f90sources.

The AMReX build system relies on the vpath functionality of
make. In a makefile, the vpath variable holds search path used
to locate source files. When a particular file is needed, the
directories listed in vpath are searched until the file is
found. The first instance of a file is used. We exploit this feature
by always putting the build directory first in vpath (this is
done in GMakerules.mak). This means that if a source file is
placed in the build directory, that copy will override any other
version in the source path.

In MAESTRO, the vpath variable is set using the macros defined
in GMaestro.mak. A user does not need to set this variable
explicitly. Additional source locations are added in the manner
described below (see § \ `2.1 <#sec:make:otherfiles>`__).

Dependencies
------------

There is no need to explicitly define the dependencies between the
source files for Fortran modules. The scripts in
AMREX_HOME/Tools/F_scripts/ are run at the start of the build
process and parse all the source files and make an explicit list of
the dependency pairs. The execution of these scripts is triggered
by including makefiles of the form .depends. On a fresh build these
will not exist. When GNU make cannot find an included makefile it will
first attempt to build it using any relevant targets before issuing an
error. The targets for the .depends files contain the recipe for
executing the dependency scripts. Once these makefiles are built by the
scripts GNU make will then read the dependencies for the current build
from them. This process is defined in GMakerules.mak.

A few files use explicit ‘include’ statements to include Fortran
source in other source files. Any include file locations should be
added to Fmincludes variable in the problem’s GNUmakefile.
This does not occur frequently. For the case of the helmholtz
equation of state, this is done automatically in GMaestro.mak.

Files Created at Compile-time
-----------------------------

Several files are created at build-time:

-  probin.f90:

   This is the module that controls runtime parameters. This is
   created by the script
   write_probin.py in ${AMREX_HOME}/Tools/F_scripts/. The
   makefile logic for creating it is in GMaestro.mak. At compile
   time, the problem, main MAESTRO/, and any microphysics
   directories (set from the EOS_DIR, CONDUCTIVITY_DIR, and NETWORK_DIR parameters in the GNUmakefile
   are searched for \_parameter files. These files
   are parsed and the probin.f90 file is output containing the
   runtime parameters, their default values, and the logic for reading
   and interpreting the inputs file.

-  build_info.f90:

   This is a module that contains basic information about the build
   environment (compiler used, build machine, build directory, compiler
   flags, etc.). This is created by the script makebuildinfo.py
   in ${AMREX_HOME}/Tools/F_scripts/ from
   GMaestro.mak by passing in a number of makefile variables. This is
   rewritten everytime the code is compiled. The primary use of this
   module is writing out the job_info file in the plotfiles.

-  (network.f90):

   This is generated at compile time *only* for the
   general_null network. The general_null network allows the
   user to define a list of non-reacting species builds the
   network.f90 based on this list. The makefile logic for building
   the network.f90 is in the GPackage.mak in
   Microphysics/networks/general_null. The script write_network.py
   in that directory does the actual parsing of the species file and
   outputs the network.f90.

MAESTRO Problem Options
=======================

.. _sec:make:otherfiles:

Problem-specific Files
----------------------

If a problem has a unique file that is needed as part of the build,
then that file should be added to a GPackage.mak file in the
problem’s directory. Since, by default, problems don’t have a
GPackage.mak, the build system needs to be told to look in the
problem directory for unique sources. This is accomplished by adding
the problem’s directory to the EXTRA_DIR variable in the
problem’s GNUmakefile.

Note that this is not necessary if you place a custom version of
a source file in the problem’s directory. Since that file is already
listed in the GPackage.mak in its original location, the build
system will know that it needs to be built. Since the vpath
variable puts the problem’s directory at the start of the search
path, the version of the file in the problem’s directory will be
found first.

Defining EOS, Network, and Conductivity Routines
------------------------------------------------

Each MAESTRO problem needs to define an equation of state, a
reaction network, and a routine to compute the conductivities (for
thermal diffusion). This is true even if the problem is not doing
reactions of thermal diffusion. These definitions are specified
in the problem’s GNUmakefile.

-  EOS_DIR:

   This variable points to the directory (by default, relative to
   Microphysics/EOS/) of the equation of state used for the build.
   Choices that work with MAESTRO are:

   -  helmholtz

   -  gamma_law_general

   -  multigamma

   To use an EOS contained in a different location, set the variable
   EOS_TOP_DIR to point to the directory above the alternate EOS
   directory.

-  CONDUCTIVITY_DIR:

   This variable points to the conductivity routine used for the build
   (by default, relative to Microphysics/conductivity/). Choices
   that work with MAESTRO are:

   -  constant

   -  timmes_stellar

   If diffusion is not being used for the problem, this should be set
   to constant. To use an alternate conductivity
   routine, set the variable CONDUCTIVITY_TOP_DIR to point
   to the directory above the alternate conductivity directory.

-  NETWORK_DIR:

   This variable points to the reaction network used for the build (by
   default, relative to Microphysics/networks/). Several options
   are present in Microphysics/networks/. A network is required even
   if you are not doing reactions, since the network defines the
   species that are advected and interpreted by the equation of state.

   A special choice, Microphysics/networks/general_null is a general
   network that simply defines the properties of one or more species.
   This requires an inputs file, specified by
   GENERAL_NET_INPUTS. This inputs file is read at compile-time and
   used to build the network.f90 file that is compiled into the
   source.

   To use an alternate reaction network, set the variable
   NETWORK_TOP_DIR to point to the directory above the alternate
   network.

Core MAESTRO modules
--------------------

Several modules are included in all MAESTRO builds by default.
From AMReX, we alway include:

-  ${AMREX_HOME}/Src/F_BaseLib

-  ${AMREX_HOME}/Src/LinearSolvers/F_MG

From Util, we always include

-  Util/model_parser

-  Util/random

The microphysics, as described above is also included. For the
networks, we include a file called NETWORK_REQUIRES into
GMaestro.mak that tells us whether to also include Util/VODE
(if NEED_VODE := t). It is assumed in this case that we need
BLAS and LINPACK, so these are compiled in from Util/BLAS and
Util/LINPACK.

You can instead link in a system-wide optimized BLAS library by setting
SYSTEM_BLAS := tGNUmakefile:BLAS in the GNUmakefile. This will
use the library specified in GNUmakefile:BLAS_LIBRARY BLAS_LIBRARY
(defaulting to -lopenblas)
to the link line, and assumes that the library is in your path. Note
that for some systems, you should have the static BLAS libraries
available.

From MAESTRO/, we add

-  MAESTRO/constants

-  MAESTRO/Source

(although see the unit tests section below regarding MAESTRO/Source.)

Finally, any extra directories listed in the EXTRA_DIR
variable are included.

For each of these included directories, GMaestro.mak adds the
list of source files defined in their GPackage.mak to the list
of files to be compiled. It also adds each of these directories to
the vpath as a directory for the build process to search in for
source files.

Unit Tests
----------

Sometimes we only want to use a few of the standard MAESTRO routines, for example in a unit test where we are testing only a small
part of the MAESTRO algorithm indepenedently. In this case, we
don’t want to comple all of the files in MAESTRO/Source. If we
set UNIT_TEST := t in our problem’s GNUmakefile, then the
GPackage.mak in MAESTRO/Source is not read, so those files
are not automatically put into the list of files to compile. Instead,
the problem should create its own GPackage.mak listing only the
subset of files that are to be compiled. MAESTRO/Source is put
into the vpath search path for sources, so those files will
still be found as needed.

AMReX-only Tests
----------------

An even more restrictive setting than UNIT_TEST := t is invoked
by setting AMREX_ONLY := t. This is like the unit test flag,
but does not include MAESTRO/Source in the vpath search
path for sources. So this is intended for cases where we don’t want
to use any MAESTRO source files. Typically, this is used in the
small unit tests that live under the various microphysics solvers. If
a probin.f90 is built for these tests, it will not include all
the MAESTRO-specific parameters, but will include any parameters from
the various microphysics routines.

Special Targets
===============

Debugging
---------

(print-\*)
~~~~~~~~~~

To see the contents of any variable in the build system, you can build
the special target print-\ *varname*, where is the name of the variable. For example, to see what the
Fortran compiler flags are, you would do:

::

    make print-FFLAGS

This would give (for gfortran, for example):

::

    FFLAGS is -Jt/Linux.gfortran/m -I t/Linux.gfortran/m -O2 -fno-range-check

This functionality is useful for debugging the makefiles.

file_locations
~~~~~~~~~~~~~~

Source files are found by searching through the make
vpath. The first instance of the file found in the vpath
is used in the build. To see which files are used and their locations,
do:

::

    make file_locations

This will also show any files that aren’t found. Some are expected
(e.g., build_info.f90 and probin.f90 are created at
compile time), but other files that are not found could indicate
an incomplete vpath.

clean and realclean
-------------------

Typing ‘make clean’ deleted the object and module files for the
current build (i.e., the current choice of MPI, NDEBUG,
COMP, and OMP). This also removes any of the compile-time
generated source files. Any other builds are left unchanged.

Typing ‘make realclean’ deletes the object and module files for
all builds—i.e., the entire t/ directory is removed.

.. _ch:makefiles:special:

Special Debugging Modes
=======================

AMReX has several options that produce executables that can help
track down memory issues, uninitialized variables, NaNs, etc.

-  NDEBUG

   GNUmakefile:NDEBUG To generate an executable
   with debugging information included in the executable (e.g., to be
   interpreted by the debugger, gdb), compile with NDEBUG
   := . This will usually add -g to the compile line and
   also lower the optimization. For gfortran it will add several
   options to catch uninitialize variables, bounds errors, etc.

-  TEST

   Setting TEST := tGNUmakefile:TEST will
   enable routines in AMReX initialize multifabs and arrays
   allowed via bl_allocate to signalliing NaNs. This behavior
   is the same as NDEBUG :=, but TEST := t uses the
   same compiler optimizations as a normal build.

   This can be useful with compiler flags that trap floating point
   exceptions (FPEs), but checks on floating point exceptions can also
   be enabled through runtime parameters passed to AMReX’s
   backtrace functionlity:

   -  boxlib_fpe_invalid: enabling FPE trapping for
      invalid operations (e.g. 0 \* inf, sqrt(-1))

   -  boxlib_fpe_zero: enable FPE trapping for
      divide-by-zero

   -  boxlib_fpe_overflow: enable FPE trapping for
      overflow

-  backtracing

   When exception trapping is enabled (either via AMReX or the
   compiler), the code will abort, and the backtrace information will
   be output to a file Backtrace.N, where N is the
   processor number. AMReX will also initialize multifabs with
   signaliing NaNs to help uncover any floating point issues.

   This is also useful to diagnose deadlocks in parallel regions.
   If the code is hanging, doing “control-C” will be intercepted
   and the code will generate a backtrace which will identify
   where in the code there was a deadlock.

   Behind the scenes, AMReX implements this capability via the
   Linux/Unix feenableexcept function (this is in
   backtrace_c.cpp in AMReX).

-  FSANITIZER

   For gfortran, gcc, g++, setting FSANITIZER :=
   tGNUmakefile:FSANITIZER will enable the
   address sanitizer support built into GCC. This is enabled through
   integration with https://github.com/google/sanitizers in GCC.

   Note: you will need to have the libraries libasan and
   libubsan installed on your machine to use this functionality.

Extending the Build System
==========================

Adding a Compiler
-----------------

Properties for different compilers are already defined in
${AMREX_HOME}/Tools/F_mk/comps/. Each compiler is given its
own file. The appropriate file is included into GMakedefs.mak
by looking at the COMP variable and the operating system. These
compiler files define the compiler flags for both optimized and debug
compiling. Additionally, the variable FCOMP_VERSION should be
defined there, based on the output from the compiler, to provide the
compiler version for output into the job_info file at runtime.

Parallel (MPI) Builds
---------------------

When building with MPI, the build system needs to know about the
location of the MPI libraries. If your local MPI has the mpif90
and mpicc wrappers installed and working, then MAESTRO will
attempt to use these. Otherwise, you will need to edit
GMakeMPI.mak and add a section specific to your machine with the
compiler and library location. It is best to simply copy an existing
similar portion of the makefile and adjust it to your system. Most
national supercomputing facilities are already supported, and parallel
builds on them should work out of the box.
