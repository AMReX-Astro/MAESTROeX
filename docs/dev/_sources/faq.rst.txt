**************************
Frequently Asked Questions
**************************

Compiling
=========

#. *What version of the Fortran standard is needed?*

   Some features of Fortran 2003 are used, in particular, the
   ISO_C_BINDING feature of Fortran 2003 is needed to define a long
   integer type for some MPI operations in Fortran. This is supported
   by most Fortran compilers even if they don’t support the entire
   Fortran 2003 standard.

   We also rely on the Fortran 95 standard that specifies that any
   local allocated arrays are automatically deallocated when a
   subroutine ends. Fortran 90 does not do this. Most
   MAESTROeX routines rely on this Fortran 95 feature.


#. *The code doesn’t compile, but complains right away that there
   is “No rule to make target ‘AMReX_constants_mod.o’, needed by ‘tmp_build_dir/d/2d.gnu.MPI/f90.depends’”*

   The environment variable ``AMREX_HOME`` needs to be the full path
   to the ``amrex/`` directory. You cannot use ‘:math:`\sim`’ as a shortcut
   for your home directory.


#. *make issues an error like:*

   ::

      $ make
      Loading /home/zingale/development/amrex//Tools/GNUMake/comps/gnu.mak...
      Loading /home/zingale/development/amrex//Tools/GNUMake/sites/Make.unknown...
      /home/zingale/development/amrex//Tools/GNUMake/Make.rules:476: tmp_build_dir/d/3d.gnu.MPI.EXE/f90.depends: No such file or directory
      make: *** No rule to make target `meth_params.F90', needed by `tmp_build_dir/d/3d.gnu.MPI.EXE/f90.depends'.  Stop.

   You need to use GNU make version 3.82 or later.

Running
=======

#. *How do we turn off all the initial projections to look at the
   initial velocity field as specified in initdata, instead of as
   modified by the velocity constraint?*

   ::

           maestro.max_step  = 1
           maestro.init_iter = 0
           maestro.init_divu_iter = 0
           maestro.do_initial_projection = false

#. *MAESTROeX crashes because the multigrid algorithm fails to
   converge—how do I get around this?*

   Setting general convergence criteria for multigrid is as much
   art as science.
   First, it is important to determine if the multigrd solver is
   close to convergence and just dancing around near the desired
   tolerance, but never reaching it, or if it is no where near
   convergence. For the latter, it may be that the multigrid
   solver was fed bad data and the problem arose in one of the earlier
   steps. To get more detail information from the multigrid solver,
   set mg_verbose to a positive integer from 1-4 (the higher
   the number the more information you receive.

   If the multigrid solver is failing during one of the initial
   “divu” iterations, it may be because the velocity is initially
   zero, so there is no velocity magnitude to use as a reference for
   convergence, and that (:math:`S - \bar{S}`) is very small (or zero). In
   this case, it is usually a good idea to perturb the initial state
   slightly, so the righthand side is non-zero.

   The tolerances used for the various multigrid solves in the code
   can be overridden on a problem-by-problem basis by setting the relevant
   parameters in the problem's inputs file (see the `solver tolerances` section in the § :ref:`runtime parameters tables<sec:runtime-parameters-tables>`).


#. *Why do the initial projection and “divu” iters sometimes
   have a harder time converging than the multigrid solves in the main algorithm?*

   The initial projection and “divu” solve sets the density to :math:`1`
   (see § :ref:`sec:Initialization`), so the coefficients in the
   elliptic solve are :math:`O(\beta_0) \sim O(\rho)`. But in the main
   algorithm, the coefficients are :math:`O(\beta_0/\rho) \sim O(1)`. Since
   :math:`\rho` can vary a lot, the variation in the coefficients in the
   initial projection and “divu” solve present a harded linear system
   to solve.


#. *How can I obtain profiling information for my run?*

   The code is already instrumented with timers. Simply compile with
   ``TINY_PROFILE=TRUE`` in the ``GNUmakefile``, or equivalently do
   ``make TINY_PROFILE=TRUE``. A summary of the timings will
   be output to ``stdout`` at the end of the run.

   With the GNU compilers, you can enabling profiling with ``gprof``
   by compiling with

   ::

         USE_GPROF=TRUE

   in your GNUmakefile.

   When you run, a file named ``gmon.out`` will be produced. This can
   be processed with ``gprof`` by running:

   ::

         gprof exec-name

   where ``exec-name`` is the name of the executable. More detailed
   line-by-line information can be obtained by passing the ``-l``
   argument to ``gprof``.


#. *How can I force MAESTROeX to output?*

   To generate a checkpoint file, in the output directory do:

   .. prompt:: bash

      touch dump_and_continue

   For a plotfile:

   .. prompt:: bash

      touch plot_and_continue

   or a small plotfile:

   .. prompt::

      touch small_plot_and_continue

   At the end of a timestep, the code will check if these files exist
   and if so do an output and then remove the file.


#. *How can I check the compilation parameters of a MAESTROeX executable?*

   The build information (including git hashes, modules, EoS, network, etc.) can be displayed by running the executable as

   ::

       ./Maestro.exe --describe

Debugging
=========

#. *How can we dump out a variable to a plotfile from any point in the
   code?*

   ::

           #include <AMReX_VisMF.H>

           VisMF::Write(uold[0],"a_uold");
           VisMF::Write(umac[0][0],"a_umacx");

   This plotfile is visualized using Amrvis using the flag ``-mf``.

#. *How can I print out a MultiFab’s contents from within the code?*

   There is a print subroutine in ``MaestroDebug.cpp`` file. This can
   be simply called as

   ::

         PrintMF(a);


   where ``a`` is a MultiFab (multi-level).

#. *How can I debug a parallel (MPI) job with gdb?*

   If you only need to use a few processors, the following command will work:

   ::

       mpiexec -n 4 xterm -e gdb ./Maestro2d.gnu.ex

   where the executable needs to be created with the ``-g`` flag to
   the compiler. This will pop up multiple xterms with gdb running
   in each. You need to then issue:

   ::

       run inputs

   where inputs is the desired inputs file *in each* xterm.

#. *How can I get more information about floating point exceptions?*

   AMReX can intercept floating point exceptions and provide a helpful
   backtrace file that shows you where they were generated.

#. *How can I get information about potential bugs before running the code?*

   We run `clang-tidy <https://clang.llvm.org/extra/clang-tidy/>`_ on all pull requests using a `GitHub action <https://github.com/AMReX-Astro/cpp-linter-action>`_. ``clang-tidy`` analyzes the source code, produces warnings for potential bugs and offers suggestions for performance improvements.

   ``clang-tidy`` can also be run locally. This requires the ``clang-tidy`` and ``bear`` packages (installed using e.g. ``sudo apt install bear clang-tidy`` on Ubuntu), and the python script
   ``run-clang-tidy.py`` (which can be downloaded from `here <https://github.com/AMReX-Astro/cpp-linter-action/blob/main/run-clang-tidy.py>`_). The analysis is performed by first compiling a problem using the ``bear`` package, then running the python script to analyze the source files. From within a problem directory, run

   .. code-block:: bash

      bear make -j 20 USE_OMP=FALSE USE_MPI=FALSE DEBUG=TRUE

      python3 run-clang-tidy.py -header-filter='MAESTROeX' -ignore-files='amrex|Microphysics' -j 20 > clang-tidy-report.txt

   The compiler flags can be modified to suit the problem to be analyzed, but the ``DEBUG`` flag must be set to ``TRUE``. The ``header-filter`` option for the python script tells the script to only analyze header files containing the given regex pattern, and the ``ignore-files`` flag tells it to ignore any source files containing the given regex pattern. The ``-j`` option tells the script to run a given number of processes in parallel. The output is then redirected to a text file.

I/O
===

#. *How can I tell from a plotfile what runtime parameters were
   used for its run? or when it was created?*

   In each plotfile directory, there is a file called ``job_info``
   (e.g. ``plt00000/job_info``) that lists the build directory and
   date, as well as the value of every runtime parameter for the run.

#. *How can I force the code to output a plotfile / checkpoint
   file at the next step?*

   In the output directory (where the code is running) do ``touch
   .dump_plotfile``. This will create an empty file called
   ``.dump_plotfile.`` At the end of each step, if the code finds
   that file, it will output a plotfile. Simply delete the file to
   restore the code to its normal plotfile behavior.

   Similarly, creating the file ``.dump_checkpoint`` will force the
   output of a checkpoint file.

Algorithm
=========

#. *Why is MAESTROeX so “hard” to use (e.g. as compared to a
   compressible code)?*

   There are several complexities to the algorithm that don’t have
   straightforward compressible counterparts. These mainly involve the
   role of the base state and the constraint equation.

   Care must be taken to setup an initial model/initial base state that
   respects the thermodynamics in MAESTROeX and is in hydrostatic equilibrium.
   Best results are attained when the model is processed with the MAESTROeX EOS and reset into HSE, as is done in the initial_model routines.
   Because MAESTROeX builds off of the base state, any flaws in that initial
   state will influence the subsequent behavior of the algorithm.

   The constraint equation brings another complexity not seen in compressible
   codes—information is instantly communicated
   across the grid. In compressible codes you can track down a problem by
   watching where it starts from and watching it move one cell per dt. In
   MAESTROeX things can go wrong in multiple places without it being obvious
   where the root problem is.

#. *In the final projection in the algorithm, we project* :math:`U^{n+1}` *,
   using a time-centered* :math:`\beta_0` *, a time-centered* :math:`\rho_0` *, but
   an* “:math:`n+1`” *-centered* :math:`S` *. Why then is the
   resulting* :math:`\phi` *(which then defines* :math:`\pi` *) is
   at* “:math:`n+1/2`” *?*

   The short answer to this question is that you should think of this
   as really projecting :math:`(U^{n+1} - U^n)` and the right hand side as having
   :math:`(S^{n+1} - S^n)`. This is because the pressure enters the dynamic equations as
   :math:`(U^{n+1} - U^n) = \ldots + \frac{1}{\rho^{n+1/2}} \nabla \pi^{n+1/2}`.
   (We approximate :math:`\pi^{n+1/2}` by :math:`\pi^{n-1/2}` then do the projection to fix the
   :math:`\pi` as well as the :math:`U`.)

   So everything is in fact time-centered.

#. *Why is* :math:`\gammabar` *computed as the average of the full
   state* :math:`\Gamma_1` *instead of computed from the base state density and
   pressure via the equation of state?*

   The primary reason is that there is no base state composition. The
   base state density is simply the average of the full state density,
   and the base state pressure is the pressure required for hydrostatic
   equilibrium. There is no thermodynamic relationship enforced between
   these base state quantities.

#. *Can I run a full star in 2-d axisymmetric geometry?*

   No. This is a design decision. There is no support for axisymmetric
   coordinates in MAESTROeX. Spherical problems must be run in 3-d.

#. *Why did we switch all the equations over to the* :math:`\tilde{\Ub}` *form
   instead of just working with* :math:`\Ub` *?*

   This is basically a numerical discretization issue. Whenever the base
   state aligns with the grid, you should be able to show that you get
   exactly the same answer each way.

   When you do a spherical star on a 3d Cartesian grid, though, the :math:`w_0`
   is defined on the radial mesh and the :math:`\tilde{\Ub}` on the Cartesian
   mesh, and the :math:`w_0` part never experiences the Cartesian projection,
   for example. So there are differences in exactly how the :math:`w_0` component
   appears (projected on the Cartesian mesh vs. interpolated from the
   radial mesh)—we made the decision at the time to separate the
   components for that reason.

#. *Why does “checkerboarding” appear in the velocity field,
   especially in regions where the flow is stagnant?*

   Checkerboarding can arise from the projection—it doesn’t see that
   mode (because it is an approximate projection) so it is unable to
   remove it. This allows the pattern to slowly build up. There are
   filtering techniques that can be used to remove these modes, but
   they are not implemented in MAESTROeX.

Analysis
========

#. *I want to open a plotfile, derive a new quantity from
   the data stored there, and write out a new plotfile with this derived
   data. How do I do this?*

   One implementation of this can be found in
   ``amrex/Tools/Postprocessing/C_Src/PtwisePltTransform.cpp``. This reads in
   the plotfile data using the the ``AMReX_DataServices`` class, performs a
   transformation on the data based on a set of components specified in the
   command line, and outputs the solution to a new plotfile.
