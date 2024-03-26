.. _sec:gpu:

***
GPU
***

In this chapter, we will present the GPU support in MAESTROeX,
including necessary build parameters, how to offload a routine
to GPU, and some basic profiling and debugging options.
Note that currently MAESTROeX only supports NVIDIA GPUs.

Requirements
============

Since MAESTROeX uses primarily CUDA C++ and CUDA Fortran,
it requires a recent version of CUDA (e.g., >= 9) and the device
must have compute capability of >= 6.

.. _sec:gpubuild:

Building GPU Support
====================

To build MAESTROeX with GPU support, add the following argument
to the GNUmakefile:

::

      USE_CUDA := TRUE

We also need to set ``USE_OMP = FALSE`` because OpenMP is currently
not compatible with building with CUDA.
``USE_CUDA = TRUE`` and ``USE_OMP = TRUE`` will fail to compile.
However, you may use MPI with CUDA for additional parallelization.

Only the IBM and PGI compilers support CUDA Fortran, so the compiler should be set as:

::

      COMP := pgi

The integrator used by the Microphysics library must be set to ``VODE90``:

::

    INTEGRATOR_DIR   := VODE90

Depending on which system you are running on, it may be necessary to specify
the CUDA Capability using the ``CUDA_ARCH`` flag. The CUDA Capability will
depend on the specific GPU hardware you are running on. On a Linux system, the
capability of your device can typically be found by compiling and running the ``deviceQuery``
script found in the CUDA samples directory:
``/usr/local/cuda/samples/1_Utilities/deviceQuery`` (its exact location may
vary depending on where CUDA is installed on your system). The default value of
this flag is 70, corresponding to a capability of 7.x. For a device with
capability 6.x, the flag should be set to:

::

    CUDA_ARCH := 60

.. _sec:gpuporting:

Offloading a routine to GPU
===========================

In order to offload a routine to the GPU, we insert a ``#pragma gpu``
statement on the line before the call. This tells the preprocessor to
generate a CUDA device version of the function, with all the
additional CUDA GPU management. In addition to this, there are a few
other modifications required to ensure the function operates correctly
on the GPU. We summarize these below.

C++
---

- Make sure that the box ``lo`` and ``hi`` are the first two arguments
  and use ``AMREX_INT_ANYD`` as the macro wrapping ``bx.loVect()`` and
  ``bx.hiVect()``

- Likewise, inplace of ``AMREX_ZFILL()``, use ``AMREX_REAL_ANYD``

- Scalars must be passed by value to Fortran functions (not by
  reference)

- Make sure that you change the variable definitions in the Fortran
  code to include value so that the Fortran knows the variables are
  being passed by value rather than by reference. If you don’t do
  this, the code is liable to segfault

- FArrayBox functions like ``setVal()`` and ``saxpy()`` are not on the
  GPU, so you should explicitly write out these operations in C++.
  The MultiFab counterparts are on the GPU.

- Use ``BL_TO_FORTRAN_ANYD()`` to wrap MultiFab arguments (and in the header file wrap with ``BL_FORT_FAB_ARG_3D``)

To illustrate these modifications, consider the function ``mk_sponge``. To call the function on the CPU, we would write

.. code-block:: c++

    for (MFIter mfi(sponge_mf, true); mfi.isValid(); ++mfi) {

        const Box& tileBox = mfi.tilebox();

        mk_sponge(AMREX_ARLIM_3D(tileBox.loVect()), AMREX_ARLIM_3D(tileBox.hiVect()),
                  BL_TO_FORTRAN_3D(sponge_mf[mfi]),
          AMREX_ZFILL(dx), &dt);
     }

Implementing the changes described above to offload this to GPU, this becomes

.. code-block:: c++

    for (MFIter mfi(sponge_mf, true); mfi.isValid(); ++mfi) {

        const Box& tileBox = mfi.tilebox();

    #pragma gpu box(tileBox)
        mk_sponge(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                  BL_TO_FORTRAN_ANYD(sponge_mf[mfi]),
          AMREX_REAL_ANYD(dx), dt);
     }

Fortran
-------

- The routine must only operate on a single zone in ``lo:hi``.

  - Even for temporary arrays, you cannot write to an ``i+1`` zone
    (e.g. when doing limiting) -- this will cause a race condition.
    If necessary, do extra computation to avoid the temporary arrays.

- Mark the routine with ``!$gpu``

- If a module defines its own variables, these variables need to be
  ``allocatable`` (even if they are scalars) and marked as
  ``attributes(managed)``. Additional routines may be needed to
  allocate these variables before they’re used and deallocate them
  when they’re no longer needed.

  - Examples of this can be seen for the sponge parameters. These are
    marked as ``allocatable`` and ``attributes(managed)`` in
    ``sponge.F90``, and therefore must be allocated (by ``init_sponge``).

- Temporary variables must be defined outside of function calls. E.g. if a
  function call contains ``foo(x(a:b)/y)``, you need to define a new variable
  ``z = x(a:b)/y`` then pass this into the function as ``foo(z)``.

  - If you don’t do this, you may see the error ``Array reshaping is
    not supported for device subprogram calls``

- If importing a function from another module, make sure to put the
  import within the function/subroutine, and put ``! function`` at the
  end of the line, e.g.

.. code-block:: fortran

   use my_module, only: my_func ! function

- Individual functions should be imported individually (so not ``use
  my_module, only: func1, func2 ! function``) and there must be a
  space either side of the ``!``

- Make sure the fortran file is ``.F90`` rather than ``.f90`` (and
  remember to update the ``Make.xx`` file to reflect this). If you
  don’t do this, you will see the error ``Label field of continuation
  line is not blank``

  - This is required as we use the convention that ``.F90`` files are
    processed by the preprocessor, and ``.f90`` files are not. The
    preprocessor will therefore only generate the required device
    function if the file has the correct extension.

We can see some of the above modifications by looking at the
subroutine ``estdt`` in ``compute_dt.F90``:

.. code-block:: fortran

   subroutine estdt(lev, dt, umax, lo, hi, dx, &
                    scal,  s_lo, s_hi, nc_s, &
            u,     u_lo, u_hi, nc_u, &
            force, f_lo, f_hi, nc_f, &
            divu,  d_lo, d_hi, &
            dSdt,  t_lo, t_hi, &
            w0_cart, w_lo, w_hi, &
            p0_cart, p_lo, p_hi, &
            gamma1bar_cart, g_lo, g_hi) bind (C,name="estdt")

      use amrex_constants_module, only: HALF
      use amrex_fort_module, only: amrex_min ! function
      use amrex_fort_module, only: amrex_max ! function

      ! input parameters
      integer  , value, intent(in   ) :: lev
      double precision, intent(inout) :: dt, umax
      integer         , intent(in   ) :: lo(3), hi(3)
      ...

      ! local variables
      double precision :: spdx, spdy, spdz, spdr, rho_min
      double precision :: fx, fy, fz, dt_temp
      double precision :: eps,denom,gradp0
      double precision :: a, b, c
      integer          :: i,j,k

      !$gpu

      rho_min = 1.d-20
      ...

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
           spdx = max(spdx ,abs(u(i,j,k,1)))

      ...

   end subroutine estdt

- Here, we can see that ``amrex_min`` and ``amrex_max`` functions from the
  ``amrex_fort_module`` are marked separately as ``! function``, which tells
  the preprocessor to generate a device version of this function.

- The scalar ``lev`` is passed in by value.

- The ``!$gpu`` directive has been inserted after the definition of
  all the variables passed into the routine and all the local
  variables, but before the main body of the function.

- The routine only operates on values in a single zone of ``lo:hi``.


.. To be documented
.. ----------------
..
.. when do we need to mark stuff as attributes(managed)?


.. _sec:gpuprofile:

Profiling with GPUs
===================

NVIDIA's profiler, ``nvprof``, is recommended when profiling for GPUs.
It returns data on how long each kernel launch lasted on the GPU,
the number of threads and registers used, the occupancy of the GPU
and provides recommendations for improving the code.  For more information on how to
use ``nvprof``, see NVIDIA's User's Guide.

If a quicker profiling method is preferred, AMReX's timers can be used
to report some generic timings that may be useful in categorizing an application.
To yield a consistent timing of a routine, a timer will need to be wrapped
around an ``MFIter`` loop that encompasses the entire set of GPU launches
contained within. For example:

.. code-block:: c++

    BL_PROFILE_VAR("A_NAME", blp);     // Profiling start
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        // code that runs on the GPU
    }
    BL_PROFILE_STOP(blp);              // Profiling stop

For now, this is the best way to profile GPU codes using the compiler flag ``TINY_PROFILE = TRUE``.
If you require further profiling detail, use ``nvprof``.

.. _sec:gpudebug:

Basic GPU Debugging
===================

- Turn off GPU offloading for some part of the code with

.. code-block:: c++

    Gpu::setLaunchRegion(0);
    ... ;
    Gpu::setLaunchRegion(1);

- To test if your kernels have launched, run

.. code-block:: sh

   nvprof ./Maestro2d.xxx

- Run under ``nvprof -o profile%p.nvvp ./Maestro2d.xxx`` for
  a small problem and examine page faults using NVIDIA's visual profiler, `nvvp`

- Run under ``cuda-memcheck``

- Run under ``cuda-gdb``

- Run with ``CUDA_LAUNCH_BLOCKING=1``.  This means that only one
  kernel will run at a time.  This can help identify if there are race
  conditions.

