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

MAESTROeX has only been tested with NVIDIA/CUDA.  In theory AMD/HIP
should work.

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

