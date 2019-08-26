#ifndef _MaestroBCThreads_H_
#define _MaestroBCThreads_H_

#include <AMReX_BLFort.H>
#include <Maestro.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <AMReX_Arena.H>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

// Define some routines to allow boundary conditions
// to be passed to offloaded Fortran subroutines
#ifdef AMREX_USE_CUDA
    static void set_bc_launch_config()
    {
        amrex::Gpu::Device::setNumThreadsMin(Maestro::minBCThreads(0), Maestro::minBCThreads(1), Maestro::minBCThreads(2));
    }

    static void clean_bc_launch_config()
    {
        amrex::Gpu::Device::setNumThreadsMin(1, 1, 1);
    }

    // Return a pointer to bc valid for use in Fortran. For the CPU this is a no-op.

    static int* prepare_bc(const int* bc, const int nvar)
    {
        int* bc_f = (int*) amrex::The_Arena()->alloc(AMREX_SPACEDIM * 2 * nvar * sizeof(int));
        AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(bc_f, bc, AMREX_SPACEDIM * 2 * nvar * sizeof(int), cudaMemcpyHostToDevice, amrex::Gpu::Device::cudaStream()));
        return bc_f;
    }

    static void clean_bc(int* bc_f)
    {
        amrex::Gpu::Device::streamSynchronize();
        amrex::The_Arena()->free(bc_f);
    }
#endif

#ifdef __cplusplus
}
#endif

#endif