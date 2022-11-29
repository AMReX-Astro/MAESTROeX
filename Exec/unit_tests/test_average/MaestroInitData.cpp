
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    // set scalars to zero
    ParallelFor(tileBox, Nscal,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    scal(i, j, k, n) = 0.0;
                });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        // set velocity to zero
        ParallelFor(tileBox, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        vel_arr(i, j, k, n) = 0.0;
                    });

        ParallelFor(tileBox, Nscal,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        scal_arr(i, j, k, n) = 0.0;
                    });
    }
}
