
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeSponge (Vector<MultiFab>& sponge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeSponge()",MakeSponge);

#ifdef AMREX_USE_CUDA
    // turn on GPU
    Cuda::setLaunchRegion(true);
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& sponge_mf = sponge[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(sponge_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "tileBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
            mk_sponge(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                      BL_TO_FORTRAN_ANYD(sponge_mf[mfi]),AMREX_REAL_ANYD(dx),dt);
        }

    }

    // average fine data onto coarser cells
    AverageDown(sponge,0,1);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    Cuda::setLaunchRegion(false);
#endif

}
