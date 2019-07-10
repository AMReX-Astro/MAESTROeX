
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeGamma1bar (const Vector<MultiFab>& scal,
                        RealVector& gamma1bar,
                        const RealVector& p0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGamma1bar()",MakeGamma1bar);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    Vector<MultiFab> gamma1(finest_level+1);
    Vector<MultiFab> p0_cart(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        gamma1[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            p0_cart[lev].setVal(0.);
        }

        Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& gamma1_mf = gamma1[lev];
        const MultiFab& scal_mf = scal[lev];
        const MultiFab& p0_mf = p0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(gamma1_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            if (spherical == 0) {
#pragma gpu box(tileBox)
                make_gamma(AMREX_INT_ANYD(tileBox.loVect()),
                           AMREX_INT_ANYD(tileBox.hiVect()),
                           lev,
                           BL_TO_FORTRAN_ANYD(gamma1_mf[mfi]),
                           BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                           p0.dataPtr());
            } else {

#pragma gpu box(tileBox)
                make_gamma_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
                                BL_TO_FORTRAN_ANYD(gamma1_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(scal_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(p0_mf[mfi]));

            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(gamma1,0,1);

    // call average to create gamma1bar
    Average(gamma1,gamma1bar,0);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}
