
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeGamma1bar (const Vector<MultiFab>& scal,
                        Vector<Real>& gamma1bar,
                        const Vector<Real>& p0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGamma1bar()",MakeGamma1bar);

    Vector<MultiFab> gamma1(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        gamma1[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& gamma1_mf = gamma1[lev];
        const MultiFab&   scal_mf =   scal[lev];
        const MultiFab&   cc_to_r = cell_cc_to_r[lev];

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
                make_gamma(&lev, ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                           BL_TO_FORTRAN_3D(gamma1_mf[mfi]),
                           BL_TO_FORTRAN_FAB(scal_mf[mfi]),
                           p0.dataPtr());
            } else {
                const Real* dx = geom[lev].CellSize();

                make_gamma_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                                BL_TO_FORTRAN_3D(gamma1_mf[mfi]),
                                BL_TO_FORTRAN_FAB(scal_mf[mfi]),
                                p0.dataPtr(), dx,
                                r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
                                BL_TO_FORTRAN_3D(cc_to_r[mfi]));
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(gamma1,0,1);

    // call average to create gamma1bar
    Average(gamma1,gamma1bar,0);
}
