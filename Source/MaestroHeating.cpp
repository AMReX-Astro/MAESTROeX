
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute heating term, rho_Hext
void Maestro::MakeHeating(Vector<MultiFab>& rho_Hext,
                          const Vector<MultiFab>& scal) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeHeating()", MakeHeating);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        MultiFab& rho_Hext_mf = rho_Hext[lev];
        const MultiFab& scal_mf = scal[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_heating(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                         BL_TO_FORTRAN_3D(rho_Hext_mf[mfi]),
                         BL_TO_FORTRAN_FAB(scal_mf[mfi]), dx, &t_old);
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
