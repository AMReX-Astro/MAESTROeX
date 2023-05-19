
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext
void Maestro::MakeHeating(Vector<MultiFab>& rho_Hext,
                          const Vector<MultiFab>& scal) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeHeating()", MakeHeating);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region

            //const Box& tileBox = mfi.tilebox();
            //const auto dx = geom[lev].CellSizeArray();

            //const Array4<const Real> scal_arr = scal[lev].array(mfi);
            //const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            Abort(
                "Make a local copy of MaestroHeating in your problem "
                "directory.");
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
