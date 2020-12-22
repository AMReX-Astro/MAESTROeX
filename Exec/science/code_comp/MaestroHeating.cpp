
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
        for (MFIter mfi(rho_Hext[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();
            const auto prob_lo = geom[lev].ProbLoArray();

            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                rho_Hext_arr(i, j, k) = 0.0;

                Real r = (AMREX_SPACEDIM == 2) ? j : k;
                Real z = prob_lo[AMREX_SPACEDIM - 1] +
                         (Real(r) + 0.5) * dx[AMREX_SPACEDIM - 1];

                if (z < 1.125 * 4.e8) {
                    Real fheat = std::sin(8.0 * M_PI * (z / 4.e8 - 1.0));

                    rho_Hext_arr(i, j, k) = heating_factor * fheat;
                }
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
