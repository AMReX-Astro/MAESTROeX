
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

            const auto lo = tileBox.loVect3d();
            const auto hi = tileBox.hiVect3d();

            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                rho_Hext_arr(i, j, k) = 0.0;

                const Real amp = 1.e17;
                Real r0 = ((Real(hi[2]) + 0.5) * dx[2] -
                           (Real(lo[2]) + 0.5) * dx[2]) /
                          2.0;
                Real sig = 0.5 * r0;

                Real x = (Real(i) + 0.5) * dx[0] + prob_lo[0];
                Real y = (Real(j) + 0.5) * dx[1] + prob_lo[1];
                Real z = (Real(k) + 0.5) * dx[2] + prob_lo[2];

                rho_Hext_arr(i, j, k) =
                    scal_arr(i, j, k, Rho) * amp *
                    std::exp(-(x - r0) * (x - r0) / (sig * sig)) *
                    std::exp(-(y - r0) * (y - r0) / (sig * sig)) *
                    std::exp(-(z - r0) * (z - r0) / (sig * sig));
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
