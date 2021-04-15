
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
            const auto prob_hi = geom[lev].ProbHiArray();

            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            const auto heat_flux_loc = constant_heat_flux;

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Real r = (AMREX_SPACEDIM == 2) ? j : k;
                // Real z = prob_lo[AMREX_SPACEDIM - 1] +
                //          (Real(r) + 0.5) * dx[AMREX_SPACEDIM - 1];

                rho_Hext_arr(i, j, k) = 0.0;

                // if (z > 0.99 * prob_hi[AMREX_SPACEDIM - 1]) {
                //     rho_Hext_arr(i, j, k) = scal_arr(i, j, k, Rho) * heat_flux_loc;
                // }

                const Real L_x = prob_hi[0];

                Real x = (Real(i) + 0.5) * dx[0] + prob_lo[0];
                Real y = (Real(j) + 0.5) * dx[1] + prob_lo[1];

                Real er = std::exp(-(prob_hi[1] - y) * (prob_hi[1] - y) / 1.e11);

                rho_Hext_arr(i, j, k) =
                    er *
                    (1.0 +
                        .00625 * std::sin(2.0 * M_PI * x / L_x)
                        +
                        .01875 * std::sin((6.0 * M_PI * x / L_x) + M_PI / 3.0)
                        +
                        .01250 * std::sin((8.0 * M_PI * x / L_x) + M_PI / 5.0)
                            ) *
                    constant_heat_flux;
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
