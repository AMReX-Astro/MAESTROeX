
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext
void Maestro::MakeHeating(Vector<MultiFab>& rho_Hext,
                          const Vector<MultiFab>& scal) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeHeating()", MakeHeating);

    const auto t_local = t_old;

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

                const Real L_x = 2.5e8;
#if (AMREX_SPACEDIM == 3)
                const Real L_y = 2.5e8;
#endif
                if (t_local <= 200.0) {
                    const Real r_layer = 1.25e8;

                    Real x = (Real(i) + 0.5) * dx[0] + prob_lo[0];
                    Real y = (Real(j) + 0.5) * dx[1] + prob_lo[1];
#if AMREX_SPACEDIM == 3
                    Real z = (Real(k) + 0.5) * dx[2] + prob_lo[2];
#else
                    Real z = 0.0;
#endif

                    Real r = (AMREX_SPACEDIM == 2) ? y : z;

                    Real er = std::exp(-(r - r_layer) * (r - r_layer) / 1.e14);

                    rho_Hext_arr(i, j, k) =
                        er *
                        (1.0 +
                         .00625 * std::sin(2.0 * M_PI * x / L_x)
#if (AMREX_SPACEDIM == 3)
                             * std::sin(2.0 * M_PI * y / L_y)
#endif
                         +
                         .01875 * std::sin((6.0 * M_PI * x / L_x) + M_PI / 3.0)
#if (AMREX_SPACEDIM == 3)
                             * std::sin((6.0 * M_PI * y / L_y) + M_PI / 3.0)
#endif
                         +
                         .01250 * std::sin((8.0 * M_PI * x / L_x) + M_PI / 5.0)
#if (AMREX_SPACEDIM == 3)
                             * std::sin((8.0 * M_PI * y / L_y) + M_PI / 5.0)
#endif
                             ) *
                        2.5e16;

                    rho_Hext_arr(i, j, k) *= scal_arr(i, j, k, Rho);
                }
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
