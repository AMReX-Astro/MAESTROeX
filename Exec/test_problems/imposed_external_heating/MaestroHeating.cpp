
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
                const Real H = 10.0;
                const Real a = 2.0;
                const Real b = 2.0;
                const Real Ts = 1.0;
                Real c = H / M_PI;
                Real EE = std::exp(-t_local / Ts);

                Real x = (Real(i) + 0.5) * dx[0] + prob_lo[0];
                Real y = (Real(j) + 0.5) * dx[1] + prob_lo[1];

                Real BOT = EE * std::cos(y / c) + b;
                Real FF = std::exp(-2.0 * t_local / Ts) * std::sin(y / c) *
                          std::sin(y / c) / (BOT * BOT);
                FF += EE * std::cos(y / c) / BOT;
                rho_Hext_arr(i, j, k) = scal_arr(i, j, k, Rho) * FF / Ts;
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
