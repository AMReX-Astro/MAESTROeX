
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    const auto tileBox = mfi.tilebox();

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    const auto ymid = prob_lo[1] + 0.5 * (prob_hi[1] - prob_lo[1]);
    const auto L_x = (prob_hi[0] - prob_lo[0]);
    const auto yshr1 = prob_lo[1] + 0.25 * (prob_hi[1] - prob_lo[1]);
    const auto yshr2 = prob_lo[1] + 0.75 * (prob_hi[1] - prob_lo[1]);

    const auto yt_loc = yt;
    const auto delx_loc = delx;

    const auto s0_arr = s0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int r = AMREX_SPACEDIM == 2 ? j : k;
        Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
        Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

        // set the scalars using s0
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;

        // x velocity as a function of y
        if (y < ymid) {
            vel(i, j, k, 0) = std::tanh((y - yshr1) / yt_loc);
        } else {
            vel(i, j, k, 0) = std::tanh((yshr2 - y) / yt_loc);
        }

        // y velocity as a function of x
        vel(i, j, k, 1) = delx_loc * std::sin(2.0 * M_PI * x / L_x);
    });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitDataSphr()", InitDataSphr);

    Abort("Error: InitLevelDataSphr not implemented");
}
