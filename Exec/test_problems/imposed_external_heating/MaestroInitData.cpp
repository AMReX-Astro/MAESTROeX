
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    if (AMREX_SPACEDIM != 2) {
        Abort("Error: InitLevelData only implemented for 2d.");
    }

    // set velocity the correct value c sin(y/c)/(cos(y/c)+b)
    const Real H = 10.0;
    const Real pi_ = 3.141592653589793238462643383279;
    const Real c = H / pi_;
    const Real a = 2;
    const Real b = 2;
    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
        vel(i, j, k, n) = 0.0;
        if (n == 2) {
            vel(i, j, k, n) = c * sin(y / c) / (cos(y / c) + b);
        }
    });

    const auto s0_arr = s0_init.const_array();
    const auto p0_arr = p0_init.const_array();

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
    });

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        const auto r = j;
        const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

        const Real rho0 = s0_arr(lev, r, Rho);

        Real rho_local = 0.0;

        rho_local += rho0;

        eos_t eos_state;

        eos_state.rho = rho_local;
        eos_state.p = p0_arr(lev, r);
        eos_state.T = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 1.0;  // single fluid
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        scal(i, j, k, Rho) = eos_state.rho;
        scal(i, j, k, RhoH) = eos_state.rho * eos_state.h;
        scal(i, j, k, Temp) = eos_state.T;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) =
                eos_state.rho * eos_state.xn[comp];
        }
    });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("InitLevelDataSphr is not implemented.");
}
