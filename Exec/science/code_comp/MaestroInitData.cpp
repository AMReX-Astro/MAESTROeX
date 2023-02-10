
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();
    const auto nrf = base_geom.nr_fine;

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    const auto rho0_loc = rho_0;

    const auto s0_arr = s0_init.const_array();
    const auto p0_arr = p0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int r = AMREX_SPACEDIM == 2 ? j : k;

        const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
#if AMREX_SPACEDIM == 3
        const Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2];
#else
        const Real z = 0.0;
#endif

#if (AMREX_SPACEDIM == 2)
        Real fheat =
            (y < 1.125 * 4.e8) ? std::sin(8.0 * M_PI * (y / 4.e8 - 1.0)) : 0.0;
#else
        Real fheat = (z < 1.125 * 4.e8) ? std::sin(8.0 * M_PI * (z / 4.e8 - 1.0)) : 0.0;
#endif

        Real rhopert =
            5.e-5 * rho0_loc * fheat *
            (std::sin(3.0 * M_PI * x / 4.e8) + std::cos(M_PI * x / 4.e8)) *
            (std::sin(3.0 * M_PI * y / 4.e8) - std::cos(M_PI * y / 4.e8));

        eos_t eos_state;

        eos_state.rho = s0_arr(lev, r, Rho) + rhopert;
        eos_state.p = p0_arr(lev, r);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] =
                s0_arr(lev, r, FirstSpec + comp) / eos_state.rho;
        }

        eos(eos_input_rp, eos_state);

        // set the scalars using eos_state
        scal(i, j, k, Rho) = eos_state.rho;
        scal(i, j, k, RhoH) = eos_state.rho * eos_state.h;
        scal(i, j, k, Temp) = eos_state.T;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) =
                eos_state.xn[comp] * eos_state.rho;
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
    });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("InitLevelDataSphr not implemented");
}
