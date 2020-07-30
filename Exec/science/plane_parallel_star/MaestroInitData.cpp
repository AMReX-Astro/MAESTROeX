
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

#if (AMREX_SPACEDIM == 3)
    Abort("This should only be run for 2d problems.");
#endif

    const auto tileBox = mfi.tilebox();
    const int max_lev = base_geom.max_radial_level + 1;
    const auto nrf = base_geom.nr_fine;

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    const auto s0_arr = s0_init.const_array();
    const auto p0_arr = p0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        // initialize rho as sum of partial densities rho*X_i
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;

        // initialize (rho h) and T using the EOS
        eos_t eos_state;
        eos_state.rho = scal(i, j, k, Rho);
        eos_state.p = p0_arr(lev, r);
        eos_state.T = scal(i, j, k, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] =
                scal(i, j, k, FirstSpec + comp) / eos_state.rho;
        }

        eos(eos_input_rp, eos_state);

        scal(i, j, k, RhoH) = eos_state.rho * eos_state.h;
        scal(i, j, k, Temp) = eos_state.T;
    });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("spherical not implemented");
}
