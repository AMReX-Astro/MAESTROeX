// Initialize and perturb the temperature

#include <Maestro.H>
#include <random>

using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

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
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
#if NAUX_NET > 0
        for (auto comp = 0; comp < NumAux; ++comp) {
            scal(i, j, k, FirstAux + comp) = s0_arr(lev, r, FirstAux + comp);
        }
#endif
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;

    });

    if (perturb_model) {

        // Generate random seed (or get from input), and RNG
        int seed = 0;
        if (problem_rp::pert_seed != -1) {
            seed = problem_rp::pert_seed;
        } else {
            std::random_device r;
            seed = r();
        }
        std::default_random_engine generator(seed);

        // Gaussian distribution
        std::normal_distribution<double> noise(1.0, 0.001); // mean of 1, std dev 10^-3

        const auto lo = amrex::lbound(tileBox);
        const auto hi = amrex::ubound(tileBox);

        // Loop through data and perturb
        // Because for and not parallelfor, need to run the first step not in serial (i.e. not with srun)
        for (int k = lo.z; k<= hi.z; k++) {
            for (int j = lo.y; j<= hi.y; j++) {
                for (int i = lo.x; i<= hi.x; i++) {
                    int r = AMREX_SPACEDIM == 2 ? j : k;  // j if 2D, k if 3D

                    Real t0 = s0_arr(lev, r, Temp);
                    Real temp = t0 * noise(generator);
                    Real dens = s0_arr(lev, r, Rho); // no change

                    // Create new eos object based on these modified values
                    eos_t eos_state;
                    eos_state.T = temp;
                    eos_state.p = p0_arr(lev, r);
                    eos_state.rho = dens;
                    for (auto comp = 0; comp < NumSpec; comp++) {
                        eos_state.xn[comp] =
                            s0_arr(lev, r, FirstSpec + comp) / s0_arr(lev, r, Rho);
                    }

                    auto eos_input_flag = eos_input_tp; // temperature & pressure eos
                    eos(eos_input_flag, eos_state);
                    scal(i, j, k, Rho) = eos_state.rho;
                    scal(i, j, k, RhoH) = eos_state.rho * eos_state.h; // re-compute enthalpy
                    scal(i, j, k, Temp) = eos_state.T;
                    for (auto comp=0; comp < NumSpec; comp++) {
                        scal(i, j, k, FirstSpec + comp) =
                            eos_state.rho * eos_state.xn[comp];
                    }
                }
            }
        }
    }
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {

    amrex::ignore_unused(lev);
    amrex::ignore_unused(time);
    amrex::ignore_unused(scal);
    amrex::ignore_unused(vel);

    Abort("InitLevelDataSphr not implemented");
}
