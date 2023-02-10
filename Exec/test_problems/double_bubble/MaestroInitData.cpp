
#include <Maestro.H>
using namespace amrex;

// prototype for pertubation function to be called on the
// device (if USE_CUDA=TRUE)
AMREX_GPU_DEVICE
void Perturb(const Real p0, const Real* s0, Real* perturbations,
             const GpuArray<Real, AMREX_SPACEDIM> prob_lo,
             const GpuArray<Real, AMREX_SPACEDIM> prob_hi, const Real x,
             const Real y, const Real y_pert_center, const Real pert_width,
             const bool single);

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    const auto tileBox = mfi.tilebox();
    const auto max_lev = base_geom.max_radial_level + 1;
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
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
    });

    // add an optional perturbation
    if (perturb_model) {
        const auto prob_lo = geom[lev].ProbLoArray();
        const auto prob_hi = geom[lev].ProbHiArray();
        const auto dx = geom[lev].CellSizeArray();

        const auto y_pert_center_loc = y_pert_center;
        const auto pert_width_loc = pert_width;
        const auto single_loc = single;

        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
            Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
#if AMREX_SPACEDIM == 3
            Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2];
#else
            Real z = 0.0;
#endif

            Real perturbations[Nscal];
            Real s0[Nscal];

            for (auto n = 0; n < Nscal; ++n) {
                s0[n] = s0_arr(lev, r, n);
            }

            Perturb(p0_arr(lev, r), s0, perturbations, prob_lo, prob_hi, x, y,
                    y_pert_center_loc, pert_width_loc, single_loc);

            scal(i, j, k, Rho) = perturbations[Rho];
            scal(i, j, k, RhoH) = perturbations[RhoH];
            scal(i, j, k, Temp) = perturbations[Temp];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) =
                    perturbations[FirstSpec + comp];
            }
        });
    }
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitDataSphr()", InitDataSphr);

#ifndef AMREX_USE_GPU
    Abort("Error: InitLevelDataSphr not implemented");
#endif
}

void Perturb(const Real p0, const Real* s0, Real* perturbations,
             const GpuArray<Real, AMREX_SPACEDIM> prob_lo,
             const GpuArray<Real, AMREX_SPACEDIM> prob_hi, const Real x,
             const Real y, const Real y_pert_center, const Real pert_width,
             const bool single) {
    // apply an optional perturbation to the initial temperature field
    // to see some bubbles

    eos_t eos_state;

#if (AMREX_SPACEDIM == 2)

    if (!single) {
        Real x1 = prob_lo[0] + (prob_hi[0] - prob_lo[0]) / 3.0;
        Real x2 = prob_lo[0] + 2.0 * (prob_hi[0] - prob_lo[0]) / 3.0;

        Real y1 = y_pert_center;
        Real y2 = y_pert_center;

        Real r1 =
            std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1)) / pert_width;
        Real r2 =
            std::sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2)) / pert_width;

        if (r1 < 2.0) {
            eos_state.rho =
                s0[Rho] * (1.0 - (pert_factor * (1.0 + std::tanh(2.0 - r1))));
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = 0.0;
            }
            eos_state.xn[1] = 1.0;
        } else if (r2 < 2.0) {
            eos_state.rho =
                s0[Rho] * (1.0 - (pert_factor * (1.0 + std::tanh(2.0 - r2))));
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = 0.0;
            }
            eos_state.xn[2] = 1.0;
        } else {
            eos_state.rho = s0[Rho];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = s0[FirstSpec + comp] / s0[Rho];
            }
        }
    } else {
        Real x1 = prob_lo[0] + 0.5 * (prob_hi[0] - prob_lo[0]);

        Real y1 = y_pert_center;

        Real r1 =
            std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1)) / pert_width;

        if (r1 < 2.0) {
            eos_state.rho =
                s0[Rho] * (1.0 - (pert_factor * (1.0 + std::tanh(2.0 - r1))));
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = 0.0;
            }
            eos_state.xn[1] = 1.0;
        } else {
            eos_state.rho = s0[Rho];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = s0[FirstSpec + comp] / s0[Rho];
            }
        }
    }

    // Use the EOS to make this temperature perturbation occur at constant
    // pressure
    eos_state.T = 10000.0;  // guess
    eos_state.p = p0;

    eos(eos_input_rp, eos_state);

    perturbations[Rho] = eos_state.rho;
    perturbations[RhoH] = eos_state.rho * eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        perturbations[FirstSpec + comp] = eos_state.rho * eos_state.xn[comp];
    }

    perturbations[Temp] = eos_state.T;

#else
#ifndef AMREX_USE_GPU
    Abort("Error: Perturb not implemented for 3d");
#endif
#endif
}
