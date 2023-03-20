
#include <Maestro.H>
using namespace amrex;

// prototype for pertubation function to be called on the
// device (if USE_CUDA=TRUE)
AMREX_GPU_DEVICE
void Perturb(const Real p0, const Real* s0, Real* scal_pert, Real* vel_pert,
             const GpuArray<Real, AMREX_SPACEDIM> prob_lo,
             const GpuArray<Real, AMREX_SPACEDIM> prob_hi, const Real x,
             const Real y, const Real* alpha, const Real* phi, const Real rho_1,
             const Real rho_2, const Real vel_amplitude, const Real vel_width,
             const int nmodes);

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

    if (perturb_model) {
        RealVector alpha(nmodes);
        RealVector phi(nmodes);

        for (auto n = 1; n < nmodes; ++n) {
            alpha[n] = 2.0 * amrex::Random() - 1.0;
            phi[n] = 2.0 * M_PI * amrex::Random();
        }

        const auto prob_lo = geom[lev].ProbLoArray();
        const auto prob_hi = geom[lev].ProbHiArray();
        const auto dx = geom[lev].CellSizeArray();

        Real* AMREX_RESTRICT alpha_p = alpha.dataPtr();
        Real* AMREX_RESTRICT phi_p = phi.dataPtr();

        const auto rho_1_loc = rho_1;
        const auto rho_2_loc = rho_2;
        const auto vel_amplitude_loc = vel_amplitude;
        const auto vel_width_loc = vel_width;
        const auto nmodes_loc = nmodes;

        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
            Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
#if AMREX_SPACEDIM == 3
            Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2];
#endif
            Real scal_pert[Nscal];
            Real vel_pert[AMREX_SPACEDIM];
            Real s0[Nscal];

            for (auto n = 0; n < Nscal; ++n) {
                s0[n] = s0_arr(lev, r, n);
            }

            Perturb(p0_arr(lev, r), s0, scal_pert, vel_pert, prob_lo, prob_hi,
                    x, y, alpha_p, phi_p, rho_1_loc, rho_2_loc,
                    vel_amplitude_loc, vel_width_loc, nmodes_loc);

            scal(i, j, k, Rho) = scal_pert[Rho];
            scal(i, j, k, RhoH) = scal_pert[RhoH];
            scal(i, j, k, Temp) = scal_pert[Temp];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) = scal_pert[FirstSpec + comp];
            }
        });
    }
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitDataSphr()", InitDataSphr);

    Abort("Error: InitLevelDataSphr not implemented");
}

void Perturb(const Real p0, const Real* s0, Real* scal_pert, Real* vel_pert,
             const GpuArray<Real, AMREX_SPACEDIM> prob_lo,
             const GpuArray<Real, AMREX_SPACEDIM> prob_hi, const Real x,
             const Real y, const Real* alpha, const Real* phi, const Real rho_1,
             const Real rho_2, const Real vel_amplitude, const Real vel_width,
             const int nmodes) {
    // apply an optional perturbation to the initial temperature field
    // to see some bubbles

    const Real L_x = prob_hi[0] - prob_lo[0];
    const Real pertheight = 0.01 * 0.5 *
                                (std::cos(2.0 * M_PI * x / L_x) +
                                 std::cos(2.0 * M_PI * (L_x - x) / L_x)) +
                            0.5;

    scal_pert[Rho] = rho_1 + ((rho_2 - rho_1) / 2.0) *
                                 (1.0 + std::tanh((y - pertheight) / 0.005));

    for (auto comp = 0; comp < NumSpec; ++comp) {
        scal_pert[FirstSpec + comp] =
            scal_pert[Rho] * s0[FirstSpec + comp] / s0[Rho];
    }

    for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
        vel_pert[n] = 0.0;
    }

#if (AMREX_SPACEDIM == 3)
    const Real y0 = 0.5 * (prob_lo[1] + prob_hi[1]);
    Real pert = 0.0;

    if (nmodes == 1) {
        pert += vel_amplitude * 0.5 *
                (std::cos(2.0 * M_PI * x / L_x) +
                 std::cos(2.0 * M_PI * (L_x - x) / L_x));
    } else {
        for (auto n = 0; n < nmodes; ++n) {
            pert += vel_amplitude * alpha[n] *
                    std::cos(2.0 * M_PI * x / L_x + phi[n]);
        }
    }
    vel_pert[1] =
        std::exp(-(y - y0) * (y - y0) / (vel_width * vel_width)) * pert;
#endif
}
