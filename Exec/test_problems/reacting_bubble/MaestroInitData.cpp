
#include <Maestro.H>
#include <actual_network.H>
#include <network_properties.H>
using namespace amrex;

// prototype for pertubation function to be called on the
// device (if USE_CUDA=TRUE)
AMREX_GPU_DEVICE
void Perturb(const Real p0_init, const Real* s0, Real* perturbations,
             const Real x, const Real y, const Real z);

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

    // add an optional perturbation
    if (perturb_model) {
        const auto prob_lo = geom[lev].ProbLoArray();
        const auto dx = geom[lev].CellSizeArray();

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

            Perturb(p0_arr(lev, r), s0, perturbations, x, y, z);

            scal(i, j, k, Rho) = perturbations[Rho];
            scal(i, j, k, RhoH) = perturbations[RhoH];
            scal(i, j, k, Temp) = perturbations[Temp];

            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) =
                    perturbations[FirstSpec + comp];
            }
#if NAUX_NET > 0
            for (auto comp = 0; comp < NumAux; ++comp) {
                scal(i, j, k, FirstAux + comp) = perturbations[FirstAux + comp];
            }
#endif
        });
    }
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        // set velocity to zero
        ParallelFor(tileBox, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        vel_arr(i, j, k, n) = 0.0;
                    });

        ParallelFor(tileBox, Nscal,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        scal_arr(i, j, k, n) = 0.0;
                    });
    }

    // if we are spherical, we want to make sure that p0 is good, since that is
    // what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    // then initialize h from rho, X, and p0.
    Vector<MultiFab> p0_cart(finest_level);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal
    Vector<MultiFab> temp_mf(finest_level);

    for (auto i = 0; i <= finest_level; ++i) {
        p0_cart[i].define(grids[lev], dmap[lev], 1, ng_s);
        temp_mf[i].define(grids[lev], dmap[lev], 1, ng_s);
    }

    BaseState<Real> temp_vec(base_geom.max_radial_level + 1, base_geom.nr_fine);
    auto temp_arr = temp_vec.array();

    const auto s0_init_arr = s0_init.const_array();

    // initialize temperature
    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            temp_arr(l, r) = s0_init_arr(l, r, Temp);
        }
    }

    Put1dArrayOnCart(temp_vec, temp_mf, false, false, bcs_s, Temp);
    MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, scal.nGrow());

    // initialize p0_cart
    Put1dArrayOnCart(p0_init, p0_cart, false, false, bcs_f, 0);

    // initialize species
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            for (auto r = 0; r < base_geom.nr_fine; ++r) {
                temp_arr(l, r) = s0_init_arr(l, r, FirstSpec + comp);
            }
        }
        Put1dArrayOnCart(temp_vec, temp_mf, false, false, bcs_s,
                         FirstSpec + comp);
        MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, scal.nGrow());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> scal_arr = scal.array(mfi);
        const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

        // initialize rho as sum of partial densities rho*X_i
        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal_arr(i, j, k, Rho) += scal_arr(i, j, k, FirstSpec + comp);
            }

            // initialize the aux variables
#ifdef AUX_THERMO
            for (auto comp = 0; comp < NumAux; ++comp) {
                scal_arr(i, j, k, FirstAux + comp) = 0.0;
            }

            for (auto comp = 0; comp < NumSpec; ++comp) {
                // set the aux quantities
                scal_arr(i, j, k, FirstAux + iye) +=
                    scal_arr(i, j, k, FirstSpec + comp) * zion[comp] *
                    aion_inv[comp];
                scal_arr(i, j, k, FirstAux + iabar) +=
                    scal_arr(i, j, k, FirstSpec + comp) * aion_inv[comp];
                scal_arr(i, j, k, FirstAux + ibea) +=
                    scal_arr(i, j, k, FirstSpec + comp) *
                    aprox19::bion(comp + 1) * aion_inv[comp];
            }

            scal_arr(i, j, k, FirstAux + iabar) =
                1.0_rt / scal_arr(i, j, k, FirstAux + iabar);

            for (auto comp = 0; comp < NumAux; ++comp) {
                scal_arr(i, j, k, FirstAux + comp) *= scal_arr(i, j, k, Rho);
            }
#endif

            // initialize (rho h) and T using the EOS
            eos_t eos_state;
            eos_state.T = scal_arr(i, j, k, Temp);
            eos_state.p = p0_arr(i, j, k);
            eos_state.rho = scal_arr(i, j, k, Rho);
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] =
                    scal_arr(i, j, k, FirstSpec + comp) / eos_state.rho;
            }
#if NAUX_NET > 0
            for (auto comp = 0; comp < NumAux; ++comp) {
                eos_state.aux[comp] =
                    scal_arr(i, j, k, FirstAux + comp) / eos_state.rho;
            }
#endif

            eos(eos_input_rp, eos_state);

            scal_arr(i, j, k, RhoH) = eos_state.rho * eos_state.h;
            scal_arr(i, j, k, Temp) = eos_state.T;
        });

        if (perturb_model) {
            const auto prob_lo = geom[lev].ProbLoArray();
            const auto dx = geom[lev].CellSizeArray();

            // add an optional perturbation
            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
                Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
                Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2];

                Real perturbations[Nscal];
                Real s0[Nscal];

                for (auto n = 0; n < Nscal; ++n) {
                    s0[n] = scal_arr(i, j, k, n);
                }

                Perturb(p0_arr(i, j, k), s0, perturbations, x, y, z);

                scal_arr(i, j, k, Rho) = perturbations[Rho];
                scal_arr(i, j, k, RhoH) = perturbations[RhoH];
                scal_arr(i, j, k, Temp) = perturbations[Temp];
#ifdef AUX_THERMO
                // initialize the aux quantities
                for (auto comp = 0; comp < NumAux; ++comp) {
                    scal_arr(i, j, k, FirstAux + comp) = 0.0;
                }
#endif

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    scal_arr(i, j, k, FirstSpec + comp) =
                        perturbations[FirstSpec + comp];
#ifdef AUX_THERMO
                    // set the aux quantities
                    scal_arr(i, j, k, FirstAux + iye) +=
                        scal_arr(i, j, k, FirstSpec + comp) * zion[comp] *
                        aion_inv[comp];
                    scal_arr(i, j, k, FirstAux + iabar) +=
                        scal_arr(i, j, k, FirstSpec + comp) * aion_inv[comp];
                    scal_arr(i, j, k, FirstAux + ibea) +=
                        scal_arr(i, j, k, FirstSpec + comp) *
                        aprox19::bion(comp + 1) * aion_inv[comp];
                }

                scal_arr(i, j, k, FirstAux + iabar) =
                    1.0_rt / scal_arr(i, j, k, FirstAux + iabar);

                for (auto comp = 0; comp < NumAux; ++comp) {
                    scal_arr(i, j, k, FirstAux + comp) *=
                        scal_arr(i, j, k, Rho);
#endif
                }
            });
        }
    }
}

void Perturb(const Real p0_init, const Real* s0, Real* perturbations,
             const Real x, const Real y, const Real z) {
    Real t0 = s0[Temp];

#if (AMREX_SPACEDIM == 2)

    Real x1 = 5.0e7;
    Real y1 = 6.5e7;
    Real r1 = std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1)) /
              (2.5e6 * pert_rad_factor);

    Real x2 = 1.2e8;
    Real y2 = 8.5e7;
    Real r2 = std::sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2)) /
              (2.5e6 * pert_rad_factor);

    Real x3 = 2.0e8;
    Real y3 = 7.5e7;
    Real r3 = std::sqrt((x - x3) * (x - x3) + (y - y3) * (y - y3)) /
              (2.5e6 * pert_rad_factor);

    // this is a tiny bubble for inputs_2d_smalldomain
    Real x4 = 0.5;
    Real y4 = 85218750.25;
    Real r4 = std::sqrt((x - x4) * (x - x4) + (y - y4) * (y - y4)) /
              (2.5e-2 * pert_rad_factor);

    Real temp = 0.0;

    if (do_small_domain) {
        temp = t0 *
               (1.0 + pert_temp_factor * (0.150 * (1.0 + std::tanh(2.0 - r1)) +
                                          0.300 * (1.0 + std::tanh(2.0 - r2)) +
                                          0.225 * (1.0 + std::tanh(2.0 - r3)) +
                                          0.300 * (1.0 + std::tanh(2.0 - r4))));
    } else {
        temp = t0 *
               (1.0 + pert_temp_factor * (0.150 * (1.0 + std::tanh(2.0 - r1)) +
                                          0.300 * (1.0 + std::tanh(2.0 - r2)) +
                                          0.225 * (1.0 + std::tanh(2.0 - r3))));
    }

#else

    Real temp = 0.0;

    if (!maestro::spherical) {
        Real x0 = 1.8e7;
        Real y0 = 1.8e7;
        Real z0 = 8.5e7;

        Real r0 = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) +
                            (z - z0) * (z - z0)) /
                  2.5e6;

        temp = t0 * (1.0 + 2.0 * (0.15 * (1.0 + std::tanh((2.0 - r0)))));
    } else {
        // center of the star is a 2.5e8
        Real x0 = 2.5e8;
        Real y0 = 2.5e8;
        Real z0 = 2.5e8;

        Real r0 = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) +
                            (z - z0) * (z - z0)) /
                  2.5e6;

        // note extra factor of 0.5 for spherical case
        temp = t0 * (1.0 + 2.0 * (0.15 * 0.5 * (1.0 + std::tanh((2.0 - r0)))));
    }

#endif

    eos_t eos_state;

    // Use the EOS to make this temperature perturbation occur at constant
    // pressure
    eos_state.T = temp;
    eos_state.p = p0_init;
    eos_state.rho = s0[Rho];
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = s0[FirstSpec + comp] / s0[Rho];
    }
#if NAUX_NET > 0
    for (auto comp = 0; comp < NumAux; ++comp) {
        eos_state.aux[comp] = s0[FirstAux + comp] / s0[Rho];
    }
#endif

    eos(eos_input_tp, eos_state);

    perturbations[Rho] = eos_state.rho;
    perturbations[RhoH] = eos_state.rho * eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        perturbations[FirstSpec + comp] = eos_state.rho * eos_state.xn[comp];
    }
#if NAUX_NET > 0
    for (auto comp = 0; comp < NumAux; ++comp) {
        perturbations[FirstAux + comp] = eos_state.rho * eos_state.aux[comp];
    }
#endif

    perturbations[Temp] = temp;
}
