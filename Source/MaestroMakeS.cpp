
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute S at cell-centers
void Maestro::Make_S_cc(
    Vector<MultiFab>& S_cc, Vector<MultiFab>& delta_gamma1_term,
    Vector<MultiFab>& delta_gamma1, const Vector<MultiFab>& scal,
    const Vector<MultiFab>& u, const Vector<MultiFab>& rho_omegadot,
    const Vector<MultiFab>& rho_Hnuc, const Vector<MultiFab>& rho_Hext,
    const Vector<MultiFab>& thermal, const BaseState<Real>& p0_s,
    const BaseState<Real>& gamma1bar, BaseState<Real>& delta_gamma1_termbar) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Make_S_cc()", Make_S_cc);

    // put 1d base state quantities on cartestian grid for spherical case
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gradp0_cart(finest_level + 1);
    Vector<MultiFab> psi_cart(finest_level + 1);

    auto nr_fine = base_geom.nr_fine;
    const auto& r_cc_loc = base_geom.r_cc_loc;

    // calculate gradp0
    BaseState<Real> gradp0(base_geom.max_radial_level + 1, nr_fine);
    auto gradp0_arr = gradp0.array();
    const auto p0 = p0_s.const_array();

    if (spherical) {
        if (use_delta_gamma1_term) {
            Real dr_loc = r_cc_loc(0, 1) - r_cc_loc(0, 0);
            gradp0_arr(0, 0) = (p0(0, 1) - p0(0, 0)) / dr_loc;

            dr_loc = r_cc_loc(0, nr_fine - 1) - r_cc_loc(0, nr_fine - 2);
            gradp0_arr(0, nr_fine - 1) =
                (p0(0, nr_fine - 1) - p0(0, nr_fine - 2)) / dr_loc;

            for (int r = 1; r < nr_fine - 1; r++) {
                dr_loc = r_cc_loc(0, r + 1) - r_cc_loc(0, r - 1);
                gradp0_arr(0, r) = (p0(0, r + 1) - p0(0, r - 1)) / dr_loc;
            }
        }

        for (int lev = 0; lev <= finest_level; ++lev) {
            gradp0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            gradp0_cart[lev].setVal(0.);
        }
    } else {
        if (use_delta_gamma1_term) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                const auto dx = geom[lev].CellSizeArray();

                // bottom and top edge cases for planar
                gradp0_arr(lev, 0) =
                    (p0(lev, 1) - p0(lev, 0)) / dx[AMREX_SPACEDIM - 1];
                gradp0_arr(lev, nr_fine - 1) =
                    (p0(lev, nr_fine - 1) - p0(lev, nr_fine - 2)) /
                    dx[AMREX_SPACEDIM - 1];

                for (int r = 1; r < nr_fine - 1; r++) {
                    gradp0_arr(lev, r) = (p0(lev, r + 1) - p0(lev, r - 1)) /
                                         (2.0 * dx[AMREX_SPACEDIM - 1]);
                }
            }
        }

        for (int lev = 0; lev <= finest_level; ++lev) {
            gradp0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            gradp0_cart[lev].setVal(0.);
        }
    }

    if (use_delta_gamma1_term) {
        Put1dArrayOnCart(gradp0, gradp0_cart, false, false, bcs_f, 0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        psi_cart[lev].define(grids[lev], dmap[lev], 1, 0);

        gamma1bar_cart[lev].setVal(0.);
        p0_cart[lev].setVal(0.);
        psi_cart[lev].setVal(0.);
    }

    if (use_delta_gamma1_term) {
        Put1dArrayOnCart(gamma1bar, gamma1bar_cart, false, false, bcs_f, 0);
        Put1dArrayOnCart(p0_s, p0_cart, false, false, bcs_f, 0);
        Put1dArrayOnCart(psi, psi_cart, false, false, bcs_f, 0);
    }

    const auto use_omegadot_terms_in_S_loc = use_omegadot_terms_in_S;
    const auto use_delta_gamma1_term_loc = use_delta_gamma1_term;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> S_cc_arr = S_cc[lev].array(mfi);
            const Array4<Real> delta_gamma1_term_arr =
                delta_gamma1_term[lev].array(mfi);
            const Array4<Real> delta_gamma1_arr = delta_gamma1[lev].array(mfi);
            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<const Real> rho_odot_arr =
                rho_omegadot[lev].array(mfi);
            const Array4<const Real> rho_Hnuc_arr = rho_Hnuc[lev].array(mfi);
            const Array4<const Real> rho_Hext_arr = rho_Hext[lev].array(mfi);
            const Array4<const Real> thermal_arr = thermal[lev].array(mfi);
            const Array4<const Real> u_arr = u[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> gradp0_cart_arr =
                gradp0_cart[lev].array(mfi);
            const Array4<const Real> gamma1bar_arr =
                gamma1bar_cart[lev].array(mfi);

            if (spherical) {
#if (AMREX_SPACEDIM == 3)
                const Array4<const Real> normal_arr = normal[lev].array(mfi);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    eos_state.rho = scal_arr(i, j, k, Rho);
                    eos_state.T = scal_arr(i, j, k, Temp);
                    for (auto n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] =
                            scal_arr(i, j, k, FirstSpec + n) / eos_state.rho;
                    }
#if NAUX_NET > 0
                    for (auto n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] =
                            scal_arr(i, j, k, FirstAux + n) / eos_state.rho;
                    }
#endif

                    // dens, temp, and xmass are inputs
                    eos(eos_input_rt, eos_state);

                    auto eos_xderivs = composition_derivatives(eos_state);

                    Real sigma =
                        eos_state.dpdT /
                        (eos_state.rho * eos_state.cp * eos_state.dpdr);

                    Real xi_term = 0.0;
                    Real pres_term = 0.0;

                    if (use_omegadot_terms_in_S_loc) {
                        for (auto comp = 0; comp < NumSpec; ++comp) {
                            xi_term -= eos_xderivs.dhdX[comp] *
                                       rho_odot_arr(i, j, k, comp) /
                                       eos_state.rho;

                            pres_term += eos_xderivs.dpdX[comp] *
                                         rho_odot_arr(i, j, k, comp) /
                                         eos_state.rho;
                        }
                    }

                    S_cc_arr(i, j, k) =
                        (sigma / eos_state.rho) *
                            (rho_Hext_arr(i, j, k) + rho_Hnuc_arr(i, j, k) +
                             thermal_arr(i, j, k)) +
                        sigma * xi_term +
                        pres_term / (eos_state.rho * eos_state.dpdr);

                    if (use_delta_gamma1_term_loc) {
                        delta_gamma1_arr(i, j, k) =
                            eos_state.gam1 - gamma1bar_arr(i, j, k);

                        Real U_dot_er = 0.0;
                        for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                            U_dot_er +=
                                u_arr(i, j, k, n) * normal_arr(i, j, k, n);
                        }

                        delta_gamma1_term_arr(i, j, k) =
                            (eos_state.gam1 - gamma1bar_arr(i, j, k)) *
                            u_arr(i, j, k, AMREX_SPACEDIM - 1) *
                            gradp0_cart_arr(i, j, k) /
                            (gamma1bar_arr(i, j, k) * gamma1bar_arr(i, j, k) *
                             p0_arr(i, j, k));
                    } else {
                        delta_gamma1_term_arr(i, j, k) = 0.0;
                        delta_gamma1_arr(i, j, k) = 0.0;
                    }
                });
#endif
            } else {
                const auto anelastic_cutoff_density_coord_lev =
                    base_geom.anelastic_cutoff_density_coord(lev);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    eos_state.rho = scal_arr(i, j, k, Rho);
                    eos_state.T = scal_arr(i, j, k, Temp);
                    for (auto n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] =
                            scal_arr(i, j, k, FirstSpec + n) / eos_state.rho;
                    }
#if NAUX_NET > 0
                    for (auto n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] =
                            scal_arr(i, j, k, FirstAux + n) / eos_state.rho;
                    }
#endif

                    // dens, temp, and xmass are inputs
                    eos(eos_input_rt, eos_state);

                    auto eos_xderivs = composition_derivatives(eos_state);

                    Real sigma =
                        eos_state.dpdT /
                        (eos_state.rho * eos_state.cp * eos_state.dpdr);

                    Real xi_term = 0.0;
                    Real pres_term = 0.0;

                    if (use_omegadot_terms_in_S_loc) {
                        for (auto comp = 0; comp < NumSpec; ++comp) {
                            xi_term -= eos_xderivs.dhdX[comp] *
                                       rho_odot_arr(i, j, k, comp) /
                                       eos_state.rho;

                            pres_term += eos_xderivs.dpdX[comp] *
                                         rho_odot_arr(i, j, k, comp) /
                                         eos_state.rho;
                        }
                    }

                    S_cc_arr(i, j, k) =
                        (sigma / eos_state.rho) *
                            (rho_Hext_arr(i, j, k) + rho_Hnuc_arr(i, j, k) +
                             thermal_arr(i, j, k)) +
                        sigma * xi_term +
                        pres_term / (eos_state.rho * eos_state.dpdr);

                    int r = AMREX_SPACEDIM == 2 ? j : k;

                    if (use_delta_gamma1_term_loc &&
                        r < anelastic_cutoff_density_coord_lev) {
                        delta_gamma1_arr(i, j, k) =
                            eos_state.gam1 - gamma1bar_arr(i, j, k);

                        delta_gamma1_term_arr(i, j, k) =
                            (eos_state.gam1 - gamma1bar_arr(i, j, k)) *
                            u_arr(i, j, k, AMREX_SPACEDIM - 1) *
                            gradp0_cart_arr(i, j, k) /
                            (gamma1bar_arr(i, j, k) * gamma1bar_arr(i, j, k) *
                             p0_arr(i, j, k));
                    } else {
                        delta_gamma1_term_arr(i, j, k) = 0.0;
                        delta_gamma1_arr(i, j, k) = 0.0;
                    }
                });
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(S_cc, 0, 1);
    AverageDown(delta_gamma1_term, 0, 1);

    if (use_delta_gamma1_term) {
        // horizontal average of delta_gamma1_term
        Average(delta_gamma1_term, delta_gamma1_termbar, 0);

        for (int lev = 0; lev <= finest_level; ++lev) {
            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(delta_gamma1_term[lev], TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();

                const Array4<Real> delta_gamma1_term_arr =
                    delta_gamma1_term[lev].array(mfi);
                const Array4<const Real> delta_gamma1_arr =
                    delta_gamma1[lev].array(mfi);
                const Array4<const Real> gamma1bar_arr =
                    gamma1bar_cart[lev].array(mfi);
                const Array4<const Real> psi_arr = psi_cart[lev].array(mfi);
                const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    delta_gamma1_term_arr(i, j, k) +=
                        delta_gamma1_arr(i, j, k) * psi_arr(i, j, k) /
                        (gamma1bar_arr(i, j, k) * gamma1bar_arr(i, j, k) *
                         p0_arr(i, j, k));
                });
            }
        }
    }
}

// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void Maestro::MakeRHCCforNodalProj(Vector<MultiFab>& rhcc,
                                   const Vector<MultiFab>& S_cc,
                                   const BaseState<Real>& Sbar,
                                   const BaseState<Real>& beta0,
                                   const Vector<MultiFab>& delta_gamma1_term) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRHCCforNodalProj()", MakeRHCCforNodalProj);

    Vector<MultiFab> Sbar_cart(finest_level + 1);
    Vector<MultiFab> beta0_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        Sbar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        Sbar_cart[lev].setVal(0.);
        beta0_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(Sbar, Sbar_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(beta0, beta0_cart, false, false, bcs_f, 0);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // fill rhcc

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_cc[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> rhcc_arr = rhcc[lev].array(mfi);
            const Array4<const Real> S_cc_arr = S_cc[lev].array(mfi);
            const Array4<const Real> Sbar_arr = Sbar_cart[lev].array(mfi);
            const Array4<const Real> beta0_arr = beta0_cart[lev].array(mfi);
            const Array4<const Real> delta_gamma1_term_arr =
                delta_gamma1_term[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                rhcc_arr(i, j, k) = beta0_arr(i, j, k) *
                                    (S_cc_arr(i, j, k) - Sbar_arr(i, j, k) +
                                     delta_gamma1_term_arr(i, j, k));
            });
        }
    }

    // averge down and fill ghost cells using first-order extrapolation
    AverageDown(rhcc, 0, 1);
    FillPatch(t_old, rhcc, rhcc, rhcc, 0, 0, 1, 0, bcs_f);
}

void Maestro::CorrectRHCCforNodalProj(Vector<MultiFab>& rhcc,
                                      const BaseState<Real>& rho0,
                                      const BaseState<Real>& beta0,
                                      const BaseState<Real>& gamma1bar,
                                      const BaseState<Real>& p0,
                                      const Vector<MultiFab>& delta_p_term) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CorrectRHCCforNodalProj()",
                   CorrectRHCCforNodalProj);

    // Local variables
    Vector<MultiFab> correction_cc(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        correction_cc[lev].define(grids[lev], dmap[lev], 1, 1);
        correction_cc[lev].setVal(0.);
    }

    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> beta0_cart(finest_level + 1);
    Vector<MultiFab> rho0_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        gamma1bar_cart[lev].setVal(0.);
        p0_cart[lev].setVal(0.);
        beta0_cart[lev].setVal(0.);
        rho0_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(gamma1bar, gamma1bar_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(beta0, beta0_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(rho0, rho0_cart, false, false, bcs_s, Rho);

    const Real dt_loc = dt;
    const Real dpdt_factor_loc = dpdt_factor;
    const Real base_cutoff_density_loc = base_cutoff_density;
    auto base_cutoff_density_coord_loc = base_geom.base_cutoff_density_coord;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(correction_cc[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> correction_arr = correction_cc[lev].array(mfi);
            const Array4<const Real> delta_p_arr = delta_p_term[lev].array(mfi);
            const Array4<const Real> beta0_arr = beta0_cart[lev].array(mfi);
            const Array4<const Real> gamma1bar_arr =
                gamma1bar_cart[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<const Real> rho0_arr = rho0_cart[lev].array(mfi);

            if (spherical) {
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Real correction_factor =
                        rho0_arr(i, j, k) > base_cutoff_density_loc
                            ? beta0_arr(i, j, k) * dpdt_factor_loc /
                                  (gamma1bar_arr(i, j, k) * p0_arr(i, j, k)) /
                                  dt_loc
                            : 0.0;

                    correction_arr(i, j, k) =
                        correction_factor * delta_p_arr(i, j, k);
                });
            } else {
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    int r = (AMREX_SPACEDIM == 2) ? j : k;
                    Real correction_factor =
                        (r < base_cutoff_density_coord_loc(lev))
                            ? beta0_arr(i, j, k) * dpdt_factor_loc /
                                  (gamma1bar_arr(i, j, k) * p0_arr(i, j, k)) /
                                  dt_loc
                            : 0.0;

                    correction_arr(i, j, k) =
                        correction_factor * delta_p_arr(i, j, k);
                });
            }
        }
    }

    // average down and fill ghost cells using first-order extrapolation
    AverageDown(correction_cc, 0, 1);
    FillPatch(t_old, correction_cc, correction_cc, correction_cc, 0, 0, 1, 0,
              bcs_f);

    // add correction term
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab::Add(rhcc[lev], correction_cc[lev], 0, 0, 1, 1);
    }
}

// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void Maestro::MakeRHCCforMacProj(
    Vector<MultiFab>& rhcc, const BaseState<Real>& rho0,
    const Vector<MultiFab>& S_cc, const BaseState<Real>& Sbar,
    const BaseState<Real>& beta0, const Vector<MultiFab>& delta_gamma1_term,
    const BaseState<Real>& gamma1bar, const BaseState<Real>& p0,
    const Vector<MultiFab>& delta_p_term, Vector<MultiFab>& delta_chi,
    const bool is_predictor) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRHCCforMacProj()", MakeRHCCforMacProj);

    // put 1d base state quantities on cartestian grid for spherical case
    Vector<MultiFab> Sbar_cart(finest_level + 1);
    Vector<MultiFab> beta0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> rho0_cart(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        Sbar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
        rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);

        Sbar_cart[lev].setVal(0.);
        beta0_cart[lev].setVal(0.);
        gamma1bar_cart[lev].setVal(0.);
        p0_cart[lev].setVal(0.);
        rho0_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(Sbar, Sbar_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(beta0, beta0_cart, false, false, bcs_f, 0);

    if (dpdt_factor > 0.0) {
        Put1dArrayOnCart(gamma1bar, gamma1bar_cart, false, false, bcs_f, 0);
        Put1dArrayOnCart(p0, p0_cart, false, false, bcs_f, 0);
        Put1dArrayOnCart(rho0, rho0_cart, false, false, bcs_s, Rho);
    }

    const Real base_cutoff_density_loc = base_cutoff_density;
    const Real dt_loc = dt;
    const Real dpdt_factor_loc = dpdt_factor;
    auto base_cutoff_density_coord_loc = base_geom.base_cutoff_density_coord;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // fill rhcc

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_cc[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> rhcc_arr = rhcc[lev].array(mfi);
            const Array4<const Real> S_cc_arr = S_cc[lev].array(mfi);
            const Array4<const Real> Sbar_arr = Sbar_cart[lev].array(mfi);
            const Array4<const Real> beta0_arr = beta0_cart[lev].array(mfi);
            const Array4<const Real> rho0_arr = rho0_cart[lev].array(mfi);
            const Array4<const Real> delta_gamma1_arr =
                delta_gamma1_term[lev].array(mfi);
            const Array4<const Real> delta_p_arr = delta_p_term[lev].array(mfi);
            const Array4<const Real> gamma1bar_arr =
                gamma1bar_cart[lev].array(mfi);
            const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
            const Array4<Real> delta_chi_arr = delta_chi[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                rhcc_arr(i, j, k) = beta0_arr(i, j, k) *
                                    (S_cc_arr(i, j, k) - Sbar_arr(i, j, k) +
                                     delta_gamma1_arr(i, j, k));
            });

            if (dpdt_factor > 0.0) {
                if (is_predictor) {
                    ParallelFor(tileBox,
                                [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                    delta_chi_arr(i, j, k) = 0.0;
                                });
                }

                if (spherical) {
                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            if (rho0_arr(i, j, k) > base_cutoff_density_loc) {
                                delta_chi_arr(i, j, k) +=
                                    dpdt_factor_loc * delta_p_arr(i, j, k) /
                                    (dt_loc * gamma1bar_arr(i, j, k) *
                                     p0_arr(i, j, k));
                                rhcc_arr(i, j, k) +=
                                    beta0_arr(i, j, k) * delta_chi_arr(i, j, k);
                            }
                        });
                } else {
                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            int r = (AMREX_SPACEDIM == 2) ? j : k;
                            if (r < base_cutoff_density_coord_loc(lev)) {
                                delta_chi_arr(i, j, k) +=
                                    dpdt_factor_loc * delta_p_arr(i, j, k) /
                                    (dt_loc * gamma1bar_arr(i, j, k) *
                                     p0_arr(i, j, k));
                                rhcc_arr(i, j, k) +=
                                    beta0_arr(i, j, k) * delta_chi_arr(i, j, k);
                            }
                        });
                }
            }
        }
    }
}
