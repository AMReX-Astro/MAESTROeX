#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::Makew0(const BaseState<Real>& w0_old, BaseState<Real>& w0_force,
                     const BaseState<Real>& Sbar_in,
                     const BaseState<Real>& rho0_old_in,
                     const BaseState<Real>& rho0_new_in,
                     const BaseState<Real>& p0_old_in,
                     const BaseState<Real>& p0_new_in,
                     const BaseState<Real>& gamma1bar_old_in,
                     const BaseState<Real>& gamma1bar_new_in,
                     const BaseState<Real>& p0_minus_peosbar, const Real dt_in,
                     const Real dtold_in, const bool is_predictor) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0()", Makew0);

    w0_force.setVal(0.);

    if (!spherical) {
        if (do_planar_invsq_grav || do_2d_planar_octant) {
            Makew0PlanarVarg(w0_old, w0_force, Sbar_in, rho0_old_in,
                             rho0_new_in, p0_old_in, p0_new_in,
                             gamma1bar_old_in, gamma1bar_new_in,
                             p0_minus_peosbar, dt_in, dtold_in);
        } else {
            Makew0Planar(w0_old, w0_force, Sbar_in,
                         p0_old_in, p0_new_in, gamma1bar_old_in,
                         gamma1bar_new_in, p0_minus_peosbar, dt_in, dtold_in,
                         is_predictor);
        }
    } else {
        if (use_exact_base_state) {
            Makew0SphrIrreg(w0_old, w0_force, Sbar_in, rho0_old_in, rho0_new_in,
                            p0_old_in, p0_new_in, gamma1bar_old_in,
                            gamma1bar_new_in, p0_minus_peosbar, dt_in,
                            dtold_in);
        } else {
            Makew0Sphr(w0_old, w0_force, Sbar_in, rho0_old_in, rho0_new_in,
                       p0_old_in, p0_new_in, gamma1bar_old_in, gamma1bar_new_in,
                       p0_minus_peosbar, dt_in, dtold_in);
        }
    }

    if (maestro_verbose >= 2) {
        const auto w0_arr = w0.const_array();
        for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
            Real max_w0 = 0.0;
            for (auto r = base_geom.r_start_coord(n, 1);
                 r <= base_geom.r_end_coord(n, 1) + 1; ++r) {
                max_w0 = amrex::max(max_w0, amrex::Math::abs(w0_arr(n, r)));
            }
            Print() << "... max CFL of w0: " << max_w0 * dt_in / base_geom.dr(n)
                    << std::endl;
        }
        Print() << std::endl;
    }
}

void Maestro::Makew0Planar(
    const BaseState<Real>& w0_old, BaseState<Real>& w0_force,
    const BaseState<Real>& Sbar_in, const BaseState<Real>& p0_old_in,
    const BaseState<Real>& p0_new_in, const BaseState<Real>& gamma1bar_old_in,
    const BaseState<Real>& gamma1bar_new_in,
    const BaseState<Real>& p0_minus_peosbar, const Real dt_in,
    const Real dtold_in, const bool is_predictor) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0Planar()", Makew0Planar);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Multilevel Outline
    //
    // Compute w0 at level 1 only
    // Initialize new w0 at bottom of coarse base array to 0.0.
    // do n=1,base_geom.finest_radial_level
    //   Compute w0 on edges at level n
    //   Obtain the starting value of w0 from the coarser grid
    //   if n>1, compare the difference between w0 at top of level n to the
    //           corresponding point on level n-1
    //   do i=n-1,1,-1
    //     Restrict w0 from level n to level i
    //     Offset the w0 on level i above the top of level n
    //   }
    // }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w0.setVal(0.0);

    // local variables
    BaseState<Real> psi_planar_state(base_geom.nr_fine);
    auto psi_planar = psi_planar_state.array();

    const auto etarho_cc_arr = etarho_cc.const_array();
    auto w0_arr = w0.array();
    const auto w0_old_arr = w0_old.const_array();
    auto w0_force_arr = w0_force.array();
    const auto Sbar_arr = Sbar_in.const_array();
    const auto p0_old_arr = p0_old_in.const_array();
    const auto p0_new_arr = p0_new_in.const_array();
    const auto gamma1bar_old_arr = gamma1bar_old_in.const_array();
    const auto gamma1bar_new_arr = gamma1bar_new_in.const_array();
    const auto p0_minus_peosbar_arr = p0_minus_peosbar.const_array();

    const Real dt_loc = dt;
    const Real grav_const_loc = grav_const;
    const Real dpdt_factor_loc = dpdt_factor;

    // pressure correction variables
    BaseState<Real> int1_over_gamma1bar_p0(base_geom.nr_fine+1);
    auto int1_over_gamma1bar_p0_planar = int1_over_gamma1bar_p0.array();

    // Compute w0 on edges at level n
    for (auto n = 0; n <= base_geom.max_radial_level; ++n) {
        psi_planar_state.setVal(0.0);
        const int base_cutoff_density_coord =
            base_geom.base_cutoff_density_coord(n);

        const Real dr_lev = base_geom.dr(n);

        for (auto j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
            if (n == 0) {
                // Initialize new w0 at bottom of coarse base array to 0.0.
                w0_arr(0, 0) = 0.0;
                int1_over_gamma1bar_p0_planar(0) = 0.0;
            } else {
                // Obtain the starting value of w0 from the coarser grid
                w0_arr(n, base_geom.r_start_coord(n, j)) =
                    w0_arr(n - 1, base_geom.r_start_coord(n, j) / 2);
            }

            // compute psi for level n
            int lo = base_geom.r_start_coord(n, j);
            int hi = base_geom.r_end_coord(n, j);
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                int r = k + lo;
                if (r < base_cutoff_density_coord) {
                    psi_planar(r) =
                        etarho_cc_arr(n, r) * amrex::Math::abs(grav_const_loc);
                }
            });
            Gpu::synchronize();

            for (auto r = base_geom.r_start_coord(n, j) + 1;
                 r <= base_geom.r_end_coord(n, j) + 1; ++r) {
                Real gamma1bar_p0_avg =
                    (gamma1bar_old_arr(n, r - 1) +
                     gamma1bar_new_arr(n, r - 1)) *
                    (p0_old_arr(n, r - 1) + p0_new_arr(n, r - 1)) / 4.0;

                Real delta_chi_w0 = 0.0;

                if (r < base_cutoff_density_coord) {
                    if (is_predictor) {
                        delta_chi_w0 = dpdt_factor_loc *
                                       p0_minus_peosbar_arr(n, r - 1) /
                                       (gamma1bar_old_arr(n, r - 1) *
                                        p0_old_arr(n, r - 1) * dt_loc);
                    } else {
                        delta_chi_w0 += dpdt_factor_loc *
                                        p0_minus_peosbar_arr(n, r - 1) /
                                        (gamma1bar_new_arr(n, r - 1) *
                                         p0_new_arr(n, r - 1) * dt_loc);
                    }
                }
                if (n == 0) {
                    int1_over_gamma1bar_p0_planar(r) =
                        int1_over_gamma1bar_p0_planar(r - 1) +
                        (1.0 / gamma1bar_p0_avg) * dr_lev;
                }
                w0_arr(n, r) = w0_arr(n, r - 1) + Sbar_arr(n, r - 1) * dr_lev -
                               psi_planar[r - 1] / gamma1bar_p0_avg * dr_lev -
                               delta_chi_w0 * dr_lev;
            }
            // add the pressure correction for a closed box for n == 0
            if (n == 0 && add_pb) {
                const int k = base_geom.r_end_coord(n, j) + 1;
                p0bdot = w0_arr(n, k) / int1_over_gamma1bar_p0_planar(k);
                // set p0b for use in EnforceHSE
                p0b = p0bdot * dt;
                for (auto r = base_geom.r_start_coord(n, j) + 1;
                     r <= base_geom.r_end_coord(n, j) + 1; ++r) {
                    w0_arr(n, r) -= p0bdot * int1_over_gamma1bar_p0_planar(r);
                }
            }
            if (n > 0) {
                // Compare the difference between w0 at top of level n to
                // the corresponding point on level n-1
                Real offset =
                    w0_arr(n, base_geom.r_end_coord(n, j) + 1) -
                    w0_arr(n - 1, (base_geom.r_end_coord(n, j) + 1) / 2);

                for (auto i = n - 1; i >= 0; --i) {
                    auto refrat = (int)amrex::Math::round(pow(2, n - i));

                    // Restrict w0 from level n to level i
                    for (auto r = base_geom.r_start_coord(n, j);
                         r <= base_geom.r_end_coord(n, j) + 1; ++r) {
                        if (r % refrat == 0) {
                            w0_arr(n, r / refrat) = w0_arr(n, r);
                        }
                    }

                    // Offset the w0 on level i above the top of level n
                    lo = (base_geom.r_end_coord(n, j) + 1) / refrat + 1;
                    hi = base_geom.nr(i);
                    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                        int r = k + lo;
                        w0_arr(i, r) += offset;
                    });
                    Gpu::synchronize();
                }
            }
        }
    }

    // zero w0 where there is no corresponding full state array
    for (auto n = 1; n <= base_geom.max_radial_level; ++n) {
        for (auto j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
            if (j == base_geom.numdisjointchunks(n)) {
                const int lo = base_geom.r_end_coord(n, j) + 2;
                const int hi = base_geom.nr(n);
                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                    int r = k + lo;
                    w0_arr(n, r) = 0.0;
                });
            } else {
                const int lo = base_geom.r_end_coord(n, j) + 2;
                const int hi = base_geom.r_start_coord(n, j + 1) - 1;
                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                    int r = k + lo;
                    w0_arr(n, r) = 0.0;
                });
            }
            Gpu::synchronize();
        }
    }

    RestrictBase(w0, false);
    FillGhostBase(w0, false);

    for (auto n = 0; n <= base_geom.max_radial_level; ++n) {
        for (auto j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            const Real dt_avg = 0.5 * (dt_in + dtold_in);
            const Real dr_lev = base_geom.dr(n);

            const int lo = base_geom.r_start_coord(n, j);
            const int hi = base_geom.r_end_coord(n, j);
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                int r = k + lo;

                Real w0_old_cen =
                    0.5 * (w0_old_arr(n, r) + w0_old_arr(n, r + 1));
                Real w0_new_cen = 0.5 * (w0_arr(n, r) + w0_arr(n, r + 1));
                Real w0_avg =
                    0.5 * (dt_in * w0_old_cen + dtold_in * w0_new_cen) / dt_avg;
                Real div_avg =
                    0.5 *
                    (dt_in * (w0_old_arr(n, r + 1) - w0_old_arr(n, r)) +
                     dtold_in * (w0_arr(n, r + 1) - w0_arr(n, r))) /
                    dt_avg;
                w0_force_arr(n, r) = (w0_new_cen - w0_old_cen) / dt_avg +
                                     w0_avg * div_avg / dr_lev;
            });
            Gpu::synchronize();
        }
    }

    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}

void Maestro::Makew0PlanarVarg(
    const BaseState<Real>& w0_old, BaseState<Real>& w0_force,
    const BaseState<Real>& Sbar_in, const BaseState<Real>& rho0_old_in,
    const BaseState<Real>& rho0_new_in, const BaseState<Real>& p0_old_in,
    const BaseState<Real>& p0_new_in, const BaseState<Real>& gamma1bar_old_in,
    const BaseState<Real>& gamma1bar_new_in,
    const BaseState<Real>& p0_minus_peosbar, const Real dt_in,
    const Real dtold_in) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0PlanarVarg()", Makew0PlanarVarg);

    const auto fine_base_density_cutoff_coord =
        base_geom.base_cutoff_density_coord(base_geom.finest_radial_level);

    const int nr_finest = base_geom.nr(base_geom.finest_radial_level);
    const Real dr_finest = base_geom.dr(base_geom.finest_radial_level);
    const Real dpdt_factor_loc = dpdt_factor;

    auto w0_arr = w0.array();
    auto w0_force_arr = w0_force.array();
    const auto w0_old_arr = w0_old.const_array();
    const auto r_edge_loc = base_geom.r_edge_loc;

    // The planar 1/r**2 gravity constraint equation is solved
    // by calling the tridiagonal solver, just like spherical.
    // This is accomplished by putting all the requisite data
    // on the finest basestate grid, solving for w0, and then
    // restricting w0 back down to the coarse grid.

    // 1) allocate the finely-gridded temporary basestate arrays
    BaseState<Real> w0_fine(nr_finest + 1);
    BaseState<Real> w0bar_fine(nr_finest + 1);
    BaseState<Real> deltaw0_fine(nr_finest + 1);
    BaseState<Real> p0_old_fine(nr_finest);
    BaseState<Real> p0_new_fine(nr_finest);
    BaseState<Real> p0_nph_fine(nr_finest);
    BaseState<Real> rho0_old_fine(nr_finest);
    BaseState<Real> rho0_new_fine(nr_finest);
    BaseState<Real> rho0_nph_fine(nr_finest);
    BaseState<Real> gamma1bar_old_fine(nr_finest);
    BaseState<Real> gamma1bar_new_fine(nr_finest);
    BaseState<Real> gamma1bar_nph_fine(nr_finest);
    BaseState<Real> p0_minus_peosbar_fine(nr_finest);
    BaseState<Real> etarho_cc_fine(nr_finest);
    BaseState<Real> Sbar_in_fine(nr_finest);
    BaseState<Real> grav_edge_fine(nr_finest + 1);

    // 2) copy the data into the temp, uniformly-gridded basestate arrays.
    ProlongBasetoUniform(p0_old_in, p0_old_fine);
    ProlongBasetoUniform(p0_new_in, p0_new_fine);
    ProlongBasetoUniform(rho0_old_in, rho0_old_fine);
    ProlongBasetoUniform(rho0_new_in, rho0_new_fine);
    ProlongBasetoUniform(gamma1bar_old_in, gamma1bar_old_fine);
    ProlongBasetoUniform(gamma1bar_new_in, gamma1bar_new_fine);
    ProlongBasetoUniform(p0_minus_peosbar, p0_minus_peosbar_fine);
    ProlongBasetoUniform(etarho_cc, etarho_cc_fine);
    ProlongBasetoUniform(Sbar_in, Sbar_in_fine);

    auto p0_nph_fine_arr = p0_nph_fine.array();
    auto w0_fine_arr = w0_fine.array();
    auto deltaw0_fine_arr = deltaw0_fine.array();
    auto w0bar_fine_arr = w0bar_fine.array();
    auto gamma1bar_nph_fine_arr = gamma1bar_nph_fine.array();
    auto p0_minus_peosbar_fine_arr = p0_minus_peosbar_fine.array();
    auto etarho_cc_fine_arr = etarho_cc_fine.array();
    auto Sbar_in_fine_arr = Sbar_in_fine.array();
    auto grav_edge_fine_arr = grav_edge_fine.array();
    auto finest_radial_level = base_geom.finest_radial_level;

    // create time-centered base-state quantities
    p0_nph_fine.copy(0.5 * (p0_old_fine + p0_new_fine));
    rho0_nph_fine.copy(0.5 * (rho0_old_fine + rho0_new_fine));
    gamma1bar_nph_fine.copy(0.5 * (gamma1bar_old_fine + gamma1bar_new_fine));

    // 3) solve to w0bar -- here we just take into account the Sbar and
    //    volume discrepancy terms
    // lower boundary condition

    w0bar_fine_arr(0) = 0.0;

    int lo = 1;
    int hi = nr_finest;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        Real gamma1bar_p0_avg =
            gamma1bar_nph_fine_arr(r - 1) * p0_nph_fine_arr(r - 1);

        Real volume_discrepancy =
            (r - 1 < fine_base_density_cutoff_coord)
                ? dpdt_factor_loc * p0_minus_peosbar_fine_arr(r - 1) / dt_in
                : 0.0;

        w0bar_fine_arr(r) = w0bar_fine_arr(r - 1) +
                            Sbar_in_fine_arr(r - 1) * dr_finest -
                            (volume_discrepancy / gamma1bar_p0_avg) * dr_finest;
    });
    Gpu::synchronize();

    // 4) get the edge-centered gravity on the uniformly-gridded
    // basestate arrays
    Abort("make_w0.f90: need to write make_grav_edge_uniform");
    //    call make_grav_edge_uniform(grav_edge_fine, rho0_nph_fine)

    // 5) solve for delta w0
    deltaw0_fine.setVal(0.0);

    // this takes the form of a tri-diagonal matrix:
    // A_j (dw_0)_{j-3/2} +
    // B_j (dw_0)_{j-1/2} +
    // C_j (dw_0)_{j+1/2} = F_j

    BaseState<Real> A_s(nr_finest + 1);
    BaseState<Real> B_s(nr_finest + 1);
    BaseState<Real> C_s(nr_finest + 1);
    BaseState<Real> u_s(nr_finest + 1);
    BaseState<Real> F_s(nr_finest + 1);

    A_s.setVal(0.0);
    B_s.setVal(0.0);
    C_s.setVal(0.0);
    F_s.setVal(0.0);
    u_s.setVal(0.0);

    auto A = A_s.array();
    auto B = B_s.array();
    auto C = C_s.array();
    auto u = u_s.array();
    auto F = A_s.array();

    lo = 1;
    hi = fine_base_density_cutoff_coord;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        A(r) = gamma1bar_nph_fine_arr(r - 1) * p0_nph_fine_arr(r - 1);
        A(r) /= dr_finest * dr_finest;

        Real dpdr = (p0_nph_fine_arr(r) - p0_nph_fine_arr(r - 1)) / dr_finest;

        B(r) = -(gamma1bar_nph_fine_arr(r - 1) * p0_nph_fine_arr(r - 1) +
                 gamma1bar_nph_fine_arr(r) * p0_nph_fine_arr(r)) /
               (dr_finest * dr_finest);
        B(r) -= 2.0 * dpdr / (r_edge_loc(finest_radial_level, r));

        C(r) = gamma1bar_nph_fine_arr(r) * p0_nph_fine_arr(r);
        C(r) /= dr_finest * dr_finest;

        F(r) = 2.0 * dpdr * w0bar_fine_arr(r) /
                   r_edge_loc(finest_radial_level, r) -
               grav_edge_fine_arr(r) *
                   (etarho_cc_fine_arr(r) - etarho_cc_fine_arr(r - 1)) /
                   dr_finest;
    });
    Gpu::synchronize();

    // Lower boundary
    A(0) = 0.0;
    B(0) = 1.0;
    C(0) = 0.0;
    F(0) = 0.0;

    // Upper boundary
    A(fine_base_density_cutoff_coord + 1) = -1.0;
    B(fine_base_density_cutoff_coord + 1) = 1.0;
    C(fine_base_density_cutoff_coord + 1) = 0.0;
    F(fine_base_density_cutoff_coord + 1) = 0.0;

    // need to synchronize gpu values with updated host values
    Gpu::synchronize();

    // Call the tridiagonal solver
    Tridiag(A, B, C, F, u, fine_base_density_cutoff_coord + 2);

    lo = 1;
    hi = fine_base_density_cutoff_coord + 1;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        deltaw0_fine_arr(r) = u(r);
    });
    Gpu::synchronize();

    lo = fine_base_density_cutoff_coord + 2;
    hi = nr_finest;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        deltaw0_fine_arr(r) =
            deltaw0_fine_arr(fine_base_density_cutoff_coord + 1);
    });
    Gpu::synchronize();

    // 6) compute w0 = w0bar + deltaw0
    ParallelFor(nr_finest + 1, [=] AMREX_GPU_DEVICE(int r) {
        w0_fine_arr(r) = w0bar_fine_arr(r) + deltaw0_fine_arr(r);
        w0_arr(finest_radial_level, r) = w0_fine_arr(r);
    });
    Gpu::synchronize();

    // 7) fill the multilevel w0 array from the uniformly-gridded w0 we
    // just solved for.  Here, we make the coarse edge underneath equal
    // to the fine edge value.
    for (auto n = base_geom.finest_radial_level; n >= 1; --n) {
        for (auto r = 0; r <= base_geom.nr(n); n += 2) {
            w0_arr(n - 1, r / 2) = w0_arr(n, r);
        }
    }

    // 8) zero w0 where there is no corresponding full state array
    for (auto n = 1; n <= base_geom.finest_radial_level; ++n) {
        for (auto j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
            if (j == base_geom.numdisjointchunks(n)) {
                lo = base_geom.r_end_coord(n, j) + 2;
                hi = base_geom.nr(n);
                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                    int r = k + lo;
                    w0_arr(n, r) = 0.0;
                });
            } else {
                lo = base_geom.r_end_coord(n, j) + 2;
                hi = base_geom.r_start_coord(n, j + 1);
                ParallelFor(hi - lo, [=] AMREX_GPU_DEVICE(int k) {
                    int r = k + lo;
                    w0_arr(n, r) = 0.0;
                });
            }
            Gpu::synchronize();
        }
    }

    RestrictBase(w0, false);
    FillGhostBase(w0, false);

    // compute the forcing terms
    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        for (auto j = 1; j <= base_geom.numdisjointchunks(n); ++j) {
            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            const Real dt_avg = 0.5 * (dt_in + dtold_in);
            const Real dr_lev = base_geom.dr(n);

            lo = base_geom.r_start_coord(n, j);
            hi = base_geom.r_end_coord(n, j);
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int k) {
                int r = k + lo;
                Real w0_old_cen =
                    0.5 * (w0_old_arr(n, r) + w0_old_arr(n, r + 1));
                Real w0_new_cen = 0.5 * (w0_arr(n, r) + w0_arr(n, r + 1));
                Real w0_avg =
                    0.5 * (dt_in * w0_old_cen + dtold_in * w0_new_cen) / dt_avg;
                Real div_avg =
                    0.5 *
                    (dt_in * (w0_old_arr(n, r + 1) - w0_old_arr(n, r)) +
                     dtold_in * (w0_arr(n, r + 1) - w0_arr(n, r))) /
                    dt_avg;
                w0_force_arr(n, r) = (w0_new_cen - w0_old_cen) / dt_avg +
                                     w0_avg * div_avg / dr_lev;
            });
            Gpu::synchronize();
        }
    }

    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}

void Maestro::Makew0Sphr(
    const BaseState<Real>& w0_old, BaseState<Real>& w0_force,
    const BaseState<Real>& Sbar_in, const BaseState<Real>& rho0_old_in,
    const BaseState<Real>& rho0_new_in, const BaseState<Real>& p0_old_in,
    const BaseState<Real>& p0_new_in, const BaseState<Real>& gamma1bar_old_in,
    const BaseState<Real>& gamma1bar_new_in,
    const BaseState<Real>& p0_minus_peosbar, const Real dt_in,
    const Real dtold_in) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0Sphr()", Makew0Sphr);

    // local variables
    const int max_lev = base_geom.max_radial_level + 1;
    BaseState<Real> gamma1bar_nph_s(base_geom.nr_fine);
    BaseState<Real> p0_nph_s(base_geom.nr_fine);
    BaseState<Real> A_s(base_geom.nr_fine + 1);
    BaseState<Real> B_s(base_geom.nr_fine + 1);
    BaseState<Real> C_s(base_geom.nr_fine + 1);
    BaseState<Real> u_s(base_geom.nr_fine + 1);
    BaseState<Real> F_s(base_geom.nr_fine + 1);
    BaseState<Real> w0_from_Sbar_s(base_geom.nr_fine + 1);
    BaseState<Real> rho0_nph_s(max_lev, base_geom.nr_fine);
    BaseState<Real> grav_edge_s(max_lev, base_geom.nr_fine + 1);

    auto gamma1bar_nph = gamma1bar_nph_s.array();
    auto p0_nph = p0_nph_s.array();
    auto A = A_s.array();
    auto B = B_s.array();
    auto C = C_s.array();
    auto u = u_s.array();
    auto F = F_s.array();
    auto w0_from_Sbar = w0_from_Sbar_s.array();
    auto rho0_nph = rho0_nph_s.array();
    auto grav_edge = grav_edge_s.array();

    const auto Sbar_arr = Sbar_in.const_array();
    const auto p0_old_arr = p0_old_in.const_array();
    const auto p0_new_arr = p0_new_in.const_array();
    const auto p0_minus_peosbar_arr = p0_minus_peosbar.const_array();
    const auto gamma1bar_old_arr = gamma1bar_old_in.const_array();
    const auto gamma1bar_new_arr = gamma1bar_new_in.const_array();
    const auto rho0_old_arr = rho0_old_in.const_array();
    const auto rho0_new_arr = rho0_new_in.const_array();
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    const auto etarho_cc_arr = etarho_cc.const_array();
    const auto etarho_ec_arr = etarho_ec.const_array();
    auto w0_arr = w0.array();
    const auto w0_old_arr = w0_old.const_array();
    auto w0_force_arr = w0_force.array();

    const auto base_cutoff_density_coord =
        base_geom.base_cutoff_density_coord(0);

    const Real dr0 = base_geom.dr(0);
    const Real dpdt_factor_loc = dpdt_factor;

    // create time-centered base-state quantities
    ParallelFor(base_geom.nr_fine, [=] AMREX_GPU_DEVICE(int r) {
        p0_nph(r) = 0.5 * (p0_old_arr(0, r) + p0_new_arr(0, r));
        rho0_nph(0, r) = 0.5 * (rho0_old_arr(0, r) + rho0_new_arr(0, r));
        gamma1bar_nph(r) =
            0.5 * (gamma1bar_old_arr(0, r) + gamma1bar_new_arr(0, r));
    });
    Gpu::synchronize();

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.
    w0_from_Sbar(0) = 0.0;

    for (auto r = 1; r <= base_geom.nr_fine; ++r) {
        Real volume_discrepancy =
            rho0_old_arr(0, r - 1) > base_cutoff_density
                ? dpdt_factor_loc * p0_minus_peosbar_arr(0, r - 1) / dt_in
                : 0.0;

        w0_from_Sbar(r) = w0_from_Sbar(r - 1) + dr0 * Sbar_arr(0, r - 1) *
                                                    r_cc_loc(0, r - 1) *
                                                    r_cc_loc(0, r - 1);
        if (volume_discrepancy != 0.0) {
            w0_from_Sbar(r) -= dr0 * volume_discrepancy * r_cc_loc(0, r - 1) *
                               r_cc_loc(0, r - 1) /
                               (gamma1bar_nph(r - 1) * p0_nph(r - 1));
        }
    }

    int lo = 1;
    int hi = base_geom.nr_fine;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        w0_from_Sbar(r) /= (r_edge_loc(0, r) * r_edge_loc(0, r));
    });
    Gpu::synchronize();

    // make the edge-centered gravity
    MakeGravEdge(grav_edge_s, rho0_nph_s);

    // NOTE:  now we solve for the remainder, (r^2 * delta w0)
    // this takes the form of a tri-diagonal matrix:
    // A_j (r^2 dw_0)_{j-3/2} +
    // B_j (r^2 dw_0)_{j-1/2} +
    // C_j (r^2 dw_0)_{j+1/2} = F_j
    A_s.setVal(0.0);
    B_s.setVal(0.0);
    C_s.setVal(0.0);
    F_s.setVal(0.0);
    u_s.setVal(0.0);

    // Note that we are solving for (r^2 delta w0), not just w0.

    int max_cutoff =
        amrex::min(base_cutoff_density_coord, base_geom.nr_fine - 1);

    lo = 1;
    hi = max_cutoff;
    amrex::ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) noexcept {
        int r = j + lo;
        A(r) = gamma1bar_nph(r - 1) * p0_nph(r - 1) /
               (r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1));
        A(r) /= dr0 * dr0;

        B(r) = -(gamma1bar_nph(r - 1) * p0_nph(r - 1) /
                     (r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1)) +
                 gamma1bar_nph(r) * p0_nph(r) /
                     (r_cc_loc(0, r) * r_cc_loc(0, r))) /
               (dr0 * dr0);

        Real dpdr = (p0_nph(r) - p0_nph(r - 1)) / dr0;

        B(r) -= 4.0 * dpdr /
                (r_edge_loc(0, r) * r_edge_loc(0, r) * r_edge_loc(0, r));

        C(r) = gamma1bar_nph(r) * p0_nph(r) / (r_cc_loc(0, r) * r_cc_loc(0, r));
        C(r) /= dr0 * dr0;

        F(r) = 4.0 * dpdr * w0_from_Sbar(r) / r_edge_loc(0, r) -
               grav_edge(0, r) *
                   (r_cc_loc(0, r) * r_cc_loc(0, r) * etarho_cc_arr(0, r) -
                    r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1) *
                        etarho_cc_arr(0, r - 1)) /
                   (dr0 * r_edge_loc(0, r) * r_edge_loc(0, r)) -
               4.0 * M_PI * Gconst * 0.5 *
                   (rho0_nph(0, r) + rho0_nph(0, r - 1)) * etarho_ec_arr(0, r);
    });
    Gpu::synchronize();

    // Lower boundary
    A(0) = 0.0;
    B(0) = 1.0;
    C(0) = 0.0;
    F(0) = 0.0;

    // Upper boundary
    A(max_cutoff + 1) = -1.0;
    B(max_cutoff + 1) = 1.0;
    C(max_cutoff + 1) = 0.0;
    F(max_cutoff + 1) = 0.0;

    // need to synchronize gpu values with updated host values
    Gpu::synchronize();

    // Call the tridiagonal solver
    Tridiag(A, B, C, F, u, max_cutoff + 2);

    w0_arr(0, 0) = w0_from_Sbar(0);

    lo = 1;
    hi = max_cutoff + 1;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        w0_arr(0, r) =
            u(r) / (r_edge_loc(0, r) * r_edge_loc(0, r)) + w0_from_Sbar(r);
    });
    Gpu::synchronize();

    lo = max_cutoff + 2;
    hi = base_geom.nr_fine;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        w0_arr(0, r) = w0_arr(0, max_cutoff + 1) *
                       r_edge_loc(0, max_cutoff + 1) *
                       r_edge_loc(0, max_cutoff + 1) /
                       (r_edge_loc(0, r) * r_edge_loc(0, r));
    });
    Gpu::synchronize();

    // Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    const Real dt_avg = 0.5 * (dt_in + dtold_in);

    ParallelFor(base_geom.nr_fine, [=] AMREX_GPU_DEVICE(int r) {
        Real w0_old_cen = 0.5 * (w0_old_arr(0, r) + w0_old_arr(0, r + 1));
        Real w0_new_cen = 0.5 * (w0_arr(0, r) + w0_arr(0, r + 1));
        Real w0_avg =
            0.5 * (dt_in * w0_old_cen + dtold_in * w0_new_cen) / dt_avg;
        Real div_avg = 0.5 *
                       (dt_in * (w0_old_arr(0, r + 1) - w0_old_arr(0, r)) +
                        dtold_in * (w0_arr(0, r + 1) - w0_arr(0, r))) /
                       dt_avg;
        w0_force_arr(0, r) =
            (w0_new_cen - w0_old_cen) / dt_avg + w0_avg * div_avg / dr0;
    });
    Gpu::synchronize();
}

void Maestro::Makew0SphrIrreg(
    const BaseState<Real>& w0_old, BaseState<Real>& w0_force,
    const BaseState<Real>& Sbar_in, const BaseState<Real>& rho0_old_in,
    const BaseState<Real>& rho0_new_in, const BaseState<Real>& p0_old_in,
    const BaseState<Real>& p0_new_in, const BaseState<Real>& gamma1bar_old_in,
    const BaseState<Real>& gamma1bar_new_in,
    const BaseState<Real>& p0_minus_peosbar, const Real dt_in,
    const Real dtold_in) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0SphrIrreg()", Makew0SphrIrreg);

    // local variables
    const int max_lev = base_geom.max_radial_level + 1;
    BaseState<Real> gamma1bar_nph_s(base_geom.nr_fine);
    BaseState<Real> p0_nph_s(base_geom.nr_fine);
    BaseState<Real> A_s(base_geom.nr_fine + 1);
    BaseState<Real> B_s(base_geom.nr_fine + 1);
    BaseState<Real> C_s(base_geom.nr_fine + 1);
    BaseState<Real> u_s(base_geom.nr_fine + 1);
    BaseState<Real> F_s(base_geom.nr_fine + 1);
    BaseState<Real> w0_from_Sbar_s(base_geom.nr_fine + 1);
    BaseState<Real> rho0_nph_s(max_lev, base_geom.nr_fine);
    BaseState<Real> grav_edge_s(max_lev, base_geom.nr_fine + 1);

    auto gamma1bar_nph = gamma1bar_nph_s.array();
    auto p0_nph = p0_nph_s.array();
    auto A = A_s.array();
    auto B = B_s.array();
    auto C = C_s.array();
    auto u = u_s.array();
    auto F = F_s.array();
    auto w0_from_Sbar = w0_from_Sbar_s.array();
    auto rho0_nph = rho0_nph_s.array();
    auto grav_edge = grav_edge_s.array();

    const auto Sbar_arr = Sbar_in.const_array();
    const auto p0_old_arr = p0_old_in.const_array();
    const auto p0_new_arr = p0_new_in.const_array();
    const auto p0_minus_peosbar_arr = p0_minus_peosbar.const_array();
    const auto gamma1bar_old_arr = gamma1bar_old_in.const_array();
    const auto gamma1bar_new_arr = gamma1bar_new_in.const_array();
    const auto rho0_old_arr = rho0_old_in.const_array();
    const auto rho0_new_arr = rho0_new_in.const_array();
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    const auto etarho_cc_arr = etarho_cc.const_array();
    const auto etarho_ec_arr = etarho_ec.const_array();
    auto w0_arr = w0.array();
    const auto w0_old_arr = w0_old.const_array();
    auto w0_force_arr = w0_force.array();

    const auto base_cutoff_density_coord =
        base_geom.base_cutoff_density_coord(0);
    const Real dpdt_factor_loc = dpdt_factor;

    // create time-centered base-state quantities
    ParallelFor(base_geom.nr_fine, [=] AMREX_GPU_DEVICE(int r) {
        p0_nph(r) = 0.5 * (p0_old_arr(0, r) + p0_new_arr(0, r));
        rho0_nph(r) = 0.5 * (rho0_old_arr(0, r) + rho0_new_arr(0, r));
        gamma1bar_nph(r) =
            0.5 * (gamma1bar_old_arr(0, r) + gamma1bar_new_arr(0, r));
    });
    Gpu::synchronize();

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.
    w0_from_Sbar(0) = 0.0;

    for (auto r = 1; r <= base_geom.nr_fine; ++r) {
        Real volume_discrepancy =
            rho0_old_arr(0, r - 1) > base_cutoff_density
                ? dpdt_factor_loc * p0_minus_peosbar_arr(0, r - 1) / dt_in
                : 0.0;

        Real dr1 = r_edge_loc(0, r) - r_edge_loc(0, r - 1);
        w0_from_Sbar(r) =
            w0_from_Sbar(r - 1) +
            dr1 * Sbar_arr(0, r - 1) * r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1) -
            dr1 * volume_discrepancy * r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1) /
                (gamma1bar_nph(r - 1) * p0_nph(r - 1));
    }

    for (auto r = 1; r <= base_geom.nr_fine; ++r) {
        w0_from_Sbar(r) /= (r_edge_loc(0, r) * r_edge_loc(0, r));
    }

    // make the edge-centered gravity
    MakeGravEdge(grav_edge_s, rho0_nph_s);

    // NOTE:  now we solve for the remainder, (r^2 * delta w0)
    // this takes the form of a tri-diagonal matrix:
    // A_j (r^2 dw_0)_{j-3/2} +
    // B_j (r^2 dw_0)_{j-1/2} +
    // C_j (r^2 dw_0)_{j+1/2} = F_j
    A_s.setVal(0.0);
    B_s.setVal(0.0);
    C_s.setVal(0.0);
    u_s.setVal(0.0);
    F_s.setVal(0.0);

    // Note that we are solving for (r^2 delta w0), not just w0.
    int max_cutoff = base_cutoff_density_coord;

    int lo = 1;
    int hi = max_cutoff;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        Real dr1 = r_edge_loc(0, r) - r_edge_loc(0, r - 1);
        Real dr2 = r_edge_loc(0, r + 1) - r_edge_loc(0, r);
        Real dr3 = r_cc_loc(0, r) - r_cc_loc(0, r - 1);

        A(r) = gamma1bar_nph(r - 1) * p0_nph(r - 1) /
               (r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1));
        A(r) /= dr1 * dr3;

        B(r) = -(gamma1bar_nph(r - 1) * p0_nph(r - 1) /
                     (r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1) * dr1) +
                 gamma1bar_nph(r) * p0_nph(r) /
                     (r_cc_loc(0, r) * r_cc_loc(0, r) * dr2)) /
               dr3;

        Real dpdr = (p0_nph(r) - p0_nph(r - 1)) / dr3;

        B(r) -= 4.0 * dpdr /
                (r_edge_loc(0, r) * r_edge_loc(0, r) * r_edge_loc(0, r));

        C(r) = gamma1bar_nph(r) * p0_nph(r) / (r_cc_loc(0, r) * r_cc_loc(0, r));
        C(r) /= dr2 * dr3;

        F(r) = 4.0 * dpdr * w0_from_Sbar(r) / r_edge_loc(0, r) -
               grav_edge(0, r) *
                   (r_cc_loc(0, r) * r_cc_loc(0, r) * etarho_cc_arr(0, r) -
                    r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1) *
                        etarho_cc_arr(0, r - 1)) /
                   (dr3 * r_edge_loc(0, r) * r_edge_loc(0, r)) -
               4.0 * M_PI * Gconst * 0.5 *
                   (rho0_nph(0, r) + rho0_nph(0, r - 1)) * etarho_ec_arr(0, r);
    });
    Gpu::synchronize();

    // Lower boundary
    A(0) = 0.0;
    B(0) = 1.0;
    C(0) = 0.0;
    F(0) = 0.0;

    // Upper boundary
    A(max_cutoff + 1) = -1.0;
    B(max_cutoff + 1) = 1.0;
    C(max_cutoff + 1) = 0.0;
    F(max_cutoff + 1) = 0.0;

    // need to synchronize gpu values with updated host values
    Gpu::synchronize();

    // Call the tridiagonal solver
    Tridiag(A, B, C, F, u, max_cutoff + 2);

    w0_arr(0, 0) = w0_from_Sbar(0);

    lo = 1;
    hi = max_cutoff + 1;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        w0_arr(0, r) =
            u(r) / (r_edge_loc(0, r) * r_edge_loc(0, r)) + w0_from_Sbar(r);
    });
    Gpu::synchronize();

    lo = max_cutoff + 2;
    hi = base_geom.nr_fine;
    ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
        int r = j + lo;
        w0_arr(0, r) = w0_arr(0, max_cutoff + 1) *
                       r_edge_loc(0, max_cutoff + 1) *
                       r_edge_loc(0, max_cutoff + 1) /
                       (r_edge_loc(0, r) * r_edge_loc(0, r));
    });
    Gpu::synchronize();

    // Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    const Real dt_avg = 0.5 * (dt_in + dtold_in);

    ParallelFor(base_geom.nr_fine, [=] AMREX_GPU_DEVICE(int r) {
        Real dr1 = r_edge_loc(0, r) - r_edge_loc(0, r - 1);
        Real w0_old_cen = 0.5 * (w0_old_arr(0, r) + w0_old_arr(0, r + 1));
        Real w0_new_cen = 0.5 * (w0_arr(0, r) + w0_arr(0, r + 1));
        Real w0_avg =
            0.5 * (dt_in * w0_old_cen + dtold_in * w0_new_cen) / dt_avg;
        Real div_avg = 0.5 *
                       (dt_in * (w0_old_arr(0, r + 1) - w0_old_arr(0, r)) +
                        dtold_in * (w0_arr(0, r + 1) - w0_arr(0, r))) /
                       dt_avg;
        w0_force_arr(0, r) =
            (w0_new_cen - w0_old_cen) / dt_avg + w0_avg * div_avg / dr1;
    });
    Gpu::synchronize();
}

void Maestro::Tridiag(const BaseStateArray<Real>& a,
                      const BaseStateArray<Real>& b,
                      const BaseStateArray<Real>& c,
                      const BaseStateArray<Real>& r,
                      const BaseStateArray<Real>& u, const int n) {
    BaseState<Real> gam_s(n);
    auto gam = gam_s.array();

    if (b(0) == 0) {
        Abort("tridiag: CANT HAVE B(0) = 0.0");
    }

    Real bet = b(0);
    u(0) = r(0) / bet;

    for (auto j = 1; j < n; j++) {
        gam(j) = c(j - 1) / bet;
        bet = b(j) - a(j) * gam(j);
        if (bet == 0) {
            Abort("tridiag: TRIDIAG FAILED");
        }
        u(j) = (r(j) - a(j) * u(j - 1)) / bet;
    }

    for (auto j = n - 2; j >= 0; --j) {
        u(j) -= gam(j + 1) * u(j + 1);
    }
}

void Maestro::ProlongBasetoUniform(const BaseState<Real>& base_ml_s,
                                   BaseState<Real>& base_fine_s)

{
    // the mask array will keep track of whether we've filled in data
    // in a corresponding radial bin.  .false. indicates that we've
    // already output there.
    IntVector imask_fine(base_geom.nr_fine);
    std::fill(imask_fine.begin(), imask_fine.end(), 1);

    // r1 is the factor between the current level grid spacing and the
    // FINEST level
    int r1 = 1;

    const auto base_ml = base_ml_s.const_array();
    auto base_fine = base_fine_s.array();

    for (auto n = base_geom.finest_radial_level; n >= 0; --n) {
        for (auto j = 1; j < base_geom.numdisjointchunks(n); ++j) {
            for (auto r = base_geom.r_start_coord(n, j);
                 r <= base_geom.r_end_coord(n, j); ++r) {
                // sum up mask to see if there are any elements set to true
                if (std::accumulate(imask_fine.begin() + r * r1 - 1,
                                    imask_fine.begin() + (r + 1) * r1 - 1,
                                    0) > 0) {
                    for (auto i = r * r1 - 1; i < (r + 1) * r1 - 1; ++r) {
                        base_fine(i) = base_ml(n, r);
                        imask_fine[i] = 0;
                    }
                }
            }
        }
        // update r1 for the next coarsest level -- assume a jump by
        // factor of 2
        r1 *= 2;
    }

    // check to make sure that no mask values are still true
    if (std::accumulate(imask_fine.begin(), imask_fine.end(), 0) > 0) {
        Abort("ERROR: unfilled cells in prolong_base_to_uniform");
    }
}
