
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Evolve()", Evolve);

    Print() << "Calling Evolve()" << std::endl;

    for (istep = start_step; istep <= max_step && t_old < stop_time; ++istep) {
        dtold = dt;

        // compute time step
        // if this is the first time step we already have a dt from either FirstDt()
        // or EstDt called during the divu_iters
        if (istep > 1) {
            EstDt();

            if (verbose > 0) {
                Print() << "Call to estdt at beginning of step " << istep
                        << " gives dt =" << dt << std::endl;
            }

            // fixme - add nuclear_dt_scalefac timestep limiter

            if (dt > max_dt_growth * dtold) {
                dt = max_dt_growth * dtold;
                if (verbose > 0) {
                    Print() << "dt_growth factor limits the new dt = " << dt
                            << std::endl;
                }
            }

            if (dt > max_dt) {
                if (verbose > 0) {
                    Print()
                        << "max_dt limits the new dt = " << max_dt << std::endl;
                }
                dt = max_dt;
            }

            if (stop_time >= 0. && t_old + dt > stop_time) {
                dt = amrex::min(dt, stop_time - t_old);
                Print() << "Stop time limits dt = " << dt << std::endl;
            }

            t_new = t_old + dt;
        }

        // advance the solution by dt

        AdvanceTimeStep(false);

        t_old = t_new;

        // move new state into old state by swapping pointers
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::swap(S_cc_old[lev], S_cc_new[lev]);
        }

        std::swap(rho0_old, rho0_new);
        rhoh0_old.swap(rhoh0_new);
        p0_nm1.swap(p0_old);
        p0_old.swap(p0_new);

        gamma1bar_old.swap(gamma1bar_new);
        grav_cell_old.swap(grav_cell_new);
    }

    // Now need to check the HSE-ness
    Real max_hse_error = -1.e30;

    const auto rho0_old_arr = rho0_old.const_array();
    const auto p0_old_arr = p0_old.const_array();

    const auto& dr = base_geom.dr;
    const Real starting_rad = (spherical) ? 0.0 : geom[0].ProbLo(0);

    for (int n = 0; n <= base_geom.max_radial_level; ++n) {
        Real mencl = 0.0;

        if (use_exact_base_state) {
            Real dr_irreg = base_geom.r_edge_loc(n, 1) -
                            base_geom.r_edge_loc(n, 0);  // edge-to-edge

            if (spherical || do_2d_planar_octant) {
                mencl = 4.0 / 3.0 * M_PI * dr_irreg * dr_irreg * dr_irreg *
                        rho0_old_arr(n, 0);
            }
        } else {
            if (spherical || do_2d_planar_octant) {
                mencl = 4.0 / 3.0 * M_PI * dr(n) * dr(n) * dr(n) *
                        rho0_old_arr(n, 0);
            }
        }

        for (auto r = 1; r < base_geom.nr(n); ++r) {
            Real rloc = use_exact_base_state
                            ? base_geom.r_cc_loc(n, r)
                            : starting_rad + (Real(r) + 0.5) * dr(n);

            if (rloc < base_geom.base_cutoff_density_coord(n)) {
                Real r_r = starting_rad;
                r_r += use_exact_base_state ? base_geom.r_edge_loc(n, r + 1)
                                            : Real(r + 1) * dr(n);
                Real r_l = starting_rad;
                r_l += use_exact_base_state ? base_geom.r_edge_loc(n, r)
                                            : Real(r) * dr(n);

                Real dr_local = r_r - r_l;

                Real g = 0.0;

                if (spherical || do_2d_planar_octant) {
                    g = -Gconst * mencl / (r_l * r_l);
                    mencl += 4.0 / 3.0 * M_PI * dr_local *
                             (r_l * r_l + r_l * r_r + r_r * r_r) *
                             rho0_old_arr(n, r);
                } else {
                    if (!do_planar_invsq_grav) {
                        g = grav_const;
                    } else {
                        g = -Gconst * planar_invsq_mass / (r_l * r_l);
                    }
                }

                Real dpdr = 0.0;
                Real rhog = 0.0;
                if (use_exact_base_state) {
                    Real dr_irreg =
                        base_geom.r_cc_loc(n, r) - base_geom.r_cc_loc(n, r - 1);
                    dpdr = (p0_old_arr(n, r) - p0_old_arr(n, r - 1)) / dr_irreg;

                    Real rfrac = (base_geom.r_edge_loc(n, r) -
                                  base_geom.r_cc_loc(n, r - 1)) /
                                 dr_irreg;
                    rhog = ((1.0 - rfrac) * rho0_old_arr(n, r) +
                            rfrac * rho0_old_arr(n, r - 1)) *
                           g;
                } else {
                    dpdr = (p0_old_arr(n, r) - p0_old_arr(n, r - 1)) / dr(n);
                    rhog =
                        0.5 * (rho0_old_arr(n, r) + rho0_old_arr(n, r - 1)) * g;
                }

                max_hse_error =
                    amrex::max(max_hse_error, amrex::Math::abs(dpdr - rhog) /
                                                  amrex::Math::abs(dpdr));
            }
        }
    }

    Print() << "Maximum HSE error = " << max_hse_error << std::endl;
}
