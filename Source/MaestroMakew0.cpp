#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::Makew0(const RealVector& w0_old, 
                RealVector& w0_force, 
                const RealVector& Sbar_in, 
                const RealVector& rho0_old_in,
                const RealVector& rho0_new_in,
                const RealVector& p0_old_in,
                const RealVector& p0_new_in,
                const RealVector& gamma1bar_old_in,
                const RealVector& gamma1bar_new_in,
                const RealVector& p0_minus_peosbar, 
                RealVector& delta_chi_w0, 
                const Real dt_in, const Real dtold_in, 
                const bool is_predictor) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0()", Makew0);

    std::fill(w0_force.begin(), w0_force.end(), 0.);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    if (!spherical) {
        if (do_planar_invsq_grav || do_2d_planar_octant) {
            Makew0PlanarVarg(w0_old, w0_force, Sbar_in, 
                             rho0_old_in, rho0_new_in,
                             p0_old_in, p0_new_in, 
                             gamma1bar_old_in, gamma1bar_new_in, 
                             p0_minus_peosbar, 
                             delta_chi_w0, dt_in, dtold_in);
        } else {
            Makew0Planar(w0_old, w0_force, Sbar_in, 
                         rho0_old_in, rho0_new_in,
                         p0_old_in, p0_new_in, 
                         gamma1bar_old_in, gamma1bar_new_in, 
                         p0_minus_peosbar, 
                         delta_chi_w0, dt_in, dtold_in,
                         is_predictor);
        }
    } else {
        if (use_exact_base_state) {
            Makew0SphrIrreg(w0_old, w0_force, Sbar_in, 
                            rho0_old_in, rho0_new_in,
                            p0_old_in, p0_new_in, 
                            gamma1bar_old_in, gamma1bar_new_in, 
                            p0_minus_peosbar, 
                            delta_chi_w0, dt_in, dtold_in);
        } else {
            Makew0Sphr(w0_old, w0_force, Sbar_in, 
                       rho0_old_in, rho0_new_in,
                       p0_old_in, p0_new_in, 
                       gamma1bar_old_in, gamma1bar_new_in, 
                       p0_minus_peosbar, 
                       delta_chi_w0, dt_in, dtold_in);
        }
    }

    if (maestro_verbose >= 2) {
        for (auto n = 0; n <= finest_radial_level; ++n) {
            Real max_w0 = 0.0;
            for (auto r = r_start_coord[n]; r <= r_end_coord[n]+1; ++r) {
                max_w0 = max(max_w0, fabs(w0[n+max_lev*r]));
            }
            Print() << "... max CFL of w0: " << max_w0 * dt_in / dr[n] << std::endl;
        }
        Print() << std::endl;
    }
}

void 
Maestro::Makew0Planar(const RealVector& w0_old, 
                      RealVector& w0_force, 
                      const RealVector& Sbar_in, 
                      const RealVector& rho0_old_in,
                      const RealVector& rho0_new_in,
                      const RealVector& p0_old_in,
                      const RealVector& p0_new_in,
                      const RealVector& gamma1bar_old_in,
                      const RealVector& gamma1bar_new_in,
                      const RealVector& p0_minus_peosbar,  
                      RealVector& delta_chi_w0, 
                      const Real dt_in, const Real dtold_in, 
                      const bool is_predictor) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0Planar()", Makew0Planar);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Multilevel Outline
    //
    // Compute w0 at level 1 only
    // Initialize new w0 at bottom of coarse base array to 0.0.
    // do n=1,finest_radial_level
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

    std::fill(w0.begin(), w0.end(), 0.);

    const int max_lev = max_radial_level+1;

    // local variables 
    RealVector psi_planar_vec(nr_fine);

    Real * AMREX_RESTRICT psi_planar = psi_planar_vec.dataPtr();
    const Real * AMREX_RESTRICT etarho_cc_p = etarho_cc.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    const Real * AMREX_RESTRICT w0_old_p = w0_old.dataPtr();
    Real * AMREX_RESTRICT w0_force_p = w0_force.dataPtr();

    const Real dt_loc = dt;
    const Real grav_const_loc = grav_const;
    const Real dpdt_factor_loc = dpdt_factor;

    // Compute w0 on edges at level n
    for (auto n = 0; n <= max_radial_level; ++n) {

        std::fill(psi_planar_vec.begin(), psi_planar_vec.end(), 0.);
        int base_cutoff_density_coord_loc = 0;
        get_base_cutoff_density_coord(n, &base_cutoff_density_coord_loc);

        const Real dr_lev = dr[n];

        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {

            if (n == 0) {
                // Initialize new w0 at bottom of coarse base array to 0.0.
                w0[0] = 0.0;
            } else {
                // Obtain the starting value of w0 from the coarser grid
                w0[n+max_lev*r_start_coord[n+max_lev*j]] = w0[n-1+max_lev*r_start_coord[n+max_lev*j]/2];
            }

            // compute psi for level n
            int lo = r_start_coord[n+max_lev*j]; 
            int hi = r_end_coord[n+max_lev*j];
            AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                int r = k + lo;
                if (r < base_cutoff_density_coord_loc) {
                    psi_planar[r] = etarho_cc_p[n+max_lev*r] * fabs(grav_const_loc);
                }
            });

            for (auto r = r_start_coord[n+max_lev*j]+1; 
                r <= r_end_coord[n+max_lev*j]+1; ++r) {

                Real gamma1bar_p0_avg = (gamma1bar_old_in[n+max_lev*(r-1)]
                    + gamma1bar_new_in[n+max_lev*(r-1)]) *
                    (p0_old_in[n+max_lev*(r-1)] + 
                    p0_new_in[n+max_lev*(r-1)])/4.0;

                if (r < base_cutoff_density_coord_loc) {
                    if (is_predictor) {
                        delta_chi_w0[n+max_lev*(r-1)] = dpdt_factor_loc * 
                            p0_minus_peosbar[n+max_lev*(r-1)] / 
                            (gamma1bar_old_in[n+max_lev*(r-1)]*
                            p0_old_in[n+max_lev*(r-1)]*dt_loc);
                    } else {
                        delta_chi_w0[n+max_lev*(r-1)] += dpdt_factor_loc *
                            p0_minus_peosbar[n+max_lev*(r-1)] / 
                            (gamma1bar_new_in[n+max_lev*(r-1)]*
                            p0_new_in[n+max_lev*(r-1)]*dt_loc);
                    }
                } else {
                    delta_chi_w0[n+max_lev*(r-1)] = 0.0;
                }

                w0[n+max_lev*r] = w0[n+max_lev*(r-1)]
                    + Sbar_in[n+max_lev*(r-1)] * dr_lev
                    - psi_planar_vec[r-1] / gamma1bar_p0_avg * dr_lev
                    - delta_chi_w0[n+max_lev*(r-1)] * dr_lev;
            }

            if (n > 0) {
                // Compare the difference between w0 at top of level n to
                // the corresponding point on level n-1
                Real offset = w0[n+max_lev*(r_end_coord[n+max_lev*j]+1)]
                    - w0[n-1+max_lev*(r_end_coord[n+max_lev*j]+1)/2];

                for (auto i = n-1; i >= 0; --i) {

                    int refrat = pow(2, n-i);

                    // Restrict w0 from level n to level i
                    for (auto r = r_start_coord[n + max_lev*j]; r <= r_end_coord[n + max_lev*j]+1; ++r) {
                        if (r % refrat == 0) {
                            w0[n+max_lev*r/refrat] = w0[n+max_lev*r];
                        }
                    }

                    // Offset the w0 on level i above the top of level n
                    // for (auto r = (r_end_coord[n+max_lev*j]+1)/refrat+1; r <= nr[i]; ++r) {
                    lo = (r_end_coord[n+max_lev*j]+1)/refrat+1; 
                    hi = nr[i];
                    AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                        int r = k + lo;
                        w0_p[i+max_lev*r] += offset;
                    });
                }
            }
        }
    }

    // zero w0 where there is no corresponding full state array
    for (auto n = 1; n <= max_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
            if (j == numdisjointchunks[n]) {
                // for (auto r = r_end_coord[n+max_lev*j]+2; 
                //      r <= nr[n]; ++r) {
                const int lo = r_end_coord[n+max_lev*j]+2; 
                const int hi = nr[n];
                AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                    int r = k + lo;
                    w0_p[n+max_lev*r] = 0.0;
                });
            } else {
                // for (auto r = r_end_coord[n+max_lev*j]+2; 
                //      r <= r_start_coord[n+max_lev*(j+1)]-1; ++r) {
                const int lo = r_end_coord[n+max_lev*j]+2; 
                const int hi = r_start_coord[n+max_lev*(j+1)]-1;
                AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                    int r = k + lo;
                    w0_p[n+max_lev*r] = 0.0;
                });
            }
        }
    }

    RestrictBase(w0, false);
    FillGhostBase(w0, false);

    for (auto n = 0; n <= max_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {

            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            const Real dt_avg = 0.5 * (dt_in + dtold_in);
            const Real dr_lev = dr[n];

            // for (auto r = r_start_coord[n + max_lev*j]; 
            //      r <= r_end_coord[n + max_lev*j]; ++r) {
            const int lo = r_start_coord[n+max_lev*j]; 
            const int hi = r_end_coord[n+max_lev*j];
            AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                int r = k + lo;

                Real w0_old_cen = 0.5 * (w0_old_p[n+max_lev*r] + 
                    w0_old_p[n+max_lev*(r+1)]);
                Real w0_new_cen = 0.5 * (w0_p[n+max_lev*r] + 
                    w0_p[n+max_lev*(r+1)]);
                Real w0_avg = 0.5 * (dt_in * w0_old_cen + dtold_in *  w0_new_cen) / dt_avg;
                Real div_avg = 0.5 * (dt_in *(w0_old_p[n+max_lev*(r+1)]
                    - w0_old_p[n+max_lev*r]) + dtold_in * (w0_p[n+max_lev*(r+1)] 
                    - w0_p[n+max_lev*r])) / dt_avg;
                w0_force_p[n+max_lev*r] = (w0_new_cen - w0_old_cen)/dt_avg 
                    + w0_avg*div_avg/dr_lev;
            });
        }
    }

    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}

void 
Maestro::Makew0PlanarVarg(const RealVector& w0_old, 
                          RealVector& w0_force, 
                          const RealVector& Sbar_in, 
                          const RealVector& rho0_old_in,
                          const RealVector& rho0_new_in,
                          const RealVector& p0_old_in,
                          const RealVector& p0_new_in,
                          const RealVector& gamma1bar_old_in,
                          const RealVector& gamma1bar_new_in,
                          const RealVector& p0_minus_peosbar,  
                          RealVector& delta_chi_w0, 
                          const Real dt_in, const Real dtold_in) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0PlanarVarg()",Makew0PlanarVarg);

    get_finest_radial_level(&finest_radial_level);
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    int fine_base_density_cutoff_coord = 0;
    get_base_cutoff_density_coord(finest_radial_level, &fine_base_density_cutoff_coord);

    const int max_lev = max_radial_level+1;
    const int nr_finest = nr[finest_radial_level];
    const Real dr_finest = dr[finest_radial_level];
    const Real dpdt_factor_loc = dpdt_factor;

    Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    Real * AMREX_RESTRICT w0_force_p = w0_force.dataPtr();
    const Real * AMREX_RESTRICT w0_old_p = w0_old.dataPtr();
    const Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();

    // The planar 1/r**2 gravity constraint equation is solved
    // by calling the tridiagonal solver, just like spherical.
    // This is accomplished by putting all the requisite data
    // on the finest basestate grid, solving for w0, and then
    // restricting w0 back down to the coarse grid.

    // 1) allocate the finely-gridded temporary basestate arrays
    RealVector w0_fine_vec(nr_finest+1);
    RealVector w0bar_fine_vec(nr_finest+1);
    RealVector deltaw0_fine_vec(nr_finest+1);
    RealVector p0_old_fine_vec(nr_finest);
    RealVector p0_new_fine_vec(nr_finest);
    RealVector p0_nph_fine_vec(nr_finest);
    RealVector rho0_old_fine_vec(nr_finest);
    RealVector rho0_new_fine_vec(nr_finest);
    RealVector rho0_nph_fine_vec(nr_finest);
    RealVector gamma1bar_old_fine_vec(nr_finest);
    RealVector gamma1bar_new_fine_vec(nr_finest);
    RealVector gamma1bar_nph_fine_vec(nr_finest);
    RealVector p0_minus_peosbar_fine_vec(nr_finest);
    RealVector etarho_cc_fine_vec(nr_finest);
    RealVector Sbar_in_fine_vec(nr_finest);
    RealVector grav_edge_fine_vec(nr_finest+1);

    // 2) copy the data into the temp, uniformly-gridded basestate arrays.
    ProlongBasetoUniform(p0_old,p0_old_fine_vec);
    ProlongBasetoUniform(p0_new,p0_new_fine_vec);
    ProlongBasetoUniform(rho0_old,rho0_old_fine_vec);
    ProlongBasetoUniform(rho0_new,rho0_new_fine_vec);
    ProlongBasetoUniform(gamma1bar_old,gamma1bar_old_fine_vec);
    ProlongBasetoUniform(gamma1bar_new,gamma1bar_new_fine_vec);
    ProlongBasetoUniform(p0_minus_peosbar,p0_minus_peosbar_fine_vec);
    ProlongBasetoUniform(etarho_cc,etarho_cc_fine_vec);
    ProlongBasetoUniform(Sbar_in,Sbar_in_fine_vec);

    Real * AMREX_RESTRICT w0_fine = w0_fine_vec.dataPtr();
    Real * AMREX_RESTRICT w0bar_fine = w0bar_fine_vec.dataPtr();
    Real * AMREX_RESTRICT deltaw0_fine = deltaw0_fine_vec.dataPtr();
    Real * AMREX_RESTRICT p0_old_fine = p0_old_fine_vec.dataPtr();
    Real * AMREX_RESTRICT p0_new_fine = p0_new_fine_vec.dataPtr();
    Real * AMREX_RESTRICT p0_nph_fine = p0_nph_fine_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_old_fine = rho0_old_fine_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_new_fine = rho0_new_fine_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_nph_fine = rho0_nph_fine_vec.dataPtr();
    Real * AMREX_RESTRICT gamma1bar_old_fine = gamma1bar_old_fine_vec.dataPtr();
    Real * AMREX_RESTRICT gamma1bar_new_fine = gamma1bar_new_fine_vec.dataPtr();
    Real * AMREX_RESTRICT gamma1bar_nph_fine = gamma1bar_nph_fine_vec.dataPtr();
    Real * AMREX_RESTRICT p0_minus_peosbar_fine = p0_minus_peosbar_fine_vec.dataPtr();
    Real * AMREX_RESTRICT etarho_cc_fine = etarho_cc_fine_vec.dataPtr();
    Real * AMREX_RESTRICT Sbar_in_fine = Sbar_in_fine_vec.dataPtr();
    Real * AMREX_RESTRICT grav_edge_fine = grav_edge_fine_vec.dataPtr();

    // create time-centered base-state quantities
    // for (auto r = 0; r < nr_finest; ++r) {
    AMREX_PARALLEL_FOR_1D(nr_finest, r, {
        p0_nph_fine[r] = 0.5*(p0_old_fine[r] + p0_new_fine[r]);
        rho0_nph_fine[r] = 0.5*(rho0_old_fine[r] + rho0_new_fine[r]);
        gamma1bar_nph_fine[r] = 0.5*(gamma1bar_old_fine[r] + gamma1bar_new_fine[r]);
    });


    // 3) solve to w0bar -- here we just take into account the Sbar and
    //    volume discrepancy terms
    // lower boundary condition
    w0bar_fine_vec[0] = 0.0;

    // for (auto r = 1; r <= nr_finest; ++r) {
    int lo = 1; 
    int hi = nr_finest;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        Real gamma1bar_p0_avg = gamma1bar_nph_fine[r-1] * p0_nph_fine[r-1];

        Real volume_discrepancy = (r-1 < fine_base_density_cutoff_coord) ? 
            dpdt_factor_loc * p0_minus_peosbar_fine[r-1]/dt_in : 0.0;

        w0bar_fine[r] =  w0bar_fine[r-1] + 
            Sbar_in_fine[r-1] * dr_finest 
            - (volume_discrepancy / gamma1bar_p0_avg ) * dr_finest;
    });

    // 4) get the edge-centered gravity on the uniformly-gridded
    // basestate arrays
    Abort("make_w0.f90: need to write make_grav_edge_uniform");
    //    call make_grav_edge_uniform(grav_edge_fine, rho0_nph_fine)


    // 5) solve for delta w0
    std::fill(deltaw0_fine_vec.begin(), deltaw0_fine_vec.end(), 0.);

    // this takes the form of a tri-diagonal matrix:
    // A_j (dw_0)_{j-3/2} +
    // B_j (dw_0)_{j-1/2} +
    // C_j (dw_0)_{j+1/2} = F_j

    RealVector A_vec(nr_finest+1);
    RealVector B_vec(nr_finest+1);
    RealVector C_vec(nr_finest+1);
    RealVector u_vec(nr_finest+1);
    RealVector F_vec(nr_finest+1);

    std::fill(A_vec.begin(), A_vec.end(), 0.);
    std::fill(B_vec.begin(), B_vec.end(), 0.);
    std::fill(C_vec.begin(), C_vec.end(), 0.);
    std::fill(F_vec.begin(), F_vec.end(), 0.);
    std::fill(u_vec.begin(), u_vec.end(), 0.);

    Real * AMREX_RESTRICT A = A_vec.dataPtr();
    Real * AMREX_RESTRICT B = B_vec.dataPtr();
    Real * AMREX_RESTRICT C = C_vec.dataPtr();
    Real * AMREX_RESTRICT F = F_vec.dataPtr();
    Real * AMREX_RESTRICT u = u_vec.dataPtr();

    // for (auto r = 1; r <= fine_base_density_cutoff_coord; ++r) {
    lo = 1; 
    hi = fine_base_density_cutoff_coord;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        A[r] = gamma1bar_nph_fine[r-1] * p0_nph_fine[r-1];
        A[r] /= dr_finest*dr_finest;

        Real dpdr = (p0_nph_fine[r]-p0_nph_fine[r-1])/dr_finest;

        B[r] = -(gamma1bar_nph_fine[r-1] * p0_nph_fine[r-1] + 
            gamma1bar_nph_fine[r] * p0_nph_fine[r]) 
            / (dr_finest*dr_finest);
        B[r] -= 2.0 * dpdr / (r_edge_loc_p[finest_radial_level+max_lev*r]);

        C[r] = gamma1bar_nph_fine[r] * p0_nph_fine[r];
        C[r] /= dr_finest*dr_finest;

        F[r] = 2.0 * dpdr * w0bar_fine[r] / 
            r_edge_loc_p[finest_radial_level+max_lev*r] -
            grav_edge_fine[r] * (etarho_cc_fine[r] - etarho_cc_fine[r-1]) / 
            dr_finest;
    });

    // Lower boundary
    A_vec[0] = 0.0;
    B_vec[0] = 1.0;
    C_vec[0] = 0.0;
    F_vec[0] = 0.0;

    // Upper boundary
    A_vec[fine_base_density_cutoff_coord+1] = -1.0;
    B_vec[fine_base_density_cutoff_coord+1] = 1.0;
    C_vec[fine_base_density_cutoff_coord+1] = 0.0;
    F_vec[fine_base_density_cutoff_coord+1] = 0.0;

    // need to synchronize gpu values with updated host values
    Gpu::synchronize();
    
    // Call the tridiagonal solver
    Tridiag(A_vec, B_vec, C_vec, F_vec, u_vec, fine_base_density_cutoff_coord+2);

    // for (auto r = 1; r <= fine_base_density_cutoff_coord+1; ++r) {
    lo = 1; 
    hi = fine_base_density_cutoff_coord+1;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        deltaw0_fine[r] = u[r];
    });

    // for (auto r = fine_base_density_cutoff_coord+2; r <= nr_finest; ++r) {
    lo = fine_base_density_cutoff_coord+2; 
    hi = nr_finest;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        deltaw0_fine[r] = deltaw0_fine[fine_base_density_cutoff_coord+1];
    });

    // 6) compute w0 = w0bar + deltaw0
    // for (auto r = 0; r < w0_fine.size(); ++r) {
    AMREX_PARALLEL_FOR_1D(w0_fine_vec.size(), r, {
        w0_fine[r] = w0bar_fine[r] + deltaw0_fine[r];
        w0_p[finest_radial_level+max_lev*r] = w0_fine[r];
    });

    // 7) fill the multilevel w0 array from the uniformly-gridded w0 we
    // just solved for.  Here, we make the coarse edge underneath equal
    // to the fine edge value.
    for (auto n = finest_radial_level; n >= 1; --n) {
        for (auto r = 0; r <= nr[n]; n+=2) {
            w0[n-1+max_lev*r/2] = w0[n+max_lev*r];
        }
    }

    // 8) zero w0 where there is no corresponding full state array
    for (auto n = 1; n <= finest_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
            if (j == numdisjointchunks[n]) {
                // for (auto r = r_end_coord[n+max_lev*j]+2; r <= nr[n]; ++r) {
                lo = r_end_coord[n+max_lev*j]+2; 
                hi = nr[n];
                AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                    int r = k + lo;
                    w0_p[n+max_lev*r] = 0.0;
                });
            } else {
                // for (auto r = r_end_coord[n+max_lev*j]+2; r < r_start_coord[n+max_lev*(j+1)]; ++r) {
                lo = r_end_coord[n+max_lev*j]+2; 
                hi = r_start_coord[n+max_lev*(j+1)];
                AMREX_PARALLEL_FOR_1D(hi-lo, k, {
                    int r = k + lo;
                    w0_p[n+max_lev*r] = 0.0;
                });
            }
        }
    }

    RestrictBase(w0, false);
    FillGhostBase(w0, false);

    // compute the forcing terms
    for (auto n = 0; n <= finest_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {

            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            const Real dt_avg = 0.5 * (dt_in + dtold_in);
            const Real dr_lev = dr[n];

            // for (auto r = r_start_coord[n+max_lev*j]; r <=r_end_coord[n+max_lev*j]; ++r) {
            lo = r_start_coord[n+max_lev*j]; 
            hi = r_end_coord[n+max_lev*j];
            AMREX_PARALLEL_FOR_1D(hi-lo+1, k, {
                int r = k + lo;
                Real w0_old_cen = 0.5 * (w0_old_p[n+max_lev*r] + w0_old_p[n+max_lev*(r+1)]);
                Real w0_new_cen = 0.5 * (w0_p[n+max_lev*r] + w0_p[n+max_lev*(r+1)]);
                Real w0_avg = 0.5 * (dt_in * w0_old_cen + dtold_in *  w0_new_cen) / dt_avg;
                Real div_avg = 0.5 * (dt_in * (w0_old_p[n+max_lev*(r+1)]-w0_old_p[n+max_lev*r]) + 
                    dtold_in * (w0_p[n+max_lev*(r+1)]-w0_p[n+max_lev*r])) / dt_avg;
                w0_force_p[n+max_lev*r] = (w0_new_cen-w0_old_cen)/dt_avg + w0_avg*div_avg/dr_lev;
            });
        }
    }

    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}

void 
Maestro::Makew0Sphr(const RealVector& w0_old, 
                    RealVector& w0_force, 
                    const RealVector& Sbar_in, 
                    const RealVector& rho0_old_in,
                    const RealVector& rho0_new_in,
                    const RealVector& p0_old_in,
                    const RealVector& p0_new_in,
                    const RealVector& gamma1bar_old_in,
                    const RealVector& gamma1bar_new_in,
                    const RealVector& p0_minus_peosbar,  
                    RealVector& delta_chi_w0, 
                    const Real dt_in, const Real dtold_in) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0Sphr()",Makew0Sphr);

    // local variables 
    const int max_lev = max_radial_level + 1;
    RealVector gamma1bar_nph_vec(nr_fine);
    RealVector p0_nph_vec(nr_fine);
    RealVector A_vec(nr_fine+1);
    RealVector B_vec(nr_fine+1);
    RealVector C_vec(nr_fine+1);
    RealVector u_vec(nr_fine+1);
    RealVector F_vec(nr_fine+1);
    RealVector w0_from_Sbar_vec(nr_fine+1);
    RealVector rho0_nph_vec(max_lev*nr_fine);
    RealVector grav_edge_vec(max_lev*(nr_fine+1));

    Real * AMREX_RESTRICT gamma1bar_nph = gamma1bar_nph_vec.dataPtr();
    Real * AMREX_RESTRICT p0_nph = p0_nph_vec.dataPtr();
    Real * AMREX_RESTRICT A = A_vec.dataPtr();
    Real * AMREX_RESTRICT B = B_vec.dataPtr();
    Real * AMREX_RESTRICT C = C_vec.dataPtr();
    Real * AMREX_RESTRICT u = u_vec.dataPtr();
    Real * AMREX_RESTRICT F = F_vec.dataPtr();
    Real * AMREX_RESTRICT w0_from_Sbar = w0_from_Sbar_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_nph = rho0_nph_vec.dataPtr();
    Real * AMREX_RESTRICT grav_edge = grav_edge_vec.dataPtr();

    const Real * AMREX_RESTRICT p0_old_p = p0_old_in.dataPtr();
    const Real * AMREX_RESTRICT p0_new_p = p0_new_in.dataPtr();
    const Real * AMREX_RESTRICT rho0_old_p = rho0_old_in.dataPtr();
    const Real * AMREX_RESTRICT rho0_new_p = rho0_new_in.dataPtr();
    const Real * AMREX_RESTRICT gamma1bar_old_p = gamma1bar_old_in.dataPtr();
    const Real * AMREX_RESTRICT gamma1bar_new_p = gamma1bar_new_in.dataPtr();
    const Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    const Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    const Real * AMREX_RESTRICT etarho_cc_p = etarho_cc.dataPtr();
    const Real * AMREX_RESTRICT etarho_ec_p = etarho_ec.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    const Real * AMREX_RESTRICT w0_old_p = w0_old.dataPtr();
    Real * AMREX_RESTRICT w0_force_p = w0_force.dataPtr();

    Real base_cutoff_dens = 0.0;
    get_base_cutoff_density(&base_cutoff_dens);
    int base_cutoff_density_coord_loc = 0;
    get_base_cutoff_density_coord(0, &base_cutoff_density_coord_loc);

    const Real dr0 = dr[0];
    const Real dpdt_factor_loc = dpdt_factor;

    // create time-centered base-state quantities
    // for (auto r = 0; r < nr_fine; ++r) {
    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        p0_nph[r] = 0.5*(p0_old_p[max_lev*r] + p0_new_p[max_lev*r]);
        rho0_nph[max_lev*r] = 0.5*(rho0_old_p[max_lev*r] + rho0_new_p[max_lev*r]);
        gamma1bar_nph[r] = 0.5*(gamma1bar_old_p[max_lev*r] + gamma1bar_new_p[max_lev*r]);
    });

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.
    w0_from_Sbar_vec[0] = 0.0;

    for (auto r = 1; r <= nr_fine; ++r) {
        // Print() << "r = " << r << " denominator = " << (gamma1bar_nph_vec[r-1]*p0_nph_vec[r-1]) << std::endl;
        Real volume_discrepancy = rho0_old_in[max_lev*(r-1)] > base_cutoff_dens ? 
            dpdt_factor_loc * p0_minus_peosbar[max_lev*(r-1)]/dt_in : 0.0;

        w0_from_Sbar_vec[r] = w0_from_Sbar_vec[r-1] + 
            dr0 * Sbar_in[max_lev*(r-1)] * r_cc_loc[max_lev*(r-1)]*r_cc_loc[max_lev*(r-1)];
        if (volume_discrepancy != 0.0) {
            w0_from_Sbar_vec[r] -= dr0 * volume_discrepancy * r_cc_loc[max_lev*(r-1)]*r_cc_loc[max_lev*(r-1)] 
            / (gamma1bar_nph_vec[r-1]*p0_nph_vec[r-1]);
        }
    }

    // for (auto r = 1; r <= nr_fine; ++r) {
    int lo = 1; 
    int hi = nr_fine;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        w0_from_Sbar[r] /= (r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]);
    });

    // make the edge-centered gravity
    MakeGravEdge(grav_edge_vec, rho0_nph_vec);

    // NOTE:  now we solve for the remainder, (r^2 * delta w0)
    // this takes the form of a tri-diagonal matrix:
    // A_j (r^2 dw_0)_{j-3/2} +
    // B_j (r^2 dw_0)_{j-1/2} +
    // C_j (r^2 dw_0)_{j+1/2} = F_j
    std::fill(A_vec.begin(), A_vec.end(), 0.);
    std::fill(B_vec.begin(), B_vec.end(), 0.);
    std::fill(C_vec.begin(), C_vec.end(), 0.);
    std::fill(F_vec.begin(), F_vec.end(), 0.);
    std::fill(u_vec.begin(), u_vec.end(), 0.);

    // Note that we are solving for (r^2 delta w0), not just w0.

    int max_cutoff = min(base_cutoff_density_coord_loc, nr_fine-1);
    
    // for (auto r = 1; r <= max_cutoff; ++r) {
    lo = 1; 
    hi = max_cutoff;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        A[r] = gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc_p[max_lev*(r-1)]*r_cc_loc_p[max_lev*(r-1)]);
        A[r] /= dr0*dr0;

        B[r] = -( gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc_p[max_lev*(r-1)]*r_cc_loc_p[max_lev*(r-1)])
                + gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r]) ) 
                / (dr0*dr0);

        Real dpdr = (p0_nph[r] - p0_nph[r-1]) / dr0;

        B[r] -= 4.0 * dpdr / (r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]);

        C[r] = gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r]);
        C[r] /= dr0*dr0;

        F[r] = 4.0 * dpdr * w0_from_Sbar[r] / r_edge_loc_p[max_lev*r] - 
                grav_edge[max_lev*r] * (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r] * etarho_cc_p[max_lev*r] - 
                r_cc_loc_p[max_lev*(r-1)]*r_cc_loc_p[max_lev*(r-1)] * etarho_cc_p[max_lev*(r-1)]) / 
                (dr0 * r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]) - 
                4.0 * M_PI * Gconst * 0.5 * 
                (rho0_nph[max_lev*r] + rho0_nph[max_lev*(r-1)]) * etarho_ec_p[max_lev*r];
    });

    // Lower boundary
    A_vec[0] = 0.0;
    B_vec[0] = 1.0;
    C_vec[0] = 0.0;
    F_vec[0] = 0.0;

    // Upper boundary
    A_vec[max_cutoff+1] = -1.0;
    B_vec[max_cutoff+1] = 1.0;
    C_vec[max_cutoff+1] = 0.0;
    F_vec[max_cutoff+1] = 0.0;

    // need to synchronize gpu values with updated host values
    Gpu::synchronize();
    
    // Call the tridiagonal solver
    Tridiag(A_vec, B_vec, C_vec, F_vec, u_vec, max_cutoff+2);

    w0[0] = w0_from_Sbar_vec[0];

    // for (auto r = 1; r <= max_cutoff+1; ++r) {
    lo = 1; 
    hi = max_cutoff+1;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        w0_p[max_lev*r] = u[r] / (r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]) + w0_from_Sbar[r];
    });

    // for (auto r = max_cutoff+2; r <= nr_fine; ++r) {
    lo = max_cutoff+2; 
    hi = nr_fine;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        w0_p[max_lev*r] = w0_p[max_cutoff+1] * r_edge_loc_p[max_lev*(max_cutoff+1)]*r_edge_loc_p[max_lev*(max_cutoff+1)]/(r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]);
    });

    // Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    const Real dt_avg = 0.5 * (dt_in + dtold_in);

    // for (auto r = 0; r < nr_fine; ++r) {
    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        Real w0_old_cen = 0.5 * (w0_old_p[max_lev*r] + w0_old_p[max_lev*r+1]);
        Real w0_new_cen = 0.5 * (w0_p[max_lev*r] + w0_p[max_lev*(r+1)]);
        Real w0_avg = 0.5 * (dt_in *  w0_old_cen + dtold_in *  w0_new_cen) / dt_avg;
        Real div_avg = 0.5 * (dt_in * (w0_old_p[max_lev*(r+1)]-w0_old_p[max_lev*r]) + dtold_in * (w0_p[max_lev*(r+1)]-w0_p[max_lev*r])) / dt_avg;
        w0_force_p[max_lev*r] = (w0_new_cen-w0_old_cen) / dt_avg + w0_avg * div_avg / dr0;
    });
}

void 
Maestro::Makew0SphrIrreg(const RealVector& w0_old, 
                        RealVector& w0_force, 
                        const RealVector& Sbar_in, 
                        const RealVector& rho0_old_in,
                        const RealVector& rho0_new_in,
                        const RealVector& p0_old_in,
                        const RealVector& p0_new_in,
                        const RealVector& gamma1bar_old_in,
                        const RealVector& gamma1bar_new_in,
                        const RealVector& p0_minus_peosbar,  
                        RealVector& delta_chi_w0, 
                        const Real dt_in, const Real dtold_in) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0SphrIrreg()",Makew0SphrIrreg);

    // local variables 
    const int max_lev = max_radial_level+1;
    RealVector gamma1bar_nph_vec(nr_fine);
    RealVector p0_nph_vec(nr_fine);
    RealVector A_vec(nr_fine+1);
    RealVector B_vec(nr_fine+1);
    RealVector C_vec(nr_fine+1);
    RealVector u_vec(nr_fine+1);
    RealVector F_vec(nr_fine+1);
    RealVector w0_from_Sbar_vec(nr_fine+1);
    RealVector rho0_nph_vec(max_lev*nr_fine);
    RealVector grav_edge_vec(max_lev*(nr_fine+1));

    Real * AMREX_RESTRICT gamma1bar_nph = gamma1bar_nph_vec.dataPtr();
    Real * AMREX_RESTRICT p0_nph = p0_nph_vec.dataPtr();
    Real * AMREX_RESTRICT A = A_vec.dataPtr();
    Real * AMREX_RESTRICT B = B_vec.dataPtr();
    Real * AMREX_RESTRICT C = C_vec.dataPtr();
    Real * AMREX_RESTRICT u = u_vec.dataPtr();
    Real * AMREX_RESTRICT F = F_vec.dataPtr();
    Real * AMREX_RESTRICT w0_from_Sbar = w0_from_Sbar_vec.dataPtr();
    Real * AMREX_RESTRICT rho0_nph = rho0_nph_vec.dataPtr();
    Real * AMREX_RESTRICT grav_edge = grav_edge_vec.dataPtr();

    const Real * AMREX_RESTRICT p0_old_p = p0_old_in.dataPtr();
    const Real * AMREX_RESTRICT p0_new_p = p0_new_in.dataPtr();
    const Real * AMREX_RESTRICT rho0_old_p = rho0_old_in.dataPtr();
    const Real * AMREX_RESTRICT rho0_new_p = rho0_new_in.dataPtr();
    const Real * AMREX_RESTRICT gamma1bar_old_p = gamma1bar_old_in.dataPtr();
    const Real * AMREX_RESTRICT gamma1bar_new_p = gamma1bar_new_in.dataPtr();
    const Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    const Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    const Real * AMREX_RESTRICT etarho_cc_p = etarho_cc.dataPtr();
    const Real * AMREX_RESTRICT etarho_ec_p = etarho_ec.dataPtr();
    Real * AMREX_RESTRICT w0_p = w0.dataPtr();
    const Real * AMREX_RESTRICT w0_old_p = w0_old.dataPtr();
    Real * AMREX_RESTRICT w0_force_p = w0_force.dataPtr();

    Real base_cutoff_dens = 0.0;
    get_base_cutoff_density(&base_cutoff_dens);
    int base_cutoff_density_coord_loc = 0;
    get_base_cutoff_density_coord(0, &base_cutoff_density_coord_loc);
    const Real dpdt_factor_loc = dpdt_factor;

    // create time-centered base-state quantities
    // for (auto r = 0; r < nr_fine; ++r) {
    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        p0_nph[r] = 0.5*(p0_old_p[max_lev*r] + p0_new_p[max_lev*r]);
        rho0_nph[r] = 0.5*(rho0_old_p[max_lev*r] + rho0_new_p[max_lev*r]);
        gamma1bar_nph[r] = 0.5*(gamma1bar_old_p[max_lev*r] + gamma1bar_new_p[max_lev*r]);
    });

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.
    w0_from_Sbar_vec[0] = 0.0;

    for (auto r = 1; r <= nr_fine; ++r) {
        Real volume_discrepancy = rho0_old_in[max_lev*(r-1)] > base_cutoff_dens ? 
            dpdt_factor_loc * p0_minus_peosbar[max_lev*(r-1)]/dt_in : 0.0;

        Real dr1 = r_edge_loc[max_lev*r] - r_edge_loc[max_lev*(r-1)];
        w0_from_Sbar_vec[r] = w0_from_Sbar_vec[r-1] + 
            dr1 * Sbar_in[max_lev*(r-1)] * r_cc_loc[max_lev*(r-1)]*r_cc_loc[max_lev*(r-1)] - 
            dr1* volume_discrepancy * r_cc_loc[max_lev*(r-1)]*r_cc_loc[max_lev*(r-1)] 
            / (gamma1bar_nph[r-1]*p0_nph[r-1]);
    }

    for (auto r = 1; r <= nr_fine; ++r) {
        w0_from_Sbar_vec[r] /= (r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r]);
    }

    // make the edge-centered gravity
    MakeGravEdge(grav_edge_vec, rho0_nph_vec);

    // NOTE:  now we solve for the remainder, (r^2 * delta w0)
    // this takes the form of a tri-diagonal matrix:
    // A_j (r^2 dw_0)_{j-3/2} +
    // B_j (r^2 dw_0)_{j-1/2} +
    // C_j (r^2 dw_0)_{j+1/2} = F_j
    std::fill(A_vec.begin(), A_vec.end(), 0.);
    std::fill(B_vec.begin(), B_vec.end(), 0.);
    std::fill(C_vec.begin(), C_vec.end(), 0.);
    std::fill(F_vec.begin(), F_vec.end(), 0.);
    std::fill(u_vec.begin(), u_vec.end(), 0.);

    // Note that we are solving for (r^2 delta w0), not just w0.
    int max_cutoff = base_cutoff_density_coord_loc;
    
    // for (auto r = 1; r <= max_cutoff; ++r) {
    int lo = 1; 
    int hi = max_cutoff;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        Real dr1 = r_edge_loc_p[max_lev*r] - r_edge_loc_p[max_lev*(r-1)];
        Real dr2 = r_edge_loc_p[max_lev*(r+1)] - r_edge_loc_p[max_lev*r];
        Real dr3 = r_cc_loc_p[max_lev*r] - r_cc_loc_p[max_lev*(r-1)];

        A[r] = gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc_p[max_lev*(r-1)]*r_cc_loc_p[max_lev*(r-1)]);
        A[r] /= dr1*dr3;

        B[r] = -( gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc_p[max_lev*(r-1)]*r_cc_loc_p[max_lev*(r-1)]*dr1) 
                + gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r]*dr2) ) 
                / dr3;

        Real dpdr = (p0_nph[r] - p0_nph[r-1]) / dr3;

        B[r] -= 4.0 * dpdr / (r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]);

        C[r] = gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r]);
        C[r] /= dr2*dr3;

        F[r] = 4.0 * dpdr * w0_from_Sbar[r] / r_edge_loc_p[max_lev*r] - 
                grav_edge[max_lev*r] * (r_cc_loc_p[max_lev*r]*r_cc_loc_p[max_lev*r] * etarho_cc_p[max_lev*r] - 
                r_cc_loc_p[max_lev*(r-1)]*r_cc_loc_p[max_lev*(r-1)] * etarho_cc_p[max_lev*(r-1)]) / 
                (dr3 * r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]) - 
                4.0 * M_PI * Gconst * 0.5 * 
                (rho0_nph[max_lev*r] + rho0_nph[max_lev*(r-1)]) * etarho_ec_p[max_lev*r];
    });

    // Lower boundary
    A_vec[0] = 0.0;
    B_vec[0] = 1.0;
    C_vec[0] = 0.0;
    F_vec[0] = 0.0;

    // Upper boundary
    A_vec[max_cutoff+1] = -1.0;
    B_vec[max_cutoff+1] = 1.0;
    C_vec[max_cutoff+1] = 0.0;
    F_vec[max_cutoff+1] = 0.0;
    
    // need to synchronize gpu values with updated host values
    Gpu::synchronize();
    
    // Call the tridiagonal solver
    Tridiag(A_vec, B_vec, C_vec, F_vec, u_vec, max_cutoff+2);

    w0_p[0] = w0_from_Sbar_vec[0];

    // for (auto r = 1; r <= max_cutoff+1; ++r) {
    lo = 1; 
    hi = max_cutoff+1;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        w0_p[max_lev*r] = u[r] / (r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]) + w0_from_Sbar[r];
    });

    // for (auto r = max_cutoff+2; r <= nr_fine; ++r) {
    lo = max_cutoff+2; 
    hi = nr_fine;
    AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
        int r = j + lo;
        w0_p[max_lev*r] = w0_p[max_lev*(max_cutoff+1)] * r_edge_loc_p[max_lev*(max_cutoff+1)]*r_edge_loc_p[max_lev*(max_cutoff+1)]/(r_edge_loc_p[max_lev*r]*r_edge_loc_p[max_lev*r]);
    });

    // Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    const Real dt_avg = 0.5 * (dt_in + dtold_in);

    // for (auto r = 0; r < nr_fine; ++r) {
    AMREX_PARALLEL_FOR_1D(nr_fine, r, {
        Real dr1 = r_edge_loc_p[max_lev*r] - r_edge_loc_p[max_lev*(r-1)];
        Real w0_old_cen = 0.5 * (w0_old_p[max_lev*r] + w0_old_p[max_lev*(r+1)]);
        Real w0_new_cen = 0.5 * (w0_p[max_lev*r] + w0_p[max_lev*(r+1)]);
        Real w0_avg = 0.5 * (dt_in *  w0_old_cen + dtold_in *  w0_new_cen) / dt_avg;
        Real div_avg = 0.5 * (dt_in * (w0_old_p[max_lev*(r+1)]-w0_old_p[max_lev*r]) + dtold_in * (w0_p[max_lev*(r+1)]-w0_p[max_lev*r])) / dt_avg;
        w0_force_p[max_lev*r] = (w0_new_cen-w0_old_cen) / dt_avg + w0_avg * div_avg / dr1;
    });
}

void
Maestro::Tridiag(const RealVector& a, const RealVector& b, 
                 const RealVector& c, const RealVector& r, 
                 RealVector& u, const int n)
{
    RealVector gam(n);

    if (b[0] == 0) Abort("tridiag: CANT HAVE B[0] = 0.0");

    Real bet = b[0];
    u[0] = r[0] / bet;

    for (auto j = 1; j < n; j++) {
        gam[j] = c[j-1] / bet;
        bet = b[j] - a[j] * gam[j];
        if (bet == 0) Abort("tridiag: TRIDIAG FAILED");
        u[j] = (r[j] - a[j] * u[j-1]) / bet;
    }

    for (auto j = n-2; j >= 0; --j) {
        u[j] -= gam[j+1] * u[j+1];
    }
}

void
Maestro::ProlongBasetoUniform(const RealVector& base_ml, 
                              RealVector& base_fine)

{
    // the mask array will keep track of whether we've filled in data
    // in a corresponding radial bin.  .false. indicates that we've
    // already output there.
    IntVector imask_fine(nr_fine);
    std::fill(imask_fine.begin(), imask_fine.end(), 1);

    // r1 is the factor between the current level grid spacing and the
    // FINEST level
    int r1 = 1;

    get_finest_radial_level(&finest_radial_level);
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    const int max_lev = max_radial_level+1;

    for (auto n = finest_radial_level; n >= 0; --n) {
        for (auto j = 1; j < numdisjointchunks[n]; ++j) {
            for (auto r = r_start_coord[n+max_lev*j]; r <= r_end_coord[n+max_lev*j]; ++r) {
                // sum up mask to see if there are any elements set to true 
                if (std::accumulate(imask_fine.begin()+r*r1-1, imask_fine.begin()+(r+1)*r1-1, 0) > 0) {
                    for (auto i = r*r1-1; i < (r+1)*r1-1; ++r) {
                        base_fine[i] = base_ml[n+max_lev*r];
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
