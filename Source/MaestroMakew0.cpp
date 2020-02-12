#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::Makew0(RealVector& w0_in, 
                const RealVector& w0_old, 
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
    BL_PROFILE_VAR("Maestro::Makew0()",Makew0);

    std::fill(w0_force.begin(), w0_force.end(), 0.);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());

    if (!spherical) {
        if (do_planar_invsq_grav || do_2d_planar_octant) {
            Makew0PlanarVarg(w0_in, w0_old, w0_force, Sbar_in, 
                         rho0_old_in, rho0_new_in,
                         p0_old_in, p0_new_in, 
                         gamma1bar_old_in, gamma1bar_new_in, 
                         p0_minus_peosbar, 
                         delta_chi_w0, dt_in, dtold_in);
        } else {
            Makew0Planar(w0_in, w0_old, w0_force, Sbar_in, 
                         rho0_old_in, rho0_new_in,
                         p0_old_in, p0_new_in, 
                         gamma1bar_old_in, gamma1bar_new_in, 
                         p0_minus_peosbar, 
                         delta_chi_w0, dt_in, dtold_in,
                         is_predictor);
        }
    } else {
        if (use_exact_base_state) {
            Makew0SphrIrreg(w0_in, w0_old, w0_force, Sbar_in, 
                            rho0_old_in, rho0_new_in,
                            p0_old_in, p0_new_in, 
                            gamma1bar_old_in, gamma1bar_new_in, 
                            p0_minus_peosbar, 
                            delta_chi_w0, dt_in, dtold_in);
        } else {
            Makew0Sphr(w0_in, w0_old, w0_force, Sbar_in, 
                      rho0_old_in, rho0_new_in,
                      p0_old_in, p0_new_in, 
                      gamma1bar_old_in, gamma1bar_new_in, 
                      p0_minus_peosbar, 
                      delta_chi_w0, dt_in, dtold_in);
        }
    }

    if (maestro_verbose >= 2) {
        for (auto n = 0; n <= finest_radial_level; ++n) {
        //    do n=0,finest_radial_level
            Real max_w0 = 0.0;
            for (auto r = r_start_coord[n]; r <= r_end_coord[n]+1; ++r) {
            // do r=r_start_coord(n,1),r_end_coord(n,1)+1
                max_w0 = max(max_w0, fabs(w0_in[n+max_lev*r]));
            }
            Print() << "... max CFL of w0: " << max_w0 * dt_in / dr[n] << std::endl;
        }
        Print() << std::endl;
    }
}

void 
Maestro::Makew0Planar(RealVector& w0_in, 
                    const RealVector& w0_old, 
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
    BL_PROFILE_VAR("Maestro::Makew0Planar()",Makew0Planar);

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

    std::fill(w0_in.begin(), w0_in.end(), 0.);

    const int max_lev = max_radial_level+1;

    // local variables 
    RealVector psi_planar(nr_fine);

    // Compute w0 on edges at level n
    for (auto n = 0; n <= max_radial_level; ++n) {

        std::fill(psi_planar.begin(), psi_planar.end(), 0.);
        int base_cutoff_density_coord = 0;
        get_base_cutoff_density_coord(n, &base_cutoff_density_coord);

        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {

            // Print() << "r_start_coord = " << r_start_coord[n+max_lev*j] << ' ' << r_end_coord[n+max_lev*j] <<std::endl;

            // Print() << "base_cutoff_density_coord[n] = " << base_cutoff_density_coord[n] << std::endl;

            if (n == 0) {
                // Initialize new w0 at bottom of coarse base array to 0.0.
                w0_in[0] = 0.0;
            } else {
                // Obtain the starting value of w0 from the coarser grid
                w0_in[n+max_lev*r_start_coord[n+max_lev*j]] = w0_in[n-1+max_lev*r_start_coord[n+max_lev*j]/2];
            }

            // compute psi for level n
            //   do r = r_start_coord[n+max_lev*j], r_end_coord[n+max_lev*j]
            for (auto r = r_start_coord[n+max_lev*j]; 
                r <= r_end_coord[n+max_lev*j]; ++r) {
                if (r < base_cutoff_density_coord) {
                    psi_planar[r] = etarho_cc[n+max_lev*r] * fabs(grav_const);
                }
            }

            for (auto r = r_start_coord[n+max_lev*j]+1; 
                r <= r_end_coord[n+max_lev*j]+1; ++r) {

                Real gamma1bar_p0_avg = (gamma1bar_old_in[n+max_lev*(r-1)]
                    + gamma1bar_new_in[n+max_lev*(r-1)]) *
                    (p0_old_in[n+max_lev*(r-1)] + 
                    p0_new_in[n+max_lev*(r-1)])/4.0;

                if (r < base_cutoff_density_coord) {
                    if (is_predictor) {
                        delta_chi_w0[n+max_lev*(r-1)] = dpdt_factor * 
                            p0_minus_peosbar[n+max_lev*(r-1)] / 
                            (gamma1bar_old_in[n+max_lev*(r-1)]*
                            p0_old_in[n+max_lev*(r-1)]*dt);
                    } else {
                        delta_chi_w0[n+max_lev*(r-1)] += dpdt_factor *
                            p0_minus_peosbar[n+max_lev*(r-1)] / 
                            (gamma1bar_new_in[n+max_lev*(r-1)]*
                            p0_new_in[n+max_lev*(r-1)]*dt);
                    }
                } else {
                    delta_chi_w0[n+max_lev*(r-1)] = 0.0;
                }

                // if (gamma1bar_p0_avg == 0) {
                    // Print() << "gamma1bar_p0_avg = " << gamma1bar_p0_avg << " gamma1bar_old, new = " << gamma1bar_old_in[n+max_lev*(r-1)] << ' ' << gamma1bar_new_in[n+max_lev*(r-1)] <<std::endl;
                // }

                w0_in[n+max_lev*r] = w0_in[n+max_lev*(r-1)]
                    + Sbar_in[n+max_lev*(r-1)] * dr[n] 
                    - psi_planar[r-1] / gamma1bar_p0_avg * dr[n] 
                    - delta_chi_w0[n+max_lev*(r-1)] * dr[n];
            }

            if (n > 0) {
                // Compare the difference between w0 at top of level n to
                // the corresponding point on level n-1
                Real offset = w0_in[n+max_lev*(r_end_coord[n+max_lev*j]+1)]
                    - w0_in[n-1+max_lev*(r_end_coord[n+max_lev*j]+1)/2];

                //  do i=n-1,0,-1
                for (auto i = n-1; i >= 0; --i) {

                    int refrat = pow(2, n-i);

                    // Restrict w0 from level n to level i
                    // do r=r_start_coord[n+max_lev*j],r_end_coord[n+max_lev*j]+1
                    for (auto r = r_start_coord[n + max_lev*j]; r <= r_end_coord[n + max_lev*j]+1; ++r) {
                        if (r % refrat == 0) {
                            w0_in[n+max_lev*r/refrat] = w0_in[n+max_lev*r];
                        }
                    }

                    // Offset the w0 on level i above the top of level n
                    for (auto r = (r_end_coord[n+max_lev*j]+1)/refrat+1; r <= nr[i]; ++r) {
                    // do r=(r_end_coord[n+max_lev*j]+1)/refrat+1,nr(i)
                        w0_in[i+max_lev*r] += offset;
                    }
                }
            }
        }
    }

    // zero w0 where there is no corresponding full state array
    for (auto n = 1; n <= max_radial_level; ++n) {
    //    do n=1,max_radial_level
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
            if (j == numdisjointchunks[n]) {
            // do r=r_end_coord[n+max_lev*j]+2,nr(n)
                for (auto r = r_end_coord[n+max_lev*j]+2; 
                     r <= nr[n]; ++r) {
                    w0_in[n+max_lev*r] = 0.0;
                }
            } else {
                // do r=r_end_coord[n+max_lev*j]+2,r_start_coord(n,j+1)-1
                for (auto r = r_end_coord[n+max_lev*j]+2; 
                     r <= r_start_coord[n+max_lev*(j+1)]-1; ++r) {
                    w0_in[n+max_lev*r] = 0.0;
                }
            }
        }
    }
    // call restrict_base(w0,0)
    RestrictBase(w0_in, false);
    // call fill_ghost_base(w0,0)
    FillGhostBase(w0_in, false);

    for (auto n = 0; n <= max_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {

            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            Real dt_avg = 0.5 * (dt_in + dtold_in);
            //   do r=r_start_coord[n+max_lev*j],r_end_coord[n+max_lev*j]
            for (auto r = r_start_coord[n + max_lev*j]; 
                 r <= r_end_coord[n + max_lev*j]; ++r) {

                Real w0_old_cen = 0.5 * (w0_old[n+max_lev*r] + 
                    w0_old[n+max_lev*(r+1)]);
                Real w0_new_cen = 0.5 * (w0_in[n+max_lev*r] + 
                    w0_in[n+max_lev*(r+1)]);
                Real w0_avg = 0.5 * (dt * w0_old_cen 
                    + dtold *  w0_new_cen) / dt_avg;
                Real div_avg = 0.5 * (dt *(w0_old[n+max_lev*(r+1)]
                    - w0_old[n+max_lev*r]) + dtold * (w0_in[n+max_lev*(r+1)] 
                    - w0_in[n+max_lev*r])) / dt_avg;
                w0_force[n+max_lev*r] = (w0_new_cen
                    - w0_old_cen)/dt_avg 
                    + w0_avg*div_avg/dr[n];
            }
        }
    }

    // call restrict_base(w0_force,1)
    // call fill_ghost_base(w0_force,1)
    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}

void 
Maestro::Makew0PlanarVarg(RealVector& w0_in, 
                        const RealVector& w0_old, 
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

    // The planar 1/r**2 gravity constraint equation is solved
    // by calling the tridiagonal solver, just like spherical.
    // This is accomplished by putting all the requisite data
    // on the finest basestate grid, solving for w0, and then
    // restricting w0 back down to the coarse grid.

    // 1) allocate the finely-gridded temporary basestate arrays
    RealVector w0_fine(nr[finest_radial_level]+1);
    RealVector w0bar_fine(nr[finest_radial_level]+1);
    RealVector deltaw0_fine(nr[finest_radial_level]+1);
    RealVector p0_old_fine(nr[finest_radial_level]);
    RealVector p0_new_fine(nr[finest_radial_level]);
    RealVector p0_nph_fine(nr[finest_radial_level]);
    RealVector rho0_old_fine(nr[finest_radial_level]);
    RealVector rho0_new_fine(nr[finest_radial_level]);
    RealVector rho0_nph_fine(nr[finest_radial_level]);
    RealVector gamma1bar_old_fine(nr[finest_radial_level]);
    RealVector gamma1bar_new_fine(nr[finest_radial_level]);
    RealVector gamma1bar_nph_fine(nr[finest_radial_level]);
    RealVector p0_minus_peosbar_fine(nr[finest_radial_level]);
    RealVector etarho_cc_fine(nr[finest_radial_level]);
    RealVector Sbar_in_fine(nr[finest_radial_level]);
    RealVector grav_edge_fine(nr[finest_radial_level]+1);

    // 2) copy the data into the temp, uniformly-gridded basestate arrays.
    ProlongBasetoUniform(p0_old,p0_old_fine);
    ProlongBasetoUniform(p0_new,p0_new_fine);
    ProlongBasetoUniform(rho0_old,rho0_old_fine);
    ProlongBasetoUniform(rho0_new,rho0_new_fine);
    ProlongBasetoUniform(gamma1bar_old,gamma1bar_old_fine);
    ProlongBasetoUniform(gamma1bar_new,gamma1bar_new_fine);
    ProlongBasetoUniform(p0_minus_peosbar,p0_minus_peosbar_fine);
    ProlongBasetoUniform(etarho_cc,etarho_cc_fine);
    ProlongBasetoUniform(Sbar_in,Sbar_in_fine);

    // create time-centered base-state quantities
    for (auto r = 0; r < nr[finest_radial_level]; ++r) {
        p0_nph_fine[r] = 0.5*(p0_old_fine[r] + p0_new_fine[r]);
        rho0_nph_fine[r] = 0.5*(rho0_old_fine[r] + rho0_new_fine[r]);
        gamma1bar_nph_fine[r] = 0.5*(gamma1bar_old_fine[r] + gamma1bar_new_fine[r]);
    }


    // 3) solve to w0bar -- here we just take into account the Sbar and
    //    volume discrepancy terms
    // lower boundary condition
    w0bar_fine[0] = 0.0;

    // do r=1,nr(finest_radial_level)
    for (auto r = 1; r <= nr[finest_radial_level]; ++r) {
        Real gamma1bar_p0_avg = gamma1bar_nph_fine[r-1] * p0_nph_fine[r-1];

        Real volume_discrepancy = (r-1 < fine_base_density_cutoff_coord) ? dpdt_factor * p0_minus_peosbar_fine[r-1]/dt : 0.0;

        w0bar_fine[r] =  w0bar_fine[r-1] + 
            Sbar_in_fine[r-1] * dr[finest_radial_level] 
            - (volume_discrepancy / gamma1bar_p0_avg ) * dr[finest_radial_level];

    }

    // 4) get the edge-centered gravity on the uniformly-gridded
    // basestate arrays
    Abort("make_w0.f90: need to write make_grav_edge_uniform");
    //    call make_grav_edge_uniform(grav_edge_fine, rho0_nph_fine)


    // 5) solve for delta w0
    std::fill(deltaw0_fine.begin(), deltaw0_fine.end(), 0.);

    // this takes the form of a tri-diagonal matrix:
    // A_j (dw_0)_{j-3/2} +
    // B_j (dw_0)_{j-1/2} +
    // C_j (dw_0)_{j+1/2} = F_j

    RealVector A(nr[finest_radial_level]+1);
    RealVector B(nr[finest_radial_level]+1);
    RealVector C(nr[finest_radial_level]+1);
    RealVector u(nr[finest_radial_level]+1);
    RealVector F(nr[finest_radial_level]+1);

    std::fill(A.begin(), A.end(), 0.);
    std::fill(B.begin(), B.end(), 0.);
    std::fill(C.begin(), C.end(), 0.);
    std::fill(F.begin(), F.end(), 0.);
    std::fill(u.begin(), u.end(), 0.);

    for (auto r = 1; r <= fine_base_density_cutoff_coord; ++r) {
        A[r] = gamma1bar_nph_fine[r-1] * p0_nph_fine[r-1];
        A[r] /= dr[finest_radial_level]*dr[finest_radial_level];

        Real dpdr = (p0_nph_fine[r]-p0_nph_fine[r-1])/dr[finest_radial_level];

        B[r] = -(gamma1bar_nph_fine[r-1] * p0_nph_fine[r-1] + 
            gamma1bar_nph_fine[r] * p0_nph_fine[r]) 
            / (dr[finest_radial_level]*dr[finest_radial_level]);
        B[r] -= 2.0 * dpdr / (r_edge_loc[finest_radial_level+max_lev*r]);

        C[r] = gamma1bar_nph_fine[r] * p0_nph_fine[r];
        C[r] /= dr[finest_radial_level]*dr[finest_radial_level];

        F[r] = 2.0 * dpdr * w0bar_fine[r] / 
            r_edge_loc[finest_radial_level+max_lev*r] -
            grav_edge_fine[r] * (etarho_cc_fine[r] - etarho_cc_fine[r-1]) / 
            dr[finest_radial_level];
    }

    // Lower boundary
    A[0] = 0.0;
    B[0] = 1.0;
    C[0] = 0.0;
    F[0] = 0.0;

    // Upper boundary
    A[fine_base_density_cutoff_coord+1] = -1.0;
    B[fine_base_density_cutoff_coord+1] = 1.0;
    C[fine_base_density_cutoff_coord+1] = 0.0;
    F[fine_base_density_cutoff_coord+1] = 0.0;

    // Call the tridiagonal solver
    Tridiag(A, B, C, F, u, fine_base_density_cutoff_coord+2);

    for (auto r = 1; r <= fine_base_density_cutoff_coord+1; ++r) {
        deltaw0_fine[r] = u[r];
    }

    for (auto r = fine_base_density_cutoff_coord+2; r <= nr[finest_radial_level]; ++r) {
        deltaw0_fine[r] = deltaw0_fine[fine_base_density_cutoff_coord+1];
    }

    // 6) compute w0 = w0bar + deltaw0
    for (auto r = 0; r < w0_fine.size(); ++r) {
        w0_fine[r] = w0bar_fine[r] + deltaw0_fine[r];
        w0_in[finest_radial_level+max_lev*r] = w0_fine[r];
    }

    // 7) fill the multilevel w0 array from the uniformly-gridded w0 we
    // just solved for.  Here, we make the coarse edge underneath equal
    // to the fine edge value.
    for (auto n = finest_radial_level; n >= 1; --n) {
    // do n = finest_radial_level, 1, -1
        for (auto r = 0; r <= nr[n]; n+=2) {
        //    do r = 0, nr(n), 2
            w0_in[n-1+max_lev*r/2] = w0_in[n+max_lev*r];
        }
    }

    // 8) zero w0 where there is no corresponding full state array
    // do n=1,finest_radial_level
    for (auto n = 1; n <= finest_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
        //    do j=1,numdisjointchunks(n)
            if (j == numdisjointchunks[n]) {
                //  do r=r_end_coord(n,j)+2,nr(n)
                for (auto r = r_end_coord[n+max_lev*j]+2; r <= nr[n]; ++r) {
                    w0_in[n+max_lev*r] = 0.0;
                }
            } else {
                for (auto r = r_end_coord[n+max_lev*j]+2; r < r_start_coord[n+max_lev*(j+1)]; ++r) {
                    w0_in[n+max_lev*r] = 0.0;
                }
            }
        }
    }

    // call restrict_base(w0,0)
    // call fill_ghost_base(w0,0)
    RestrictBase(w0_in, false);
    FillGhostBase(w0_in, false);

    // compute the forcing terms
    for (auto n = 0; n <= finest_radial_level; ++n) {
        for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
        // do n=0,finest_radial_level
        //    do j=1,numdisjointchunks(n)

            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            Real dt_avg = 0.5 * (dt + dtold);
            //   do r=r_start_coord(n,j),r_end_coord(n,j)
            for (auto r = r_start_coord[n+max_lev*j]; r <=r_end_coord[n+max_lev*j]; ++r) {
                Real w0_old_cen = 0.5 * (w0_old[n+max_lev*r] + w0_old[n+max_lev*(r+1)]);
                Real w0_new_cen = 0.5 * (w0_in[n+max_lev*r] + w0_in[n+max_lev*(r+1)]);
                Real w0_avg = 0.5 * (dt * w0_old_cen + dtold *  w0_new_cen) / dt_avg;
                Real div_avg = 0.5 * (dt * (w0_old[n+max_lev*(r+1)]-w0_old[n+max_lev*r]) + 
                    dtold * (w0_in[n+max_lev*(r+1)]-w0_in[n+max_lev*r])) / dt_avg;
                w0_force[n+max_lev*r] = (w0_new_cen-w0_old_cen)/dt_avg + w0_avg*div_avg/dr[n];
            }
        }
    }

    // call restrict_base(w0_force,1)
    // call fill_ghost_base(w0_force,1)
    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}

void 
Maestro::Makew0Sphr(RealVector& w0_in, 
                    const RealVector& w0_old, 
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
    RealVector gamma1bar_nph(nr_fine);
    RealVector p0_nph(nr_fine);
    RealVector A(nr_fine+1);
    RealVector B(nr_fine+1);
    RealVector C(nr_fine+1);
    RealVector u(nr_fine+1);
    RealVector F(nr_fine+1);
    RealVector w0_from_Sbar(nr_fine+1);
    RealVector rho0_nph(nr_fine);
    RealVector grav_edge(nr_fine+1);

    get_base_cutoff_density(&base_cutoff_density);
    int base_cutoff_density_coord = 0;
    get_base_cutoff_density_coord(0, &base_cutoff_density_coord);

    Print() << "dr = " << dr[0] << std::endl;

    const int max_lev = max_radial_level+1;

    // create time-centered base-state quantities
    for (auto r = 0; r < nr_fine; ++r) {
        p0_nph[r] = 0.5*(p0_old_in[r] + p0_new_in[r]);
        rho0_nph[r] = 0.5*(rho0_old_in[r] + rho0_new_in[r]);
        gamma1bar_nph[r] = 0.5*(gamma1bar_old_in[r] + gamma1bar_new_in[r]);
    }

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.

    w0_from_Sbar[0] = 0.0;

    for (auto r = 1; r <= nr_fine; ++r) {

        Real volume_discrepancy = rho0_old_in[r-1] > base_cutoff_density ? 
            dpdt_factor * p0_minus_peosbar[r-1]/dt : 0.0;

        w0_from_Sbar[r] = w0_from_Sbar[r-1] + 
            dr[0] * Sbar_in[r-1] * r_cc_loc[r-1]*r_cc_loc[r-1] - 
            dr[0]* volume_discrepancy * r_cc_loc[r-1]*r_cc_loc[r-1] 
            / (gamma1bar_nph[r-1]*p0_nph[r-1]);

    }

    for (auto r = 1; r <= nr_fine; ++r) {
        w0_from_Sbar[r] /= (r_edge_loc[r]*r_edge_loc[r]);
    }

    // make the edge-centered gravity
    MakeGravEdge(grav_edge, rho0_nph);

    // NOTE:  now we solve for the remainder, (r^2 * delta w0)
    // this takes the form of a tri-diagonal matrix:
    // A_j (r^2 dw_0)_{j-3/2} +
    // B_j (r^2 dw_0)_{j-1/2} +
    // C_j (r^2 dw_0)_{j+1/2} = F_j
    std::fill(A.begin(), A.end(), 0.);
    std::fill(B.begin(), B.end(), 0.);
    std::fill(C.begin(), C.end(), 0.);
    std::fill(F.begin(), F.end(), 0.);
    std::fill(u.begin(), u.end(), 0.);

    // Note that we are solving for (r^2 delta w0), not just w0.

    int max_cutoff = min(base_cutoff_density_coord, nr_fine-1);
    
    for (auto r = 1; r <= max_cutoff; ++r) {
        A[r] = gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc[r-1]*r_cc_loc[r-1]);
        A[r] /= dr[0]*dr[0];

        B[r] = -( gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc[r-1]*r_cc_loc[r-1])
                + gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc[r]*r_cc_loc[r]) ) 
                / (dr[0]*dr[0]);

        Real dpdr = (p0_nph[r] - p0_nph[r-1]) / dr[0];

        B[r] -= 4.0 * dpdr / (r_edge_loc[r]*r_edge_loc[r]*r_edge_loc[r]);

        C[r] = gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc[r]*r_cc_loc[r]);
        C[r] /= dr[0]*dr[0];

        F[r] = 4.0 * dpdr * w0_from_Sbar[r] / r_edge_loc[r] - 
                grav_edge[r] * (r_cc_loc[r]*r_cc_loc[r] * etarho_cc[r] - 
                r_cc_loc[r-1]*r_cc_loc[r-1] * etarho_cc[r-1]) / 
                (dr[0] * r_edge_loc[r]*r_edge_loc[r]) - 
                4.0 * M_PI * Gconst * 0.5 * 
                (rho0_nph[r] + rho0_nph[r-1]) * etarho_ec[r];
    }

    // Lower boundary
    A[0] = 0.0;
    B[0] = 1.0;
    C[0] = 0.0;
    F[0] = 0.0;

    // Upper boundary
    A[max_cutoff+1] = -1.0;
    B[max_cutoff+1] = 1.0;
    C[max_cutoff+1] = 0.0;
    F[max_cutoff+1] = 0.0;

    // Call the tridiagonal solver
    Tridiag(A, B, C, F, u, max_cutoff+2);

    w0_in[0] = w0_from_Sbar[0];

    // do r=1,max_cutoff+1
    for (auto r = 1; r <= max_cutoff+1; ++r) {
        w0_in[r] = u[r] / (r_edge_loc[r]*r_edge_loc[r]) + w0_from_Sbar[r];
    }

    // do r=max_cutoff+2,nr_fine
    for (auto r = max_cutoff+2; r <= nr_fine; ++r) {
        w0_in[r] = w0_in[max_cutoff+1] * r_edge_loc[max_cutoff+1]*r_edge_loc[max_cutoff+1]/(r_edge_loc[r]*r_edge_loc[r]);
    }

    // Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    Real dt_avg = 0.5 * (dt_in + dtold_in);

    // do r = 0,nr_fine-1
    for (auto r = 0; r < nr_fine; ++r) {
        Real w0_old_cen = 0.5 * (w0_old[r] + w0_old[r+1]);
        Real w0_new_cen = 0.5 * (w0_in[r] + w0_in[r+1]);
        Real w0_avg = 0.5 * (dt_in *  w0_old_cen + dtold_in *  w0_new_cen) / dt_avg;
        Real div_avg = 0.5 * (dt_in * (w0_old[r+1]-w0_old[r]) + dtold_in * (w0_in[r+1]-w0_in[r])) / dt_avg;
        w0_force[r] = (w0_new_cen-w0_old_cen) / dt_avg + w0_avg * div_avg / dr[0];
    }
}

void 
Maestro::Makew0SphrIrreg(RealVector& w0_in, 
                        const RealVector& w0_old, 
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
    RealVector gamma1bar_nph(nr_fine);
    RealVector p0_nph(nr_fine);
    RealVector A(nr_fine+1);
    RealVector B(nr_fine+1);
    RealVector C(nr_fine+1);
    RealVector u(nr_fine+1);
    RealVector F(nr_fine+1);
    RealVector w0_from_Sbar(nr_fine+1);
    RealVector rho0_nph(nr_fine);
    RealVector grav_edge(nr_fine+1);

    get_base_cutoff_density(&base_cutoff_density);
    int base_cutoff_density_coord = 0;
    get_base_cutoff_density_coord(0, &base_cutoff_density_coord);

    const int max_lev = max_radial_level+1;

    // create time-centered base-state quantities
    for (auto r = 0; r < nr_fine; ++r) {
        p0_nph[r] = 0.5*(p0_old_in[r] + p0_new_in[r]);
        rho0_nph[r] = 0.5*(rho0_old_in[r] + rho0_new_in[r]);
        gamma1bar_nph[r] = 0.5*(gamma1bar_old_in[r] + gamma1bar_new_in[r]);
    }

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.

    w0_from_Sbar[0] = 0.0;

    for (auto r = 1; r <= nr_fine; ++r) {

        Real volume_discrepancy = rho0_old_in[r-1] > base_cutoff_density ? 
            dpdt_factor * p0_minus_peosbar[r-1]/dt : 0.0;

        Real dr1 = r_edge_loc[max_lev*r] - r_edge_loc[max_lev*(r-1)];
        w0_from_Sbar[r] = w0_from_Sbar[r-1] + 
            dr1 * Sbar_in[r-1] * r_cc_loc[r-1]*r_cc_loc[r-1] - 
            dr1* volume_discrepancy * r_cc_loc[r-1]*r_cc_loc[r-1] 
            / (gamma1bar_nph[r-1]*p0_nph[r-1]);
    }

    for (auto r = 1; r <= nr_fine; ++r) {
        w0_from_Sbar[r] /= (r_edge_loc[r]*r_edge_loc[r]);
    }

    // make the edge-centered gravity
    MakeGravEdge(grav_edge, rho0_nph);

    // NOTE:  now we solve for the remainder, (r^2 * delta w0)
    // this takes the form of a tri-diagonal matrix:
    // A_j (r^2 dw_0)_{j-3/2} +
    // B_j (r^2 dw_0)_{j-1/2} +
    // C_j (r^2 dw_0)_{j+1/2} = F_j
    std::fill(A.begin(), A.end(), 0.);
    std::fill(B.begin(), B.end(), 0.);
    std::fill(C.begin(), C.end(), 0.);
    std::fill(F.begin(), F.end(), 0.);
    std::fill(u.begin(), u.end(), 0.);

    // Note that we are solving for (r^2 delta w0), not just w0.

    int max_cutoff = base_cutoff_density_coord;
    
    for (auto r = 1; r <= max_cutoff; ++r) {
        Real dr1 = r_edge_loc[max_lev*r] - r_edge_loc[max_lev*(r-1)];
        Real dr2 = r_edge_loc[max_lev*(r+1)] - r_edge_loc[max_lev*r];
        Real dr3 = r_cc_loc[max_lev*r] - r_cc_loc[max_lev*(r-1)];

        A[r] = gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc[r-1]*r_cc_loc[r-1]);
        A[r] /= dr1*dr3;

        B[r] = -( gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc[r-1]*r_cc_loc[r-1]*dr1) 
                + gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc[r]*r_cc_loc[r]*dr2) ) 
                / dr3;

        Real dpdr = (p0_nph[r] - p0_nph[r-1]) / dr3;

        B[r] -= 4.0 * dpdr / (r_edge_loc[r]*r_edge_loc[r]*r_edge_loc[r]);

        C[r] = gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc[r]*r_cc_loc[r]);
        C[r] /= dr2*dr3;

        F[r] = 4.0 * dpdr * w0_from_Sbar[r] / r_edge_loc[r] - 
                grav_edge[r] * (r_cc_loc[r]*r_cc_loc[r] * etarho_cc[r] - 
                r_cc_loc[r-1]*r_cc_loc[r-1] * etarho_cc[r-1]) / 
                (dr3 * r_edge_loc[r]*r_edge_loc[r]) - 
                4.0 * M_PI * Gconst * 0.5 * 
                (rho0_nph[r] + rho0_nph[r-1]) * etarho_ec[r];
    }

    // Lower boundary
    A[0] = 0.0;
    B[0] = 1.0;
    C[0] = 0.0;
    F[0] = 0.0;

    // Upper boundary
    A[max_cutoff+1] = -1.0;
    B[max_cutoff+1] = 1.0;
    C[max_cutoff+1] = 0.0;
    F[max_cutoff+1] = 0.0;

    // Call the tridiagonal solver
    Tridiag(A, B, C, F, u, max_cutoff+2);

    w0_in[0] = w0_from_Sbar[0];

    // do r=1,max_cutoff+1
    for (auto r = 1; r <= max_cutoff+1; ++r) {
        w0_in[r] = u[r] / (r_edge_loc[r]*r_edge_loc[r]) + w0_from_Sbar[r];
    }

    // do r=max_cutoff+2,nr_fine
    for (auto r = max_cutoff+2; r <= nr_fine; ++r) {
        w0_in[r] = w0_in[max_cutoff+1] * r_edge_loc[max_cutoff+1]*r_edge_loc[max_cutoff+1]/(r_edge_loc[r]*r_edge_loc[r]);
    }

    // Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    Real dt_avg = 0.5 * (dt_in + dtold_in);

    // do r = 0,nr_fine-1
    for (auto r = 0; r < nr_fine; ++r) {
        Real dr1 = r_edge_loc[max_lev*r] - r_edge_loc[max_lev*(r-1)];
        Real w0_old_cen = 0.5 * (w0_old[r] + w0_old[r+1]);
        Real w0_new_cen = 0.5 * (w0_in[r] + w0_in[r+1]);
        Real w0_avg = 0.5 * (dt_in *  w0_old_cen + dtold_in *  w0_new_cen) / dt_avg;
        Real div_avg = 0.5 * (dt_in * (w0_old[r+1]-w0_old[r]) + dtold_in * (w0_in[r+1]-w0_in[r])) / dt_avg;
        w0_force[r] = (w0_new_cen-w0_old_cen) / dt_avg + w0_avg * div_avg / dr1;
    }
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