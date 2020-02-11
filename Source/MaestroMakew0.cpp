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

    if (!spherical) {
        if (do_planar_invsq_grav || do_2d_planar_octant) {
            Makew0PlanarVarg();
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
            Makew0SphrIrreg();
        } else {
            Makew0Sphr(w0_in, w0_old, w0_force, Sbar_in, 
                        rho0_old_in, rho0_new_in,
                        p0_old_in, p0_new_in, 
                        gamma1bar_old_in, gamma1bar_new_in, 
                        p0_minus_peosbar, 
                        delta_chi_w0, dt_in, dtold_in,
                        is_predictor);
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
    RealVector w0_old_cen(max_lev*nr_fine);
    RealVector w0_new_cen(max_lev*nr_fine);
    RealVector psi_planar(nr_fine);

    // Compute w0 on edges at level n
    for (auto n = 0; n <= max_radial_level; ++n) {

        std::fill(psi_planar.begin(), psi_planar.end(), 0.);

        for (auto j = 0; j < numdisjointchunks[n]; ++j) {

            if (n == 0) {
                // Initialize new w0 at bottom of coarse base array to 0.0.
                w0_in[0] = 0.0;
            } else {
                // Obtain the starting value of w0 from the coarser grid
                w0_in[n+max_lev*r_start_coord[n+max_lev*j]] = w0_in[n-1+max_lev*r_start_coord[n+max_lev*j]/2];
            }

            // compute psi for level n
            //   do r = r_start_coord[n+max_lev*j], r_end_coord[n+max_lev*j]
            for (auto r = r_start_coord[n + max_lev*j]; 
                j <= r_end_coord[n + max_lev*j]; ++j) {
                if (r < base_cutoff_density_coord[n]) {
                    psi_planar[r] = etarho_cc[n + max_lev*r] * fabs(grav_const);
                }
            }

            for (auto r = r_start_coord[n + max_lev*j]+1; 
                r <= r_end_coord[n + max_lev*j]+1; ++r) {

                Real gamma1bar_p0_avg = (gamma1bar_old_in[n+max_lev*(r-1)]
                    + gamma1bar_new_in[n+max_lev*(r-1)]) *
                    (p0_old_in[n+max_lev*(r-1)] + 
                    p0_new_in[n+max_lev*(r-1)])/4.0;

                if (r < base_cutoff_density_coord[n]) {
                    if (is_predictor == 1) {
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
        for (auto j = 0; j < numdisjointchunks[n]; ++j) {
            if (j == numdisjointchunks[n]-1) {
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
        for (auto j = 0; j < numdisjointchunks[n]; ++j) {

            // Compute the forcing term in the base state velocity
            // equation, - 1/rho0 grad pi0
            Real dt_avg = 0.5 * (dt_in + dtold_in);
            //   do r=r_start_coord[n+max_lev*j],r_end_coord[n+max_lev*j]
            for (auto r = r_start_coord[n + max_lev*j]; 
                 r <= r_end_coord[n + max_lev*j]; ++r) {

                // Print() << "j , r, n = " << j << ' ' << r << ' ' << n <<std::endl;
                w0_old_cen[n+max_lev*r] = 0.5 * (w0_old[n+max_lev*r] + 
                    w0_old[n+max_lev*(r+1)]);
                w0_new_cen[n+max_lev*r] = 0.5 * (w0_in[n+max_lev*r] + 
                    w0_in[n+max_lev*(r+1)]);
                Real w0_avg = 0.5 * (dt * w0_old_cen[n+max_lev*r] 
                    + dtold *  w0_new_cen[n+max_lev*r]) / dt_avg;
                Real div_avg = 0.5 * (dt *(w0_old[n+max_lev*(r+1)]
                    - w0_old[n+max_lev*r]) + dtold * (w0_in[n+max_lev*(r+1)] 
                    - w0_in[n+max_lev*r])) / dt_avg;
                w0_force[n+max_lev*r] = (w0_new_cen[n+max_lev*r] 
                    - w0_old_cen[n+max_lev*r])/dt_avg 
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
                    const Real dt_in, const Real dtold_in, 
                    const bool is_predictor) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0Sphr()",Makew0Sphr);

    // local variables 
    RealVector w0_old_cen(nr_fine);
    RealVector w0_new_cen(nr_fine);
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

    // create time-centered base-state quantities
    for (auto r = 0; r < nr_fine; ++r) {
    // do r=0,nr_fine-1
        p0_nph[r] = 0.5*(p0_old_in[r] + p0_new_in[r]);
        rho0_nph[r] = 0.5*(rho0_old_in[r] + rho0_new_in[r]);
        gamma1bar_nph[r] = 0.5*(gamma1bar_old_in[r] + gamma1bar_new_in[r]);
    }

    // NOTE: We first solve for the w0 resulting only from Sbar,
    //      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    //      (r^2 Sbar).  Then we will solve for the update, delta w0.

    // w0_from_Sbar = 0.0
    std::fill(w0_from_Sbar.begin(), w0_from_Sbar.end(), 0.);

    for (auto r = 1; r <= nr_fine; ++r) {

        Real volume_discrepancy = rho0_old_in[r-1] > base_cutoff_density ? dpdt_factor * p0_minus_peosbar[r-1]/dt : 0.0;

        w0_from_Sbar[r] = w0_from_Sbar[r-1] + dr[0] * Sbar_in[r-1] * r_cc_loc[r-1]*r_cc_loc[r-1] - 
                dr[0]* volume_discrepancy * r_cc_loc[r-1]*r_cc_loc[r-1] / (gamma1bar_nph[r-1]*p0_nph[r-1]);

    }

    for (auto r = 1; r <= nr_fine; ++r) {
        w0_from_Sbar[r] = w0_from_Sbar[r] / (r_edge_loc[r]*r_edge_loc[r]);
    }

    // make the edge-centered gravity
    // call make_grav_edge(grav_edge,rho0_nph,r_edge_loc)
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

    int max_cutoff = min(base_cutoff_density_coord[0], nr_fine-1);
    
    for (auto r = 1; r <= max_cutoff; ++r) {
        A[r] = gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc[r-1]*r_cc_loc[r-1]);
        A[r] /= dr[0]*dr[0];

        B[r] = -( gamma1bar_nph[r-1] * p0_nph[r-1] / (r_cc_loc[r-1]*r_cc_loc[r-1])
                + gamma1bar_nph[r] * p0_nph[r] / (r_cc_loc[r]*r_cc_loc[r]) ) / (dr[0]*dr[0]);

        Real dpdr = (p0_nph[r]-p0_nph[r-1])/dr[0];

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
        w0_old_cen[r] = 0.5 * (w0_old[r] + w0_old[r+1]);
        w0_new_cen[r] = 0.5 * (w0_in[r] + w0_in[r+1]);
        Real w0_avg = 0.5 * (dt_in *  w0_old_cen[r] + dtold_in *  w0_new_cen[r]) / dt_avg;
        Real div_avg = 0.5 * (dt_in * (w0_old[r+1]-w0_old[r]) + dtold_in * (w0_in[r+1]-w0_in[r])) / dt_avg;
        w0_force[r] = (w0_new_cen[r]-w0_old_cen[r]) / dt_avg + w0_avg * div_avg / dr[0];
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
        gam[j] = c[j-1]/bet;
        bet = b[j] - a[j]*gam[j];
        if (bet == 0) Abort("tridiag: TRIDIAG FAILED");
        u[j] = (r[j]-a[j]*u[j-1])/bet;
    }

    for (auto j = n-2; j >= 0; --j) {
        u[j] -= gam[j+1]*u[j+1];
    }
}