#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::Makew0(const RealVector& w0_old, RealVector& w0_force, 
                const RealVector& Sbar_in, const RealVector& p0_minus_peosbar, 
                RealVector& delta_chi_w0, const bool is_predictor) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0()",Makew0);

    std::fill(w0_force.begin(), w0_force.end(), 0.);

    if (!spherical) {
        if (do_planar_invsq_grav || do_2d_planar_octant) {
            Makew0PlanarVarg();
        } else {
            Makew0Planar(w0_old, w0_force, Sbar_in, p0_minus_peosbar, 
                         delta_chi_w0, is_predictor);
        }
    } else {
        if (use_exact_base_state) {
            Makew0SphrIrreg();
        } else {
            Makew0Sphr();
        }
    }

    if (maestro_verbose >= 2) {
        for (auto n = 0; n <= finest_radial_level; ++n) {
        //    do n=0,finest_radial_level
            Real max_w0 = 0.0;
            for (auto r = r_start_coord[n]; r <= r_end_coord[n]+1; ++r) {
            // do r=r_start_coord(n,1),r_end_coord(n,1)+1
                max_w0 = max(max_w0, fabs(w0[n+(max_radial_level+1)*r]));
            }
            Print() << "... max CFL of w0: " << max_w0 * dt / dr[n] << std::endl;
        }
        Print() << std::endl;
    }
}

void 
Maestro::Makew0Planar(const RealVector& w0_old, RealVector& w0_force, 
                      const RealVector& Sbar_in, 
                      const RealVector& p0_minus_peosbar, 
                      RealVector& delta_chi_w0, const bool is_predictor) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Makew0Planar()",Makew0Planar);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ! Multilevel Outline
    // !
    // ! Compute w0 at level 1 only
    // ! Initialize new w0 at bottom of coarse base array to zero.
    // ! do n=1,finest_radial_level
    // !   Compute w0 on edges at level n
    // !   Obtain the starting value of w0 from the coarser grid
    // !   if n>1, compare the difference between w0 at top of level n to the
    // !           corresponding point on level n-1
    // !   do i=n-1,1,-1
    // !     Restrict w0 from level n to level i
    // !     Offset the w0 on level i above the top of level n
    // !   }
    // ! }
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::fill(w0.begin(), w0.end(), 0.);

    // local variables 
    RealVector w0_old_cen((max_radial_level+1)*nr_fine);
    RealVector w0_new_cen((max_radial_level+1)*nr_fine);
    RealVector psi_planar(nr_fine);

    // ! Compute w0 on edges at level n
    for (auto n = 0; n <= max_radial_level; ++n) {

        std::fill(psi_planar.begin(), psi_planar.end(), 0.);

        for (auto j = 0; j < numdisjointchunks[n]; ++j) {

            if (n == 0) {
                //  ! Initialize new w0 at bottom of coarse base array to zero.
                w0[0] = 0.0;
            } else {
                //  ! Obtain the starting value of w0 from the coarser grid
                w0[n+(max_radial_level+1)*r_start_coord[n+(max_radial_level+1)*j]] = w0[n-1+(max_radial_level+1)*r_start_coord[n+(max_radial_level+1)*j]/2];
            }

            //   ! compute psi for level n
            //   do r = r_start_coord[n+(max_radial_level+1)*j], r_end_coord[n+(max_radial_level+1)*j]
            for (auto r = r_start_coord[n + (max_radial_level+1)*j]; 
                j <= r_end_coord[n + (max_radial_level+1)*j]; ++j) {
                if (r < base_cutoff_density_coord[n]) {
                    psi_planar[r] = etarho_cc[n + (max_radial_level+1)*r] * fabs(grav_const);
                }
            }

            for (auto r = r_start_coord[n + (max_radial_level+1)*j]+1; 
                r <= r_end_coord[n + (max_radial_level+1)*j]+1; ++r) {

                Real gamma1bar_p0_avg = (gamma1bar_old[n+(max_radial_level+1)*(r-1)]
                    + gamma1bar_new[n+(max_radial_level+1)*(r-1)]) *
                    (p0_old[n+(max_radial_level+1)*(r-1)] + 
                    p0_new[n+(max_radial_level+1)*(r-1)])/4.0;

                if (r < base_cutoff_density_coord[n]) {
                    if (is_predictor == 1) {
                        delta_chi_w0[n+(max_radial_level+1)*(r-1)] = dpdt_factor * 
                            p0_minus_peosbar[n+(max_radial_level+1)*(r-1)] / 
                            (gamma1bar_old[n+(max_radial_level+1)*(r-1)]*
                            p0_old[n+(max_radial_level+1)*(r-1)]*dt);
                    } else {
                        delta_chi_w0[n+(max_radial_level+1)*(r-1)] += dpdt_factor *
                            p0_minus_peosbar[n+(max_radial_level+1)*(r-1)] / 
                            (gamma1bar_new[n+(max_radial_level+1)*(r-1)]*
                            p0_new[n+(max_radial_level+1)*(r-1)]*dt);
                    }
                } else {
                    delta_chi_w0[n+(max_radial_level+1)*(r-1)] = 0.0;
                }

                w0[n+(max_radial_level+1)*r] = w0[n+(max_radial_level+1)*(r-1)]
                    + Sbar_in[n+(max_radial_level+1)*(r-1)] * dr[n] 
                    - psi_planar[r-1] / gamma1bar_p0_avg * dr[n] 
                    - delta_chi_w0[n+(max_radial_level+1)*(r-1)] * dr[n];
            }

            if (n > 0) {
                // Compare the difference between w0 at top of level n to
                // the corresponding point on level n-1
                Real offset = w0[n+(max_radial_level+1)*(r_end_coord[n+(max_radial_level+1)*j]+1)]
                    - w0[n-1+(max_radial_level+1)*(r_end_coord[n+(max_radial_level+1)*j]+1)/2];

                //  do i=n-1,0,-1
                for (auto i = n-1; i >= 0; --i) {

                    int refrat = pow(2, n-i);

                    // Restrict w0 from level n to level i
                    // do r=r_start_coord[n+(max_radial_level+1)*j],r_end_coord[n+(max_radial_level+1)*j]+1
                    for (auto r = r_start_coord[n + (max_radial_level+1)*j]; r <= r_end_coord[n + (max_radial_level+1)*j]+1; ++r) {
                        if (r % refrat == 0) {
                            w0[n+(max_radial_level+1)*r/refrat] = w0[n+(max_radial_level+1)*r];
                        }
                    }

                    // Offset the w0 on level i above the top of level n
                    for (auto r = (r_end_coord[n+(max_radial_level+1)*j]+1)/refrat+1; r <= nr[i]; ++r) {
                    // do r=(r_end_coord[n+(max_radial_level+1)*j]+1)/refrat+1,nr(i)
                        w0[i+(max_radial_level+1)*r] += offset;
                    }
                }
            }
        }
    }

    //    ! zero w0 where there is no corresponding full state array
    for (auto n = 1; n <= max_radial_level; ++n) {
    //    do n=1,max_radial_level
        for (auto j = 0; j < numdisjointchunks[n]; ++j) {
            if (j == numdisjointchunks[n]-1) {
            // do r=r_end_coord[n+(max_radial_level+1)*j]+2,nr(n)
                for (auto r = r_end_coord[n+(max_radial_level+1)*j]+2; 
                     r <= nr[n]; ++r) {
                    w0[n+(max_radial_level+1)*r] = 0.0;
                }
            } else {
                // do r=r_end_coord[n+(max_radial_level+1)*j]+2,r_start_coord(n,j+1)-1
                for (auto r = r_end_coord[n+(max_radial_level+1)*j]+2; 
                     r <= r_start_coord[n+(max_radial_level+1)*(j+1)]-1; ++r) {
                    w0[n+(max_radial_level+1)*r] = 0.0;
                }
            }
        }
    }

    // call restrict_base(w0,0)
    RestrictBase(w0, false);
    // call fill_ghost_base(w0,0)
    FillGhostBase(w0, false);

    for (auto n = 0; n <= max_radial_level; ++n) {
        for (auto j = 0; j < numdisjointchunks[n]; ++j) {

            //   ! Compute the forcing term in the base state velocity
            //   ! equation, - 1/rho0 grad pi0
            Real dt_avg = 0.5 * (dt + dtold);
            //   do r=r_start_coord[n+(max_radial_level+1)*j],r_end_coord[n+(max_radial_level+1)*j]
            for (auto r = r_start_coord[n + (max_radial_level+1)*j]; 
                 r <= r_end_coord[n + (max_radial_level+1)*j]; ++r) {
                w0_old_cen[n+(max_radial_level+1)*r] = 0.5 * (
                    w0_old[n+(max_radial_level+1)*r] + 
                    w0_old[n+(max_radial_level+1)*(r+1)]);
                w0_new_cen[n+(max_radial_level+1)*r] = 0.5 * (
                    w0[n+(max_radial_level+1)*r] + 
                    w0[n+(max_radial_level+1)*(r+1)]);
                Real w0_avg = 0.5 * (dt * w0_old_cen[n+(max_radial_level+1)*r] 
                    + dtold *  w0_new_cen[n+(max_radial_level+1)*r]) / dt_avg;
                Real div_avg = 0.5 * (dt *(w0_old[n+(max_radial_level+1)*(r+1)]
                    - w0_old[n+(max_radial_level+1)*r]) + 
                    dtold * (w0[n+(max_radial_level+1)*(r+1)] 
                    - w0[n+(max_radial_level+1)*r])) / dt_avg;
                w0_force[n+(max_radial_level+1)*r] = (
                        w0_new_cen[n+(max_radial_level+1)*r] 
                        - w0_old_cen[n+(max_radial_level+1)*r])/dt_avg 
                        + w0_avg*div_avg/dr[n];
            }
        }
    }

    // call restrict_base(w0_force,1)
    // call fill_ghost_base(w0_force,1)
    RestrictBase(w0_force, true);
    FillGhostBase(w0_force, true);
}