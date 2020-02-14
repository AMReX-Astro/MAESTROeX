#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::MakeBeta0(RealVector& beta, 
                   const RealVector& rho0,
                   const RealVector& p0,
                   const RealVector& gamma1bar,
                   const RealVector& grav_cell,
                   const bool is_irreg) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeBeta0()",MakeBeta0);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());
    get_finest_radial_level(&finest_radial_level);

    RealVector beta0_edge((finest_radial_level+1)*(nr_fine+1));

    std::fill(beta0.begin(), beta0.end(), 0.);

    if (beta0_type == 1) {
        ///////////////////////////////////////////////////////////////////////
        // Compute beta0 on the edges and average to the center      
        //
        // Multilevel Outline:
        //
        // First, compute beta0 on edges and centers at level 0 only
        // Obtain the starting value from rho0 at the bottom of the domain.
        // do n=1,finest_radial_level
        //   Compute beta0 on edges and centers at level n
        //   Obtain the starting value of beta0_edge_lo from the coarser grid
        //   if n>0, compare the difference between beta0 at the top of level n to the
        //           corresponding point on level n-1
        //   do i=n-1,0,-1
        //     Offset the centered beta on level i above this point so the total integral 
        //      is consistent
        //     Redo the anelastic cutoff part
        //   }
        // }
        // call restrict_base and fill_ghost_base
        //////////////////////////////////////////////////////////////////////

        //    do n=0,finest_radial_level
        for (auto n = 0; n <= finest_radial_level; ++n) {
            //   do j=1,numdisjointchunks[n]
            for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
                // Compute beta0 on edges and centers at level n
                if (n == 0) {
                    beta0_edge[0] = rho0[0];
                } else {
                    // Obtain the starting value of beta0_edge_lo from the coarser grid
                    beta0_edge(n,r_start_coord[n+max_lev*j]) = beta0_edge(n-1,r_start_coord[n+max_lev*j]/2)
                }

                // do r=r_start_coord[n+max_lev*j],r_end_coord[n+max_lev*j]
                for (auto r = r_start_coord[n+max_lev*j]; r <= r_end_coord[n+max_lev*j]; ++r) {
                    Real lambda = 0.0;
                    Real mu = 0.0;
                    Real nu = 0.0;

                    if (r < anelastic_cutoff_density_coord[n]) {

                        Real drp = is_irreg ? 
                            r_edge_loc[n+max_lev*(r+1)] - r_edge_loc[n+max_lev*r] : dr[n];
                        Real drm = dr[n];
                        if (is_irreg) {
                            if (r > 0) {
                                drm = r_edge_loc[n+max_lev*r] - r_edge_loc[n+max_lev*(r-1)];
                            } else {
                                drm = drp;
                            }
                        }

                        if (r == 0 || r == nr[n]-1) {
                            // lambda = 0.0;
                            // mu = 0.0;
                            // nu = 0.0;
                        } else {
                            Real drc = is_irreg ? 
                                r_cc_loc[n+max_lev*(r+1)] - r_cc_loc[n+max_lev*(r-1)] : drc;

                            // piecewise linear reconstruction of rho0,
                            // gamma1bar, and p0 -- see paper III, appendix C
                            Real del = 0.5 * (rho0[n+max_lev*(r+1)] - rho0[n+max_lev*(r-1)])/drc;
                            Real dpls = 2.0 * (rho0[n+max_lev*(r+1)] - rho0[n+max_lev*r])/drp;
                            Real dmin = 2.0 * (rho0[n+max_lev*r] - rho0[n+max_lev*(r-1)])/drm;
                            Real slim = min(fabs(dpls), fabs(dmin));
                            slim = dpls*dmin > 0.0 ? slim : 0.0;
                            Real sflag  = copysign(1.0, del);
                            lambda = sflag*min(slim, fabs(del));
                            
                            del = 0.5* (gamma1bar[n+max_lev*(r+1)] - gamma1bar[n+max_lev*(r-1)])/drc;
                            dpls  = 2.0 * (gamma1bar[n+max_lev*(r+1)] - gamma1bar[n+max_lev*r])/drp;
                            dmin  = 2.0 * (gamma1bar[n+max_lev*r] - gamma1bar[n+max_lev*(r-1)])/drm;
                            slim  = min(fabs(dpls), fabs(dmin));
                            slim  = dpls*dmin > 0.0 ? slim : 0.0;
                            sflag = copysign(1.0, del);
                            mu = sflag*min(slim, fabs(del));
                            
                            del = 0.5* (p0[n+max_lev*(r+1)] - p0[n+max_lev*(r-1)])/drc;
                            dpls  = 2.0 * (p0[n+max_lev*(r+1)] - p0[n+max_lev*r])/drp;
                            dmin  = 2.0 * (p0[n+max_lev*r] - p0[n+max_lev*(r-1)])/drm;
                            slim  = min(fabs(dpls), fabs(dmin));
                            slim  = dpls*dmin > 0.0 ? slim : 0.0;
                            sflag = copysign(1.0, del);
                            nu = sflag*min(slim, fabs(del));
                        }

                        if (is_irreg) {
                            // edge-to-cell-center spacings 
                            drp = 2.0 * (r_edge_loc[n+max_lev*(r+1)] - r_cc_loc[n+max_lev*r]);
                            drm = 2.0 * (r_cc_loc[n+max_lev*r] - r_edge_loc[n+max_lev*r]);
                        }

                        Real integral = 0.0;

                        if (nu == 0.0 || mu == 0.0 ||
                            (nu*gamma1bar[n+max_lev*r] - mu*p0[n+max_lev*r]) == 0.0 ||
                            ((gamma1bar[n+max_lev*r] + 0.5*mu*drp)/
                            (gamma1bar[n+max_lev*r] - 0.5*mu*drm)) <= 0.0 ||
                            ((p0[n+max_lev*r] + 0.5*nu*drp)/
                            (p0[n+max_lev*r] - 0.5*nu*drm)) <= 0.0) {
                            
                            // just do piecewise constant integration
                            integral = fabs(grav_cell[n+max_lev*r])*rho0[n+max_lev*r]*0.5*(drp+drm)
                                / (p0[n+max_lev*r]*gamma1bar[n+max_lev*r]);
                            
                        } else {
                            if (use_linear_grav_in_beta0 && !is_irreg) {
                                // also do piecewise linear reconstruction of
                                // gravity -- not documented in publication yet.
                                Real del = 0.5* (grav_cell[n+max_lev*(r+1)] - grav_cell[n+max_lev*(r-1)])/dr[n];
                                Real dpls  = 2.0 * (grav_cell[n+max_lev*(r+1)] - grav_cell[n+max_lev*r])/dr[n];
                                Real dmin  = 2.0 * (grav_cell[n+max_lev*r] - grav_cell[n+max_lev*(r-1)])/dr[n];
                                Real slim  = min(fabs(dpls), fabs(dmin));
                                slim  = dpls*dmin > 0.0 ? slim : 0.0;
                                Real sflag = copysign(1.0, del);
                                Real kappa = sflag*min(slim, fabs(del));
                                
                                Real denom = nu*gamma1bar[n+max_lev*r] - mu*p0[n+max_lev*r];
                                Real coeff1 = (lambda*gamma1bar[n+max_lev*r] - mu*rho0[n+max_lev*r]) *
                                    (kappa *gamma1bar[n+max_lev*r] + mu*fabs(grav_cell[n+max_lev*r])) /
                                    (mu*mu*denom);
                                Real coeff2 = (lambda*p0[n+max_lev*r] - nu*rho0[n+max_lev*r])*
                                    (-kappa*p0[n+max_lev*r] - nu*fabs(grav_cell[n+max_lev*r])) /
                                    (nu*nu*denom);
                                Real coeff3 = kappa*lambda / (mu*nu);
                                
                                integral = 
                                    coeff1*log( (gamma1bar[n+max_lev*r] + 0.5*mu*dr[n])/
                                                (gamma1bar[n+max_lev*r] - 0.5*mu*dr[n]) ) +
                                    coeff2*log( (p0[n+max_lev*r] + 0.5*nu*dr[n])/
                                                (p0[n+max_lev*r] - 0.5*nu*dr[n]) ) -
                                    coeff3*dr[n];

                            } else {
                                // paper III, equation C2
                                Real denom = nu*gamma1bar[n+max_lev*r] - mu*p0[n+max_lev*r];
                                Real coeff1 = lambda*gamma1bar[n+max_lev*r]/mu - rho0[n+max_lev*r];
                                Real coeff2 = lambda*p0[n+max_lev*r]/nu - rho0[n+max_lev*r];

                                integral = (fabs(grav_cell[n+max_lev*r])/denom)*
                                    (coeff1*log( (gamma1bar[n+max_lev*r] + 0.5*mu*drp)/
                                                (gamma1bar[n+max_lev*r] - 0.5*mu*drm)) -
                                    coeff2*log( (p0[n+max_lev*r] + 0.5*nu*drp)/
                                                (p0[n+max_lev*r] - 0.5*nu*drm)) );
                            }
                        }

                        beta0_edge[n+max_lev*(r+1)] = beta0_edge[n+max_lev*r] * exp(-integral);
                        beta0[n+max_lev*r] = 0.5*(beta0_edge[n+max_lev*r] + 
                            beta0_edge[n+max_lev*(r+1)]);

                    } else {// r >= anelastic_cutoff_density

                        beta0[n+max_lev*r] = beta0[n+max_lev*(r-1)] * (rho0[n+max_lev*r]/rho0[n+max_lev*(r-1)]);
                        beta0_edge[n+max_lev*(r+1)] = 2.0*beta0[n+max_lev*r] - beta0_edge[n+max_lev*r];
                    }
                }

                if (n  >  0) {
                    // Compare the difference between beta0 at the top of level n to the 
                    // corresponding point on level n-1
                    Real offset = beta0_edge[n+max_lev*(r_end_coord[n+max_lev*j]+1)]
                        - beta0_edge[n-1+max_lev*(r_end_coord[n+max_lev*j]+1)/2];

                    // do i=n-1,0,-1
                    for (auto i = n-1; i >= 0; --i) {

                        int refrat = pow(2, n-i);

                        // Offset the centered beta on level i above this point so the total 
                        // integral is consistent
                        // do r=r_end_coord[n+max_lev*j]/refrat+1,nr(i)
                        for (auto r = r_end_coord[n+max_lev*j]/refrat+1; r <= nr[i]; ++r) {
                            beta0[i+max_lev*r] += offset;
                        }

                        // Redo the anelastic cutoff part
                        // do r=anelastic_cutoff_density_coord(i),nr(i)
                        for (auto r = anelastic_cutoff_density_coord[i]; r <= nr[i]; ++r) {
                            if (rho0[i+max_lev*(r-1)] /= 0.0) {
                                beta0[i+max_lev*r] = beta0[i+max_lev*(r-1)] * 
                                    (rho0[i+max_lev*r]/rho0[i+max_lev*(r-1)]);
                            }
                        }

                        // This next piece of coded is needed for the case when the anelastic 
                        // cutoff coordinate lives on level n.  We first average beta0 from 
                        // level i+1 to level i in the region between the anelastic cutoff and 
                        // the top of grid n.  Then recompute beta0 at level i above the top 
                        // of grid n.
                        if (r_end_coord[n+max_lev*j] >= anelastic_cutoff_density_coord[n]) {

                            // do r=anelastic_cutoff_density_coord(i),(r_end_coord[n+max_lev*j]+1)/refrat-1
                            for (auto r = anelastic_cutoff_density_coord[i]; 
                                 r <= (r_end_coord[n+max_lev*j]+1)/refrat-1; ++r) {
                                beta0[i+max_lev*r] = 0.5*(beta0[i+1+max_lev*2*r] + 
                                    beta0[i+1+max_lev*(2*r+1)]);
                            }

                            // do r=(r_end_coord[n+max_lev*j]+1)/refrat,nr(i)
                            for (auto r = (r_end_coord[n+max_lev*j]+1)/refrat; 
                                 r <= nr[i]; ++r) {
                                if (rho0[i+max_lev*(r-1)] /= 0.0) {
                                    beta0[i+max_lev*r] = beta0[i+max_lev*(r-1)] * 
                                        (rho0[i+max_lev*r]/rho0[i+max_lev*(r-1)]);
                                }
                            }
                        }
                    } // end loop over i=n-1,0,-1
                } // } (n  >  0)
            } // end loop over disjoint chunks
        } // end loop over levels

        // 0.0 the beta0 where there is no corresponding full state array
        //    do n=1,finest_radial_level
        //       do j=1,numdisjointchunks[n]
        for (auto n = 1; n <= finest_radial_level; ++n) {
            for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
                if (j == numdisjointchunks[n]) {
                    // do r=r_end_coord[n+max_lev*j]+1,nr[n]-1
                    for (auto r = r_end_coord[n+max_lev*j]+1; r < nr[n]; ++r) {
                        beta0[n+max_lev*r] = 0.0;
                    }
                } else {
                    // do r=r_end_coord[n+max_lev*j]+1,r_start_coord(n,j+1)-1
                    for (auto r = r_end_coord[n+max_lev*j]+1; r < r_start_coord[n+max_lev*(j+1)]; ++r) {
                        beta0[n+max_lev*r] = 0.0;
                    }
                }
            }
        }
    } else if (beta0_type == 2) {
        // beta_0 = rho_0
        //    do n=0,finest_radial_level
        //       do j=1,numdisjointchunks[n]
        for (auto n = 0; n <= finest_radial_level; ++n) {
            for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
                // do r=r_start_coord[n+max_lev*j],r_end_coord[n+max_lev*j]
                for (auto r = r_start_coord[n+max_lev*j]; r <= r_end_coord[n+max_lev*j]; ++r) {
                    beta0[n+max_lev*r] = rho0[n+max_lev*r];
                }
            }
        }
    } else if (beta0_type == 3) {
        // beta_0 = 1.0
        //    do n=0,finest_radial_level
        //       do j=1,numdisjointchunks[n]
        //          do r=r_start_coord[n+max_lev*j],r_end_coord[n+max_lev*j]
        for (auto n = 0; n <= finest_radial_level; ++n) {
            for (auto j = 1; j <= numdisjointchunks[n]; ++j) {
                for (auto r = r_start_coord[n+max_lev*j]; r <= r_end_coord[n+max_lev*j]; ++r) {
                    beta0[n+max_lev*r] = 1.0;
                }
            }
        }
    }

    RestrictBase(beta0, true);
    FillGhostBase(beta0, true);
}
