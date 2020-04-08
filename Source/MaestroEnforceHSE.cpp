#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::EnforceHSE(const RealVector& rho0, 
                    RealVector& p0,
                    const RealVector& grav_cell) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EnforceHSE()", EnforceHSE);

    const int max_lev = max_radial_level+1;

    RealVector grav_edge((finest_radial_level+1)*(nr_fine+1));
    RealVector p0old((finest_radial_level+1)*(nr_fine+1));

    Real offset = 0.0;

    MakeGravEdge(grav_edge, rho0);

    // create a copy of the input pressure to help us with initial
    // conditions
    for (auto i = 0; i < p0.size(); ++i) {
        p0old[i] = p0[i];
    }

    // zero the new pressure so we don't leave a non-zero pressure in
    // fine radial regions that no longer have a corresponding full
    // state
    std::fill(p0.begin(), p0.end(), 0.);

    // integrate all of level 1 first
    // use the old pressure at r=0 as a reference point
    p0[0] = p0old[0];

    // now integrate upwards from the bottom later, we will offset the
    // entire pressure so we have effectively integrated from the "top"
    if (use_exact_base_state && spherical) {
        for (auto r = 1; r <= min(r_end_coord(0,1),base_cutoff_density_coord(0)); ++r) {
            // uneven grid spacing
            Real dr1 = r_edge_loc_b(0,r)-r_cc_loc_b(0,r-1);
            Real dr2 = r_cc_loc_b(0,r)-r_edge_loc_b(0,r);
            p0[max_lev*r] = p0[max_lev*(r-1)] + (dr1*rho0[max_lev*(r-1)] + 
                dr2*rho0[max_lev*r])*grav_edge[max_lev*r];
        }
    } else {
        for (auto r = 1; r <= min(r_end_coord(0,1),base_cutoff_density_coord(0)); r++) {
            // assume even grid spacing
            p0[max_lev*r] = p0[max_lev*(r-1)] + 0.5*dr(0)*
                (rho0[max_lev*(r-1)] + rho0[max_lev*r])*grav_edge[max_lev*r];
        }
    }
    for (auto r = base_cutoff_density_coord(0)+1; r <= r_end_coord(0,1); ++r) {
        p0[max_lev*r] = p0[max_lev*(r-1)];
    }

    if (!spherical) {

        for (auto n = 1; n <= finest_radial_level; ++n) {
            for (auto i = 1; i <= numdisjointchunks(n); ++i) {

                // get pressure in the bottom cell of this disjointchunk
                if (r_start_coord(n,i) == 0) {
                    // if we are at the bottom of the domain, use the old
                    // pressure as reference
                    p0[n] = p0old[n];
                } else if (r_start_coord(n,i) <= base_cutoff_density_coord(n)) {
                    // we integrate upwards starting from the nearest coarse
                    // cell at a lower physical height

                    if (do_planar_invsq_grav || do_2d_planar_octant) {
                        // we have variable gravity
                        p0[n+max_lev*r_start_coord(n,i)] = p0[n-1+max_lev*(r_start_coord(n,i)/2-1)]
                                + (dr(n)/4.0)* 
                                (2.0*rho0[n+max_lev*r_start_coord(n,i)]/3.0 + 
                                4.0*rho0[n-1+max_lev*(r_start_coord(n,i)/2-1)]/3.0)* 
                                (grav_edge[n+max_lev*r_start_coord(n,i)] + 
                                grav_cell[n-1+max_lev*(r_start_coord(n,i)/2-1)])  
                                + (dr(n)/8.0)* 
                                (5.0*rho0[n+max_lev*r_start_coord(n,i)]/3.0 + 
                                1.0*rho0[n-1+max_lev*(r_start_coord(n,i)/2-1)]/3.0)* 
                                (grav_edge[n+max_lev*r_start_coord(n,i)] + 
                                grav_cell[n+max_lev*r_start_coord(n,i)]);
                    } else {
                        // assuming constant g here
                        p0[n+max_lev*r_start_coord(n,i)] = p0[n-1+max_lev*(r_start_coord(n,i)/2-1)] 
                                + (3.0*grav_cell[1]*dr(n)/4.0)* 
                                (rho0[n-1+max_lev*(r_start_coord(n,i)/2-1)]+rho0[n+max_lev*r_start_coord(n,i)]);
                    }
                } else {
                    // copy pressure from below
                    p0[n+max_lev*r_start_coord(n,i)] = p0[n-1+max_lev*(r_start_coord(n,i)/2-1)];
                }

                // integrate upwards as normal
                for (auto r = r_start_coord(n,i)+1; 
                     r <= min(r_end_coord(n,i),base_cutoff_density_coord(n)); ++r) {
                    p0[n+max_lev*r] = p0[n+max_lev*(r-1)] + 0.5 * dr(n) * 
                        (rho0[n+max_lev*r]+rho0[n+max_lev*(r-1)])*grav_edge[n+max_lev*r];
                }
                for (auto r = base_cutoff_density_coord(n)+1; 
                     r <= r_end_coord(n,i); ++r) {
                    p0[n+max_lev*r] = p0[n+max_lev*(r-1)];
                }

                // now we need to look at the first coarser cell above this
                // disjoint chunk and look at the mismatch between the
                // integration performed over the finer cells vs. the coarse
                // cells.  Then we need to offset the coarse cell above this
                // point to sync up.

                // first, compute the value of the pressure in the coarse
                // cell above the disjointchunk.
                if (r_end_coord(n,i) == nr(n)-1) {
                    // for (auto nothing - we are at the top of the domain
                    offset = 0.0;
                } else if (r_end_coord(n,i) <= base_cutoff_density_coord(n)) {
                    // use fine -> coarse stencil in notes
                    if (do_planar_invsq_grav || do_2d_planar_octant) {
                        // we have variable gravity
                        Real temp = p0[n+max_lev*r_end_coord(n,i)] 
                                + (dr(n)/4.0)* 
                                (2.0*rho0[n+max_lev*r_end_coord(n,i)]/3.0 + 
                                4.0*rho0[n-1+max_lev*(r_end_coord(n,i)+1)/2]/3.0)* 
                                (grav_edge[n-1+max_lev*(r_end_coord(n,i)+1)/2] + 
                                grav_cell[n-1+max_lev*(r_end_coord(n,i)+1)/2]) 
                                + (dr(n)/8.0)* 
                                (5.0*rho0[n+max_lev*r_end_coord(n,i)]/3.0 + 
                                1.0*rho0[n-1+max_lev*(r_end_coord(n,i)+1)/2]/3.0)* 
                                (grav_cell[n+max_lev*r_end_coord(n,i)] + 
                                grav_edge[n-1+max_lev*(r_end_coord(n,i)+1)/2]);
                        offset = p0[n-1+max_lev*(r_end_coord(n,i)+1)/2] - temp;
                    } else {
                        // assuming constant g here
                        Real temp = p0[n+max_lev*r_end_coord(n,i)] + (3.0*grav_cell[1]*dr(n)/4.0)* 
                                (rho0[n+max_lev*r_end_coord(n,i)]+rho0[n-1+max_lev*(r_end_coord(n,i)+1)/2]);
                        offset = p0[n-1+max_lev*(r_end_coord(n,i)+1)/2] - temp;
                    }
                } else {
                    // copy pressure from below
                    Real temp = p0[n+max_lev*r_end_coord(n,i)];
                    offset = p0[n-1+max_lev*(r_end_coord(n,i)+1)/2] - temp;
                }

                // if we are not at the top of the domain, we need to
                // subtract the offset for all values at and above this point
                if (r_end_coord(n,i) != nr(n)-1) {
                    for (auto l = n-1; l >= 0; --l) {
                        for (int r = round( (r_end_coord[n+max_lev*i]+1)/pow(2, n-l) ); 
                             r <= nr[l]-1; ++r) {
                            p0[l+max_lev*r] -= offset;
                        }
                    }
                }
            } // end loop over disjoint chunks
        } // end loop over levels
    } // spherical

    // now compare pressure in the last cell and offset to make sure we
    // are integrating "from the top"
    // we use the coarsest level as the reference point
    offset = p0[max_lev*(nr(0)-1)] - p0old[max_lev*(nr(0)-1)];

    // offset level 0
    for (auto r = 0; r < nr_fine; ++r) {
        p0[max_lev*r] -= offset;
    }

    // offset remaining levels
    for (auto n = 1; n <= finest_radial_level; ++n) {
        for (auto i = 1; i <= numdisjointchunks(n); ++i) {
            for (auto r = r_start_coord(n,i); r <= r_end_coord(n,i); ++r) {
                p0[n+max_lev*r] -= offset;
            }
        }
    }

    // zero p0 where there is no corresponding full state array
    for (auto n = 1; n <= finest_radial_level; ++n) {
        for (auto i = 1; i <= numdisjointchunks(n); ++i) {
            if (i == numdisjointchunks(n)) {
                for (auto r = r_end_coord(n,i)+1; r < nr(n); ++r) {
                    p0[n+max_lev*r] = 0.0;
                }
            } else {
                for (auto r = r_end_coord(n,i)+1; r < r_start_coord(n,i+1); ++r) {
                    p0[n+max_lev*r] = 0.0;
                }
            }
        }
    }

    RestrictBase(p0, true);
    FillGhostBase(p0, true);
}
