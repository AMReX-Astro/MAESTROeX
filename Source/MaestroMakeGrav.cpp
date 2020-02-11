#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::MakeGravCell(RealVector& grav_cell, 
                      const RealVector& rho0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravCell()",MakeGravCell);

    const int max_lev = max_radial_level+1;

    if (!spherical) {
        if (do_planar_invsq_grav)  {
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            //   do n=0,finest_radial_level
            for (auto n = 0; n <= finest_radial_level; ++n) {
                for (auto r = 0; r < nr[n]; ++r) {
                    // do r = 0, nr(n)-1
                    grav_cell[n+max_lev*r] = -Gconst*planar_invsq_mass / (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r]);
                }
            }
        } else if (do_2d_planar_octant) {
            //   compute gravity as in the spherical case
            RealVector m((finest_radial_level+1)*nr_fine);
            //   allocate(m(0:finest_radial_level,0:nr_fine-1))

            const int n = 0;
            m[n] = 4.0/3.0*M_PI*rho0[n]*r_cc_loc[n]*r_cc_loc[n]*r_cc_loc[n];
            grav_cell[n] = -Gconst * m[n] / (r_cc_loc[n]*r_cc_loc[n]);

            // do r=1,nr(n)-1
            for (auto r = 1; r < nr[n]; ++r) {

                // the mass is defined at the cell-centers, so to compute
                // the mass at the current center, we need to add the
                // contribution of the upper half of the zone below us and
                // the lower half of the current zone.

                // don't add any contributions from outside the star --
                // i.e.  rho < base_cutoff_density
                Real term1 = 0.0;
                if (rho0[n+max_lev*(r-1)] > base_cutoff_density) {
                    term1 = 4.0/3.0*M_PI*rho0[n+max_lev*(r-1)] *
                        (r_edge_loc[n+max_lev*r] - r_cc_loc[n+max_lev*(r-1)]) *
                        (r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r] +
                        r_edge_loc[n+max_lev*r]*r_cc_loc[n+max_lev*(r-1)] +
                        r_cc_loc[n+max_lev*(r-1)]*r_cc_loc[n+max_lev*(r-1)]);
                } 

                Real term2 = 0.0;
                if (rho0[n+max_lev*r] > base_cutoff_density) {
                    term2 = 4.0/3.0*M_PI*rho0[n+max_lev*r]*
                        (r_cc_loc[n+max_lev*r] - r_edge_loc[n+max_lev*r]) *
                        (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r] +
                        r_cc_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r] +
                        r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r]);
                } 

                m[n+max_lev*r] = m[n+max_lev*(r-1)] + term1 + term2;

                grav_cell[n+max_lev*r] = -Gconst * m[n+max_lev*r] / (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r]);
            }

            for (auto n = 1; n <= finest_radial_level; ++n) {
            //   do n=1,finest_radial_level
                for (auto i = 1; i <= numdisjointchunks[n]; ++i) {
                //  do i=1,numdisjointchunks(n)

                    if (r_start_coord[n+max_lev*i] == 0) {
                        m[n] = 4.0/3.0*M_PI*rho0[n]*r_cc_loc[n]*r_cc_loc[n]*r_cc_loc[n];
                        grav_cell[n] = -Gconst * m[n] / (r_cc_loc[n]*r_cc_loc[n]);
                    } else {
                        int r = r_start_coord[n+max_lev*i];
                        m[n+max_lev*r] = m[n-1+max_lev*(r/2-1)];

                        // the mass is defined at the cell-centers, so to compute
                        // the mass at the current center, we need to add the
                        // contribution of the upper half of the zone below us and
                        // the lower half of the current zone.

                        // don't add any contributions from outside the star --
                        // i.e.  rho < base_cutoff_density
                        Real term1 = 0.0;
                        if (rho0[n-1+max_lev*(r/2-1)] > base_cutoff_density) {
                            term1 = 4.0/3.0*M_PI*rho0[n-1+max_lev*(r/2-1)] *
                                (r_edge_loc[n-1+max_lev*r/2] - r_cc_loc[n-1+max_lev*(r/2-1)]) *
                                (r_edge_loc[n-1+max_lev*r/2]*r_edge_loc[n-1+max_lev*r/2] +
                                r_edge_loc[n-1+max_lev*r/2]*r_cc_loc[n-1+max_lev*(r/2-1)] +
                                r_cc_loc[n-1+max_lev*(r/2-1)]*r_cc_loc[n-1+max_lev*(r/2-1)]);
                        } 

                        Real term2 = 0.0;
                        if (rho0[n+max_lev*r] > base_cutoff_density) {
                            term2 = 4.0/3.0*M_PI*rho0[n+max_lev*r]*
                                (r_cc_loc[n+max_lev*r] - r_edge_loc[n+max_lev*r]) *
                                (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r] +
                                r_cc_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r] +
                                r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r]);
                        } 

                        m[n+max_lev*r] += term1 + term2;
                        grav_cell[n+max_lev*r] = -Gconst * m[n+max_lev*r] / (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r]);
                    }

                    // do r=r_start_coord[n+max_lev*i]+1,r_end_coord[n+max_lev*i]
                    for (auto r = r_start_coord[n+max_lev*i]+1;
                         r <= r_end_coord[n+max_lev*i]; ++r) {

                        // the mass is defined at the cell-centers, so to compute
                        // the mass at the current center, we need to add the
                        // contribution of the upper half of the zone below us and
                        // the lower half of the current zone.

                        // don't add any contributions from outside the star --
                        // i.e.  rho < base_cutoff_density
                        Real term1 = 0.0;
                        if (rho0[n+max_lev*(r-1)] > base_cutoff_density) {
                            term1 = 4.0/3.0*M_PI*rho0[n+max_lev*(r-1)] *
                                (r_edge_loc[n+max_lev*r] - r_cc_loc[n+max_lev*(r-1)]) *
                                (r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r] +
                                r_edge_loc[n+max_lev*r]*r_cc_loc[n+max_lev*(r-1)] +
                                r_cc_loc[n+max_lev*(r-1)]*r_cc_loc[n+max_lev*(r-1)]);
                        } 

                        Real term2 = 0.0;
                        if (rho0[n+max_lev*r] > base_cutoff_density) {
                            term2 = 4.0/3.0*M_PI*rho0[n+max_lev*r]*
                                (r_cc_loc[n+max_lev*r] - r_edge_loc[n+max_lev*r]) *
                                (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r] +
                                r_cc_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r] +
                                r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r]);
                        } 

                        m[n+max_lev*r] = m[n+max_lev*(r-1)] + term1 + term2;

                        grav_cell[n+max_lev*r] = -Gconst * m[n+max_lev*r] / (r_cc_loc[n+max_lev*r]*r_cc_loc[n+max_lev*r]);
                    }
                }
            }

            //   call restrict_base(grav_cell,1)
            //   call fill_ghost_base(grav_cell,1)
            RestrictBase(grav_cell, true);
            FillGhostBase(grav_cell, true);
        } else {
            // constant gravity
            std::fill(grav_cell.begin(), grav_cell.end(), grav_const);
        }
    } else { // spherical = 1

        //    allocate(m(0:0,0:nr_fine-1))
        RealVector m(nr_fine);

        m[0] = 4.0/3.0*M_PI*rho0[0]*r_cc_loc[0]*r_cc_loc[0]*r_cc_loc[0];
        grav_cell[0] = -Gconst * m[0] / (r_cc_loc[0]*r_cc_loc[0]);

        // do r=1,nr_fine-1
        for (auto r = 1; r < nr_fine; ++r) {

            // the mass is defined at the cell-centers, so to compute
            // the mass at the current center, we need to add the
            // contribution of the upper half of the zone below us and
            // the lower half of the current zone.

            // don't add any contributions from outside the star --
            // i.e.  rho < base_cutoff_density
            Real term1 = 0.0;
            if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                term1 = 4.0/3.0*M_PI*rho0[max_lev*(r-1)] *
                    (r_edge_loc[max_lev*r] - r_cc_loc[max_lev*(r-1)]) *
                    (r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r] +
                    r_edge_loc[max_lev*r]*r_cc_loc[max_lev*(r-1)] +
                    r_cc_loc[max_lev*(r-1)]*r_cc_loc[max_lev*(r-1)]);
            } 

            Real term2 = 0.0;
            if (rho0[max_lev*r] > base_cutoff_density) {
                term2 = 4.0/3.0*M_PI*rho0[max_lev*r]*
                    (r_cc_loc[max_lev*r] - r_edge_loc[max_lev*r]) *
                    (r_cc_loc[max_lev*r]*r_cc_loc[max_lev*r] +
                    r_cc_loc[max_lev*r]*r_edge_loc[max_lev*r] +
                    r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r]);
            } 

            m[max_lev*r] = m[max_lev*(r-1)] + term1 + term2;

            grav_cell[max_lev*r] = -Gconst * m[max_lev*r] / (r_cc_loc[max_lev*r]*r_cc_loc[max_lev*r]);
        }
    }
}

void
Maestro::MakeGravEdge(RealVector& grav_edge, 
                      const RealVector& rho0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravEdge()",MakeGravEdge);

    RealVector m;
    const int max_lev = max_radial_level+1;

    if (!spherical) {
        if (do_planar_invsq_grav)  {
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            //   
            for (auto n = 0; n <= finest_radial_level; ++n) {
                for (auto r = 0; r < nr[n]; ++r) {
                // do r = 0, nr(n)-1
                    grav_edge[n+max_lev*r] = -Gconst*planar_invsq_mass / (r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r]);
                }
            }
        } else if (do_2d_planar_octant) {
            // compute gravity as in spherical geometry

            // allocate(m(0:finest_radial_level,0:nr_fine))
            m.resize((finest_radial_level+1)*(nr_fine+1));
            m.shrink_to_fit();

            grav_edge[0] = 0.0;
            m[0] = 0.0;

            //   do r=1,nr(0)-1
            for (auto r = 0; r < nr[0]; ++r) {

                // only add to the enclosed mass if the density is
                // > base_cutoff_density
                if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                    m[max_lev*r] = m[max_lev*(r-1)] + 4.0/3.0*M_PI *
                        (r_edge_loc[max_lev*r] - r_edge_loc[max_lev*(r-1)]) *
                        (r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r] +
                        r_edge_loc[max_lev*r]*r_edge_loc[max_lev*(r-1)] +
                        r_edge_loc[max_lev*(r-1)]*r_edge_loc[max_lev*(r-1)]) * rho0[max_lev*(r-1)];
                } else {
                    m[max_lev*r] = m[max_lev*(r-1)];
                }

                grav_edge[max_lev*r] = -Gconst * m[max_lev*r] / (r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r]);
            }

            for (auto n = 1; n <= finest_radial_level; ++n) {
            //   do n=1,finest_radial_level
                for (auto i = 1; i <= numdisjointchunks[n]; ++i) {
                //  do i=1,numdisjointchunks(n)

                    if (r_start_coord[n+max_lev*i] == 0) {
                        m[n] = 0.0;
                    } else {
                        m[n+max_lev*r_start_coord[n+max_lev*i]] = m[n-1*max_lev*r_start_coord[n+max_lev*i]/2];
                        grav_edge[n+max_lev*r_start_coord[n+max_lev*i]] = grav_edge[n-1+max_lev*r_start_coord[n+max_lev*i]/2];
                    }

                    // do r=r_start_coord[n+max_lev*i]+1,r_end_coord[n+max_lev*i]+1
                    for (auto r = r_start_coord[n+max_lev*i]+1; 
                         r <= r_end_coord[n+max_lev*i]+1; ++r) {

                        // only add to the enclosed mass if the density is
                        // > base_cutoff_density
                        if (rho0[n+max_lev*(r-1)] > base_cutoff_density) {
                            m[n+max_lev*r] = m[n+max_lev*(r-1)] + 4.0/3.0*M_PI *
                                (r_edge_loc[n+max_lev*r] - r_edge_loc[n+max_lev*(r-1)]) *
                                (r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r] +
                                r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*(r-1)] +
                                r_edge_loc[n+max_lev*(r-1)]*r_edge_loc[n+max_lev*(r-1)]) * rho0[n+max_lev*(r-1)];
                        } else {
                            m[n+max_lev*r] = m[n+max_lev*(r-1)];
                        }

                        grav_edge[n+max_lev*r] = -Gconst * m[n+max_lev*r] / (r_edge_loc[n+max_lev*r]*r_edge_loc[n+max_lev*r]);
                    }
                }
            }
            //   call restrict_base(grav_edge,0)
            //   call fill_ghost_base(grav_edge,0)
            RestrictBase(grav_edge, false);
            FillGhostBase(grav_edge, false);
        } else {
            // constant gravity
            std::fill(grav_edge.begin(), grav_edge.end(), grav_const);
        }
        
    } else {

       grav_edge[0] = 0.0;
       Real mencl = 0.0;

        for (auto r = 1; r <= nr_fine; ++r) {
        //    do r=1,nr_fine

            // only add to the enclosed mass if the density is
            // > base_cutoff_density
            if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                mencl += 4.0/3.0 * M_PI *
                    (r_edge_loc[max_lev*r] - r_edge_loc[max_lev*(r-1)]) *
                    (r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r] +
                    r_edge_loc[max_lev*r]*r_edge_loc[max_lev*(r-1)] +
                    r_edge_loc[max_lev*(r-1)]*r_edge_loc[max_lev*(r-1)]) * rho0[max_lev*(r-1)];
            }

            grav_edge[max_lev*r] = -Gconst * mencl / (r_edge_loc[max_lev*r]*r_edge_loc[max_lev*r]);
        }
    }
}