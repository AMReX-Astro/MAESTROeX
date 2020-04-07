#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::MakeGravCell(BaseState<Real>& grav_cell, 
                      const BaseState<Real>& rho0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravCell()", MakeGravCell);

    const int max_lev = max_radial_level+1;

    if (!spherical) {
        if (do_planar_invsq_grav)  {
            const auto r_cc_loc_p = r_cc_loc_b;
            const Real planar_invsq_mass_loc = planar_invsq_mass;
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            for (auto n = 0; n <= finest_radial_level; ++n) {
                // for (auto r = 0; r < nr(n); ++r) {
                const int nr_lev = nr(n);
                AMREX_PARALLEL_FOR_1D(nr_lev, r, {
                    grav_cell(n,r) = -Gconst*planar_invsq_mass_loc / (r_cc_loc_p(n,r)*r_cc_loc_p(n,r));
                });
            }
        } else if (do_2d_planar_octant) {
            //   compute gravity as in the spherical case
            BaseState<Real> m(finest_radial_level+1, nr_fine);

            // level = 0
            m(0,0) = 4.0/3.0*M_PI*rho0(0,0)*r_cc_loc_b(0,0)*r_cc_loc_b(0,0)*r_cc_loc_b(0,0);
            grav_cell(0,0) = -Gconst * m(0,0) / (r_cc_loc_b(0,0)*r_cc_loc_b(0,0));

            int nr_lev = nr(0);

            for (auto r = 1; r < nr_lev; ++r) {

                // the mass is defined at the cell-centers, so to compute
                // the mass at the current center, we need to add the
                // contribution of the upper half of the zone below us and
                // the lower half of the current zone.

                // don't add any contributions from outside the star --
                // i.e.  rho < base_cutoff_density
                Real term1 = 0.0;
                if (rho0(0,r-1) > base_cutoff_density) {
                    term1 = 4.0/3.0*M_PI*rho0(0,r-1) *
                        (r_edge_loc_b(0,r) - r_cc_loc_b(0,r-1)) *
                        (r_edge_loc_b(0,r)*r_edge_loc_b(0,r) +
                        r_edge_loc_b(0,r)*r_cc_loc_b(0,r-1) +
                        r_cc_loc_b(0,r-1)*r_cc_loc_b(0,r-1));
                } 

                Real term2 = 0.0;
                if (rho0(0,r) > base_cutoff_density) {
                    term2 = 4.0/3.0*M_PI*rho0(0,r)*
                        (r_cc_loc_b(0,r) - r_edge_loc_b(0,r)) *
                        (r_cc_loc_b(0,r)*r_cc_loc_b(0,r) +
                        r_cc_loc_b(0,r)*r_edge_loc_b(0,r) +
                        r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
                } 

                m(0,r) = m(0,r-1) + term1 + term2;

                grav_cell(0,r) = -Gconst * m(0,r) / (r_cc_loc_b(0,r)*r_cc_loc_b(0,r));
            }

            // level > 0

            for (auto n = 1; n <= finest_radial_level; ++n) {
                for (auto i = 1; i <= numdisjointchunks(n); ++i) {

                    if (r_start_coord(n,i) == 0) {
                        m(n,0) = 4.0/3.0*M_PI*rho0(n,0)*r_cc_loc_b(n,0)*r_cc_loc_b(n,0)*r_cc_loc_b(n,0);
                        grav_cell(n,0) = -Gconst * m(n,0) / (r_cc_loc_b(n,0)*r_cc_loc_b(n,0));
                    } else {
                        int r = r_start_coord(n,i);
                        m(n,r) = m(n-1,r/2-1);

                        // the mass is defined at the cell-centers, so to compute
                        // the mass at the current center, we need to add the
                        // contribution of the upper half of the zone below us and
                        // the lower half of the current zone.

                        // don't add any contributions from outside the star --
                        // i.e.  rho < base_cutoff_density
                        Real term1 = 0.0;
                        if (rho0(n-1,r/2-1) > base_cutoff_density) {
                            term1 = 4.0/3.0*M_PI*rho0(n-1,r/2-1) *
                                (r_edge_loc_b(n-1,r/2) - r_cc_loc_b(n-1,r/2-1)) *
                                (r_edge_loc_b(n-1,r/2)*r_edge_loc_b(n-1,r/2) +
                                r_edge_loc_b(n-1,r/2)*r_cc_loc_b(n-1,r/2-1) +
                                r_cc_loc_b(n-1,r/2-1)*r_cc_loc_b(n-1,r/2-1));
                        } 

                        Real term2 = 0.0;
                        if (rho0(n,r) > base_cutoff_density) {
                            term2 = 4.0/3.0*M_PI*rho0(n,r)*
                                (r_cc_loc_b(n,r) - r_edge_loc_b(n,r)) *
                                (r_cc_loc_b(n,r)*r_cc_loc_b(n,r) +
                                r_cc_loc_b(n,r)*r_edge_loc_b(n,r) +
                                r_edge_loc_b(n,r)*r_edge_loc_b(n,r));
                        } 

                        m(n,r) += term1 + term2;
                        grav_cell(n,r) = -Gconst * m(n,r) / (r_cc_loc_b(n,r)*r_cc_loc_b(n,r));
                    }

                    for (auto r = r_start_coord(n,i)+1;
                         r <= r_end_coord(n,i); ++r) {

                        // the mass is defined at the cell-centers, so to compute
                        // the mass at the current center, we need to add the
                        // contribution of the upper half of the zone below us and
                        // the lower half of the current zone.

                        // don't add any contributions from outside the star --
                        // i.e.  rho < base_cutoff_density
                        Real term1 = 0.0;
                        if (rho0(n,r-1) > base_cutoff_density) {
                            term1 = 4.0/3.0*M_PI*rho0(n,r-1) *
                                (r_edge_loc_b(n,r) - r_cc_loc_b(n,r-1)) *
                                (r_edge_loc_b(n,r)*r_edge_loc_b(n,r) +
                                r_edge_loc_b(n,r)*r_cc_loc_b(n,r-1) +
                                r_cc_loc_b(n,r-1)*r_cc_loc_b(n,r-1));
                        } 

                        Real term2 = 0.0;
                        if (rho0(n,r) > base_cutoff_density) {
                            term2 = 4.0/3.0*M_PI*rho0(n,r)*
                                (r_cc_loc_b(n,r) - r_edge_loc_b(n,r)) *
                                (r_cc_loc_b(n,r)*r_cc_loc_b(n,r) +
                                r_cc_loc_b(n,r)*r_edge_loc_b(n,r) +
                                r_edge_loc_b(n,r)*r_edge_loc_b(n,r));
                        } 

                        m(n,r) = m(n,r-1) + term1 + term2;

                        grav_cell(n,r) = -Gconst * m(n,r) / (r_cc_loc_b(n,r)*r_cc_loc_b(n,r));
                    }
                }
            }

            RestrictBase(grav_cell, true);
            FillGhostBase(grav_cell, true);
        } else {
            // constant gravity
            grav_cell.setVal(grav_const);
        }
    } else { // spherical = 1

        BaseState<Real> m(1, nr_fine);

        m(0,0) = 4.0/3.0*M_PI*rho0(0,0)*r_cc_loc_b(0,0)*r_cc_loc_b(0,0)*r_cc_loc_b(0,0);
        grav_cell(0,0) = -Gconst * m(0,0) / (r_cc_loc_b(0,0)*r_cc_loc_b(0,0));

        for (auto r = 1; r < nr_fine; ++r) {

            // the mass is defined at the cell-centers, so to compute
            // the mass at the current center, we need to add the
            // contribution of the upper half of the zone below us and
            // the lower half of the current zone.

            // don't add any contributions from outside the star --
            // i.e.  rho < base_cutoff_density
            Real term1 = 0.0;
            if (rho0(0,r-1) > base_cutoff_density) {
                term1 = 4.0/3.0*M_PI*rho0(0,r-1) *
                    (r_edge_loc_b(0,r) - r_cc_loc_b(0,r-1)) *
                    (r_edge_loc_b(0,r)*r_edge_loc_b(0,r) +
                    r_edge_loc_b(0,r)*r_cc_loc_b(0,r-1) +
                    r_cc_loc_b(0,r-1)*r_cc_loc_b(0,r-1));
            } 

            Real term2 = 0.0;
            if (rho0(0,r) > base_cutoff_density) {
                term2 = 4.0/3.0*M_PI*rho0(0,r)*
                    (r_cc_loc_b(0,r) - r_edge_loc_b(0,r)) *
                    (r_cc_loc_b(0,r)*r_cc_loc_b(0,r) +
                    r_cc_loc_b(0,r)*r_edge_loc_b(0,r) +
                    r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
            } 

            m(0,r) = m(0,r-1) + term1 + term2;

            grav_cell(0,r) = -Gconst * m(0,r) / (r_cc_loc_b(0,r)*r_cc_loc_b(0,r));
        }
    }
}

void
Maestro::MakeGravEdge(RealVector& grav_edge, 
                      const RealVector& rho0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravEdge()",MakeGravEdge);

    const int max_lev = max_radial_level+1;

    get_base_cutoff_density(&base_cutoff_density);

    if (!spherical) {
        if (do_planar_invsq_grav)  {
            Real * AMREX_RESTRICT grav_edge_p = grav_edge.dataPtr();
            const auto r_edge_loc_p = r_edge_loc_b;
            const Real planar_invsq_mass_loc = planar_invsq_mass;
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            //   
            for (auto n = 0; n <= finest_radial_level; ++n) {
                // for (auto r = 0; r < nr(n); ++r) {
                const int nr_lev = nr(n);
                AMREX_PARALLEL_FOR_1D(nr_lev, r, {
                    grav_edge_p[n+max_lev*r] = -Gconst*planar_invsq_mass_loc / (r_edge_loc_p(n,r)*r_edge_loc_p(n,r));
                });
            }
        } else if (do_2d_planar_octant) {
            // compute gravity as in spherical geometry

            BaseState<Real> m(finest_radial_level+1, nr_fine+1);

            grav_edge[0] = 0.0;
            m(0,0) = 0.0;

            for (auto r = 0; r < nr(0); ++r) {

                // only add to the enclosed mass if the density is
                // > base_cutoff_density
                if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                    m(0,r) = m(0,r-1) + 4.0/3.0*M_PI *
                        (r_edge_loc_b(0,r) - r_edge_loc_b(0,r-1)) *
                        (r_edge_loc_b(0,r)*r_edge_loc_b(0,r) +
                        r_edge_loc_b(0,r)*r_edge_loc_b(0,r-1) +
                        r_edge_loc_b(0,r-1)*r_edge_loc_b(0,r-1)) * rho0[max_lev*(r-1)];
                } else {
                    m(0,r) = m(0,r-1);
                }

                grav_edge[max_lev*r] = -Gconst * m(0,r) / (r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
            }

            for (auto n = 1; n <= finest_radial_level; ++n) {
                for (auto i = 1; i <= numdisjointchunks(n); ++i) {

                    if (r_start_coord(n,i) == 0) {
                        m(n,0) = 0.0;
                    } else {
                        m(n,r_start_coord(n,i)) = m(n-1,r_start_coord(n,i)/2);
                        grav_edge[n+max_lev*r_start_coord(n,i)] = grav_edge[n-1+max_lev*r_start_coord(n,i)/2];
                    }

                    for (auto r = r_start_coord(n,i)+1; 
                         r <= r_end_coord(n,i)+1; ++r) {

                        // only add to the enclosed mass if the density is
                        // > base_cutoff_density
                        if (rho0[n+max_lev*(r-1)] > base_cutoff_density) {
                            m(n,r) = m(n,r-1) + 4.0/3.0*M_PI *
                                (r_edge_loc_b(n,r) - r_edge_loc_b(n,r-1)) *
                                (r_edge_loc_b(n,r)*r_edge_loc_b(n,r) +
                                r_edge_loc_b(n,r)*r_edge_loc_b(n,r-1) +
                                r_edge_loc_b(n,r-1)*r_edge_loc_b(n,r-1)) * rho0[n+max_lev*(r-1)];
                        } else {
                            m(n,r) = m(n,r-1);
                        }

                        grav_edge[n+max_lev*r] = -Gconst * m(n,r) / (r_edge_loc_b(n,r)*r_edge_loc_b(n,r));
                    }
                }
            }
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

            // only add to the enclosed mass if the density is
            // > base_cutoff_density
            if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                mencl += 4.0/3.0 * M_PI *
                    (r_edge_loc_b(0,r) - r_edge_loc_b(0,r-1)) *
                    (r_edge_loc_b(0,r) * r_edge_loc_b(0,r) +
                    r_edge_loc_b(0,r) * r_edge_loc_b(0,r-1) +
                    r_edge_loc_b(0,r-1)*r_edge_loc_b(0,r-1)) * rho0[max_lev*(r-1)];
            }

            grav_edge[max_lev*r] = -Gconst * mencl / (r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
        }
    }
}

void
Maestro::MakeGravEdge(BaseState<Real>& grav_edge, 
                      const RealVector& rho0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravEdge()",MakeGravEdge);

    const int max_lev = max_radial_level+1;

    get_base_cutoff_density(&base_cutoff_density);

    if (!spherical) {
        if (do_planar_invsq_grav)  {
            const auto r_edge_loc_p = r_edge_loc_b;
            const Real planar_invsq_mass_loc = planar_invsq_mass;
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            //   
            for (auto n = 0; n <= finest_radial_level; ++n) {
                // for (auto r = 0; r < nr(n); ++r) {
                const int nr_lev = nr(n);
                AMREX_PARALLEL_FOR_1D(nr_lev, r, {
                    grav_edge(n,r) = -Gconst*planar_invsq_mass_loc / (r_edge_loc_p(n,r)*r_edge_loc_p(n,r));
                });
            }
        } else if (do_2d_planar_octant) {
            // compute gravity as in spherical geometry

            BaseState<Real> m(finest_radial_level+1, nr_fine+1);

            grav_edge[0] = 0.0;
            m(0,0) = 0.0;

            for (auto r = 0; r < nr(0); ++r) {

                // only add to the enclosed mass if the density is
                // > base_cutoff_density
                if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                    m(0,r) = m(0,r-1) + 4.0/3.0*M_PI *
                        (r_edge_loc_b(0,r) - r_edge_loc_b(0,r-1)) *
                        (r_edge_loc_b(0,r)*r_edge_loc_b(0,r) +
                        r_edge_loc_b(0,r)*r_edge_loc_b(0,r-1) +
                        r_edge_loc_b(0,r-1)*r_edge_loc_b(0,r-1)) * rho0[max_lev*(r-1)];
                } else {
                    m(0,r) = m(0,r-1);
                }

                grav_edge(0,r) = -Gconst * m(0,r) / (r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
            }

            for (auto n = 1; n <= finest_radial_level; ++n) {
                for (auto i = 1; i <= numdisjointchunks(n); ++i) {

                    if (r_start_coord(n,i) == 0) {
                        m(n,0) = 0.0;
                    } else {
                        m(n,r_start_coord(n,i)) = m(n-1,r_start_coord(n,i)/2);
                        grav_edge(n,r_start_coord(n,i)) = grav_edge(n-1,r_start_coord(n,i)/2);
                    }

                    for (auto r = r_start_coord(n,i)+1; 
                         r <= r_end_coord(n,i)+1; ++r) {

                        // only add to the enclosed mass if the density is
                        // > base_cutoff_density
                        if (rho0[n+max_lev*(r-1)] > base_cutoff_density) {
                            m(n,r) = m(n,r-1) + 4.0/3.0*M_PI *
                                (r_edge_loc_b(n,r) - r_edge_loc_b(n,r-1)) *
                                (r_edge_loc_b(n,r)*r_edge_loc_b(n,r) +
                                r_edge_loc_b(n,r)*r_edge_loc_b(n,r-1) +
                                r_edge_loc_b(n,r-1)*r_edge_loc_b(n,r-1)) * rho0[n+max_lev*(r-1)];
                        } else {
                            m(n,r) = m(n,r-1);
                        }

                        grav_edge(n,r) = -Gconst * m(n,r) / (r_edge_loc_b(n,r)*r_edge_loc_b(n,r));
                    }
                }
            }
            RestrictBase(grav_edge, false);
            FillGhostBase(grav_edge, false);
        } else {
            // constant gravity
            grav_edge.setVal(grav_const);
        }
        
    } else {

        grav_edge[0] = 0.0;
        Real mencl = 0.0;

        for (auto r = 1; r <= nr_fine; ++r) {

            // only add to the enclosed mass if the density is
            // > base_cutoff_density
            if (rho0[max_lev*(r-1)] > base_cutoff_density) {
                mencl += 4.0/3.0 * M_PI *
                    (r_edge_loc_b(0,r) - r_edge_loc_b(0,r-1)) *
                    (r_edge_loc_b(0,r) * r_edge_loc_b(0,r) +
                    r_edge_loc_b(0,r) * r_edge_loc_b(0,r-1) +
                    r_edge_loc_b(0,r-1)*r_edge_loc_b(0,r-1)) * rho0[max_lev*(r-1)];
            }

            grav_edge(0,r) = -Gconst * mencl / (r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
        }
    }
}

void
Maestro::MakeGravEdge(BaseState<Real>& grav_edge, 
                      const BaseState<Real>& rho0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravEdge()", MakeGravEdge);

    const int max_lev = max_radial_level+1;

    get_base_cutoff_density(&base_cutoff_density);

    if (!spherical) {
        if (do_planar_invsq_grav)  {
            const auto r_edge_loc_p = r_edge_loc_b;
            const Real planar_invsq_mass_loc = planar_invsq_mass;
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            //   
            for (auto n = 0; n <= finest_radial_level; ++n) {
                // for (auto r = 0; r < nr(n); ++r) {
                const int nr_lev = nr(n);
                AMREX_PARALLEL_FOR_1D(nr_lev, r, {
                    grav_edge(n,r) = -Gconst*planar_invsq_mass_loc / (r_edge_loc_p(n,r)*r_edge_loc_p(n,r));
                });
            }
        } else if (do_2d_planar_octant) {
            // compute gravity as in spherical geometry

            BaseState<Real> m(finest_radial_level+1, nr_fine+1);

            grav_edge(0,0) = 0.0;
            m(0,0) = 0.0;

            for (auto r = 0; r < nr(0); ++r) {

                // only add to the enclosed mass if the density is
                // > base_cutoff_density
                if (rho0(0,r-1) > base_cutoff_density) {
                    m(0,r) = m(0,r-1) + 4.0/3.0*M_PI *
                        (r_edge_loc_b(0,r) - r_edge_loc_b(0,r-1)) *
                        (r_edge_loc_b(0,r)*r_edge_loc_b(0,r) +
                        r_edge_loc_b(0,r)*r_edge_loc_b(0,r-1) +
                        r_edge_loc_b(0,r-1)*r_edge_loc_b(0,r-1)) * rho0(0,r-1);
                } else {
                    m(0,r) = m(0,r-1);
                }

                grav_edge(0,r) = -Gconst * m(0,r) / (r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
            }

            for (auto n = 1; n <= finest_radial_level; ++n) {
                for (auto i = 1; i <= numdisjointchunks(n); ++i) {

                    if (r_start_coord(n,i) == 0) {
                        m(n,0) = 0.0;
                    } else {
                        m(n,r_start_coord(n,i)) = m(n-1,r_start_coord(n,i)/2);
                        grav_edge(n,r_start_coord(n,i)) = grav_edge(n-1,r_start_coord(n,i)/2);
                    }

                    for (auto r = r_start_coord(n,i)+1; 
                         r <= r_end_coord(n,i)+1; ++r) {

                        // only add to the enclosed mass if the density is
                        // > base_cutoff_density
                        if (rho0(n,r-1) > base_cutoff_density) {
                            m(n,r) = m(n,r-1) + 4.0/3.0*M_PI *
                                (r_edge_loc_b(n,r) - r_edge_loc_b(n,r-1)) *
                                (r_edge_loc_b(n,r)*r_edge_loc_b(n,r) +
                                r_edge_loc_b(n,r)*r_edge_loc_b(n,r-1) +
                                r_edge_loc_b(n,r-1)*r_edge_loc_b(n,r-1)) * rho0(n,r-1);
                        } else {
                            m(n,r) = m(n,r-1);
                        }

                        grav_edge(n,r) = -Gconst * m(n,r) / (r_edge_loc_b(n,r)*r_edge_loc_b(n,r));
                    }
                }
            }
            RestrictBase(grav_edge, false);
            FillGhostBase(grav_edge, false);
        } else {
            // constant gravity
            // std::fill(grav_edge.begin(), grav_edge.end(), grav_const);
            grav_edge.setVal(grav_const);
        }
        
    } else {

        grav_edge(0,0) = 0.0;
        Real mencl = 0.0;

        for (auto r = 1; r <= nr_fine; ++r) {

            // only add to the enclosed mass if the density is
            // > base_cutoff_density
            if (rho0(0,r-1) > base_cutoff_density) {
                mencl += 4.0/3.0 * M_PI *
                    (r_edge_loc_b(0,r) - r_edge_loc_b(0,r-1)) *
                    (r_edge_loc_b(0,r) * r_edge_loc_b(0,r) +
                    r_edge_loc_b(0,r) * r_edge_loc_b(0,r-1) +
                    r_edge_loc_b(0,r-1)*r_edge_loc_b(0,r-1)) * rho0(0,r-1);
            }

            grav_edge(0,r) = -Gconst * mencl / (r_edge_loc_b(0,r)*r_edge_loc_b(0,r));
        }
    }
}