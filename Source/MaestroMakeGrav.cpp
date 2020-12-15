#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::MakeGravCell(BaseState<Real>& grav_cell,
                           const BaseState<Real>& rho0_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravCell()", MakeGravCell);

    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    auto grav_cell_arr = grav_cell.array();
    const auto rho0 = rho0_s.const_array();

    if (!spherical) {
        if (do_planar_invsq_grav) {
            const Real planar_invsq_mass_loc = planar_invsq_mass;
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
                const int nr_lev = base_geom.nr(n);
                ParallelFor(nr_lev, [=] AMREX_GPU_DEVICE(long r) {
                    grav_cell_arr(n, r) = -Gconst * planar_invsq_mass_loc /
                                          (r_cc_loc(n, r) * r_cc_loc(n, r));
                });
                Gpu::synchronize();
            }
        } else if (do_2d_planar_octant) {
            //   compute gravity as in the spherical case
            BaseState<Real> m_state(base_geom.finest_radial_level + 1,
                                    base_geom.nr_fine);
            auto m = m_state.array();

            // level = 0
            m(0, 0) = 4.0 / 3.0 * M_PI * rho0(0, 0) * r_cc_loc(0, 0) *
                      r_cc_loc(0, 0) * r_cc_loc(0, 0);
            grav_cell_arr(0, 0) =
                -Gconst * m(0, 0) / (r_cc_loc(0, 0) * r_cc_loc(0, 0));

            for (auto r = 1; r < base_geom.nr(0); ++r) {
                // the mass is defined at the cell-centers, so to compute
                // the mass at the current center, we need to add the
                // contribution of the upper half of the zone below us and
                // the lower half of the current zone.

                // don't add any contributions from outside the star --
                // i.e.  rho < base_cutoff_density
                Real term1 = 0.0;
                if (rho0(0, r - 1) > base_cutoff_density) {
                    term1 = 4.0 / 3.0 * M_PI * rho0(0, r - 1) *
                            (r_edge_loc(0, r) - r_cc_loc(0, r - 1)) *
                            (r_edge_loc(0, r) * r_edge_loc(0, r) +
                             r_edge_loc(0, r) * r_cc_loc(0, r - 1) +
                             r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1));
                }

                Real term2 = 0.0;
                if (rho0(0, r) > base_cutoff_density) {
                    term2 = 4.0 / 3.0 * M_PI * rho0(0, r) *
                            (r_cc_loc(0, r) - r_edge_loc(0, r)) *
                            (r_cc_loc(0, r) * r_cc_loc(0, r) +
                             r_cc_loc(0, r) * r_edge_loc(0, r) +
                             r_edge_loc(0, r) * r_edge_loc(0, r));
                }

                m(0, r) = m(0, r - 1) + term1 + term2;

                grav_cell_arr(0, r) =
                    -Gconst * m(0, r) / (r_cc_loc(0, r) * r_cc_loc(0, r));
            }

            // level > 0

            for (auto n = 1; n <= base_geom.finest_radial_level; ++n) {
                for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
                    if (base_geom.r_start_coord(n, i) == 0) {
                        m(n, 0) = 4.0 / 3.0 * M_PI * rho0(n, 0) *
                                  r_cc_loc(n, 0) * r_cc_loc(n, 0) *
                                  r_cc_loc(n, 0);
                        grav_cell_arr(n, 0) = -Gconst * m(n, 0) /
                                              (r_cc_loc(n, 0) * r_cc_loc(n, 0));
                    } else {
                        int r = base_geom.r_start_coord(n, i);
                        m(n, r) = m(n - 1, r / 2 - 1);

                        // the mass is defined at the cell-centers, so to compute
                        // the mass at the current center, we need to add the
                        // contribution of the upper half of the zone below us and
                        // the lower half of the current zone.

                        // don't add any contributions from outside the star --
                        // i.e.  rho < base_cutoff_density
                        Real term1 = 0.0;
                        if (rho0(n - 1, r / 2 - 1) > base_cutoff_density) {
                            term1 = 4.0 / 3.0 * M_PI * rho0(n - 1, r / 2 - 1) *
                                    (r_edge_loc(n - 1, r / 2) -
                                     r_cc_loc(n - 1, r / 2 - 1)) *
                                    (r_edge_loc(n - 1, r / 2) *
                                         r_edge_loc(n - 1, r / 2) +
                                     r_edge_loc(n - 1, r / 2) *
                                         r_cc_loc(n - 1, r / 2 - 1) +
                                     r_cc_loc(n - 1, r / 2 - 1) *
                                         r_cc_loc(n - 1, r / 2 - 1));
                        }

                        Real term2 = 0.0;
                        if (rho0(n, r) > base_cutoff_density) {
                            term2 = 4.0 / 3.0 * M_PI * rho0(n, r) *
                                    (r_cc_loc(n, r) - r_edge_loc(n, r)) *
                                    (r_cc_loc(n, r) * r_cc_loc(n, r) +
                                     r_cc_loc(n, r) * r_edge_loc(n, r) +
                                     r_edge_loc(n, r) * r_edge_loc(n, r));
                        }

                        m(n, r) += term1 + term2;
                        grav_cell_arr(n, r) = -Gconst * m(n, r) /
                                              (r_cc_loc(n, r) * r_cc_loc(n, r));
                    }

                    for (auto r = base_geom.r_start_coord(n, i) + 1;
                         r <= base_geom.r_end_coord(n, i); ++r) {
                        // the mass is defined at the cell-centers, so to compute
                        // the mass at the current center, we need to add the
                        // contribution of the upper half of the zone below us and
                        // the lower half of the current zone.

                        // don't add any contributions from outside the star --
                        // i.e.  rho < base_cutoff_density
                        Real term1 = 0.0;
                        if (rho0(n, r - 1) > base_cutoff_density) {
                            term1 = 4.0 / 3.0 * M_PI * rho0(n, r - 1) *
                                    (r_edge_loc(n, r) - r_cc_loc(n, r - 1)) *
                                    (r_edge_loc(n, r) * r_edge_loc(n, r) +
                                     r_edge_loc(n, r) * r_cc_loc(n, r - 1) +
                                     r_cc_loc(n, r - 1) * r_cc_loc(n, r - 1));
                        }

                        Real term2 = 0.0;
                        if (rho0(n, r) > base_cutoff_density) {
                            term2 = 4.0 / 3.0 * M_PI * rho0(n, r) *
                                    (r_cc_loc(n, r) - r_edge_loc(n, r)) *
                                    (r_cc_loc(n, r) * r_cc_loc(n, r) +
                                     r_cc_loc(n, r) * r_edge_loc(n, r) +
                                     r_edge_loc(n, r) * r_edge_loc(n, r));
                        }

                        m(n, r) = m(n, r - 1) + term1 + term2;

                        grav_cell_arr(n, r) = -Gconst * m(n, r) /
                                              (r_cc_loc(n, r) * r_cc_loc(n, r));
                    }
                }
            }

            RestrictBase(grav_cell, true);
            FillGhostBase(grav_cell, true);
        } else {
            // constant gravity
            grav_cell.setVal(grav_const);
        }
    } else {  // spherical = 1

        BaseState<Real> m_state(1, base_geom.nr_fine);
        auto m = m_state.array();

        m(0, 0) = 4.0 / 3.0 * M_PI * rho0(0, 0) * r_cc_loc(0, 0) *
                  r_cc_loc(0, 0) * r_cc_loc(0, 0);
        grav_cell_arr(0, 0) =
            -Gconst * m(0, 0) / (r_cc_loc(0, 0) * r_cc_loc(0, 0));

        for (auto r = 1; r < base_geom.nr_fine; ++r) {
            // the mass is defined at the cell-centers, so to compute
            // the mass at the current center, we need to add the
            // contribution of the upper half of the zone below us and
            // the lower half of the current zone.

            // don't add any contributions from outside the star --
            // i.e.  rho < base_cutoff_density
            Real term1 = 0.0;
            if (rho0(0, r - 1) > base_cutoff_density) {
                term1 = 4.0 / 3.0 * M_PI * rho0(0, r - 1) *
                        (r_edge_loc(0, r) - r_cc_loc(0, r - 1)) *
                        (r_edge_loc(0, r) * r_edge_loc(0, r) +
                         r_edge_loc(0, r) * r_cc_loc(0, r - 1) +
                         r_cc_loc(0, r - 1) * r_cc_loc(0, r - 1));
            }

            Real term2 = 0.0;
            if (rho0(0, r) > base_cutoff_density) {
                term2 = 4.0 / 3.0 * M_PI * rho0(0, r) *
                        (r_cc_loc(0, r) - r_edge_loc(0, r)) *
                        (r_cc_loc(0, r) * r_cc_loc(0, r) +
                         r_cc_loc(0, r) * r_edge_loc(0, r) +
                         r_edge_loc(0, r) * r_edge_loc(0, r));
            }

            m(0, r) = m(0, r - 1) + term1 + term2;

            grav_cell_arr(0, r) =
                -Gconst * m(0, r) / (r_cc_loc(0, r) * r_cc_loc(0, r));
        }
    }
}

void Maestro::MakeGravEdge(BaseState<Real>& grav_edge_state,
                           const BaseState<Real>& rho0_state) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravEdge()", MakeGravEdge);

    const auto& r_edge_loc = base_geom.r_edge_loc;
    auto grav_edge = grav_edge_state.array();
    const auto rho0 = rho0_state.const_array();

    if (!spherical) {
        if (do_planar_invsq_grav) {
            const Real planar_invsq_mass_loc = planar_invsq_mass;
            // we are doing a plane-parallel geometry with a 1/r**2
            // gravitational acceleration.  The mass is assumed to be
            // at the origin.  The mass in the computational domain
            // does not contribute to the gravitational acceleration.
            //
            for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
                const int nr_lev = base_geom.nr(n);
                ParallelFor(nr_lev, [=] AMREX_GPU_DEVICE(long r) {
                    grav_edge(n, r) = -Gconst * planar_invsq_mass_loc /
                                      (r_edge_loc(n, r) * r_edge_loc(n, r));
                });
                Gpu::synchronize();
            }
        } else if (do_2d_planar_octant) {
            // compute gravity as in spherical geometry

            BaseState<Real> m_state(base_geom.finest_radial_level + 1,
                                    base_geom.nr_fine + 1);
            auto m = m_state.array();

            grav_edge(0, 0) = 0.0;
            m(0, 0) = 0.0;

            for (auto r = 0; r < base_geom.nr(0); ++r) {
                // only add to the enclosed mass if the density is
                // > base_cutoff_density
                if (rho0(0, r - 1) > base_cutoff_density) {
                    m(0, r) =
                        m(0, r - 1) +
                        4.0 / 3.0 * M_PI *
                            (r_edge_loc(0, r) - r_edge_loc(0, r - 1)) *
                            (r_edge_loc(0, r) * r_edge_loc(0, r) +
                             r_edge_loc(0, r) * r_edge_loc(0, r - 1) +
                             r_edge_loc(0, r - 1) * r_edge_loc(0, r - 1)) *
                            rho0(0, r - 1);
                } else {
                    m(0, r) = m(0, r - 1);
                }

                grav_edge(0, r) =
                    -Gconst * m(0, r) / (r_edge_loc(0, r) * r_edge_loc(0, r));
            }

            for (auto n = 1; n <= base_geom.finest_radial_level; ++n) {
                for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
                    if (base_geom.r_start_coord(n, i) == 0) {
                        m(n, 0) = 0.0;
                    } else {
                        m(n, base_geom.r_start_coord(n, i)) =
                            m(n - 1, base_geom.r_start_coord(n, i) / 2);
                        grav_edge(n, base_geom.r_start_coord(n, i)) =
                            grav_edge(n - 1, base_geom.r_start_coord(n, i) / 2);
                    }

                    for (auto r = base_geom.r_start_coord(n, i) + 1;
                         r <= base_geom.r_end_coord(n, i) + 1; ++r) {
                        // only add to the enclosed mass if the density is
                        // > base_cutoff_density
                        if (rho0(n, r - 1) > base_cutoff_density) {
                            m(n, r) =
                                m(n, r - 1) +
                                4.0 / 3.0 * M_PI *
                                    (r_edge_loc(n, r) - r_edge_loc(n, r - 1)) *
                                    (r_edge_loc(n, r) * r_edge_loc(n, r) +
                                     r_edge_loc(n, r) * r_edge_loc(n, r - 1) +
                                     r_edge_loc(n, r - 1) *
                                         r_edge_loc(n, r - 1)) *
                                    rho0(n, r - 1);
                        } else {
                            m(n, r) = m(n, r - 1);
                        }

                        grav_edge(n, r) = -Gconst * m(n, r) /
                                          (r_edge_loc(n, r) * r_edge_loc(n, r));
                    }
                }
            }
            RestrictBase(grav_edge, false);
            FillGhostBase(grav_edge, false);
        } else {
            // constant gravity
            grav_edge_state.setVal(grav_const);
        }

    } else {
        grav_edge(0, 0) = 0.0;
        Real mencl = 0.0;

        for (auto r = 1; r <= base_geom.nr_fine; ++r) {
            // only add to the enclosed mass if the density is
            // > base_cutoff_density
            if (rho0(0, r - 1) > base_cutoff_density) {
                mencl += 4.0 / 3.0 * M_PI *
                         (r_edge_loc(0, r) - r_edge_loc(0, r - 1)) *
                         (r_edge_loc(0, r) * r_edge_loc(0, r) +
                          r_edge_loc(0, r) * r_edge_loc(0, r - 1) +
                          r_edge_loc(0, r - 1) * r_edge_loc(0, r - 1)) *
                         rho0(0, r - 1);
            }

            grav_edge(0, r) =
                -Gconst * mencl / (r_edge_loc(0, r) * r_edge_loc(0, r));
        }
    }
}
