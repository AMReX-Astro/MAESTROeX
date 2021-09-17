#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::EnforceHSE(const BaseState<Real>& rho0_s, BaseState<Real>& p0_s,
                         const BaseState<Real>& grav_cell_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EnforceHSE()", EnforceHSE);

    const auto& dr = base_geom.dr;
    const auto& r_start_coord = base_geom.r_start_coord;
    const auto& r_end_coord = base_geom.r_end_coord;

    BaseState<Real> grav_edge_s(base_geom.max_radial_level + 1,
                                base_geom.nr_fine + 1);
    BaseState<Real> p0old_s(base_geom.max_radial_level + 1, base_geom.nr_fine);
    auto grav_edge = grav_edge_s.array();
    auto p0old = p0old_s.array();
    const auto rho0 = rho0_s.const_array();
    auto p0 = p0_s.array();
    const auto grav_cell = grav_cell_s.const_array();

    Real offset = 0.0;

    // const BaseState<Real> rho0_b(rho0, base_geom.max_radial_level+1, base_geom.nr_fine);
    MakeGravEdge(grav_edge_s, rho0_s);

    // create a copy of the input pressure to help us with initial
    // conditions
    p0old_s.copy(p0_s);

    // zero the new pressure so we don't leave a non-zero pressure in
    // fine radial regions that no longer have a corresponding full
    // state
    p0_s.setVal(0.0);

    // integrate all of level 1 first
    // use the old pressure at r=0 as a reference point
    p0(0, 0) = p0old(0, 0);

    // now integrate upwards from the bottom later, we will offset the
    // entire pressure so we have effectively integrated from the "top"
    if (use_exact_base_state && spherical) {
        for (auto r = 1;
             r <= amrex::min(r_end_coord(0, 1),
                             base_geom.base_cutoff_density_coord(0));
             ++r) {
            // uneven grid spacing
            Real dr1 =
                base_geom.r_edge_loc(0, r) - base_geom.r_cc_loc(0, r - 1);
            Real dr2 = base_geom.r_cc_loc(0, r) - base_geom.r_edge_loc(0, r);
            p0(0, r) =
                p0(0, r - 1) +
                (dr1 * rho0(0, r - 1) + dr2 * rho0(0, r)) * grav_edge(0, r);
        }
    } else {
        for (auto r = 1;
             r <= amrex::min(r_end_coord(0, 1),
                             base_geom.base_cutoff_density_coord(0));
             r++) {
            // assume even grid spacing
            p0(0, r) = p0(0, r - 1) + 0.5 * dr(0) *
                                          (rho0(0, r - 1) + rho0(0, r)) *
                                          grav_edge(0, r);
        }
    }
    for (auto r = base_geom.base_cutoff_density_coord(0) + 1;
         r <= base_geom.r_end_coord(0, 1); ++r) {
        p0(0, r) = p0(0, r - 1);
    }

    if (!spherical) {
        for (auto n = 1; n <= base_geom.finest_radial_level; ++n) {
            for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
                // get pressure in the bottom cell of this disjointchunk
                if (r_start_coord(n, i) == 0) {
                    // if we are at the bottom of the domain, use the old
                    // pressure as reference
                    p0(n, 0) = p0old(n, 0);
                } else if (r_start_coord(n, i) <=
                           base_geom.base_cutoff_density_coord(n)) {
                    // we integrate upwards starting from the nearest coarse
                    // cell at a lower physical height

                    if (do_planar_invsq_grav || do_2d_planar_octant) {
                        // we have variable gravity
                        p0(n, r_start_coord(n, i)) =
                            p0(n - 1, r_start_coord(n, i) / 2 - 1) +
                            (dr(n) / 4.0) *
                                (2.0 * rho0(n, r_start_coord(n, i)) / 3.0 +
                                 4.0 *
                                     rho0(n - 1, r_start_coord(n, i) / 2 - 1) /
                                     3.0) *
                                (grav_edge(n, r_start_coord(n, i)) +
                                 grav_cell(n - 1,
                                           r_start_coord(n, i) / 2 - 1)) +
                            (dr(n) / 8.0) *
                                (5.0 * rho0(n, r_start_coord(n, i)) / 3.0 +
                                 1.0 *
                                     rho0(n - 1, r_start_coord(n, i) / 2 - 1) /
                                     3.0) *
                                (grav_edge(n, r_start_coord(n, i)) +
                                 grav_cell(n, r_start_coord(n, i)));
                    } else {
                        // assuming constant g here
                        p0(n, r_start_coord(n, i)) =
                            p0(n - 1, r_start_coord(n, i) / 2 - 1) +
                            (3.0 * grav_cell(1, 0) * dr(n) / 4.0) *
                                (rho0(n - 1, r_start_coord(n, i) / 2 - 1) +
                                 rho0(n, r_start_coord(n, i)));
                    }
                } else {
                    // copy pressure from below
                    p0(n, r_start_coord(n, i)) =
                        p0(n - 1, r_start_coord(n, i) / 2 - 1);
                }

                // integrate upwards as normal
                for (auto r = r_start_coord(n, i) + 1;
                     r <= amrex::min(r_end_coord(n, i),
                                     base_geom.base_cutoff_density_coord(n));
                     ++r) {
                    p0(n, r) = p0(n, r - 1) +
                               0.5 * dr(n) * (rho0(n, r) + rho0(n, r - 1)) *
                                   grav_edge(n, r);
                }
                for (auto r = base_geom.base_cutoff_density_coord(n) + 1;
                     r <= r_end_coord(n, i); ++r) {
                    p0(n, r) = p0(n, r - 1);
                }

                // now we need to look at the first coarser cell above this
                // disjoint chunk and look at the mismatch between the
                // integration performed over the finer cells vs. the coarse
                // cells.  Then we need to offset the coarse cell above this
                // point to sync up.

                // first, compute the value of the pressure in the coarse
                // cell above the disjointchunk.
                if (r_end_coord(n, i) == base_geom.nr(n) - 1) {
                    // for (auto nothing - we are at the top of the domain
                    offset = 0.0;
                } else if (r_end_coord(n, i) <=
                           base_geom.base_cutoff_density_coord(n)) {
                    // use fine -> coarse stencil in notes
                    if (do_planar_invsq_grav || do_2d_planar_octant) {
                        // we have variable gravity
                        Real temp =
                            p0(n, r_end_coord(n, i)) +
                            (dr(n) / 4.0) *
                                (2.0 * rho0(n, r_end_coord(n, i)) / 3.0 +
                                 4.0 *
                                     rho0(n - 1, (r_end_coord(n, i) + 1) / 2) /
                                     3.0) *
                                (grav_edge(n - 1, (r_end_coord(n, i) + 1) / 2) +
                                 grav_cell(n - 1,
                                           (r_end_coord(n, i) + 1) / 2)) +
                            (dr(n) / 8.0) *
                                (5.0 * rho0(n, r_end_coord(n, i)) / 3.0 +
                                 1.0 *
                                     rho0(n - 1, (r_end_coord(n, i) + 1) / 2) /
                                     3.0) *
                                (grav_cell(n, r_end_coord(n, i)) +
                                 grav_edge(n - 1, (r_end_coord(n, i) + 1) / 2));
                        offset = p0(n - 1, (r_end_coord(n, i) + 1) / 2) - temp;
                    } else {
                        // assuming constant g here
                        Real temp =
                            p0(n, r_end_coord(n, i)) +
                            (3.0 * grav_cell(1, 0) * dr(n) / 4.0) *
                                (rho0(n, r_end_coord(n, i)) +
                                 rho0(n - 1, (r_end_coord(n, i) + 1) / 2));
                        offset = p0(n - 1, (r_end_coord(n, i) + 1) / 2) - temp;
                    }
                } else {
                    // copy pressure from below
                    Real temp = p0(n, r_end_coord(n, i));
                    offset = p0(n - 1, (r_end_coord(n, i) + 1) / 2) - temp;
                }

                // if we are not at the top of the domain, we need to
                // subtract the offset for all values at and above this point
                if (r_end_coord(n, i) != base_geom.nr(n) - 1) {
                    for (auto l = n - 1; l >= 0; --l) {
                        for (auto r = (int)amrex::Math::round(
                                 (r_end_coord(n, i) + 1) / pow(2, n - l));
                             r <= base_geom.nr(l) - 1; ++r) {
                            p0(l, r) -= offset;
                        }
                    }
                }
            }  // end loop over disjoint chunks
        }      // end loop over levels
    }          // spherical

    // if the top is closed in planar geometry, offset the pressure to
    // make it consistent with that boundary condition
    //
    // otherwise, now compare pressure in the last cell and offset to
    // make sure we are integrating "from the top"
    // we use the coarsest level as the reference point
    if (add_pb && !spherical) {
        offset = -1.0 * p0b;
    } else {
        offset = p0(0, base_geom.nr(0) - 1) - p0old(0, base_geom.nr(0) - 1);
    }

    // offset level 0
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        p0(0, r) -= offset;
    }

    // offset remaining levels
    for (auto n = 1; n <= base_geom.finest_radial_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            for (auto r = r_start_coord(n, i); r <= r_end_coord(n, i); ++r) {
                p0(n, r) -= offset;
            }
        }
    }

    // zero p0 where there is no corresponding full state array
    for (auto n = 1; n <= base_geom.finest_radial_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            if (i == base_geom.numdisjointchunks(n)) {
                for (auto r = r_end_coord(n, i) + 1; r < base_geom.nr(n); ++r) {
                    p0(n, r) = 0.0;
                }
            } else {
                for (auto r = r_end_coord(n, i) + 1;
                     r < r_start_coord(n, i + 1); ++r) {
                    p0(n, r) = 0.0;
                }
            }
        }
    }

    RestrictBase(p0, true);
    FillGhostBase(p0, true);
}
