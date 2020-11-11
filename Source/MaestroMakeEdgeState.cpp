#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::MakeEdgeState1d(BaseState<Real>& s, BaseState<Real>& sedge,
                              BaseState<Real>& force) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1d()", MakeEdgeState1d);

    if (spherical) {
        MakeEdgeState1dSphr(s, sedge, force);
    } else {
        MakeEdgeState1dPlanar(s, sedge, force);
    }
}

void Maestro::MakeEdgeState1dSphr(BaseState<Real>& s_state,
                                  BaseState<Real>& sedge_state,
                                  BaseState<Real>& force_state) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1dSphr()", MakeEdgeState1dSphr);

    const auto rel_eps_local = rel_eps;

    const Real dth = 0.5 * dt;
    const Real C = 1.25;
    const int cen = 0;
    const int lim = 1;
    const int flag = 2;
    const int fromm = 3;
    const Real dr0 = base_geom.dr(0);
    const Real dtdr = dt / dr0;
    const auto nr_fine = base_geom.nr_fine;

    const int slope_order_loc = slope_order;
    auto s = s_state.array();
    auto sedge = sedge_state.array();
    auto force = force_state.array();

    BaseState<Real> sedgel_state(nr_fine + 1);
    BaseState<Real> sedger_state(nr_fine + 1);
    auto sedgel = sedgel_state.array();
    auto sedger = sedger_state.array();

    // copy valid data into array with ghost cells
    const int ng = 3;  // number of ghost cells
    BaseState<Real> s_ghost_state(nr_fine + 2 * ng);
    auto s_ghost = s_ghost_state.array();
    for (int i = 0; i < nr_fine; ++i) {
        s_ghost(i + ng) = s(0, i);
    }

    for (int i = 0; i < ng; i++) {
        // symmetry boundary condition at center
        s_ghost(ng - 1 - i) = s(0, i);
        // first-order extrapolation at top of star
        s_ghost(ng + nr_fine + i) = s(0, nr_fine - 1);
    }

    const auto w0_arr = w0.const_array();

    if (ppm_type == 0) {
        ParallelFor(nr_fine, [=] AMREX_GPU_DEVICE(long r) {
            Real slope = 0.0;

            // this will hold values at r-1, r and r+1
            Real dsscr[4 * 3];

            if (slope_order_loc == 0) {
                slope = 0.0;
            } else if (slope_order_loc == 2) {
                // index of ghost array offset by ng
                int g = r + ng;

                Real del = 0.5 * (s_ghost(g + 1) - s_ghost(g - 1));
                Real dpls = 2.0 * (s_ghost(g + 1) - s_ghost(g));
                Real dmin = 2.0 * (s_ghost(g) - s_ghost(g - 1));
                Real slim =
                    amrex::min(amrex::Math::abs(dpls), amrex::Math::abs(dmin));
                slim = dpls * dmin > 0.0 ? slim : 0.0;
                Real sflag = amrex::Math::copysign(1.0, del);
                slope = sflag * amrex::min(slim, amrex::Math::abs(del));
            } else if (slope_order_loc == 4) {
                for (int i = 0; i < 3; ++i) {
                    // index of ghost array offset by ng
                    int g = r + ng - 1 + i;
                    // do standard limiting to compute temporary slopes
                    dsscr[i * 4 + cen] =
                        0.5 * (s_ghost(g + 1) - s_ghost(g - 1));
                    Real dpls = 2.0 * (s_ghost(g + 1) - s_ghost(g));
                    Real dmin = 2.0 * (s_ghost(g) - s_ghost(g - 1));
                    dsscr[i * 4 + lim] = amrex::min(amrex::Math::abs(dmin),
                                                    amrex::Math::abs(dpls));
                    dsscr[i * 4 + lim] =
                        dpls * dmin > 0.0 ? dsscr[i * 4 + lim] : 0.0;
                    dsscr[i * 4 + flag] =
                        amrex::Math::copysign(1.0, dsscr[i * 4 + cen]);
                    dsscr[i * 4 + fromm] =
                        dsscr[i * 4 + flag] *
                        amrex::min(dsscr[i * 4 + lim],
                                   amrex::Math::abs(dsscr[i * 4 + cen]));
                }

                // fourth-order limited slopes
                Real ds = 4.0 / 3.0 * dsscr[1 * 4 + cen] -
                          (dsscr[2 * 4 + fromm] + dsscr[fromm]) / 6.0;
                slope = dsscr[1 * 4 + flag] *
                        amrex::min(amrex::Math::abs(ds), dsscr[1 * 4 + lim]);
            }

            // compute sedgel and sedger
            Real u = 0.5 * (w0_arr(0, r) + w0_arr(0, r + 1));
            Real ubardth =
                dth * u / dr0;  // NOTE: ubardth=0 for use_exact_base_state case
            sedgel(r + 1) =
                s(0, r) + (0.5 - ubardth) * slope + dth * force(0, r);
            sedger(r) = s(0, r) - (0.5 + ubardth) * slope + dth * force(0, r);
        });
        Gpu::synchronize();

    } else if (ppm_type == 1) {
        ParallelFor(nr_fine, [=] AMREX_GPU_DEVICE(long r) {
            // interpolate s to radial edges

            // sm

            // right side
            int g = r + ng;
            // compute van Leer slopes
            Real del = 0.5 * (s_ghost(g + 1) - s_ghost(g - 1));
            Real dmin = 2.0 * (s_ghost(g) - s_ghost(g - 1));
            Real dpls = 2.0 * (s_ghost(g + 1) - s_ghost(g));
            Real dsscrr = dmin * dpls > 0.0
                              ? amrex::Math::copysign(1.0, del) *
                                    amrex::min(amrex::Math::abs(del),
                                               amrex::Math::abs(dmin),
                                               amrex::Math::abs(dpls))
                              : 0.0;

            // left side
            g = r + ng - 1;
            // compute van Leer slopes
            del = 0.5 * (s_ghost(g + 1) - s_ghost(g - 1));
            dmin = 2.0 * (s_ghost(g) - s_ghost(g - 1));
            dpls = 2.0 * (s_ghost(g + 1) - s_ghost(g));
            Real dsscrl = dmin * dpls > 0.0
                              ? amrex::Math::copysign(1.0, del) *
                                    amrex::min(amrex::Math::abs(del),
                                               amrex::Math::abs(dmin),
                                               amrex::Math::abs(dpls))
                              : 0.0;

            // sm
            g = r + ng;
            // 4th order interpolation of s to radial faces
            Real sm =
                0.5 * (s_ghost(g) + s_ghost(g - 1)) - (dsscrr - dsscrl) / 6.0;
            // make sure sm lies in between adjacent cell-centered values
            sm = amrex::max(sm, amrex::min(s_ghost(g), s_ghost(g - 1)));
            sm = amrex::min(sm, amrex::max(s_ghost(g), s_ghost(g - 1)));

            // sp

            // left side
            dsscrl = dsscrr;

            // right side
            g = r + ng + 1;
            // compute van Leer slopes
            del = 0.5 * (s_ghost(g + 1) - s_ghost(g - 1));
            dmin = 2.0 * (s_ghost(g) - s_ghost(g - 1));
            dpls = 2.0 * (s_ghost(g + 1) - s_ghost(g));
            dsscrr = dmin * dpls > 0.0 ? amrex::Math::copysign(1.0, del) *
                                             amrex::min(amrex::Math::abs(del),
                                                        amrex::Math::abs(dmin),
                                                        amrex::Math::abs(dpls))
                                       : 0.0;

            // sp
            // g = i + ng + 1;
            // 4th order interpolation of s to radial faces
            Real sp =
                0.5 * (s_ghost(g) + s_ghost(g - 1)) - (dsscrr - dsscrl) / 6.0;
            // make sure sp lies in between adjacent cell-centered values
            sp = amrex::max(sp, amrex::min(s_ghost(g), s_ghost(g - 1)));
            sp = amrex::min(sp, amrex::max(s_ghost(g), s_ghost(g - 1)));

            // modify using quadratic limiters
            if ((sp - s(0, r)) * (s(0, r) - sm) <= 0.0) {
                sp = s(0, r);
                sm = s(0, r);
            } else if (amrex::Math::abs(sp - s(0, r)) >=
                       2.0 * amrex::Math::abs(sm - s(0, r))) {
                sp = 3.0 * s(0, r) - 2.0 * sm;
            } else if (amrex::Math::abs(sm - s(0, r)) >=
                       2.0 * amrex::Math::abs(sp - s(0, r))) {
                sm = 3.0 * s(0, r) - 2.0 * sp;
            }

            // compute Ip and Im
            Real sigmap = amrex::Math::abs(w0_arr(0, r + 1)) *
                          dtdr;  // NOTE: sigmap=0 for use_exact_base_state case
            Real sigmam = amrex::Math::abs(w0_arr(0, r)) *
                          dtdr;  // NOTE: sigmam=0 for use_exact_base_state case
            Real s6 = 6.0 * s(0, r) - 3.0 * (sm + sp);

            Real Ip = w0_arr(0, r + 1) > rel_eps_local
                          ? sp - 0.5 * sigmap *
                                     (sp - sm - (1.0 - 2.0 / 3.0 * sigmap) * s6)
                          : s(0, r);

            Real Im = w0_arr(0, r) < -rel_eps_local
                          ? sm + 0.5 * sigmam *
                                     (sp - sm + (1.0 - 2.0 / 3.0 * sigmam) * s6)
                          : s(0, r);

            // // compute sedgel and sedger
            sedgel(r + 1) = Ip + dth * force(0, r);
            sedger(r) = Im + dth * force(0, r);
        });
        Gpu::synchronize();
    } else if (ppm_type == 2) {
        ParallelFor(nr_fine, [=] AMREX_GPU_DEVICE(long r) {
            // interpolate s to radial edges, store these temporary values into sedgel

            // this will hold values at r-1, r, r+1 and r+2
            Real sedget[4];

            for (int j = 0; j < 4; j++) {
                int g = r + ng - 1 + j;

                sedget[j] = (7.0 / 12.0) * (s_ghost(g - 1) + s_ghost(g)) -
                            (1.0 / 12.0) * (s_ghost(g - 2) + s_ghost(g + 1));

                // limit sedge
                if ((sedget[j] - s_ghost(g - 1)) * (s_ghost(g) - sedget[j]) <
                    0.0) {
                    Real D2 =
                        3.0 * (s_ghost(g - 1) - 2.0 * sedget[j] + s_ghost(g));
                    Real D2L =
                        s_ghost(g - 2) - 2.0 * s_ghost(g - 1) + s_ghost(g);
                    Real D2R =
                        s_ghost(g - 1) - 2.0 * s_ghost(g) + s_ghost(g + 1);
                    Real sgn = amrex::Math::copysign(1.0, D2);
                    Real D2LIM =
                        sgn * amrex::max(amrex::min(C * sgn * D2L,
                                                    amrex::min(C * sgn * D2R,
                                                               sgn * D2)),
                                         0.0);
                    sedget[j] =
                        0.5 * (s_ghost(g - 1) + s_ghost(g)) - D2LIM / 6.0;
                }
            }

            int g = r + ng;

            // use Colella 2008 limiters
            // This is a new version of the algorithm
            // to eliminate sensitivity to roundoff.
            Real alphap = sedget[2] - s_ghost(g);
            Real alpham = sedget[1] - s_ghost(g);
            bool bigp =
                amrex::Math::abs(alphap) > 2.0 * amrex::Math::abs(alpham);
            bool bigm =
                amrex::Math::abs(alpham) > 2.0 * amrex::Math::abs(alphap);
            bool extremum = false;

            if (alpham * alphap >= 0.0) {
                extremum = true;
            } else if (bigp || bigm) {
                // Possible extremum. We look at cell centered values and face
                // centered values for a change in copysign in the differences adjacent to
                // the cell. We use the pair of differences whose minimum magnitude is
                // the largest, and thus least susceptible to sensitivity to roundoff.
                Real dafacem = sedget[1] - sedget[0];
                Real dafacep = sedget[3] - sedget[2];
                Real dabarm = s_ghost(g) - s_ghost(g - 1);
                Real dabarp = s_ghost(g + 1) - s_ghost(g);
                Real dafacemin = amrex::min(amrex::Math::abs(dafacem),
                                            amrex::Math::abs(dafacep));
                Real dabarmin = amrex::min(amrex::Math::abs(dabarm),
                                           amrex::Math::abs(dabarp));
                Real dachkm = 0.0;
                Real dachkp = 0.0;
                if (dafacemin >= dabarmin) {
                    dachkm = dafacem;
                    dachkp = dafacep;
                } else {
                    dachkm = dabarm;
                    dachkp = dabarp;
                }
                extremum = (dachkm * dachkp <= 0.0);
            }

            if (extremum) {
                Real D2 = 6.0 * (alpham + alphap);
                Real D2L = s_ghost(g - 2) - 2.0 * s_ghost(g - 1) + s_ghost(g);
                Real D2R = s_ghost(g) - 2.0 * s_ghost(g + 1) + s_ghost(g + 2);
                Real D2C = s_ghost(g - 1) - 2.0 * s_ghost(g) + s_ghost(g + 1);
                Real sgn = amrex::Math::copysign(1.0, D2);
                Real D2LIM = amrex::max(
                    amrex::min(sgn * D2, amrex::min(C * sgn * D2L,
                                                    amrex::min(C * sgn * D2R,
                                                               C * sgn * D2C))),
                    0.0);
                Real D2ABS = amrex::max(amrex::Math::abs(D2), 1.e-10);
                alpham = alpham * D2LIM / D2ABS;
                alphap = alphap * D2LIM / D2ABS;
            } else {
                if (bigp) {
                    Real sgn = amrex::Math::copysign(1.0, alpham);
                    Real amax = -alphap * alphap / (4.0 * (alpham + alphap));
                    Real delam = s_ghost(g - 1) - s_ghost(g);
                    if (sgn * amax >= sgn * delam) {
                        if (sgn * (delam - alpham) >= 1.e-10) {
                            alphap = (-2.0 * delam -
                                      2.0 * sgn *
                                          sqrt(delam * delam - delam * alpham));
                        } else {
                            alphap = -2.0 * alpham;
                        }
                    }
                }
                if (bigm) {
                    Real sgn = amrex::Math::copysign(1.0, alphap);
                    Real amax = -alpham * alpham / (4.0 * (alpham + alphap));
                    Real delap = s_ghost(g + 1) - s_ghost(g);
                    if (sgn * amax >= sgn * delap) {
                        if (sgn * (delap - alphap) >= 1.e-10) {
                            alpham = (-2.0 * delap -
                                      2.0 * sgn *
                                          sqrt(delap * delap - delap * alphap));
                        } else {
                            alpham = -2.0 * alphap;
                        }
                    }
                }
            }

            Real sm = s_ghost(g) + alpham;
            Real sp = s_ghost(g) + alphap;

            // compute Ip and Im
            Real sigmap = amrex::Math::abs(w0_arr(0, r + 1)) *
                          dtdr;  // NOTE: sigmap=0 for use_exact_base_state case
            Real sigmam = amrex::Math::abs(w0_arr(0, r)) *
                          dtdr;  // NOTE: sigmam=0 for use_exact_base_state case
            Real s6 = 6.0 * s(0, r) - 3.0 * (sm + sp);

            Real Ip = w0_arr(0, r + 1) > rel_eps_local
                          ? sp - 0.5 * sigmap *
                                     (sp - sm - (1.0 - 2.0 / 3.0 * sigmap) * s6)
                          : s(0, r);

            Real Im = w0_arr(0, r) < -rel_eps_local
                          ? sm + 0.5 * sigmam *
                                     (sp - sm + (1.0 - 2.0 / 3.0 * sigmam) * s6)
                          : s(0, r);

            // // compute sedgel and sedger
            sedgel(r + 1) = Ip + dth * force(0, r);
            sedger(r) = Im + dth * force(0, r);
        });
        Gpu::synchronize();
    }

    ParallelFor(nr_fine + 1, [=] AMREX_GPU_DEVICE(long r) {
        // Fix center and edge of star by reflecting the extrapolated state.
        // An alternate way would be to compute these values using the entire algorithm,
        // but that would require more ghost cells at several stages.
        // By symmetry arguments, this would make no difference at the center of the star
        // and the accuracy at the edge of the star is not important here
        if (r == 0) {
            sedgel(r) = sedger(r);
        } else if (r == nr_fine) {
            sedger(r) = sedgel(r);
        }

        // solve Riemann problem to get final edge state
        sedge(0, r) = w0_arr(0, r) > 0.0 ? sedgel(r) : sedger(r);
        sedge(0, r) = amrex::Math::abs(w0_arr(0, r)) < rel_eps_local
                          ? 0.5 * (sedger(r) + sedgel(r))
                          : sedge(0, r);
    });
    Gpu::synchronize();
}

void Maestro::MakeEdgeState1dPlanar(BaseState<Real>& s_state,
                                    BaseState<Real>& sedge_state,
                                    BaseState<Real>& force_state) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeState1dPlanar()", MakeEdgeState1dPlanar);

    const auto rel_eps_local = rel_eps;

    const Real dth = 0.5 * dt;
    const Real C = 1.25;
    const int cen = 0;
    const int lim = 1;
    const int flag = 2;

    const int slope_order_loc = slope_order;
    auto s = s_state.array();
    auto sedge = sedge_state.array();
    auto force = force_state.array();

    BaseState<Real> sedgel_state(base_geom.max_radial_level + 1,
                                 base_geom.nr_fine + 1);
    BaseState<Real> sedger_state(base_geom.max_radial_level + 1,
                                 base_geom.nr_fine + 1);
    auto sedgel = sedgel_state.array();
    auto sedger = sedger_state.array();

    const auto w0_arr = w0.const_array();

    for (int n = 0; n <= base_geom.finest_radial_level; ++n) {
        const int nr_lev = base_geom.nr(n);
        const Real dr_lev = base_geom.dr(n);
        const Real dtdr = dt / dr_lev;

        for (int i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i);

            // error checking to make sure that there is a 2 cell buffer at the top and bottom
            // of the domain for finer levels in planar geometry.  This can be removed if
            // blocking_factor is implemented at set > 1.
            if (ppm_type == 1 || ppm_type == 2) {
                if (base_geom.r_start_coord(n, i) == 2) {
                    Abort(
                        "make_edge_state assumes blocking_factor > 1 at lo "
                        "boundary");
                } else if (base_geom.r_end_coord(n, i) == nr_lev - 3) {
                    Abort(
                        "make_edge_state assumes blocking_factor > 1 at hi "
                        "boundary");
                }
            }

            if (ppm_type == 0) {
                // compute slopes
                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(long j) {
                    Real slope = 0.0;
                    int r = j + lo;

                    Real dxscr[3];
                    for (Real& d : dxscr) {
                        d = 0.0;
                    }

                    if (slope_order_loc == 0) {
                        slope = 0.0;
                    } else if (slope_order_loc == 2) {
                        if (r == 0) {
                            // one-sided difference
                            slope = s(n, r + 1) - s(n, r);
                        } else if (r == nr_lev - 1) {
                            // one-sided difference
                            slope = s(n, r) - s(n, r - 1);
                        } else {
                            // do standard limiting on interior cells
                            Real del = 0.5 * (s(n, r + 1) - s(n, r - 1));
                            Real dpls = 2.0 * (s(n, r + 1) - s(n, r));
                            Real dmin = 2.0 * (s(n, r) - s(n, r - 1));
                            Real slim = amrex::min(amrex::Math::abs(dpls),
                                                   amrex::Math::abs(dmin));
                            slim = dpls * dmin > 0.0 ? slim : 0.0;
                            Real sflag = amrex::Math::copysign(1.0, del);
                            slope =
                                sflag * amrex::min(slim, amrex::Math::abs(del));
                        }

                    } else if (slope_order_loc == 4) {
                        // we need to calculate dxscr(fromm) for r-1 and r+1
                        Real dxscrm = 0.0;
                        Real dxscrp = 0.0;
                        // r-1
                        int rm = r - 1;

                        if (rm == 0) {
                            // one-sided difference
                            dxscrm = s(n, rm + 1) - s(n, rm);
                        } else if (rm == nr_lev - 1) {
                            // one-sided difference
                            dxscrm = s(n, rm) - s(n, rm - 1);
                        } else if (rm > 0 && rm < nr_lev - 1) {
                            // do standard limiting to compute temporary slopes
                            dxscr[cen] = 0.5 * (s(n, rm + 1) - s(n, rm - 1));
                            Real dpls = 2.0 * (s(n, rm + 1) - s(n, rm));
                            Real dmin = 2.0 * (s(n, rm) - s(n, rm - 1));
                            dxscr[lim] = amrex::min(amrex::Math::abs(dmin),
                                                    amrex::Math::abs(dpls));
                            dxscr[lim] = dpls * dmin > 0.0 ? dxscr[lim] : 0.0;
                            dxscr[flag] =
                                amrex::Math::copysign(1.0, dxscr[cen]);
                            dxscrm = dxscr[flag] *
                                     amrex::min(dxscr[lim],
                                                amrex::Math::abs(dxscr[cen]));
                        }

                        int rp = r + 1;

                        if (rp == 0) {
                            // one-sided difference
                            dxscrp = s(n, rp + 1) - s(n, rp);
                        } else if (rp == nr_lev - 1) {
                            // one-sided difference
                            dxscrp = s(n, rp) - s(n, rp - 1);
                        } else if (rp > 0 && rp < nr_lev - 1) {
                            // do standard limiting to compute temporary slopes
                            dxscr[cen] = 0.5 * (s(n, rp + 1) - s(n, rp - 1));
                            Real dpls = 2.0 * (s(n, rp + 1) - s(n, rp));
                            Real dmin = 2.0 * (s(n, rp) - s(n, rp - 1));
                            dxscr[lim] = amrex::min(amrex::Math::abs(dmin),
                                                    amrex::Math::abs(dpls));
                            dxscr[lim] = dpls * dmin > 0.0 ? dxscr[lim] : 0.0;
                            dxscr[flag] =
                                amrex::Math::copysign(1.0, dxscr[cen]);
                            dxscrp = dxscr[flag] *
                                     amrex::min(dxscr[lim],
                                                amrex::Math::abs(dxscr[cen]));
                        }

                        // now find dxscr for r
                        if (r > 0 && r < nr_lev - 1) {
                            // do standard limiting to compute temporary slopes
                            dxscr[cen] = 0.5 * (s(n, r + 1) - s(n, r - 1));
                            Real dpls = 2.0 * (s(n, r + 1) - s(n, r));
                            Real dmin = 2.0 * (s(n, r) - s(n, r - 1));
                            dxscr[lim] = amrex::min(amrex::Math::abs(dmin),
                                                    amrex::Math::abs(dpls));
                            dxscr[lim] = dpls * dmin > 0.0 ? dxscr[lim] : 0.0;
                            dxscr[flag] =
                                amrex::Math::copysign(1.0, dxscr[cen]);
                        }

                        if (r == 0) {
                            // one-sided difference
                            slope = s(n, r + 1) - s(n, r);
                        } else if (r == nr_lev - 1) {
                            // one-sided difference
                            slope = s(n, r) - s(n, r - 1);
                        } else {
                            // fourth-order limited slopes on interior
                            Real ds = 4.0 / 3.0 * dxscr[cen] -
                                      (dxscrp + dxscrm) / 6.0;
                            slope =
                                dxscr[flag] *
                                amrex::min(amrex::Math::abs(ds), dxscr[lim]);
                        }
                    }  // which slope order

                    // compute sedgel and sedger
                    Real u = 0.5 * (w0_arr(n, r) + w0_arr(n, r + 1));
                    Real ubardth = dth * u / dr_lev;
                    sedgel(n, r + 1) =
                        s(n, r) + (0.5 - ubardth) * slope + dth * force(n, r);
                    sedger(n, r) =
                        s(n, r) - (0.5 + ubardth) * slope + dth * force(n, r);
                });
                Gpu::synchronize();

            } else if (ppm_type == 1) {
                // interpolate s to radial edges, store these temporary values into sedgel

                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(long j) {
                    int r = j + lo;

                    // calculate sm

                    // compute van Leer slopes
                    // r - 1
                    Real dsvlm = 0.0;
                    int rm = r - 1;
                    if (rm == 0) {
                        // one-sided difference
                        dsvlm = s(n, rm + 1) - s(n, rm);
                    } else if (rm == nr_lev - 1) {
                        // one-sided difference
                        dsvlm = s(n, rm) - s(n, rm - 1);
                    } else if (rm > 0 && rm < nr_lev - 1) {
                        Real del = 0.5 * (s(n, rm + 1) - s(n, rm - 1));
                        Real dmin = 2.0 * (s(n, rm) - s(n, rm - 1));
                        Real dpls = 2.0 * (s(n, rm + 1) - s(n, rm));
                        dsvlm =
                            dmin * dpls > 0.0
                                ? amrex::Math::copysign(1.0, del) *
                                      amrex::min(
                                          amrex::Math::abs(del),
                                          amrex::min(amrex::Math::abs(dmin),
                                                     amrex::Math::abs(dpls)))
                                : 0.0;
                    }

                    // r
                    Real dsvl = 0.0;
                    if (r == 0) {
                        // one-sided difference
                        dsvl = s(n, r + 1) - s(n, r);
                    } else if (r == nr_lev - 1) {
                        // one-sided difference
                        dsvl = s(n, r) - s(n, r - 1);
                    } else if (r > 0 && r < nr_lev - 1) {
                        Real del = 0.5 * (s(n, r + 1) - s(n, r - 1));
                        Real dmin = 2.0 * (s(n, r) - s(n, r - 1));
                        Real dpls = 2.0 * (s(n, r + 1) - s(n, r));
                        dsvl = dmin * dpls > 0.0
                                   ? amrex::Math::copysign(1.0, del) *
                                         amrex::min(
                                             amrex::Math::abs(del),
                                             amrex::min(amrex::Math::abs(dmin),
                                                        amrex::Math::abs(dpls)))
                                   : 0.0;
                    }

                    Real sm = 0.0;
                    if (r == 0) {
                        // 2nd order interpolation to boundary face
                        sm = s(n, r) - 0.5 * dsvl;
                    } else if (r == nr_lev) {
                        // 2nd order interpolation to boundary face
                        sm = s(n, r - 1) + 0.5 * dsvl;
                    } else {
                        // 4th order interpolation of s to radial faces
                        sm = 0.5 * (s(n, r) + s(n, r - 1)) -
                             (dsvl - dsvlm) / 6.0;
                        // make sure sm lies in between adjacent cell-centered values
                        sm = amrex::max(sm, amrex::min(s(n, r), s(n, r - 1)));
                        sm = amrex::min(sm, amrex::max(s(n, r), s(n, r - 1)));
                    }

                    // calculate sp
                    // compute van Leer slopes
                    // r + 1 - 1
                    dsvlm = dsvl;

                    // r + 1
                    int rp = r + 1;
                    dsvl = 0.0;
                    if (rp == 0) {
                        // one-sided difference
                        dsvl = s(n, rp + 1) - s(n, rp);
                    } else if (rp == nr_lev - 1) {
                        // one-sided difference
                        dsvl = s(n, rp) - s(n, rp - 1);
                    } else if (rp > 0 && rp < nr_lev - 1) {
                        Real del = 0.5 * (s(n, rp + 1) - s(n, rp - 1));
                        Real dmin = 2.0 * (s(n, rp) - s(n, rp - 1));
                        Real dpls = 2.0 * (s(n, rp + 1) - s(n, rp));
                        dsvl = dmin * dpls > 0.0
                                   ? amrex::Math::copysign(1.0, del) *
                                         amrex::min(
                                             amrex::Math::abs(del),
                                             amrex::min(amrex::Math::abs(dmin),
                                                        amrex::Math::abs(dpls)))
                                   : 0.0;
                    }

                    Real sp = 0.0;
                    if (rp == 0) {
                        // 2nd order interpolation to boundary face
                        sp = s(n, rp) - 0.5 * dsvl;
                    } else if (rp == nr_lev) {
                        // 2nd order interpolation to boundary face
                        sp = s(n, rp - 1) + 0.5 * dsvl;
                    } else {
                        // 4th order interpolation of s to radial faces
                        sp = 0.5 * (s(n, rp) + s(n, rp - 1)) -
                             (dsvl - dsvlm) / 6.0;
                        // make sure sedgel lies in between adjacent cell-centered values
                        sp = amrex::max(sp, amrex::min(s(n, rp), s(n, rp - 1)));
                        sp = amrex::min(sp, amrex::max(s(n, rp), s(n, rp - 1)));
                    }

                    // modify using quadratic limiters
                    if ((sp - s(n, r)) * (s(n, r) - sm) <= 0.0) {
                        sp = s(n, r);
                        sm = s(n, r);
                    } else if (amrex::Math::abs(sp - s(n, r)) >=
                               2.0 * amrex::Math::abs(sm - s(n, r))) {
                        sp = 3.0 * s(n, r) - 2.0 * sm;
                    } else if (amrex::Math::abs(sm - s(n, r)) >=
                               2.0 * amrex::Math::abs(sp - s(n, r))) {
                        sm = 3.0 * s(n, r) - 2.0 * sp;
                    }

                    // compute Ip and Im
                    Real sigmap = amrex::Math::abs(w0_arr(n, r + 1)) * dtdr;
                    Real sigmam = amrex::Math::abs(w0_arr(n, r)) * dtdr;
                    Real s6 = 6.0 * s(n, r) - 3.0 * (sm + sp);
                    Real Ip = 0.0;
                    Real Im = 0.0;
                    if (w0_arr(n, r + 1) > rel_eps_local) {
                        Ip = sp -
                             (sigmap / 2.0) *
                                 (sp - sm - (1.0 - 2.0 / 3.0 * sigmap) * s6);
                    } else {
                        Ip = s(n, r);
                    }
                    if (w0_arr(n, r) < -rel_eps_local) {
                        Im = sm +
                             (sigmam / 2.0) *
                                 (sp - sm + (1.0 - 2.0 / 3.0 * sigmam) * s6);
                    } else {
                        Im = s(n, r);
                    }

                    // compute sedgel and sedger
                    sedgel(n, r + 1) = Ip + dth * force(n, r);
                    sedger(n, r) = Im + dth * force(n, r);
                });
                Gpu::synchronize();

            } else if (ppm_type == 2) {
                // interpolate s to radial edges

                // need a vector to store intermediate values
                BaseState<Real> sedget_s(base_geom.nr_fine + 1);
                auto sedget = sedget_s.array();

                ParallelFor(hi - lo + 2, [=] AMREX_GPU_DEVICE(long j) {
                    int r = j + lo;

                    // left side
                    Real dsvl = 0.0;
                    if (r - 1 == 0) {
                        // one-sided difference
                        dsvl = s(n, r) - s(n, r - 1);
                    } else if (r - 1 == nr_lev - 1) {
                        // one-sided difference
                        dsvl = s(n, r - 1) - s(n, r - 2);
                    } else if (r - 1 > 0 && r - 1 < nr_lev - 1) {
                        // centered difference
                        dsvl = 0.5 * (s(n, r) - s(n, r - 2));
                    }

                    // right side
                    Real dsvr = 0.0;
                    if (r == 0) {
                        // one-sided difference
                        dsvr = s(n, r + 1) - s(n, r);
                    } else if (r == nr_lev - 1) {
                        // one-sided difference
                        dsvr = s(n, r) - s(n, r - 1);
                    } else if (r > 0 && r < nr_lev - 1) {
                        // centered difference
                        dsvr = 0.5 * (s(n, r + 1) - s(n, r - 1));
                    }

                    if (r == 0) {
                        // 2nd order interpolation to boundary face
                        sedget(r) = s(n, r) - 0.5 * dsvr;
                    } else if (r == nr_lev) {
                        // 2nd order interpolation to boundary face
                        sedget(r) = s(n, r - 1) + 0.5 * dsvr;
                    } else if (r > 0 && r < nr_lev) {
                        // 4th order interpolation of s to radial faces
                        sedget(r) =
                            0.5 * (s(n, r) + s(n, r - 1)) - (dsvr - dsvl) / 6.0;
                        if (r >= 2 && r <= nr_lev - 2) {
                            // limit sedge
                            if ((sedget(r) - s(n, r - 1)) *
                                    (s(n, r) - sedget(r)) <
                                0.0) {
                                Real D2 = 3.0 * (s(n, r - 1) - 2.0 * sedget(r) +
                                                 s(n, r));
                                Real D2L =
                                    s(n, r - 2) - 2.0 * s(n, r - 1) + s(n, r);
                                Real D2R =
                                    s(n, r - 1) - 2.0 * s(n, r) + s(n, r + 1);
                                Real sgn = amrex::Math::copysign(1.0, D2);
                                Real D2LIM =
                                    sgn *
                                    amrex::max(
                                        amrex::min(C * sgn * D2L,
                                                   amrex::min(C * sgn * D2R,
                                                              sgn * D2)),
                                        0.0);
                                sedget(r) =
                                    0.5 * (s(n, r - 1) + s(n, r)) - D2LIM / 6.0;
                            }
                        }
                    }
                });
                Gpu::synchronize();

                ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(long j) {
                    int r = j + lo;

                    // use Colella 2008 limiters
                    // This is a new version of the algorithm
                    // to eliminate sensitivity to roundoff.
                    Real sm = 0.0;
                    Real sp = 0.0;

                    if (r >= 2 && r <= nr_lev - 3) {
                        Real alphap = sedget(r + 1) - s(n, r);
                        Real alpham = sedget(r) - s(n, r);
                        bool bigp = amrex::Math::abs(alphap) >
                                    2.0 * amrex::Math::abs(alpham);
                        bool bigm = amrex::Math::abs(alpham) >
                                    2.0 * amrex::Math::abs(alphap);
                        bool extremum = false;

                        if (alpham * alphap >= 0.0) {
                            extremum = true;
                        } else if (bigp || bigm) {
                            // Possible extremum. We look at cell centered values and face
                            // centered values for a change in copysign in the differences adjacent to
                            // the cell. We use the pair of differences whose minimum magnitude is
                            // the largest, and thus least susceptible to sensitivity to roundoff.
                            Real dafacem = sedget(r) - sedget(r - 1);
                            Real dafacep = sedget(r + 2) - sedget(r + 1);
                            Real dabarm = s(n, r) - s(n, r - 1);
                            Real dabarp = s(n, r + 1) - s(n, r);
                            Real dafacemin =
                                amrex::min(amrex::Math::abs(dafacem),
                                           amrex::Math::abs(dafacep));
                            Real dabarmin =
                                amrex::min(amrex::Math::abs(dabarm),
                                           amrex::Math::abs(dabarp));
                            Real dachkm = 0.0;
                            Real dachkp = 0.0;
                            if (dafacemin >= dabarmin) {
                                dachkm = dafacem;
                                dachkp = dafacep;
                            } else {
                                dachkm = dabarm;
                                dachkp = dabarp;
                            }
                            extremum = (dachkm * dachkp <= 0.0);
                        }

                        if (extremum) {
                            Real D2 = 6.0 * (alpham + alphap);
                            Real D2L =
                                s(n, r - 2) - 2.0 * s(n, r - 1) + s(n, r);
                            Real D2R =
                                s(n, r) - 2.0 * s(n, r + 1) + s(n, r + 2);
                            Real D2C =
                                s(n, r - 1) - 2.0 * s(n, r) + s(n, r + 1);
                            Real sgn = amrex::Math::copysign(1.0, D2);
                            Real D2LIM = amrex::max(
                                amrex::min(
                                    sgn * D2,
                                    amrex::min(C * sgn * D2L,
                                               amrex::min(C * sgn * D2R,
                                                          C * sgn * D2C))),
                                0.0);
                            Real D2ABS =
                                amrex::max(amrex::Math::abs(D2), 1.e-10);
                            alpham = alpham * D2LIM / D2ABS;
                            alphap = alphap * D2LIM / D2ABS;
                        } else {
                            if (bigp) {
                                Real sgn = amrex::Math::copysign(1.0, alpham);
                                Real amax = -alphap * alphap /
                                            (4.0 * (alpham + alphap));
                                Real delam = s(n, r - 1) - s(n, r);
                                if (sgn * amax >= sgn * delam) {
                                    if (sgn * (delam - alpham) >= 1.e-10) {
                                        alphap = (-2.0 * delam -
                                                  2.0 * sgn *
                                                      sqrt(delam * delam -
                                                           delam * alpham));
                                    } else {
                                        alphap = -2.0 * alpham;
                                    }
                                }
                            }
                            if (bigm) {
                                Real sgn = amrex::Math::copysign(1.0, alphap);
                                Real amax = -alpham * alpham /
                                            (4.0 * (alpham + alphap));
                                Real delap = s(n, r + 1) - s(n, r);
                                if (sgn * amax >= sgn * delap) {
                                    if (sgn * (delap - alphap) >= 1.e-10) {
                                        alpham = (-2.0 * delap -
                                                  2.0 * sgn *
                                                      sqrt(delap * delap -
                                                           delap * alphap));
                                    } else {
                                        alpham = -2.0 * alphap;
                                    }
                                }
                            }
                        }

                        sm = s(n, r) + alpham;
                        sp = s(n, r) + alphap;

                    } else {
                        sp = sedget(r + 1);
                        sm = sedget(r);
                    }  // test (r >= 2 && r <= nr-3)

                    // compute Ip and Im
                    Real sigmap = amrex::Math::abs(w0_arr(n, r + 1)) * dtdr;
                    Real sigmam = amrex::Math::abs(w0_arr(n, r)) * dtdr;
                    Real s6 = 6.0 * s(n, r) - 3.0 * (sm + sp);
                    Real Ip = 0.0;
                    Real Im = 0.0;
                    if (w0_arr(n, r + 1) > rel_eps_local) {
                        Ip = sp -
                             0.5 * sigmap *
                                 (sp - sm - (1.0 - 2.0 / 3.0 * sigmap) * s6);
                    } else {
                        Ip = s(n, r);
                    }
                    if (w0_arr(n, r) < -rel_eps_local) {
                        Im = sm +
                             0.5 * sigmam *
                                 (sp - sm + (1.0 - 2.0 / 3.0 * sigmam) * s6);
                    } else {
                        Im = s(n, r);
                    }

                    // compute sedgel and sedger
                    sedgel(n, r + 1) = Ip + dth * force(n, r);
                    sedger(n, r) = Im + dth * force(n, r);
                });
                Gpu::synchronize();
            }
        }
    }

    for (int n = 0; n <= base_geom.finest_radial_level; ++n) {
        for (int i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i);

            // sync up edge states at coarse-fine interface

            // if we are not at the finest level, copy in the sedger and sedgel states
            // from the next finer level at the c-f interface
            if (n < base_geom.finest_radial_level) {
                sedger(n, base_geom.r_start_coord(n + 1, i) / 2) =
                    sedger(n + 1, base_geom.r_start_coord(n + 1, i));
                sedgel(n, (base_geom.r_end_coord(n + 1, i) + 1) / 2) =
                    sedgel(n + 1, base_geom.r_end_coord(n + 1, i) + 1);
            }

            // if we are not at the coarsest level, copy in the sedgel and sedger states
            // from the next coarser level at the c-f interface
            if (n > 0) {
                sedgel(n, lo) = sedgel(n - 1, lo / 2);
                sedger(n, hi + 1) = sedger(n - 1, (hi + 1) / 2);
            }
        }
    }

    for (int n = 0; n <= base_geom.finest_radial_level; ++n) {
        const auto nr_lev = base_geom.nr(n);

        for (int i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i);

            // solve Riemann problem to get final edge state
            ParallelFor(hi - lo + 2, [=] AMREX_GPU_DEVICE(long j) {
                int r = j + lo;

                if (r == 0) {
                    // pick interior state at lo domain boundary
                    sedge(n, r) = sedger(n, r);
                } else if (r == nr_lev) {
                    // pick interior state at hi domain boundary
                    sedge(n, r) = sedgel(n, r);
                } else {
                    // upwind
                    sedge(n, r) =
                        w0_arr(n, r) > 0.0 ? sedgel(n, r) : sedger(n, r);
                    sedge(n, r) = amrex::Math::abs(w0_arr(n, r)) < rel_eps_local
                                      ? 0.5 * (sedger(n, r) + sedgel(n, r))
                                      : sedge(n, r);
                }
            });
            Gpu::synchronize();
        }  // loop over disjointchunks
    }      // loop over levels
}
