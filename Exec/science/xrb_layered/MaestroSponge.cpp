
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::SpongeInit(const BaseState<Real>& rho0_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SpongeInit()", SpongeInit);

    // There is a bottom sponge and a top sponge for xrb_mixed problem.
    //
    // The start of the bottom sponge is r_sp = botsponge_lo_r, and
    // the end of the bottom sponge is r_tp = botsponge_hi_r.
    //
    // The start of the top sponge is r_sp_outer = topsponge_lo_r, and
    // the end of the top sponge is r_tp_outer = topsponge_hi_r.

    const auto rho0 = rho0_s.const_array();

    Real prob_lo_r = geom[0].ProbLo(AMREX_SPACEDIM - 1);

    // Top sponge
    Real r_top = prob_lo_r + base_geom.r_edge_loc(0, 1);
    r_sp_outer = r_top;

    // we center the sponge on the anelastic cutoff
    sponge_start_density = sponge_start_factor * anelastic_cutoff_density;

    // set topsponge_lo_r = r_sp_outer;
    for (auto r = 0; r <= base_geom.r_end_coord(0, 1); ++r) {
        if (rho0(0, r) < sponge_start_density) {
            r_sp_outer = prob_lo_r + (static_cast<Real>(r) + 0.5) * base_geom.dr(0);
            break;
        }
    }

    // set topsponge_hi_r = r_tp_outer;
    r_tp_outer = r_top;
    for (auto r = 0; r <= base_geom.r_end_coord(0, 1); ++r) {
        if (rho0(0, r) < anelastic_cutoff_density) {
            r_tp_outer = prob_lo_r + (static_cast<Real>(r) + 0.5) * base_geom.dr(0);
            break;
        }
    }

    // set botsponge_lo_r = r_sp;
    for (auto r = 0; r <= base_geom.r_end_coord(0, 1); ++r) {
        if (rho0(0, r) < 6.0e7) {
            r_sp = prob_lo_r + (Real(r) + 0.5) * base_geom.dr(0);
            break;
        }
    }

    // set botsponge_hi_r = r_tp;
    for (auto r = 0; r <= base_geom.r_end_coord(0, 1); ++r) {
        if (rho0(0, r) < 5.0e7) {
            r_tp = prob_lo_r + (Real(r) + 0.5) * base_geom.dr(0);
            break;
        }
    }

    if (maestro_verbose >= 1) {
        if (xrb_use_bottom_sponge) {
            Print() << "inner sponge: r_sp      , r_tp      : " << r_sp << ", "
                    << r_tp << std::endl;
        }

        Print() << "outer sponge: r_sp_outer, r_tp_outer: " << r_sp_outer
                << ", " << r_tp_outer << std::endl;
    }
}

void Maestro::MakeSponge(Vector<MultiFab>& sponge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeSponge()", MakeSponge);

    const Real botsponge_lo_r = r_sp;
    const Real botsponge_hi_r = r_tp;
    const Real topsponge_lo_r = r_sp_outer;
    const Real topsponge_hi_r = r_tp_outer;
    const Real sponge_min_loc = sponge_min;

    if (AMREX_SPACEDIM != 2) {
        Abort("ERROR: sponge only supported for 2d in xrb_mixed");
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sponge[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();

            const Array4<Real> sponge_arr = sponge[lev].array(mfi);
            const auto prob_lo = geom[lev].ProbLoArray();

            Real smdamp = 1.0;

            const auto lo = tileBox.loVect3d()[AMREX_SPACEDIM - 1];
            const auto hi = tileBox.hiVect3d()[AMREX_SPACEDIM - 1];

            ParallelFor(tileBox, [=] (int i, int j, int k)
                        { sponge_arr(i, j, k) = 1.0; });

            if (xrb_use_bottom_sponge) {

                // do both the top and bottom sponges

                // look over the vertical direction
                for (int j = lo; j <= hi; ++j) {
                    Real y = prob_lo[1] + (static_cast<Real>(j) + 0.5) * dx[1];

                    if (y <= botsponge_lo_r) {
                        ParallelFor(tileBox, [=](int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = sponge_min_loc;
                            }
                        });

                    } else if (y <= botsponge_hi_r) {
                        smdamp = -0.5 * (1.0 - sponge_min_loc) *
                            std::cos(M_PI * (y - botsponge_lo_r) /
                                     (botsponge_hi_r - botsponge_lo_r)) +
                            0.5 * (1.0 + sponge_min_loc);

                        ParallelFor(tileBox, [=](int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = smdamp;
                            }
                        });

                    } else if (y <= topsponge_lo_r) {
                        ParallelFor(tileBox, [=] (int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = 1.0;
                            }
                        });
                    } else if (y <= topsponge_hi_r) {
                        smdamp = 0.5 * (1.0 - sponge_min_loc) *
                            std::cos(M_PI * (y - topsponge_lo_r) /
                                     (topsponge_hi_r - topsponge_lo_r)) +
                            0.5 * (1.0 + sponge_min_loc);

                        ParallelFor(tileBox, [=] (int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = smdamp;
                            }
                        });
                    } else {
                        ParallelFor(tileBox, [=] (int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = sponge_min_loc;
                            }
                        });
                    }

                }

            } else {

                // just top sponge

                for (int j = lo; j <= hi; ++j) {
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

                    if (y <= topsponge_lo_r) {
                        ParallelFor(tileBox, [=] (int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = 1.0;
                            }
                        });
                    } else if (y <= topsponge_hi_r) {
                        smdamp = 0.5 * (1.0 - sponge_min_loc) *
                            std::cos(M_PI * (y - topsponge_lo_r) /
                                     (topsponge_hi_r - topsponge_lo_r)) +
                            0.5 * (1.0 + sponge_min_loc);

                        ParallelFor(tileBox, [=] (int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = smdamp;
                            }
                        });
                    } else {
                        ParallelFor(tileBox, [=] (int ii, int jj, int kk) {
                            if (jj == j) {
                                sponge_arr(ii, jj, kk) = sponge_min_loc;
                            }
                        });
                    }
                }
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(sponge, 0, 1);
}
