
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::SpongeInit(const BaseState<Real>& rho0_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SpongeInit()", SpongeInit);

    // The sponge has a HALF * ( 1 - cos( (r - r_sp)/L)) profile, where
    // the width, L, is r_tp - r_sp.
    //
    // The center of the sponge, r_md, is set to the radius where r =
    // sponge_center_density
    //
    // The start of the sponge, r_sp, (moving outward from the center)
    // is the radius where r = sponge_start_factor * sponge_center_density
    //
    // The top of the sponge is then 2 * r_md - r_tp

    const auto rho0 = rho0_s.const_array();

    Real prob_lo_r = spherical ? 0.0 : geom[0].ProbLo(AMREX_SPACEDIM - 1);

    Real r_top = prob_lo_r;
    if (!use_exact_base_state) {
        r_top += Real(base_geom.r_end_coord(0, 1) + 1) * base_geom.dr(0);
    } else {
        r_top += base_geom.r_edge_loc(0, base_geom.nr_fine);
    }
    r_sp = r_top;

    sponge_start_density = sponge_start_factor * sponge_center_density;

    if (!use_exact_base_state) {
        // set r_sp;
        for (auto r = 0; r <= base_geom.r_end_coord(0, 1); ++r) {
            if (rho0(0, r) < sponge_start_density) {
                r_sp = prob_lo_r + (Real(r) + 0.5) * base_geom.dr(0);
                break;
            }
        }

        // set r_md
        r_md = r_top;
        for (auto r = 0; r <= base_geom.r_end_coord(0, 1); ++r) {
            if (rho0(0, r) < sponge_center_density) {
                r_md = prob_lo_r + (Real(r) + 0.5) * base_geom.dr(0);
                break;
            }
        }
    } else {
        // set r_sp;
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            if (rho0(0, r) < sponge_start_density) {
                r_sp = prob_lo_r + base_geom.r_cc_loc(0, r);
                break;
            }
        }

        // set r_md
        r_md = r_top;
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            if (rho0(0, r) < sponge_center_density) {
                r_md = prob_lo_r + base_geom.r_cc_loc(0, r);
                break;
            }
        }
    }

    // set r_tp
    r_tp = 2.0 * r_md - r_sp;

    // outer sponge parameters used for spherical problems
    if (spherical) {
        r_sp_outer = r_tp;
        if (!use_exact_base_state) {
            r_tp_outer = r_sp_outer + 4.0 * drdxfac * base_geom.dr(0);
        } else {
            r_tp_outer = r_sp_outer +
                         4.0 * 2.0 / std::sqrt(3) * base_geom.r_cc_loc(0, 0);
        }
    }

    if (maestro_verbose >= 1) {
        Print() << "inner sponge: r_sp      , r_tp      : " << r_sp << ", "
                << r_tp << std::endl;
        if (spherical) {
            Print() << "outer sponge: r_sp_outer, r_tp_outer: " << r_sp_outer
                    << ", " << r_tp_outer << std::endl;
        }
    }
}

void Maestro::MakeSponge(Vector<MultiFab>& sponge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeSponge()", MakeSponge);

    const Real dt_loc = dt;
    const Real r_sp_loc = r_sp;
    const Real r_tp_loc = r_tp;
    const Real sponge_kappa_loc = sponge_kappa;

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

            if (!spherical) {
                const auto lo = tileBox.loVect3d()[AMREX_SPACEDIM - 1];
                const auto hi = tileBox.hiVect3d()[AMREX_SPACEDIM - 1];
                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    sponge_arr(i, j, k) = 1.0;
                });

                for (int n = 0; n <= hi - lo; ++n) {
#if (AMREX_SPACEDIM == 2)
                    int j = lo + n;
                    Real r = prob_lo[1] + (Real(j) + 0.5) * dx[1];
#else
                    int k = lo + n;
                    Real r = prob_lo[2] + (Real(k) + 0.5) * dx[2];
#endif
                    if (r >= r_sp_loc) {
                        if (r < r_tp_loc) {
                            Real smdamp =
                                0.5 * (1.0 - std::cos(M_PI * (r - r_sp_loc) /
                                                      (r_tp_loc - r_sp_loc)));
                            int smdamp_idx = lo + n;

                            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(
                                                     int ii, int jj, int kk) {
#if (AMREX_SPACEDIM == 2)
                                if (jj == smdamp_idx)
#else
                                if (kk == smdamp_idx)
#endif
                                    sponge_arr(ii, jj, kk) =
                                        1.0 / (1.0 + dt_loc * smdamp *
                                                         sponge_kappa_loc);
                            });
                        }
                    }
                }
            } else {
#if (AMREX_SPACEDIM == 3)
                const Real r_sp_outer_loc = r_sp_outer;
                const Real r_tp_outer_loc = r_tp_outer;
                const auto& center_p = center;

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    sponge_arr(i, j, k) = 1.0;

                    Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                    Real r = std::sqrt(x * x + y * y + z * z);

                    // Inner sponge: damps velocities at edge of star
                    if (r >= r_sp_loc) {
                        Real smdamp = 1.0;
                        if (r < r_tp_loc) {
                            smdamp =
                                0.5 * (1.0 - std::cos(M_PI * (r - r_sp_loc) /
                                                      (r_tp_loc - r_sp_loc)));
                        }
                        sponge_arr(i, j, k) =
                            1.0 / (1.0 + dt_loc * smdamp * sponge_kappa_loc);
                    }

                    // Outer sponge: damps velocities in the corners of the domain
                    if (r >= r_sp_outer_loc) {
                        Real smdamp = 1.0;
                        if (r < r_tp_outer_loc) {
                            smdamp =
                                0.5 *
                                (1.0 -
                                 std::cos(M_PI * (r - r_sp_outer_loc) /
                                          (r_tp_outer_loc - r_sp_outer_loc)));
                        }
                        sponge_arr(i, j, k) /=
                            1.0 + dt_loc * smdamp * 10.0 * sponge_kappa_loc;
                    }
                });
#endif
            }
        }
    }

    // average fine data onto coarser cells
    AverageDown(sponge, 0, 1);
}
