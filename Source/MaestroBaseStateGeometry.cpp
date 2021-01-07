#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 3)
void Maestro::InitBaseStateMapSphr(
    const int lev, const MFIter& mfi,
    const GpuArray<Real, AMREX_SPACEDIM> dx_fine,
    const GpuArray<Real, AMREX_SPACEDIM> dx_lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseStateMapSphr()", InitBaseStateMapSphr);

    if (!spherical) {
        Abort("init_base_state_map_sphr() does not work for planar");
    }

    if (use_exact_base_state) {
        const Box& tilebox = mfi.tilebox();
        const Array4<int> cc_to_r = cell_cc_to_r[lev].array(mfi);

        const auto probLo = geom[0].ProbLoArray();
        const auto center_p = center;

        ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real x = probLo[0] + (static_cast<Real>(i) + 0.5) * dx_lev[0] -
                     center_p[0];
            Real y = probLo[1] + (static_cast<Real>(j) + 0.5) * dx_lev[1] -
                     center_p[1];
            Real z = probLo[2] + (static_cast<Real>(k) + 0.5) * dx_lev[2] -
                     center_p[2];

            Real index =
                (x * x + y * y + z * z) / (2.0 * dx_fine[0] * dx_fine[0]) -
                0.375;
            cc_to_r(i, j, k) = (int)amrex::Math::round(index);
        });
    }
}
#endif

void Maestro::ComputeCutoffCoords(const BaseState<Real>& rho0_state) {
    // compute the coordinates of the anelastic cutoff
    bool found = false;
    int which_lev = 0;
    const auto rho0 = rho0_state.const_array();

    // find the finest level containing the anelastic cutoff density,
    // and set the anelastic cutoff coord for this level
    for (auto n = base_geom.finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            if (!found) {
                int lo = base_geom.r_start_coord(n, i);
                int hi = base_geom.r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) <= anelastic_cutoff_density) {
                        base_geom.anelastic_cutoff_density_coord(n) = r;
                        which_lev = n;
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    // if the anelastic cutoff density was not found anywhere, { set
    // it to above the top of the domain on the finest level
    if (!found) {
        which_lev = base_geom.finest_radial_level;
        base_geom.anelastic_cutoff_density_coord(
            base_geom.finest_radial_level) =
            base_geom.nr(base_geom.finest_radial_level);
    }

    // set the anelastic cutoff coordinate on the finer levels
    for (auto n = which_lev + 1; n <= base_geom.finest_radial_level; ++n) {
        base_geom.anelastic_cutoff_density_coord(n) =
            2 * base_geom.anelastic_cutoff_density_coord(n - 1) + 1;
    }

    // set the anelastic cutoff coordinate on the coarser levels
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (base_geom.anelastic_cutoff_density_coord(n + 1) % 2 == 0) {
            base_geom.anelastic_cutoff_density_coord(n) =
                base_geom.anelastic_cutoff_density_coord(n + 1) / 2;
        } else {
            base_geom.anelastic_cutoff_density_coord(n) =
                base_geom.anelastic_cutoff_density_coord(n + 1) / 2 + 1;
        }
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // compute the coordinates of the base cutoff density
    found = false;

    // find the finest level containing the base cutoff density,
    // and set the base cutoff coord for this level
    for (auto n = base_geom.finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            if (!found) {
                int lo = base_geom.r_start_coord(n, i);
                int hi = base_geom.r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) <= base_cutoff_density) {
                        base_geom.base_cutoff_density_coord(n) = r;
                        which_lev = n;
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    // if the base cutoff density was not found anywhere, { set
    // it to above the top of the domain on the finest level
    if (!found) {
        which_lev = base_geom.finest_radial_level;
        base_geom.base_cutoff_density_coord(base_geom.finest_radial_level) =
            base_geom.nr(base_geom.finest_radial_level);
    }

    // set the base cutoff coordinate on the finer levels
    // do n=which_lev+1,base_geom.finest_radial_level
    for (auto n = which_lev + 1; n <= base_geom.finest_radial_level; ++n) {
        base_geom.base_cutoff_density_coord(n) =
            2 * base_geom.base_cutoff_density_coord(n - 1) + 1;
    }

    // set the base cutoff coordinate on the coarser levels
    // do n=which_lev-1,0,-1
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (base_geom.base_cutoff_density_coord(n + 1) % 2 == 0) {
            base_geom.base_cutoff_density_coord(n) =
                base_geom.base_cutoff_density_coord(n + 1) / 2;
        } else {
            base_geom.base_cutoff_density_coord(n) =
                base_geom.base_cutoff_density_coord(n + 1) / 2 + 1;
        }
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // compute the coordinates of the burning cutoff density
    found = false;

    // find the finest level containing the burning cutoff density,
    // and set the burning cutoff coord for this level
    for (auto n = base_geom.finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            if (!found) {
                int lo = base_geom.r_start_coord(n, i);
                int hi = base_geom.r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) <= burning_cutoff_density_lo) {
                        base_geom.burning_cutoff_density_lo_coord(n) = r;
                        which_lev = n;
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    // if the burning cutoff density was not found anywhere, { set
    // it to above the top of the domain on the finest level
    if (!found) {
        which_lev = base_geom.finest_radial_level;
        base_geom.burning_cutoff_density_lo_coord(
            base_geom.finest_radial_level) =
            base_geom.nr(base_geom.finest_radial_level);
    }

    // set the burning cutoff coordinate on the finer levels
    // do n=which_lev+1,base_geom.finest_radial_level
    for (auto n = which_lev + 1; n <= base_geom.finest_radial_level; ++n) {
        base_geom.burning_cutoff_density_lo_coord(n) =
            2 * base_geom.burning_cutoff_density_lo_coord(n - 1) + 1;
    }

    // set the burning cutoff coordinate on the coarser levels
    // do n=which_lev-1,0,-1
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (base_geom.burning_cutoff_density_lo_coord(n + 1) % 2 == 0) {
            base_geom.burning_cutoff_density_lo_coord(n) =
                base_geom.burning_cutoff_density_lo_coord(n + 1) / 2;
        } else {
            base_geom.burning_cutoff_density_lo_coord(n) =
                base_geom.burning_cutoff_density_lo_coord(n + 1) / 2 + 1;
        }
    }

    // compute the coordinates of the burning cutoff density upper limit
    found = false;

    // find the finest level containing the burning cutoff density,
    // and set the burning cutoff coord for this level
    for (auto n = base_geom.finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            if (!found) {
                int lo = base_geom.r_start_coord(n, i);
                int hi = base_geom.r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) >= burning_cutoff_density_hi) {
                        base_geom.burning_cutoff_density_hi_coord(n) = r;
                        which_lev = n;
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    // if the burning cutoff density was not found anywhere, { set
    // it to above the bottom of the domain
    if (!found) {
        which_lev = base_geom.finest_radial_level;
        base_geom.burning_cutoff_density_hi_coord(
            base_geom.finest_radial_level) = 0;
    }

    // set the burning cutoff coordinate on the finer levels
    for (auto n = which_lev + 1; n <= base_geom.finest_radial_level; ++n) {
        base_geom.burning_cutoff_density_hi_coord(n) = 0;
    }

    // set the burning cutoff coordinate on the coarser levels
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (base_geom.burning_cutoff_density_hi_coord(n + 1) % 2 == 0) {
            base_geom.burning_cutoff_density_hi_coord(n) = 0;
        } else {
            base_geom.burning_cutoff_density_hi_coord(n) = 0;
        }
    }
}

void Maestro::RestrictBase(BaseState<Real>& s0, const bool is_cell_centered) {
    RestrictBase(s0.array(), is_cell_centered);
}

void Maestro::RestrictBase(const BaseStateArray<Real>& s0,
                           const bool is_cell_centered) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RestrictBase()", RestrictBase);

    for (int n = base_geom.finest_radial_level; n >= 1; --n) {
        for (int i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            if (is_cell_centered) {
                // for level n, make the coarser cells underneath simply the average of the fine
                for (auto j = base_geom.r_start_coord(n, i);
                     j < base_geom.r_end_coord(n, i); j += 2) {
                    s0(n - 1, j / 2) = 0.5 * (s0(n, j) + s0(n, j + 1));
                }
            } else {
                // for level n, make the coarse edge underneath equal to the fine edge value
                for (auto j = base_geom.r_start_coord(n, i);
                     j <= base_geom.r_end_coord(n, i) + 1; j += 2) {
                    s0(n - 1, j / 2) = s0(n, j);
                }
            }
        }
    }
}

void Maestro::FillGhostBase(BaseState<Real>& s0, const bool is_cell_centered) {
    FillGhostBase(s0.array(), is_cell_centered);
}

void Maestro::FillGhostBase(const BaseStateArray<Real>& s0,
                            const bool is_cell_centered) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillGhostBase()", FillGhostBase);

    for (int n = base_geom.finest_radial_level; n >= 1; --n) {
        for (int i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i);

            if (is_cell_centered) {
                if (lo != 0) {
                    int r_crse = lo / 2 - 1;
                    Real del =
                        0.5 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse - 1));
                    Real dpls =
                        2.0 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse));
                    Real dmin =
                        2.0 * (s0(n - 1, r_crse) - s0(n - 1, r_crse - 1));
                    Real slim = amrex::min(amrex::Math::abs(dpls),
                                           amrex::Math::abs(dmin));
                    slim = dpls * dmin > 0.0 ? slim : 0.0;
                    Real slope = amrex::Math::copysign(1.0, del) *
                                 amrex::min(slim, amrex::Math::abs(del));
                    s0(n, lo - 1) = s0(n - 1, r_crse) + 0.25 * slope;
                    s0(n, lo - 2) = s0(n - 1, r_crse) - 0.25 * slope;

                    r_crse = lo / 2 - 2;
                    del = 0.5 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse - 1));
                    dpls = 2.0 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse));
                    dmin = 2.0 * (s0(n - 1, r_crse) - s0(n - 1, r_crse - 1));
                    slim = amrex::min(amrex::Math::abs(dpls),
                                      amrex::Math::abs(dmin));
                    slim = dpls * dmin > 0.0 ? slim : 0.0;
                    slope = amrex::Math::copysign(1.0, del) *
                            amrex::min(slim, amrex::Math::abs(del));
                    s0(n, lo - 3) = s0(n - 1, r_crse) + 0.25 * slope;
                    s0(n, lo - 4) = s0(n - 1, r_crse) - 0.25 * slope;
                }

                if (hi != base_geom.nr(n) - 1) {
                    int r_crse = (hi + 1) / 2;
                    Real del =
                        0.5 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse - 1));
                    Real dpls =
                        2.0 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse));
                    Real dmin =
                        2.0 * (s0(n - 1, r_crse) - s0(n - 1, r_crse - 1));
                    Real slim = amrex::min(amrex::Math::abs(dpls),
                                           amrex::Math::abs(dmin));
                    slim = dpls * dmin > 0.0 ? slim : 0.0;
                    Real slope = amrex::Math::copysign(1.0, del) *
                                 amrex::min(slim, amrex::Math::abs(del));
                    s0(n, hi + 1) = s0(n - 1, r_crse) - 0.25 * slope;
                    s0(n, hi + 2) = s0(n - 1, r_crse) + 0.25 * slope;

                    r_crse = (hi + 1) / 2 + 1;
                    del = 0.5 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse - 1));
                    dpls = 2.0 * (s0(n - 1, r_crse + 1) - s0(n - 1, r_crse));
                    dmin = 2.0 * (s0(n - 1, r_crse) - s0(n - 1, r_crse - 1));
                    slim = amrex::min(amrex::Math::abs(dpls),
                                      amrex::Math::abs(dmin));
                    slim = dpls * dmin > 0.0 ? slim : 0.0;
                    slope = amrex::Math::copysign(1.0, del) *
                            amrex::min(slim, amrex::Math::abs(del));
                    s0(n, hi + 3) = s0(n - 1, r_crse) - 0.25 * slope;
                    s0(n, hi + 4) = s0(n - 1, r_crse) + 0.25 * slope;
                }
            } else {
                if (lo != 0) {
                    // quadratic interpolation from the three closest points
                    s0(n, lo - 1) = -s0(n, lo + 1) / 3.0 + s0(n, lo) +
                                    s0(n - 1, lo / 2 - 1) / 3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0(n, lo - 2) = s0(n - 1, (lo - 2) / 2);
                }

                if (hi + 1 != base_geom.nr(n)) {
                    // quadratic interpolation from the three closest points
                    s0(n, hi + 2) = -s0(n, hi) / 3.0 + s0(n, hi + 1) +
                                    s0(n - 1, (hi + 3) / 2) / 3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0(n, hi + 3) = s0(n - 1, (hi + 3) / 2);
                }
            }
        }
    }
}
