
#include <BaseStateGeometry.H>
#include <maestro_params.H>

using namespace amrex;
using namespace maestro;

void BaseStateGeometry::Init(const int max_radial_level_in,
                             const int nr_fine_in, const Real dr_fine_in,
                             const int nr_irreg_in,
                             const Vector<Geometry>& geom, const int max_level,
                             GpuArray<Real, 3>& center) {
    // timer for profiling
    BL_PROFILE_VAR("BaseStateGeometry::Init()", Init);

    Print() << "Calling BaseStateGeometry::Init()" << std::endl;

    const auto probLo = geom[0].ProbLoArray();
    const auto probHi = geom[0].ProbHiArray();

    max_radial_level = max_radial_level_in;
    finest_radial_level = 0;  // This will be reset after regridding.
    nr_fine = nr_fine_in;
    dr_fine = dr_fine_in;
    nr_irreg = nr_irreg_in;

    dr_d.define(max_radial_level + 1);
    nr_d.define(max_radial_level + 1);
    r_cc_loc_d.define(max_radial_level + 1, nr_fine);
    r_edge_loc_d.define(max_radial_level + 1, nr_fine + 1);

    base_cutoff_density_coord_d.define(max_radial_level + 1);
    anelastic_cutoff_density_coord_d.define(max_radial_level + 1);
    burning_cutoff_density_lo_coord_d.define(max_radial_level + 1);
    burning_cutoff_density_hi_coord_d.define(max_radial_level + 1);

    // initialize BaseStateArrays
    dr.init(dr_d);
    nr.init(nr_d);
    r_cc_loc.init(r_cc_loc_d);
    r_edge_loc.init(r_edge_loc_d);

    base_cutoff_density_coord.init(base_cutoff_density_coord_d);
    anelastic_cutoff_density_coord.init(anelastic_cutoff_density_coord_d);
    burning_cutoff_density_lo_coord.init(burning_cutoff_density_lo_coord_d);
    burning_cutoff_density_hi_coord.init(burning_cutoff_density_hi_coord_d);

    // compute center(:)
    if (octant) {
        for (auto i = 0; i < 3; ++i) {
            if (!(spherical && AMREX_SPACEDIM == 3 && probLo[i] == 0.0)) {
                Abort("ERROR: octant requires spherical with prob_lo = 0.0");
            }
            center[i] = 0.0;
        }
    } else {
        for (auto i = 0; i < AMREX_SPACEDIM; ++i) {
            center[i] = 0.5 * (probLo[i] + probHi[i]);
        }
    }

    // compute nr(:) and dr(:)
    nr(max_radial_level) = nr_fine;
    dr(max_radial_level) = dr_fine;

    // computes dr, nr, r_cc_loc, r_edge_loc
    if (!spherical) {
        // cartesian case

        // compute nr(:) and dr(:) assuming refinement ratio = 2
        for (auto n = max_radial_level - 1; n >= 0; --n) {
            nr(n) = nr(n + 1) / 2;
            dr(n) = dr(n + 1) * 2.0;
        }

        // compute r_cc_loc, r_edge_loc
        for (auto n = 0; n <= max_radial_level; ++n) {
            for (auto i = 0; i < nr(n); ++i) {
                r_cc_loc(n, i) =
                    probLo[AMREX_SPACEDIM - 1] + (Real(i) + 0.5) * dr(n);
            }
            for (auto i = 0; i <= nr(n); ++i) {
                r_edge_loc(n, i) =
                    probLo[AMREX_SPACEDIM - 1] + (Real(i)) * dr(n);
            }
        }
    } else {
        // spherical case
        // compute r_cc_loc, r_edge_loc
        if (use_exact_base_state) {
            const Real* dx_fine = geom[max_level].CellSize();
            // nr_fine = nr_irreg + 1
            for (auto i = 0; i < nr_fine; ++i) {
                r_cc_loc(0, i) = sqrt(0.75 + 2.0 * Real(i)) * dx_fine[0];
            }
            r_edge_loc(0) = 0.0;
            for (auto i = 0; i < nr_fine; ++i) {
                r_edge_loc(0, i + 1) =
                    sqrt(0.75 + 2.0 * (Real(i) + 0.5)) * dx_fine[0];
            }
        } else {
            for (auto i = 0; i < nr_fine; ++i) {
                r_cc_loc(0, i) = (Real(i) + 0.5) * dr(0);
            }
            for (auto i = 0; i <= nr_fine; ++i) {
                r_edge_loc(0, i) = Real(i) * dr(0);
            }
        }
    }
}

void BaseStateGeometry::ComputeCutoffCoords(const BaseStateArray<Real>& rho0) {
    // timer for profiling
    BL_PROFILE_VAR("BaseStateGeometry::ComputeCutoffCoords",
                   ComputeCutoffCoords);

    // compute the coordinates of the anelastic cutoff
    bool found = false;
    int which_lev = 0;

    // find the finest level containing the anelastic cutoff density,
    // and set the anelastic cutoff coord for this level
    for (auto n = finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks(n); ++i) {
            if (!found) {
                int lo = r_start_coord(n, i);
                int hi = r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) <= anelastic_cutoff_density) {
                        anelastic_cutoff_density_coord(n) = r;
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
        which_lev = finest_radial_level;
        anelastic_cutoff_density_coord(finest_radial_level) =
            nr(finest_radial_level);
    }

    // set the anelastic cutoff coordinate on the finer levels
    for (auto n = which_lev + 1; n <= finest_radial_level; ++n) {
        anelastic_cutoff_density_coord(n) =
            2 * anelastic_cutoff_density_coord(n - 1) + 1;
    }

    // set the anelastic cutoff coordinate on the coarser levels
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (anelastic_cutoff_density_coord(n + 1) % 2 == 0) {
            anelastic_cutoff_density_coord(n) =
                anelastic_cutoff_density_coord(n + 1) / 2;
        } else {
            anelastic_cutoff_density_coord(n) =
                anelastic_cutoff_density_coord(n + 1) / 2 + 1;
        }
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // compute the coordinates of the base cutoff density
    found = false;

    // find the finest level containing the base cutoff density,
    // and set the base cutoff coord for this level
    for (auto n = finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks(n); ++i) {
            if (!found) {
                int lo = r_start_coord(n, i);
                int hi = r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) <= base_cutoff_density) {
                        base_cutoff_density_coord(n) = r;
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
        which_lev = finest_radial_level;
        base_cutoff_density_coord(finest_radial_level) =
            nr(finest_radial_level);
    }

    // set the base cutoff coordinate on the finer levels
    // do n=which_lev+1,finest_radial_level
    for (auto n = which_lev + 1; n <= finest_radial_level; ++n) {
        base_cutoff_density_coord(n) = 2 * base_cutoff_density_coord(n - 1) + 1;
    }

    // set the base cutoff coordinate on the coarser levels
    // do n=which_lev-1,0,-1
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (base_cutoff_density_coord(n + 1) % 2 == 0) {
            base_cutoff_density_coord(n) = base_cutoff_density_coord(n + 1) / 2;
        } else {
            base_cutoff_density_coord(n) =
                base_cutoff_density_coord(n + 1) / 2 + 1;
        }
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // compute the coordinates of the burning cutoff density
    found = false;

    // find the finest level containing the burning cutoff density,
    // and set the burning cutoff coord for this level
    for (auto n = finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks(n); ++i) {
            if (!found) {
                //  do r=r_start_coord(n,i),r_end_coord(n,i)
                int lo = r_start_coord(n, i);
                int hi = r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) <= burning_cutoff_density_lo) {
                        burning_cutoff_density_lo_coord(n) = r;
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
        which_lev = finest_radial_level;
        burning_cutoff_density_lo_coord(finest_radial_level) =
            nr(finest_radial_level);
    }

    // set the burning cutoff coordinate on the finer levels
    // do n=which_lev+1,finest_radial_level
    for (auto n = which_lev + 1; n <= finest_radial_level; ++n) {
        burning_cutoff_density_lo_coord(n) =
            2 * burning_cutoff_density_lo_coord(n - 1) + 1;
    }

    // set the burning cutoff coordinate on the coarser levels
    // do n=which_lev-1,0,-1
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (burning_cutoff_density_lo_coord(n + 1) % 2 == 0) {
            burning_cutoff_density_lo_coord(n) =
                burning_cutoff_density_lo_coord(n + 1) / 2;
        } else {
            burning_cutoff_density_lo_coord(n) =
                burning_cutoff_density_lo_coord(n + 1) / 2 + 1;
        }
    }

    // compute the coordinates of the burning cutoff density upper limit
    found = false;

    // find the finest level containing the burning cutoff density,
    // and set the burning cutoff coord for this level
    for (auto n = finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks(n); ++i) {
            if (!found) {
                //  do r=r_end_coord(n,i),r_start_coord(n,i),1
                int lo = r_start_coord(n, i);
                int hi = r_end_coord(n, i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0(n, r) >= burning_cutoff_density_hi) {
                        burning_cutoff_density_hi_coord(n) = r;
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
        which_lev = finest_radial_level;
        burning_cutoff_density_hi_coord(finest_radial_level) = 0;
    }

    // set the burning cutoff coordinate on the finer levels
    for (auto n = which_lev + 1; n <= finest_radial_level; ++n) {
        burning_cutoff_density_hi_coord(n) = 0;
    }

    // set the burning cutoff coordinate on the coarser levels
    for (auto n = which_lev - 1; n >= 0; --n) {
        if (burning_cutoff_density_hi_coord(n + 1) % 2 == 0) {
            burning_cutoff_density_hi_coord(n) = 0;
        } else {
            burning_cutoff_density_hi_coord(n) = 0;
        }
    }
}

void BaseStateGeometry::InitMultiLevel(const int finest_radial_level_in,
                                       const BaseStateArray<int>& tag_array) {
    // compute numdisjointchunks, r_start_coord, r_end_coord
    // FIXME - right now there is one chunk at each level that spans the domain

    // NOTE: in the Fortran r_start_coord and r_end_coord had the shapes
    // r_start_coord(0:finest_radial_level,1:maxchunks), so here have offset
    // second index by 1 so it can be indexed from 0.

    // timer for profiling
    BL_PROFILE_VAR("BaseStateGeometry::InitMultiLevel", InitMultilevel);

    finest_radial_level = spherical ? 0 : finest_radial_level_in;

    numdisjointchunks_d.define(finest_radial_level + 1);
    numdisjointchunks.init(numdisjointchunks_d);

    // loop through tag_array first to determine the maximum number of chunks
    // to use for allocating r_start_coord and r_end_coord
    int maxchunks = 1;

    if (!spherical) {
        for (auto n = 1; n <= finest_radial_level; ++n) {
            // initialize variables
            bool chunk_start = false;
            int nchunks = 0;

            // increment nchunks at beginning of each chunk
            // (ex. when the tagging index changes from 0 to 1)
            for (auto r = 0; r < nr(n - 1); ++r) {
                if (tag_array(n - 1, r) > 0 && !chunk_start) {
                    chunk_start = true;
                    nchunks++;
                } else if (tag_array(n - 1, r) == 0 && chunk_start) {
                    chunk_start = false;
                }
            }
            maxchunks = amrex::max(nchunks, maxchunks);
        }
    }

    // maxchunks+1 here because we will be using indices 1:maxchunks
    r_start_coord_d.define(finest_radial_level + 1, maxchunks + 1);
    r_end_coord_d.define(finest_radial_level + 1, maxchunks + 1);

    r_start_coord.init(r_start_coord_d);
    r_end_coord.init(r_end_coord_d);

    if (!spherical) {
        // coarsest grid always has 1 chunk of data
        numdisjointchunks(0) = 1;
        r_start_coord(0, 1) = 0;
        r_end_coord(0, 1) = nr(0) - 1;

        // for > 1 chunks (multilevel)
        for (auto n = 1; n <= finest_radial_level; ++n) {
            // initialize variables
            bool chunk_start = false;
            numdisjointchunks(n) = 0;

            // increment numdisjointchunks at beginning of each chunk
            // (ex. when the tagging index changes from 0 to 1)
            for (auto r = 0; r < nr(n - 1); ++r) {
                if (tag_array(n - 1, r) > 0 && !chunk_start) {
                    chunk_start = true;
                    numdisjointchunks(n)++;
                    r_start_coord(n, numdisjointchunks(n)) = 2 * r;
                } else if (tag_array(n - 1, r) == 0 && chunk_start) {
                    r_end_coord(n, numdisjointchunks(n)) = 2 * r - 1;
                    chunk_start = false;
                } else if (r == nr(n - 1) - 1 && chunk_start) {
                    // if last chunk is at the end of array
                    r_end_coord(n, numdisjointchunks(n)) = 2 * r + 1;
                }
            }
        }
    } else {
        numdisjointchunks(0) = 1;
        r_start_coord(0, 1) = 0;
        r_end_coord(0, 1) = nr(0) - 1;
    }
}
