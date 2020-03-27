#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseStateGeometry(const int max_radial_level_in, 
                               const int nr_fine_in,
                               const Real dr_fine_in,
                               const int nr_irreg_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseStateGeometry()", InitBaseStateGeometry); 

    Print() << "Calling InitBaseStateGeometry()" << std::endl;

    const auto probLo = geom[0].ProbLoArray();
    const auto probHi = geom[0].ProbHiArray();

    max_radial_level = max_radial_level_in;
    finest_radial_level = max_radial_level_in; // FIXME - we want to set this after regridding
    nr_fine = nr_fine_in;
    dr_fine = dr_fine_in;
    nr_irreg = nr_irreg_in;

    dr.resize(max_radial_level+1);
    nr.resize(max_radial_level+1);
    dr.resize(max_radial_level+1);
    nr.resize(max_radial_level+1);

    base_cutoff_density_coord.resize(max_radial_level+1);
    anelastic_cutoff_density_coord.resize(max_radial_level+1);
    burning_cutoff_density_lo_coord.resize(max_radial_level+1);
    burning_cutoff_density_hi_coord.resize(max_radial_level+1);

    const int max_lev = max_radial_level+1;

    // compute center(:)
    if (octant) {
        for (auto i = 0; i < 3; ++i) {
            if (!(spherical == 1 && AMREX_SPACEDIM == 3 && 
                    probLo[i] == 0.0 ) ) {
                Abort("ERROR: octant requires spherical with prob_lo = 0.0");
            }
            center[i] = 0.0;
        }
    } else {
        for (auto i = 0; i < AMREX_SPACEDIM; ++i) {
            center[i] = 0.5*(probLo[i] + probHi[i]);
        }
    }

    // compute nr(:) and dr(:)
    nr(max_radial_level) = nr_fine;
    dr(max_radial_level) = dr_fine;

    // computes dr, nr, r_cc_loc, r_edge_loc
    if (spherical == 0) {
        // cartesian case

        // compute nr(:) and dr(:) assuming refinement ratio = 2
        for (auto n = max_radial_level-1; n >= 0; --n) {
          nr(n) = nr(n+1) / 2;
          dr(n) = dr(n+1) * 2.0;
        }

        // compute r_cc_loc, r_edge_loc
        for (auto n = 0; n <= max_radial_level; ++n) {
            for (auto i = 0; i < nr(n); ++i) {
                r_cc_loc_b(n,i) = probLo[AMREX_SPACEDIM-1] + (Real(i)+0.5)*dr(n);

                r_cc_loc_b(n,i) = probLo[AMREX_SPACEDIM-1] + (Real(i)+0.5)*dr(n);
            }
            for (auto i = 0; i <= nr(n); ++i) {
                r_edge_loc_b(n,i) = probLo[AMREX_SPACEDIM-1] + (Real(i))*dr(n);
                r_edge_loc_b(n,i) = probLo[AMREX_SPACEDIM-1] + (Real(i))*dr(n);
            }
        }
    } else {

        // spherical case
        // compute r_cc_loc, r_edge_loc
        if (use_exact_base_state) {
            const Real* dx_fine = geom[max_level].CellSize();
            // nr_fine = nr_irreg + 1
            for (auto i = 0; i < nr_fine; ++i) {
                r_cc_loc_b(0,i) = sqrt(0.75+2.0*Real(i))*dx_fine[0];
                r_cc_loc_b(0,i) = sqrt(0.75+2.0*Real(i))*dx_fine[0];
            }
            r_edge_loc_b(0,0) = 0.0;
            r_edge_loc_b(0,0) = 0.0;
            for (auto i = 0; i < nr_fine; ++i) {
                r_edge_loc_b(0,i+1) = sqrt(0.75+2.0*(Real(i)+0.5))*dx_fine[0];
                r_edge_loc_b(0,i+1) = sqrt(0.75+2.0*(Real(i)+0.5))*dx_fine[0];
            }
        } else {
            for (auto i = 0; i < nr_fine; ++i) {
                r_cc_loc_b(0,i) = (Real(i)+0.5)*dr(0);
                r_cc_loc_b(0,i) = (Real(i)+0.5)*dr(0);
            }
            for (auto i = 0; i <= nr_fine; ++i) {
                r_edge_loc_b(0,i) = Real(i)*dr(0);
                r_edge_loc_b(0,i) = Real(i)*dr(0);
            }
        }
    }
}

#if (AMREX_SPACEDIM == 3)
void 
Maestro::InitBaseStateMapSphr(const int lev, const MFIter& mfi, 
                              const GpuArray<Real,AMREX_SPACEDIM> dx_fine, const GpuArray<Real,AMREX_SPACEDIM> dx_lev) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseStateMapSphr()", InitBaseStateMapSphr); 

    if (!spherical) {
        Abort("init_base_state_map_sphr() does not work for planar");
    }

    if (use_exact_base_state) {

        const Box& tilebox = mfi.tilebox();
        const Array4<Real> cc_to_r = cell_cc_to_r[lev].array(mfi);

        const auto probLo = geom[0].ProbLoArray();

        AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
            Real x = probLo[0] + (static_cast<Real>(i)+0.5)*dx_lev[0] - center[0];
            Real y = probLo[1] + (static_cast<Real>(j)+0.5)*dx_lev[1] - center[1];
            Real z = probLo[2] + (static_cast<Real>(k)+0.5)*dx_lev[2] - center[2];

            Real index = (x*x + y*y + z*z)/(2.0*dx_fine[0]*dx_fine[0]) - 0.375;
            cc_to_r(i,j,k) = round(index);
        })
    }
}
#endif

void 
Maestro::ComputeCutoffCoords(const RealVector& rho0)
{
    // compute the coordinates of the anelastic cutoff
    bool found = false;
    const int max_lev = max_radial_level + 1;
    int which_lev = 0;

    // find the finest level containing the anelastic cutoff density,
    // and set the anelastic cutoff coord for this level
    for (auto n=finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks_b(n); ++i) {
            if (!found) {
                int lo = r_start_coord_b(n,i);
                int hi = r_end_coord_b(n,i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0[n + max_lev*r] <= anelastic_cutoff_density) {
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
        anelastic_cutoff_density_coord(finest_radial_level) = nr(finest_radial_level);
    }

    // set the anelastic cutoff coordinate on the finer levels
    for (auto n = which_lev+1; n <= finest_radial_level; ++n) {
        anelastic_cutoff_density_coord(n) = 2*anelastic_cutoff_density_coord(n-1)+1;
    }

    // set the anelastic cutoff coordinate on the coarser levels
    for (auto n = which_lev-1; n >= 0; --n) {
        if (anelastic_cutoff_density_coord(n+1) % 2 == 0) {
            anelastic_cutoff_density_coord(n) = anelastic_cutoff_density_coord(n+1) / 2;
        } else {
            anelastic_cutoff_density_coord(n) = anelastic_cutoff_density_coord(n+1) / 2 + 1;
        }
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // compute the coordinates of the base cutoff density
    found = false;

    // find the finest level containing the base cutoff density,
    // and set the base cutoff coord for this level
    for (auto n=finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks_b(n); ++i) {
            if (!found) {
                int lo = r_start_coord_b(n,i);
                int hi = r_end_coord_b(n,i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0[n + max_lev*r] <= base_cutoff_density) {
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
        base_cutoff_density_coord(finest_radial_level) = nr(finest_radial_level);
    }

    // set the base cutoff coordinate on the finer levels
    // do n=which_lev+1,finest_radial_level
    for (auto n = which_lev+1; n <= finest_radial_level; ++n) {
        base_cutoff_density_coord(n) = 2*base_cutoff_density_coord(n-1)+1;
    }

    // set the base cutoff coordinate on the coarser levels
    // do n=which_lev-1,0,-1
    for (auto n = which_lev-1; n >= 0; --n) {
        if (base_cutoff_density_coord(n+1) % 2 == 0) {
            base_cutoff_density_coord(n) = base_cutoff_density_coord(n+1) / 2;
        } else {
            base_cutoff_density_coord(n) = base_cutoff_density_coord(n+1) / 2 + 1;
        }
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // compute the coordinates of the burning cutoff density
    found = false;

    // find the finest level containing the burning cutoff density,
    // and set the burning cutoff coord for this level
    for (auto n=finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks_b(n); ++i) {
            if (!found) {
                //  do r=r_start_coord(n,i),r_end_coord(n,i)
                int lo = r_start_coord_b(n,i);
                int hi = r_end_coord_b(n,i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0[n + max_lev*r] <= burning_cutoff_density_lo) {
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
        burning_cutoff_density_lo_coord(finest_radial_level) = nr(finest_radial_level);
    }

    // set the burning cutoff coordinate on the finer levels
    // do n=which_lev+1,finest_radial_level
    for (auto n = which_lev+1; n <= finest_radial_level; ++n) {
        burning_cutoff_density_lo_coord(n) = 2*burning_cutoff_density_lo_coord(n-1)+1;
    }

    // set the burning cutoff coordinate on the coarser levels
    // do n=which_lev-1,0,-1
    for (auto n = which_lev-1; n >= 0; --n) {
        if (burning_cutoff_density_lo_coord(n+1)%2 == 0) {
            burning_cutoff_density_lo_coord(n) = burning_cutoff_density_lo_coord(n+1) / 2;
        } else {
            burning_cutoff_density_lo_coord(n) = burning_cutoff_density_lo_coord(n+1) / 2 + 1;
        }
    }

    // compute the coordinates of the burning cutoff density upper limit
    found = false;

    // find the finest level containing the burning cutoff density,
    // and set the burning cutoff coord for this level
    for (auto n=finest_radial_level; n >= 0; --n) {
        for (auto i = 1; i <= numdisjointchunks_b(n); ++i) {
            if (!found) {
                //  do r=r_end_coord(n,i),r_start_coord(n,i),1
                int lo = r_start_coord_b(n,i);
                int hi = r_end_coord_b(n,i);
                for (auto r = lo; r <= hi; ++r) {
                    if (rho0[n + max_lev*r] >= burning_cutoff_density_hi) {
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
    for (auto n = which_lev+1; n <= finest_radial_level; ++n) {
        burning_cutoff_density_hi_coord(n) = 0;
    }

    // set the burning cutoff coordinate on the coarser levels
    for (auto n = which_lev-1; n >= 0; --n) {
        if (burning_cutoff_density_hi_coord(n+1)%2 == 0) {
            burning_cutoff_density_hi_coord(n) = 0;
        } else {
            burning_cutoff_density_hi_coord(n) = 0;
        }
    }
}

void 
Maestro::InitMultilevel(const int finest_radial_level_in) {
    // compute numdisjointchunks, r_start_coord, r_end_coord
    // FIXME - right now there is one chunk at each level that spans the domain

    // NOTE: in the Fortran r_start_coord and r_end_coord had the shapes
    // r_start_coord(0:finest_radial_level,1:maxchunks), so here have offset
    // second index by 1 so it can be indexed from 0.

    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitMultilevel()", InitMultilevel); 

    const int max_lev = max_radial_level+1;

    if (spherical) {
       finest_radial_level = 0;
    } else {
       finest_radial_level = finest_radial_level_in;
    }

    numdisjointchunks_b.resize(finest_radial_level+1);

    // loop through tag_array first to determine the maximum number of chunks
    // to use for allocating r_start_coord and r_end_coord
    int maxchunks = 1;
    for (auto n = 1; n <= finest_radial_level; ++n) {

        // initialize variables
        bool chunk_start = false;
        int nchunks = 0;

        // increment nchunks at beginning of each chunk
        // (ex. when the tagging index changes from 0 to 1)
        for (auto r = 0; r < nr(n-1); ++r) {
            if (tag_array[n-1 + max_lev*r] > 0 && !chunk_start) {
                chunk_start = true;
                nchunks++;
            } else if (tag_array[n-1 + max_lev*r] == 0 && chunk_start) {
                chunk_start = false;
            }
        }
        maxchunks = max(nchunks, maxchunks);
    }

    r_start_coord_b.resize(finest_radial_level+1, maxchunks+1);
    r_end_coord_b.resize(finest_radial_level+1, maxchunks+1);

    if (!spherical) {

        // coarsest grid always has 1 chunk of data
        numdisjointchunks_b(0) = 1;
        r_start_coord_b(0,1) = 0;
        r_end_coord_b(0,1) = nr(0)-1;

        // for > 1 chunks (multilevel)
        for (auto n = 1; n <= finest_radial_level; ++n) {
            // initialize variables
            bool chunk_start = false;
            numdisjointchunks_b(n) = 0;

            // increment numdisjointchunks at beginning of each chunk
            // (ex. when the tagging index changes from 0 to 1)
            for (auto r = 0; r < nr(n-1); ++r) {
                if (tag_array[n-1 + max_lev*r] > 0 && !chunk_start) {
                    chunk_start = true;
                    numdisjointchunks_b(n)++;
                    r_start_coord_b(n,numdisjointchunks_b(n)) = 2*r;
                } else if (tag_array[n-1+max_lev*r]==0 && chunk_start) {
                    r_end_coord_b(n, numdisjointchunks_b(n)) = 2*r-1;
                    chunk_start = false;
                } else if (r==nr(n-1)-1 && chunk_start) {
                    // if last chunk is at the end of array
                    r_end_coord_b(n, numdisjointchunks_b(n)) = 2*r-1;
                }
            }
        }
    } else {
        numdisjointchunks_b(0) = 1;
        r_start_coord_b(0,1) = 0;
        r_end_coord_b(0,1) = nr(0)-1;
    }
}

void 
Maestro::RestrictBase(RealVector& s0, bool is_cell_centered)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RestrictBase()", RestrictBase); 

    const int max_lev = max_radial_level + 1;

    for (int n = finest_radial_level; n >= 1; --n) {        
        for (int i = 1; i <= numdisjointchunks_b(n); ++i) {
            if (is_cell_centered) {
                // for level n, make the coarser cells underneath simply the average of the fine
                for (auto j = r_start_coord_b(n,i); j < r_end_coord_b(n,i); j+=2) {
                    s0[n-1 + max_lev*j/2] = 0.5 * (s0[n + max_lev*j] + s0[n + max_lev*(j+1)]);
                }
            } else {
                // for level n, make the coarse edge underneath equal to the fine edge value
                for (auto j = r_start_coord_b(n,i); j <= r_end_coord_b(n,i)+1; j+=2) {
                    s0[n-1 + max_lev*j/2] = s0[n + max_lev*j];
                }
            }
        }
    }
}

void 
Maestro::RestrictBase(BaseState<Real>& s0, bool is_cell_centered)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RestrictBase()", RestrictBase); 

    const int max_lev = max_radial_level + 1;

    for (int n = finest_radial_level; n >= 1; --n) {        
        for (int i = 1; i <= numdisjointchunks_b(n); ++i) {
            if (is_cell_centered) {
                // for level n, make the coarser cells underneath simply the average of the fine
                for (auto j = r_start_coord_b(n,i); j < r_end_coord_b(n,i); j+=2) {
                    s0(n-1,j/2) = 0.5 * (s0(n,j) + s0(n,j+1));
                }
            } else {
                // for level n, make the coarse edge underneath equal to the fine edge value
                for (auto j = r_start_coord_b(n,i); j <= r_end_coord_b(n,i)+1; j+=2) {
                    s0(n-1,j/2) = s0(n,j);
                }
            }
        }
    }
}

void 
Maestro::FillGhostBase(RealVector& s0, bool is_cell_centered)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillGhostBase()", FillGhostBase); 

    const int max_lev = max_radial_level + 1;

    for (int n = finest_radial_level; n >= 1; --n) {
        for (int i = 1; i <= numdisjointchunks_b(n); ++i) {

            const int lo = r_start_coord_b(n,i);
            const int hi = r_end_coord_b(n,i);

            if (is_cell_centered) {

                if (lo != 0) {
                    int r_crse = lo/2-1;
                    Real del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    Real dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    Real dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    Real slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    Real slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(lo-1)] = s0[n-1+max_lev*r_crse] + 0.25*slope;
                    s0[n+max_lev*(lo-2)] = s0[n-1+max_lev*r_crse] - 0.25*slope;

                    r_crse = lo/2-2;
                    del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(lo-3)] = s0[n-1+max_lev*r_crse] + 0.25*slope;
                    s0[n+max_lev*(lo-4)] = s0[n-1+max_lev*r_crse] - 0.25*slope;
                }

                if (hi != nr(n)-1) {
                    int r_crse = (hi+1)/2;
                    Real del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    Real dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    Real dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    Real slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    Real slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(hi+1)] = s0[n-1+max_lev*r_crse] - 0.25*slope;
                    s0[n+max_lev*(hi+2)] = s0[n-1+max_lev*r_crse] + 0.25*slope;

                    r_crse = (hi+1)/2+1;
                    del = 0.5*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*(r_crse-1)]);
                    dpls = 2.0*(s0[n-1+max_lev*(r_crse+1)]-s0[n-1+max_lev*r_crse]);
                    dmin = 2.0*(s0[n-1+max_lev*r_crse]-s0[n-1+max_lev*(r_crse-1)]);
                    slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0[n+max_lev*(hi+3)] = s0[n-1+max_lev*r_crse] - 0.25*slope;
                    s0[n+max_lev*(hi+4)] = s0[n-1+max_lev*r_crse] + 0.25*slope;
                }
            } else {
                if (lo != 0) {
                    // quadratic interpolation from the three closest points
                    s0[n+max_lev*(lo-1)] = -s0[n+max_lev*(lo+1)]/3.0
                        + s0[n+max_lev*lo] + s0[n-1+max_lev*(lo/2-1)]/3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0[n+max_lev*(lo-2)] = s0[n-1+max_lev*((lo-2)/2)];
                }

                if (hi+1 != nr(n)) {
                    // quadratic interpolation from the three closest points
                    s0[n+max_lev*(hi+2)] = -s0[n+max_lev*hi]/3.0
                        + s0[n+max_lev*(hi+1)] + s0[n-1+max_lev*((hi+3)/2)]/3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0[n+max_lev*(hi+3)] = s0[n-1+max_lev*((hi+3)/2)];
                }
            }
        }
    }
}

void 
Maestro::FillGhostBase(BaseState<Real>& s0, bool is_cell_centered)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillGhostBase()", FillGhostBase); 

    const int max_lev = max_radial_level + 1;

    for (int n = finest_radial_level; n >= 1; --n) {
        for (int i = 1; i <= numdisjointchunks_b(n); ++i) {

            const int lo = r_start_coord_b(n,i);
            const int hi = r_end_coord_b(n,i);

            if (is_cell_centered) {

                if (lo != 0) {
                    int r_crse = lo/2-1;
                    Real del = 0.5*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1));
                    Real dpls = 2.0*(s0(n-1,r_crse+1)-s0(n-1,r_crse));
                    Real dmin = 2.0*(s0(n-1,r_crse)-s0(n-1,r_crse-1));
                    Real slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    Real slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0(n,lo-1) = s0(n-1,r_crse) + 0.25*slope;
                    s0(n,lo-2) = s0(n-1,r_crse) - 0.25*slope;

                    r_crse = lo/2-2;
                    del = 0.5*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1));
                    dpls = 2.0*(s0(n-1,r_crse+1)-s0(n-1,r_crse));
                    dmin = 2.0*(s0(n-1,r_crse)-s0(n-1,r_crse-1));
                    slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0(n,lo-3) = s0(n-1,r_crse) + 0.25*slope;
                    s0(n,lo-4) = s0(n-1,r_crse) - 0.25*slope;
                }

                if (hi != nr(n)-1) {
                    int r_crse = (hi+1)/2;
                    Real del = 0.5*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1));
                    Real dpls = 2.0*(s0(n-1,r_crse+1)-s0(n-1,r_crse));
                    Real dmin = 2.0*(s0(n-1,r_crse)-s0(n-1,r_crse-1));
                    Real slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    Real slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0(n,hi+1) = s0(n-1,r_crse) - 0.25*slope;
                    s0(n,hi+2) = s0(n-1,r_crse) + 0.25*slope;

                    r_crse = (hi+1)/2+1;
                    del = 0.5*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1));
                    dpls = 2.0*(s0(n-1,r_crse+1)-s0(n-1,r_crse));
                    dmin = 2.0*(s0(n-1,r_crse)-s0(n-1,r_crse-1));
                    slim = min(fabs(dpls),fabs(dmin));
                    slim = dpls*dmin > 0.0 ? slim : 0.0;
                    slope = copysign(1.0,del)*min(slim,fabs(del));
                    s0(n,hi+3) = s0(n-1,r_crse) - 0.25*slope;
                    s0(n,hi+4) = s0(n-1,r_crse) + 0.25*slope;
                }
            } else {
                if (lo != 0) {
                    // quadratic interpolation from the three closest points
                    s0(n,lo-1) = -s0(n,lo+1)/3.0
                        + s0(n,lo) + s0(n-1,lo/2-1)/3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0(n,lo-2) = s0(n-1,(lo-2)/2);
                }

                if (hi+1 != nr(n)) {
                    // quadratic interpolation from the three closest points
                    s0(n,hi+2) = -s0(n,hi)/3.0
                        + s0(n,hi+1) + s0(n-1,(hi+3)/2)/3.0;
                    // copy the next ghost cell value directly in from the coarser level
                    s0(n,hi+3) = s0(n-1,(hi+3)/2);
                }
            }
        }
    }
}
