
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// Given a multifab of data (phi), average down to a base state quantity, phibar.
// If we are in plane-parallel, the averaging is at constant height.
// If we are spherical, { the averaging is done at constant radius.

void Maestro::Average (const Vector<MultiFab>& phi,
                       RealVector& phibar,
                       int comp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Average()",Average);

    const int max_lev = max_radial_level+1;
    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());
    get_finest_radial_level(&finest_radial_level);

    for (int i = 0; i < phibar.size(); ++i) {
        phibar[i] = 0.0;
    }

    if (spherical == 0) {

        // planar case

        // phibar is dimensioned to "max_radial_level" so we must mimic that for phisum
        // so we can simply swap this result with phibar
        RealVector phisum((max_radial_level+1)*nr_fine,0.0);
        phisum.shrink_to_fit();
        Real * AMREX_RESTRICT phisum_p = phisum.dataPtr();

        // this stores how many cells there are laterally at each level
        Vector<int> ncell(max_radial_level+1);

        // loop is over the existing levels (up to finest_level)
        for (int lev=0; lev<=finest_level; ++lev) {

            // Get the index space of the domain
            const Box domainBox = geom[lev].Domain();

            // compute number of cells at any given height for each level
            if (AMREX_SPACEDIM==2) {
                ncell[lev] = domainBox.bigEnd(0)+1;
            }
            else if (AMREX_SPACEDIM==3) {
                ncell[lev] = (domainBox.bigEnd(0)+1)*(domainBox.bigEnd(1)+1);
            }

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                // Get the index space of the valid region
                const Box& tilebox = mfi.tilebox();

                const Array4<const Real> phi_arr = phi[lev].array(mfi, comp);

                AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                    int r = AMREX_SPACEDIM == 2 ? j : k;
                    amrex::Gpu::Atomic::Add(&(phisum_p[lev+max_lev*r]), phi_arr(i, j, k));
                });
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),(max_radial_level+1)*nr_fine);

        // divide phisum by ncell so it stores "phibar"
        for (int lev = 0; lev < max_lev; ++lev) {
            for (auto i = 1; i <= numdisjointchunks[lev]; ++i) { 
                const int lo = r_start_coord[lev+max_lev*i];
                const int hi = r_end_coord[lev+max_lev*i];
                Real ncell_lev = ncell[lev];
                AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                    int r = j + lo;
                    phisum_p[lev+max_lev*r] /= ncell_lev;
                });
            }
        }

        RestrictBase(phisum, true);
        FillGhostBase(phisum, true);

        // swap pointers so phibar contains the computed average
        std::swap(phisum,phibar);

    } else if (spherical == 1 && use_exact_base_state) {
        // spherical case with uneven base state spacing

        // phibar is dimensioned to "max_radial_level" so we must mimic that for phisum
        // so we can simply swap this result with phibar
        RealVector phisum((max_radial_level+1)*nr_fine,0.0);

        // this stores how many cells there are at each level
        Vector<int> ncell((max_radial_level+1)*nr_fine,0);

        // loop is over the existing levels (up to finest_level)
        for (int lev=0; lev<=finest_level; ++lev) {

            // get references to the MultiFabs at level lev
            const MultiFab& phi_mf = phi[lev];
            const MultiFab& cc_to_r = cell_cc_to_r[lev];

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
            for (MFIter mfi(phi_mf); mfi.isValid(); ++mfi) {

                // Get the index space of the valid region
                const Box& validBox = mfi.validbox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                average_sphr_irreg(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                                   BL_TO_FORTRAN_N_3D(phi_mf[mfi],comp),
                                   phisum.dataPtr(), ncell.dataPtr(),
                                   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),(max_radial_level+1)*nr_fine);
        ParallelDescriptor::ReduceIntSum(ncell.dataPtr(),(max_radial_level+1)*nr_fine);

        // divide phisum by ncell so it stores "phibar"
        for (int lev = 0; lev < max_lev; ++lev) {
            for (auto r = 0; r < nr_fine; ++r) {
                if (ncell[lev+max_lev*r] > 0) {
                    phisum[lev+max_lev*r] /= ncell[lev+max_lev*r];
                } else {
                    // keep value constant if it is outside the cutoff coords
                    phisum[lev+max_lev*r] = phisum[lev+max_lev*(r-1)];
                }
            }           
        }

        RestrictBase(phisum, true);
        FillGhostBase(phisum, true);

        // swap pointers so phibar contains the computed average
        std::swap(phisum,phibar);
    } else {
        // spherical case with even base state spacing

        // For spherical, we construct a 1D array at each level, phisum, that has space
        // allocated for every possible radius that a cell-center at each level can
        // map into.  The radial locations have been precomputed and stored in radii.
        RealVector phisum((finest_level+1)*(nr_irreg+2), 0.0);
        RealVector radii((finest_level+1)*(nr_irreg+3));
        IntVector ncell((finest_level+1)*(nr_irreg+2), 0);

        Real * AMREX_RESTRICT radii_p = radii.dataPtr();
        Real * AMREX_RESTRICT phisum_p = phisum.dataPtr();
        int * AMREX_RESTRICT ncell_p = ncell.dataPtr();
        const Real * AMREX_RESTRICT center_p = center.dataPtr();

        const int fine_lev = finest_level + 1;
        const int nr_irreg_loc = nr_irreg;

        // radii contains every possible distance that a cell-center at the finest
        // level can map into
        for (int lev=0; lev<=finest_level; ++lev) {

            // Get the index space of the domain
            const auto dx = geom[lev].CellSizeArray();

            AMREX_PARALLEL_FOR_1D(nr_irreg+1, r, {
                radii_p[lev+fine_lev*(r+1)] = std::sqrt(0.75+2.0*Real(r)) * dx[0];
            });

            radii[lev+fine_lev*(nr_irreg+2)] = 1.e99;
            radii[lev] = 0.0;
        }

        // loop is over the existing levels (up to finest_level)
        for (int lev=finest_level; lev>=0; --lev) {

            // Get the grid size of the domain
            const auto dx = geom[lev].CellSizeArray();
            const auto prob_lo = geom[lev].ProbLoArray();

            // get references to the MultiFabs at level lev
            const MultiFab& phi_mf = phi[lev];

            // create mask assuming refinement ratio = 2
            int finelev = lev+1;
            if (lev == finest_level) finelev = finest_level;

            const BoxArray& fba = phi[finelev].boxArray();
            const iMultiFab& mask = makeFineMask(phi_mf, fba, IntVect(2));

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)

            for ( MFIter mfi(phi_mf); mfi.isValid(); ++mfi ) {
                // Get the index space of the valid region
                const Box& tilebox = mfi.validbox();

                const Array4<const int> mask_arr = mask.array(mfi);
                const Array4<const Real> phi_arr = phi[lev].array(mfi, comp);

                bool use_mask = !(lev==fine_lev-1);

                AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                    Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                    Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                    Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                    // make sure the cell isn't covered by finer cells
                    bool cell_valid = true;
                    if (use_mask) {
                        if (mask_arr(i,j,k) == 1) cell_valid = false;
                    }

                    if (cell_valid) {
                        // compute distance to center
                        Real radius = sqrt(x*x + y*y + z*z);

                        // figure out which radii index this point maps into
                        int index = ((radius/dx[0])*(radius/dx[0]) - 0.75) / 2.0;

                        // due to roundoff error, need to ensure that we are in the proper radial bin
                        if (index < nr_irreg_loc) {
                            if (fabs(radius-radii_p[lev+fine_lev*(index+1)]) > fabs(radius-radii_p[lev+fine_lev*(index+2)])) {
                                index++;
                            }
                        }

                        amrex::Gpu::Atomic::Add(&(phisum_p[lev+fine_lev*(index+1)]), phi_arr(i,j,k));
                        amrex::Gpu::Atomic::Add(&(ncell_p[lev+fine_lev*(index+1)]), 1);
                    }
                });
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),(finest_level+1)*(nr_irreg+2));
        ParallelDescriptor::ReduceIntSum(ncell.dataPtr(),(finest_level+1)*(nr_irreg+2));

        // normalize phisum so it actually stores the average at a radius
        for (auto n = 0; n <= finest_level; ++n) {
            for (auto r = 0; r <= nr_irreg; ++r) {
                if (ncell[n+fine_lev*(r+1)] != 0) {
                    phisum[n+fine_lev*(r+1)] /= Real(ncell[n+fine_lev*(r+1)]);
                }
            }
        }

        IntVector which_lev(nr_fine);
        IntVector max_rcoord(fine_lev);

        // compute center point for the finest level
        phisum[finest_level] = (11.0/8.0) * phisum[finest_level+fine_lev]
            - (3.0/8.0) * phisum[finest_level+fine_lev*2];
        ncell[finest_level] = 1;

        // choose which level to interpolate from
        int * AMREX_RESTRICT which_lev_p = which_lev.dataPtr();
        const int dr0 = dr[0];
        const int nrf = nr_fine;

        AMREX_PARALLEL_FOR_1D(nrf, r, {

            Real radius = (Real(r) + 0.5) * dr0;
            // Vector<int> rcoord_p(fine_lev, 0);
            int rcoord_p[MAESTRO_MAX_LEVELS];

            for (auto i = 0; i < fine_lev; ++i) {
                rcoord_p[i] = 0.0;
            }

            // for each level, find the closest coordinate
            for (auto n = 0; n < fine_lev; ++n) {
                for (auto j = rcoord_p[n]; j <= nr_irreg_loc; ++j) {
                    if (fabs(radius-radii_p[n+fine_lev*(j+1)]) < fabs(radius-radii_p[n+fine_lev*(j+2)])) {
                        rcoord_p[n] = j;
                        break;
                    }
                }
            }

            // make sure closest coordinate is in bounds
            for (auto n = 0; n < fine_lev-1; ++n) {
                rcoord_p[n] = max(rcoord_p[n],1);
            }
            for (auto n = 0; n < fine_lev; ++n) {
                rcoord_p[n] = min(rcoord_p[n],nr_irreg_loc-1);
            }

            // choose the level with the largest min over the ncell interpolation points
            which_lev_p[r] = 0;

            int min_all = min(ncell_p[fine_lev*(rcoord_p[0])], 
                ncell_p[fine_lev*(rcoord_p[0]+1)], 
                ncell_p[fine_lev*(rcoord_p[0]+2)]);

            for (auto n = 1; n < fine_lev; ++n) {
                int min_lev = min(ncell_p[n+fine_lev*(rcoord_p[n])], 
                    ncell_p[n+fine_lev*(rcoord_p[n]+1)], 
                    ncell_p[n+fine_lev*(rcoord_p[n]+2)]);

                if (min_lev > min_all) {
                    min_all = min_lev;
                    which_lev_p[r] = n;
                }
            }

            // if the min hit count at all levels is zero, we expand the search
            // to find the closest instance of where the hitcount becomes nonzero
            int j = 1;
            while (min_all == 0) {
                j++;
                for (auto n = 0; n < fine_lev; ++n) {
                    int min_lev = max(ncell_p[n+fine_lev*(max(1,rcoord_p[n]-j)+1)], 
                        ncell_p[n+fine_lev*(min(rcoord_p[n]+j,nr_irreg_loc-1)+1)]);
                    if (min_lev != 0) {
                        which_lev_p[r] = n;
                        min_all = min_lev;
                        break;
                    }
                }
            }
        });

        // squish the list at each level down to exclude points with no contribution
        for (auto n = 0; n <= finest_level; ++n) {
            int j = 0;
            for (auto r = 0; r <= nr_irreg; ++r) {
                while (ncell[n+fine_lev*(j+1)] == 0) {
                    j++;
                    if (j > nr_irreg) {
                        break;
                    }
                }
                if (j > nr_irreg) {
                    for (auto i = r; i <= nr_irreg; ++i) {
                        phisum[n+fine_lev*(i+1)] = 1.e99;
                    }
                    for (auto i = r; i <= nr_irreg+1; ++i) {
                        radii[n+fine_lev*(i+1)] = 1.e99;
                    }
                    max_rcoord[n] = r - 1;
                    break;
                }
                phisum[n+fine_lev*(r+1)] = phisum[n+fine_lev*(j+1)];
                radii [n+fine_lev*(r+1)] = radii [n+fine_lev*(j+1)];
                ncell [n+fine_lev*(r+1)] = ncell [n+fine_lev*(j+1)];
                j++;
                if (j > nr_irreg) {
                    max_rcoord[n] = r;
                    break;
                }
            }
        }

        // compute phibar
        int * AMREX_RESTRICT max_rcoord_p = max_rcoord.dataPtr();
        Real * AMREX_RESTRICT phibar_p = phibar.dataPtr();

        const Real drdxfac_loc = drdxfac;

        AMREX_PARALLEL_FOR_1D(nrf, r, {

            Real radius = (Real(r) + 0.5) * dr0;
            int stencil_coord = 0;

            // find the closest coordinate
            for (auto j = stencil_coord; j <= max_rcoord_p[which_lev_p[r]]; ++j) {
                if (fabs(radius-radii_p[which_lev_p[r]+fine_lev*(j+1)]) < 
                    fabs(radius-radii_p[which_lev_p[r]+fine_lev*(j+2)])) {
                    stencil_coord = j;
                    break;
                }
            }

            // make sure the interpolation points will be in bounds
            if (which_lev_p[r] != fine_lev-1) {
                stencil_coord = max(stencil_coord, 1);
            }
            stencil_coord = min(stencil_coord, 
                    max_rcoord_p[which_lev_p[r]]-1);

            bool limit = (r > nrf - 1 - drdxfac_loc*pow(2.0, (fine_lev-2))) ? false : true;

            phibar_p[max_lev*r] = QuadInterp(radius, 
                    radii_p[which_lev_p[r]+fine_lev*(stencil_coord)], 
                    radii_p[which_lev_p[r]+fine_lev*(stencil_coord+1)], 
                    radii_p[which_lev_p[r]+fine_lev*(stencil_coord+2)], 
                    phisum_p[which_lev_p[r]+fine_lev*(stencil_coord)], 
                    phisum_p[which_lev_p[r]+fine_lev*(stencil_coord+1)], 
                    phisum_p[which_lev_p[r]+fine_lev*(stencil_coord+2)], limit);
        });
    }
}