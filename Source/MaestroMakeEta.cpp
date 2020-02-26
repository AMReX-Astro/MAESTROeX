
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute eta_rho at edge- and cell-centers
void
Maestro::MakeEtarho (RealVector& etarho_edge,
                     RealVector& etarho_cell,
                     const Vector<MultiFab>& etarho_flux)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEtarho()",MakeEtarho);

#ifdef AMREX_USE_CUDA
    bool launched;
    if (deterministic_nodal_solve) {
        launched = !Gpu::notInLaunchRegion();
        // turn off GPU
        if (launched) Gpu::setLaunchRegion(false);
    }
#endif

    // Local variables
    const int nrf = nr_fine + 1;
    const int max_lev = max_radial_level + 1;
    RealVector etarhosum( (nr_fine+1)*(max_radial_level+1), 0.0 );
    etarhosum.shrink_to_fit();

    // this stores how many cells there are laterally at each level
    RealVector ncell(max_radial_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space of the domain
        const Box& domainBox = geom[lev].Domain();

        // compute number of cells at any given height for each level
        if (AMREX_SPACEDIM==2) {
            ncell[lev] = domainBox.bigEnd(0)+1;
        }
        else if (AMREX_SPACEDIM==3) {
            ncell[lev] = (domainBox.bigEnd(0)+1)*(domainBox.bigEnd(1)+1);
        }

        // get references to the MultiFabs at level lev
        const MultiFab& sold_mf = sold[lev];
        
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sold_mf, false); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid tile region
            const Box& validbox = mfi.validbox();

            const Array4<const Real> etarhoflux_arr = etarho_flux[lev].array(mfi);
            Real * AMREX_RESTRICT etarhosum_p = etarhosum.dataPtr();

#if (AMREX_SPACEDIM == 2)
            int zlo = 0;
            AMREX_PARALLEL_FOR_3D(validbox, i, j, k, {
                if (k == zlo) {
                    amrex::Gpu::Atomic::Add(&(etarhosum_p[j+nrf*lev]), etarhoflux_arr(i,j,k));
                }
            });

            // we only add the contribution at the top edge if we are at the top of the domain
            // this prevents double counting
            auto top_edge = false;
            for (auto i = 1; i <= numdisjointchunks[lev]; ++i) {
                if (validbox.hiVect()[1] == r_end_coord[lev+max_lev*i]) {
                    top_edge = true;
                }
            }
            if (top_edge) {
                int k = 0;
                int j = validbox.hiVect()[1]+1;
                int lo = validbox.loVect()[0];
                int hi = validbox.hiVect()[0];

                AMREX_PARALLEL_FOR_1D(hi-lo+1, n, {
                    int i = n + lo;
                    amrex::Gpu::Atomic::Add(&(etarhosum_p[j+nrf*lev]), etarhoflux_arr(i,j,k));
                });
            }
#else 
            AMREX_PARALLEL_FOR_3D(validbox, i, j, k, {
                amrex::Gpu::Atomic::Add(&(etarhosum_p[k+nrf*lev]), etarhoflux_arr(i,j,k));
            });

            // we only add the contribution at the top edge if we are at the top of the domain
            // this prevents double counting
            auto top_edge = false;

            for (auto i = 1; i <= numdisjointchunks[lev]; ++i) {
                if (validbox.hiVect()[2] == r_end_coord[lev+max_lev*i]) {
                    top_edge = true;
                }
            }
            
            if (top_edge) {
                int zhi = validbox.hiVect()[2] + 1;
                const auto& zbx = mfi.nodaltilebox(2);
                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    if (k == zhi) {
                        amrex::Gpu::Atomic::Add(&(etarhosum_p[k+nrf*lev]), etarhoflux_arr(i,j,k));
                    }
                });
            }
#endif
        }
    }

    ParallelDescriptor::ReduceRealSum(etarhosum.dataPtr(),(nr_fine+1)*(max_radial_level+1));

    get_numdisjointchunks(numdisjointchunks.dataPtr());
    get_r_start_coord(r_start_coord.dataPtr());
    get_r_end_coord(r_end_coord.dataPtr());
    get_finest_radial_level(&finest_radial_level);

    std::fill(etarho_edge.begin(), etarho_edge.end(), 0.0);
    std::fill(etarho_cell.begin(), etarho_cell.end(), 0.0);

    Real * AMREX_RESTRICT etarho_edge_p = etarho_edge.dataPtr();
    const Real * AMREX_RESTRICT etarhosum_p = etarhosum.dataPtr();
    Real * AMREX_RESTRICT etarho_cell_p = etarho_cell.dataPtr();

    for (auto n = 0; n <= finest_radial_level; ++n) {
        for (auto i = 1; i <= numdisjointchunks[n]; ++i) {
            Real ncell_lev = ncell[n];
            const int lo = r_start_coord[n+max_lev*i];
            const int hi = r_end_coord[n+max_lev*i]+1;
            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                // note swapped shaping for etarhosum
                etarho_edge_p[n+max_lev*r] = etarhosum_p[r+nrf*n] / ncell_lev;
            });
        }
    }

    // These calls shouldn't be needed since the planar algorithm doesn't use
    // these outside of this function, but this is just to be safe in case
    // things change in the future.
    RestrictBase(etarho_edge, false);
    FillGhostBase(etarho_edge, false);

    // make the cell-centered etarho_cc by averaging etarho to centers
    for (auto n = 0; n <= finest_radial_level; ++n) {
        for (auto i = 1; i <= numdisjointchunks[n]; ++i) {
            const int lo = r_start_coord[n+max_lev*i];
            const int hi = r_end_coord[n+max_lev*i]+1;
            AMREX_PARALLEL_FOR_1D(hi-lo+1, j, {
                int r = j + lo;
                etarho_cell_p[n+max_lev*r] = 0.5 * (etarho_edge_p[n+max_lev*r] + 
                    etarho_edge_p[n+max_lev*(r+1)]);
            });
        }
    }

    // These calls shouldn't be needed since the planar algorithm only uses
    // etarho_cell to make_psi, and then we fill ghost cells in make_psi, but
    // this is just to be safe in case things change in the future
    RestrictBase(etarho_cell, true);
    FillGhostBase(etarho_cell, true);

#ifdef AMREX_USE_CUDA
    if (deterministic_nodal_solve) {
        // turn GPU back on
        if (launched) Gpu::setLaunchRegion(true);
    }
#endif
}


void
Maestro::MakeEtarhoSphr (const Vector<MultiFab>& scal_old,
                         const Vector<MultiFab>& scal_new,
                         const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                         const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                         RealVector& etarho_edge,
                         RealVector& etarho_cell)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEtarhoSphr()",MakeEtarhoSphr);

    const int max_lev = max_radial_level + 1;
    const int nrf = nr_fine + 1;

    Vector<MultiFab> eta_cart(finest_level+1);
    Vector<MultiFab> rho0_nph_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        eta_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_nph_cart[lev].define(grids[lev],dmap[lev],1,0);
    }

    RealVector rho0_nph(max_lev*nr_fine);

    const Real * AMREX_RESTRICT rho0_old_p = rho0_old.dataPtr();
    const Real * AMREX_RESTRICT rho0_new_p = rho0_new.dataPtr();
    Real * AMREX_RESTRICT rho0_nph_p = rho0_nph.dataPtr();

    AMREX_PARALLEL_FOR_1D(max_lev*nr_fine, i, {
        rho0_nph_p[i] = 0.5*(rho0_old_p[i]+rho0_new_p[i]);
    });

    Put1dArrayOnCart(rho0_nph, rho0_nph_cart, 0, 0, bcs_f, 0);

#if (AMREX_SPACEDIM == 3)
    for (int lev=0; lev<=finest_level; ++lev) {

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scal_old[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> rho_old = scal_old[lev].array(mfi);
            const Array4<const Real> rho_new = scal_new[lev].array(mfi);
            const Array4<const Real> rho0_nph_cart_arr = rho0_nph_cart[lev].array(mfi);
            const Array4<const Real> umac_arr = umac[lev][0].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
            const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
            const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
            const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
            const Array4<const Real> normal_arr = normal[lev].array(mfi);
            const Array4<Real> eta_cart_arr = eta_cart[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                Real U_dot_er = 0.5*(umac_arr(i,j,k) + umac_arr(i+1,j,k) 
                    + w0macx(i,j,k) + w0macx(i+1,j,k)) * normal_arr(i,j,k,0) + 
                    0.5*(vmac(i,j,k) + vmac(i,j+1,k) 
                    + w0macy(i,j,k) + w0macy(i,j+1,k)) * normal_arr(i,j,k,1) + 
                    0.5*(wmac(i,j,k) + wmac(i,j,k+1) 
                    + w0macz(i,j,k) + w0macz(i,j,k+1)) * normal_arr(i,j,k,2);

                // construct time-centered [ rho' (U dot e_r) ]
                eta_cart_arr(i,j,k) = (0.5*(rho_old(i,j,k) + rho_new(i,j,k)) 
                    - rho0_nph_cart_arr(i,j,k)) * U_dot_er;
            });
        }         // end MFIter loop
    }     // end loop over levels

#else
    Abort("MakeEtarhoSphr: Spherical is not valid for DIM != 3");
#endif

    // average fine data onto coarser cells & fill ghost cells
    AverageDown(eta_cart,0,1);
    FillPatch(t_old, eta_cart, eta_cart, eta_cart, 0, 0, 1, 0, bcs_f);

    // compute etarho_cc as the average of eta_cart = [ rho' (U dot e_r) ]
    Average(eta_cart,etarho_cell,0);

    const Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    const Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();
    Real * AMREX_RESTRICT etarho_edge_p = etarho_edge.dataPtr();
    const Real * AMREX_RESTRICT etarho_cell_p = etarho_cell.dataPtr();

    // put eta on base state edges
    // note that in spherical the base state has no refinement
    // the 0th value of etarho = 0, since U dot . e_r must be
    // zero at the center (since e_r is not defined there)
    const bool sph_loc = spherical;
    AMREX_PARALLEL_FOR_1D(nrf, r, {
        if (r == 0) {
            etarho_edge_p[r] = 0.0;
        } else if (r == nrf-1) {
            // probably should do some better extrapolation here eventually
            etarho_edge_p[r] = etarho_cell_p[r-1];
        } else {
            if (sph_loc) {
                Real dr1 = r_cc_loc_p[r] - r_edge_loc_p[r];
                Real dr2 = r_edge_loc_p[r] - r_cc_loc_p[r-1];
                etarho_edge_p[r] = (dr2*etarho_cell_p[r] + dr1*etarho_cell_p[r-1])/(dr1+dr2);
            } else {
                etarho_edge_p[r] = 0.5*(etarho_cell_p[r] + etarho_cell_p[r-1]);
            }
        }
    });
}
