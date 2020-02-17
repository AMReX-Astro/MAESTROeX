
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
        const MultiFab& etarhoflux_mf = etarho_flux[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid tile region
            const Box& validbox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells)
            // We will also pass "bx", which specifies the "valid" tile region.
            // sum_etarho(&lev, ARLIM_3D(domainBox.loVect()), ARLIM_3D(domainBox.hiVect()),
            //            ARLIM_3D(validbox.loVect()), ARLIM_3D(validbox.hiVect()),
            //            BL_TO_FORTRAN_3D(etarhoflux_mf[mfi]),
            //            etarhosum.dataPtr());

            const Array4<const Real> etarhoflux_arr = etarho_flux[lev].array(mfi);
            Real * AMREX_RESTRICT etarhosum_p = etarhosum.dataPtr();

#ifdef AMREX_USE_CUDA
            auto launched = !Gpu::notInLaunchRegion();
            // turn off GPU
            if (launched) Gpu::setLaunchRegion(false);
#endif

#if (AMREX_SPACEDIM == 2)
            int zlo = validbox.loVect()[2];
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(validbox, i, j, k, {
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
                int k = validbox.loVect()[2];
                int j = validbox.loVect()[1]+1;
                int lo = validbox.loVect()[0];
                int hi = validbox.hiVect()[0];

                AMREX_HOST_DEVICE_PARALLEL_FOR_1D(hi-lo+1, k, {
                    amrex::Gpu::Atomic::Add(&(etarhosum_p[j+nrf*lev]), etarhoflux_arr(i,j,k));
                });
            }
#else 
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(validbox, i, j, k, {
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
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(zbx, i, j, k, {
                    if (k == zhi) {
                        amrex::Gpu::Atomic::Add(&(etarhosum_p[k+nrf*lev]), etarhoflux_arr(i,j,k));
                    }
                });
            }
#endif

#ifdef AMREX_USE_CUDA
    // turn GPU back on
    if (launched) Gpu::setLaunchRegion(true);
#endif
        }
    }

    ParallelDescriptor::ReduceRealSum(etarhosum.dataPtr(),(nr_fine+1)*(max_radial_level+1));

    // const int max_lev = max_radial_level+1;
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
                etarho_edge_p[n+max_lev*r] = etarhosum_p[r+(nr_fine+1)*n] / ncell_lev;
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
            Real ncell_lev = ncell[n];
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

    // make_etarho_planar(etarho_edge.dataPtr(), etarho_cell.dataPtr(),
    //                    etarhosum.dataPtr(), ncell.dataPtr());
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

    Vector<MultiFab> eta_cart(finest_level+1);
    Vector<MultiFab> rho0_nph_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        eta_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_nph_cart[lev].define(grids[lev],dmap[lev],1,0);
    }

    RealVector rho0_nph( (max_radial_level+1)*nr_fine );

    for (int i=0; i<rho0_nph.size(); ++i) {
        rho0_nph[i] = 0.5*(rho0_old[i]+rho0_new[i]);
    }

    Put1dArrayOnCart(rho0_nph, rho0_nph_cart, 0, 0, bcs_f, 0);

#if (AMREX_SPACEDIM == 3)

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scalold_mf = scal_old[lev];
        const MultiFab& scalnew_mf = scal_new[lev];
        const MultiFab& umac_mf   = umac[lev][0];
        const MultiFab& vmac_mf   = umac[lev][1];
        const MultiFab& wmac_mf   = umac[lev][2];
        const MultiFab& w0macx_mf = w0mac[lev][0];
        const MultiFab& w0macy_mf = w0mac[lev][1];
        const MultiFab& w0macz_mf = w0mac[lev][2];
        const MultiFab& normal_mf = normal[lev];
        MultiFab& etacart_mf = eta_cart[lev];
        const MultiFab& rho0_nph_cart_mf = rho0_nph_cart[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scalold_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

#pragma gpu box(tilebox)
            construct_eta_cart( AMREX_INT_ANYD(tilebox.loVect()), AMREX_INT_ANYD(tilebox.hiVect()),
                                scalold_mf[mfi].dataPtr(Rho),
                                AMREX_INT_ANYD(scalold_mf[mfi].loVect()), 
                                AMREX_INT_ANYD(scalold_mf[mfi].hiVect()),
                                scalnew_mf[mfi].dataPtr(Rho),
                                AMREX_INT_ANYD(scalnew_mf[mfi].loVect()), 
                                AMREX_INT_ANYD(scalnew_mf[mfi].hiVect()),
                                BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(normal_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(etacart_mf[mfi]),
                                BL_TO_FORTRAN_ANYD(rho0_nph_cart_mf[mfi]));
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

    // put eta on base state edges
    // note that in spherical the base state has no refinement
    // the 0th value of etarho = 0, since U dot . e_r must be
    // zero at the center (since e_r is not defined there)
    etarho_edge[0] = 0.0;
    if (spherical == 1) {

        double dr1, dr2;
        for (int r=1; r<nr_fine; ++r) {
            dr1 = r_cc_loc[r] - r_edge_loc[r];
            dr2 = r_edge_loc[r] - r_cc_loc[r-1];
            etarho_edge[r] = (dr2*etarho_cell[r] + dr1*etarho_cell[r-1])/(dr1+dr2);
        }

    } else {
        for (int r=1; r<nr_fine; ++r) {
            etarho_edge[r] = 0.5*(etarho_cell[r] + etarho_cell[r-1]);
        }
    }
    // probably should do some better extrapolation here eventually
    etarho_edge[nr_fine] = etarho_cell[nr_fine-1];
}
