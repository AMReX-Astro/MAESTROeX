
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
    RealVector etarhosum( (nr_fine+1)*(max_radial_level+1), 0.0 );
    etarhosum.shrink_to_fit();

    // this stores how many cells there are laterally at each level
    RealVector ncell(max_radial_level+1);

    const int max_lev = max_radial_level;
    get_r_end_coord(r_end_coord.dataPtr());

    Real * AMREX_RESTRICT etarhosum_p = etarhosum.dataPtr();
    const int* AMREX_RESTRICT r_end_coord_p = r_end_coord.dataPtr();

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
        
        const int nchunks = numdisjointchunks[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid tile region
            const Box& tilebox = mfi.tilebox();
            const Box& nbx = mfi.nodaltilebox(AMREX_SPACEDIM-1);

            const Array4<const Real> etarhoflux = etarho_flux[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                int r = AMREX_SPACEDIM == 2 ? j : k;
                etarhosum_p[lev+r*(max_lev+1)] += etarhoflux(i,j,k);
            });

            // This has large errors so going to leave it commented out for now
            // AMREX_PARALLEL_FOR_3D(nbx, i, j, k, {
            //     int r = AMREX_SPACEDIM == 2 ? j : k;
            //     for (int n=0; n<nchunks; n++) {
            //         if (r-1 == r_end_coord_p[n+lev*(nr_fine+1)]) {
            //             etarhosum_p[lev+r*(max_lev+1)] += etarhoflux(i,j,k);
            //         }
            //     }
            // });

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells)
            // We will also pass "bx", which specifies the "valid" tile region.
            sum_etarho(&lev, ARLIM_3D(domainBox.loVect()), ARLIM_3D(domainBox.hiVect()),
                       ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
                       BL_TO_FORTRAN_3D(etarhoflux_mf[mfi]),
                       etarhosum.dataPtr());
        }
    }

    ParallelDescriptor::ReduceRealSum(etarhosum.dataPtr(),(nr_fine+1)*(max_radial_level+1));

    make_etarho_planar(etarho_edge.dataPtr(), etarho_cell.dataPtr(),
                       etarhosum.dataPtr(), ncell.dataPtr());
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

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scalold_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tilebox = mfi.tilebox();

            const Array4<const Real> rho_old = scal_old[lev].array(mfi, Rho);
            const Array4<const Real> rho_new = scal_new[lev].array(mfi, Rho);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
            const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
            const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
            const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
            const Array4<const Real> normal_arr = normal[lev].array(mfi);
            const Array4<Real> etacart_arr = eta_cart[lev].array(mfi);
            const Array4<const Real> rho0_nph_arr = rho0_nph_cart[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
                Real U_dot_er = 0.5*(umacx(i,j,k) + umacx(i+1,j,k) 
                    + w0macx(i,j,k) + w0macx(i+1,j,k)) * normal_arr(i,j,k,0) + 
                    0.5*(vmac(i,j,k) + vmac(i,j+1,k) 
                    + w0macy(i,j,k) + w0macy(i,j+1,k)) * normal_arr(i,j,k,1) + 
                    0.5*(wmac(i,j,k) + wmac(i,j,k+1) 
                    + w0macz(i,j,k) + w0macz(i,j,k+1)) * normal_arr(i,j,k,2);

                // construct time-centered [ rho' (U dot e_r) ]
                etacart_arr(i,j,k) = (0.5*(rho_old(i,j,k) + rho_new(i,j,k)) - 
                    rho0_nph_arr(i,j,k)) * U_dot_er;
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
