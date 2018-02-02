
#include <Maestro.H>

using namespace amrex;

// compute unprojected mac velocities
void
Maestro::AdvancePremac (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                        const Vector<Real>& w0_force)
{

    // create a uold with filled ghost cells
    Vector<MultiFab> utilde(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        utilde[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
    }

    FillPatch(t_new, utilde, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

    // create a MultiFab to hold uold + w0
    Vector<MultiFab>      ufull(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
    }

    // create ufull = uold + w0
    Put1dArrayOnCart(w0,ufull,1,1,bcs_u,0);
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Add(ufull[lev],utilde[lev],0,0,AMREX_SPACEDIM,ng_adv);
    }

    // create a face-centered MultiFab to hold utrans
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > utrans(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        utrans[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        utrans[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
        utrans[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
    }

    // create utrans
    MakeUtrans(utilde,ufull,utrans);
    
    // create a MultiFab to hold the velocity forcing
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    int do_add_utilde_force = 1;
    MakeVelForce(vel_force,utrans,sold,rho0_old,grav_cell_old,w0_force,do_add_utilde_force);

    // add w0 to trans velocities
    Addw0 (utrans,1.);

    VelPred(utilde,ufull,utrans,umac,vel_force);
}

void
Maestro::MakeUtrans (const Vector<MultiFab>& utilde,
                     const Vector<MultiFab>& ufull,
                     Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& utrans_mf  = utrans[lev][0];
#if (AMREX_SPACEDIM >= 2)
              MultiFab& vtrans_mf  = utrans[lev][1];
#if (AMREX_SPACEDIM == 3)
              MultiFab& wtrans_mf  = utrans[lev][2];
#endif
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
            mkutrans_1d(
#elif (AMREX_SPACEDIM == 2)
            mkutrans_2d(
#elif (AMREX_SPACEDIM == 3)
            mkutrans_3d(
#endif
                        &lev, domainBox.loVect(), domainBox.hiVect(),
                        validBox.loVect(), validBox.hiVect(),
                        BL_TO_FORTRAN_FAB(utilde_mf[mfi]), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_FAB(ufull_mf[mfi]),
                        BL_TO_FORTRAN_3D(utrans_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
                        BL_TO_FORTRAN_3D(vtrans_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wtrans_mf[mfi]),
#endif
#endif
                        w0.dataPtr(), dx, &dt, bcs_u[0].data(), phys_bc.dataPtr());

        } // end MFIter loop
    } // end loop over levels
    
    // fill peroidic ghost cells
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            utrans[lev][d].FillBoundary(geom[lev].periodicity());
        }
    }

    // fill ghost cells behind physical boundaries
    FillUmacGhost(utrans);

    // FIXME need to add edge_restriction and create_umac_grown
    //
    //
}

void
Maestro::VelPred (const Vector<MultiFab>& utilde,
                  const Vector<MultiFab>& ufull,
                  const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans,
                  Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                  const Vector<MultiFab>& force)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& umac_mf    = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
        const MultiFab& utrans_mf  = utrans[lev][0];
        const MultiFab& vtrans_mf  = utrans[lev][1];
              MultiFab& vmac_mf    = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wtrans_mf  = utrans[lev][2];
              MultiFab& wmac_mf    = umac[lev][2];
#endif
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
            velpred_1d(
#elif (AMREX_SPACEDIM == 2)
            velpred_2d(
#elif (AMREX_SPACEDIM == 3)
            velpred_3d(
#endif
                        &lev, domainBox.loVect(), domainBox.hiVect(),
                        validBox.loVect(), validBox.hiVect(),
                        BL_TO_FORTRAN_FAB(utilde_mf[mfi]), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_FAB(ufull_mf[mfi]),
                        BL_TO_FORTRAN_3D(utrans_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
                        BL_TO_FORTRAN_3D(vtrans_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wtrans_mf[mfi]),
#endif
#endif
                        BL_TO_FORTRAN_3D(umac_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
                        BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
#endif
                        BL_TO_FORTRAN_FAB(force_mf[mfi]),
                        w0.dataPtr(), dx, &dt, bcs_u[0].data(), phys_bc.dataPtr());
        } // end MFIter loop
    } // end loop over levels

    // FIXME need to add edge_restriction
    //
    //

}

void
    Maestro::MakeEdgeScal (const Vector<MultiFab>& state,
			   Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
			   const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
			   const Vector<MultiFab>& force,
			   int is_vel, const Vector<BCRec>& bcs, int nbccomp, 
			   int start_scomp, int start_bccomp, int num_comp, int is_conservative)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf   = state[lev];
              MultiFab& sedgex_mf = sedge[lev][0];
        const MultiFab& umac_mf   = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
              MultiFab& sedgey_mf = sedge[lev][1];
        const MultiFab& vmac_mf   = umac[lev][1];

#if (AMREX_SPACEDIM == 3)
              MultiFab& sedgez_mf = sedge[lev][2];
        const MultiFab& wmac_mf   = umac[lev][2];

#endif
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();
            const Real* dx = geom[lev].CellSize();

	    // Be careful to pass in comp+1 for fortran indexing
            for (int scomp = start_scomp+1; scomp <= start_scomp + num_comp; ++scomp) {

                int bccomp = start_bccomp + scomp - start_scomp;

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
                make_edge_scal_1d(
#elif (AMREX_SPACEDIM == 2)
                make_edge_scal_2d(
#elif (AMREX_SPACEDIM == 3)
                make_edge_scal_3d(
#endif
                    domainBox.loVect(), domainBox.hiVect(),
                    validBox.loVect(), validBox.hiVect(),
                    BL_TO_FORTRAN_FAB(scal_mf[mfi]), scal_mf.nGrow(),
                    BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
                    BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
#endif
#endif
                    BL_TO_FORTRAN_3D(umac_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
                    BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
#endif
                    BL_TO_FORTRAN_FAB(force_mf[mfi]),
                    dx, &dt, &is_vel, bcs[0].data(),
                    &nbccomp, &scomp, &bccomp, &is_conservative);
            } // end loop over components 
        } // end MFIter loop
    } // end loop over levels

    // FIXME need to add edge_restriction
    //
    //

}

void
    Maestro::MakeRhoXFlux (const Vector<MultiFab>& state,
			   Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux, 
			   Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
			   const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
			   const Vector<Real>& r0_old,
			   const Vector<Real>& r0_edge_old,
			   const Vector<Real>& r0_new,
			   const Vector<Real>& r0_edge_new,
			   int start_comp, int num_comp)
{
    // Make sure to pass in comp+1 for fortran indexing
    const int startcomp = start_comp + 1;
    const int endcomp = startcomp + num_comp;

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf   = state[lev];
              MultiFab& sedgex_mf = sedge[lev][0];
	      MultiFab& sfluxx_mf = sflux[lev][0];
        const MultiFab& umac_mf   = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
              MultiFab& sedgey_mf = sedge[lev][1];
	      MultiFab& sfluxy_mf = sflux[lev][1];
        const MultiFab& vmac_mf   = umac[lev][1];

#if (AMREX_SPACEDIM == 3)
              MultiFab& sedgez_mf = sedge[lev][2];
	      MultiFab& sfluxz_mf = sflux[lev][2];
        const MultiFab& wmac_mf   = umac[lev][2];

#endif
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();

	    // call fortran subroutine
	    // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
	    // lo/hi coordinates (including ghost cells), and/or the # of components
	    // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
	    make_rhoX_flux_1d(
#elif (AMREX_SPACEDIM == 2)
	    make_rhoX_flux_2d(
#elif (AMREX_SPACEDIM == 3)
            make_rhoX_flux_3d(
#endif
			      &lev, domainBox.loVect(), domainBox.hiVect(),
			      validBox.loVect(), validBox.hiVect(),
			      BL_TO_FORTRAN_FAB(sfluxx_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
			      BL_TO_FORTRAN_FAB(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
			      BL_TO_FORTRAN_FAB(sfluxz_mf[mfi]),
#endif
#endif
			      BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
			      BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
			      BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
#endif
#endif
			      BL_TO_FORTRAN_3D(umac_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
			      BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
			      BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
#endif
			      r0_old.dataPtr(), r0_edge_old.dataPtr(), 
			      r0_new.dataPtr(), r0_edge_new.dataPtr(),
			      w0.dataPtr(), 
			      &startcomp, &endcomp);
	    
	} // end MFIter loop
    } // end loop over levels

    // FIXME need to add edge_restriction
    //
    //

}
