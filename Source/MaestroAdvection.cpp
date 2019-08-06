
#include <Maestro.H>

using namespace amrex;

// compute unprojected mac velocities
void
Maestro::AdvancePremac (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                        const Vector<std::array< MultiFab,AMREX_SPACEDIM > >& w0mac,
                        const RealVector& w0_force,
                        const Vector<MultiFab>& w0_force_cart)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::AdvancePremac()",AdvancePremac);

	// create a uold with filled ghost cells
	Vector<MultiFab> utilde(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		utilde[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
		utilde[lev].setVal(0.);
	}

	FillPatch(t_new, utilde, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);

	// create a MultiFab to hold uold + w0
	Vector<MultiFab>      ufull(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
		ufull[lev].setVal(0.);
	}

	// create ufull = uold + w0
	Put1dArrayOnCart(w0,ufull,1,1,bcs_u,0,1);
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
		for (int j=0; j < AMREX_SPACEDIM; j++)
			utrans[lev][j].setVal(0.);
	}

	// create utrans
	MakeUtrans(utilde,ufull,utrans,w0mac);

	// create a MultiFab to hold the velocity forcing
	Vector<MultiFab> vel_force(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		if (ppm_trace_forces == 0) {
			vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
		} else {
			// tracing needs more ghost cells
			vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
		}
		vel_force[lev].setVal(0.);

	}

	int do_add_utilde_force = 1;
	MakeVelForce(vel_force,utrans,sold,rho0_old,grav_cell_old,
	             w0_force,w0_force_cart,do_add_utilde_force);

	// add w0 to trans velocities
	Addw0 (utrans,w0mac,1.);

	VelPred(utilde,ufull,utrans,umac,w0mac,vel_force);
}

void
Maestro::MakeUtrans (const Vector<MultiFab>& utilde,
                     const Vector<MultiFab>& ufull,
                     Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans,
		     const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeUtrans()",MakeUtrans);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& utrans_mf  = utrans[lev][0];
              MultiFab& vtrans_mf  = utrans[lev][1];
#if (AMREX_SPACEDIM == 3)
              MultiFab& wtrans_mf  = utrans[lev][2];
	const MultiFab& w0macx_mf  = w0mac[lev][0];
	const MultiFab& w0macy_mf  = w0mac[lev][1];
	const MultiFab& w0macz_mf  = w0mac[lev][2];
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)

	// NOTE: don't tile, but threaded in fortran subroutine
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 2)
            mkutrans_2d(
#elif (AMREX_SPACEDIM == 3)
            mkutrans_3d(
#endif
                        &lev, domainBox.loVect(), domainBox.hiVect(),
                        tileBox.loVect(), tileBox.hiVect(),
                        BL_TO_FORTRAN_FAB(utilde_mf[mfi]), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_FAB(ufull_mf[mfi]), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_3D(utrans_mf[mfi]),
                        BL_TO_FORTRAN_3D(vtrans_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wtrans_mf[mfi]),
                        BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                        BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                        BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
#endif
                        w0.dataPtr(), dx, &dt, bcs_u[0].data(), phys_bc.dataPtr());

        } // end MFIter loop
    } // end loop over levels

    if (finest_level == 0) {
	// fill periodic ghost cells
	for (int lev=0; lev<=finest_level; ++lev) {
	    for (int d=0; d<AMREX_SPACEDIM; ++d) {
				utrans[lev][d].FillBoundary(geom[lev].periodicity());
	    }
	}

	// fill ghost cells behind physical boundaries
	FillUmacGhost(utrans);
    } else {
	// edge_restriction
	AverageDownFaces(utrans);

	// fill ghost cells for all levels
	FillPatchUedge(utrans);
    }

}

void
Maestro::VelPred (const Vector<MultiFab>& utilde,
                  const Vector<MultiFab>& ufull,
                  const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans,
                  Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
		  const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                  const Vector<MultiFab>& force)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelPred()",VelPred);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& umac_mf    = umac[lev][0];
        const MultiFab& utrans_mf  = utrans[lev][0];
        const MultiFab& vtrans_mf  = utrans[lev][1];
              MultiFab& vmac_mf    = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wtrans_mf  = utrans[lev][2];
              MultiFab& wmac_mf    = umac[lev][2];
	const MultiFab& w0macx_mf  = w0mac[lev][0];
	const MultiFab& w0macy_mf  = w0mac[lev][1];
	const MultiFab& w0macz_mf  = w0mac[lev][2];
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)

	// NOTE: don't tile, but threaded in fortran subroutine
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 2)
            velpred_2d(
#elif (AMREX_SPACEDIM == 3)
            velpred_3d(
#endif
                        &lev, domainBox.loVect(), domainBox.hiVect(),
                        tileBox.loVect(), tileBox.hiVect(),
                        BL_TO_FORTRAN_FAB(utilde_mf[mfi]), utilde_mf.nGrow(),
                        BL_TO_FORTRAN_FAB(ufull_mf[mfi]), ufull_mf.nGrow(),
                        BL_TO_FORTRAN_3D(utrans_mf[mfi]),
                        BL_TO_FORTRAN_3D(vtrans_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wtrans_mf[mfi]),
#endif
                        BL_TO_FORTRAN_3D(umac_mf[mfi]),
                        BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wmac_mf[mfi]),
			BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
			BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
			BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
#endif
                        BL_TO_FORTRAN_FAB(force_mf[mfi]), force_mf.nGrow(),
                        w0.dataPtr(), dx, &dt, bcs_u[0].data(), phys_bc.dataPtr());
        } // end MFIter loop
    } // end loop over levels

    // edge_restriction
    AverageDownFaces(umac);

}

void
Maestro::MakeEdgeScal (const Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<MultiFab>& force,
                       int is_vel, const Vector<BCRec>& bcs, int nbccomp,
                       int start_scomp, int start_bccomp, int num_comp, int is_conservative)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEdgeScal()",MakeEdgeScal);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Get the index space and grid spacing of the domain
        const Box& domainBox = geom[lev].Domain();
        const Real* dx = geom[lev].CellSize();

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf   = state[lev];
              MultiFab& sedgex_mf = sedge[lev][0];
        const MultiFab& umac_mf   = umac[lev][0];
              MultiFab& sedgey_mf = sedge[lev][1];
        const MultiFab& vmac_mf   = umac[lev][1];

#if (AMREX_SPACEDIM == 3)
              MultiFab& sedgez_mf = sedge[lev][2];
        const MultiFab& wmac_mf   = umac[lev][2];
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)

	// NOTE: don't tile, but threaded in fortran subroutine
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

	    // Be careful to pass in comp+1 for fortran indexing
            for (int scomp = start_scomp+1; scomp <= start_scomp + num_comp; ++scomp) {

                int bccomp = start_bccomp + scomp - start_scomp;

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 2)
                make_edge_scal_2d(
#elif (AMREX_SPACEDIM == 3)
                make_edge_scal_3d(
#endif
                    domainBox.loVect(), domainBox.hiVect(),
                    tileBox.loVect(), tileBox.hiVect(),
                    BL_TO_FORTRAN_FAB(scal_mf[mfi]), scal_mf.nGrow(),
                    BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
                    BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
#endif
                    BL_TO_FORTRAN_3D(umac_mf[mfi]),
                    BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                    umac_mf.nGrow(),
                    BL_TO_FORTRAN_FAB(force_mf[mfi]),
                    dx, &dt, &is_vel, bcs[0].data(),
                    &nbccomp, &scomp, &bccomp, &is_conservative);
            } // end loop over components
        } // end MFIter loop
    } // end loop over levels

    // We use edge_restriction for the output velocity if is_vel == 1
    // we do not use edge_restriction for scalars because instead we will use
    // reflux on the fluxes in make_flux.
    if (is_vel == 1) {
        if (reflux_type == 1 || reflux_type == 2) {
            AverageDownFaces(sedge);
        }
    }
}

void
Maestro::MakeRhoXFlux (const Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                       Vector<MultiFab>& etarhoflux,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                       const RealVector& r0_old,
                       const RealVector& r0_edge_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_old,
                       const RealVector& r0_new,
                       const RealVector& r0_edge_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_new,
                       const RealVector& r0_predicted_edge,
                       int start_comp, int num_comp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoXFlux()",MakeRhoXFlux);

    // Make sure to pass in comp+1 for fortran indexing
    const int startcomp = start_comp + 1;
    const int endcomp = startcomp + num_comp - 1;

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf        = state[lev];
              MultiFab& sedgex_mf      = sedge[lev][0];
	      MultiFab& sfluxx_mf      = sflux[lev][0];
	      MultiFab& etarhoflux_mf = etarhoflux[lev];
        const MultiFab& umac_mf        = umac[lev][0];

              MultiFab& sedgey_mf      = sedge[lev][1];
	      MultiFab& sfluxy_mf      = sflux[lev][1];
        const MultiFab& vmac_mf        = umac[lev][1];

#if (AMREX_SPACEDIM == 3)
              MultiFab& sedgez_mf      = sedge[lev][2];
	      MultiFab& sfluxz_mf      = sflux[lev][2];
        const MultiFab& wmac_mf        = umac[lev][2];

	// if spherical == 1
	const MultiFab& w0macx_mf = w0mac[lev][0];
	const MultiFab& w0macy_mf = w0mac[lev][1];
	const MultiFab& w0macz_mf = w0mac[lev][2];
	MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;
	rho0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
	rho0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
	rho0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);

	if (spherical == 1) {
	    MultiFab::LinComb(rho0mac_edgex,0.5,r0mac_old[lev][0],0,0.5,r0mac_new[lev][0],0,0,1,1);
	    MultiFab::LinComb(rho0mac_edgey,0.5,r0mac_old[lev][1],0,0.5,r0mac_new[lev][1],0,0,1,1);
	    MultiFab::LinComb(rho0mac_edgez,0.5,r0mac_old[lev][2],0,0.5,r0mac_new[lev][2],0,0,1,1);
	}
#endif


        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		// NOTE: not sure if this should be tiled
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

	    // call fortran subroutine
	    // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
	    // lo/hi coordinates (including ghost cells), and/or the # of components
	    // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 0) {

#if (AMREX_SPACEDIM == 2)
	        make_rhoX_flux_2d(
#elif (AMREX_SPACEDIM == 3)
                make_rhoX_flux_3d(
#endif
				  &lev, tileBox.loVect(), tileBox.hiVect(),
				  BL_TO_FORTRAN_FAB(sfluxx_mf[mfi]),
				  BL_TO_FORTRAN_FAB(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				  BL_TO_FORTRAN_FAB(sfluxz_mf[mfi]),
#endif
				  BL_TO_FORTRAN_3D(etarhoflux_mf[mfi]),
				  BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
				  BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				  BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
#endif
				  BL_TO_FORTRAN_3D(umac_mf[mfi]),
				  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
				  r0_old.dataPtr(), r0_edge_old.dataPtr(),
				  r0_new.dataPtr(), r0_edge_new.dataPtr(),
				  r0_predicted_edge.dataPtr(),
				  w0.dataPtr(),
				  &startcomp, &endcomp);
	    } else {

#if (AMREX_SPACEDIM == 3)
		// if (use_exact_base_state)
		// {
		//     // add make_rhoX_flux_3d_sphr_irreg()
		// }
		// else
		// {
		    make_rhoX_flux_3d_sphr(tileBox.loVect(), tileBox.hiVect(),
			               BL_TO_FORTRAN_FAB(sfluxx_mf[mfi]),
			               BL_TO_FORTRAN_FAB(sfluxy_mf[mfi]),
			               BL_TO_FORTRAN_FAB(sfluxz_mf[mfi]),
				       BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
				       BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
				       BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
				       BL_TO_FORTRAN_3D(umac_mf[mfi]),
				       BL_TO_FORTRAN_3D(vmac_mf[mfi]),
				       BL_TO_FORTRAN_3D(wmac_mf[mfi]),
				       BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
				       BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
				       BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
				       BL_TO_FORTRAN_3D(rho0mac_edgex[mfi]),
				       BL_TO_FORTRAN_3D(rho0mac_edgey[mfi]),
				       BL_TO_FORTRAN_3D(rho0mac_edgez[mfi]),
				       &startcomp, &endcomp);
		// }
#else
	        Abort("MakeRhoXFlux: Spherical is not valid for DIM < 3");
#endif
	    } // end spherical
	} // end MFIter loop


	// increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {

	    // Get the grid size
	    const Real* dx = geom[lev].CellSize();
	    // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
	    const Real area[3] = {dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1]};
#else
	    const Real area[2] = {dx[1], dx[0]};
#endif

	    if (flux_reg_s[lev+1])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,start_comp,start_comp,num_comp, -1.0*dt*area[i]);
		    // also include density flux
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,Rho,Rho,1, -1.0*dt*area[i]);
                }
            }
	    if (flux_reg_s[lev])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,start_comp,start_comp,num_comp, 1.0*dt*area[i]);
		    // also include density flux
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,Rho,Rho,1, 1.0*dt*area[i]);
                }
            }

	    if (spherical == 0) {
		// need edge_restrict for etarhoflux
	    }
        }

    } // end loop over levels

    // average down fluxes
    if (reflux_type == 1) {
	AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()
}

void
Maestro::MakeRhoHFlux (const Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                       const RealVector& r0_old,
                       const RealVector& r0_edge_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_old,
                       const RealVector& r0_new,
                       const RealVector& r0_edge_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_new,
                       const RealVector& rh0_old,
                       const RealVector& rh0_edge_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& rh0mac_old,
                       const RealVector& rh0_new,
                       const RealVector& rh0_edge_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& rh0mac_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& h0mac_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& h0mac_new)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoHFlux()",MakeRhoHFlux);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf   = state[lev];
              MultiFab& sedgex_mf = sedge[lev][0];
	      MultiFab& sfluxx_mf = sflux[lev][0];
        const MultiFab& umac_mf   = umac[lev][0];

              MultiFab& sedgey_mf = sedge[lev][1];
	      MultiFab& sfluxy_mf = sflux[lev][1];
        const MultiFab& vmac_mf   = umac[lev][1];

#if (AMREX_SPACEDIM == 3)
              MultiFab& sedgez_mf = sedge[lev][2];
	      MultiFab& sfluxz_mf = sflux[lev][2];
        const MultiFab& wmac_mf   = umac[lev][2];

	// if spherical == 1
	const MultiFab& w0macx_mf = w0mac[lev][0];
	const MultiFab& w0macy_mf = w0mac[lev][1];
	const MultiFab& w0macz_mf = w0mac[lev][2];
	MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;
	MultiFab h0mac_edgex, h0mac_edgey, h0mac_edgez;
	MultiFab rhoh0mac_edgex, rhoh0mac_edgey, rhoh0mac_edgez;
	rho0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);
	rho0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);
	rho0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0);
	h0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);
	h0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);
	h0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0);
	rhoh0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);
	rhoh0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);
	rhoh0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0);

	if (spherical == 1) {
	    if (use_exact_base_state) {
		MultiFab::LinComb(rhoh0mac_edgex,0.5,rh0mac_old[lev][0],0,0.5,rh0mac_new[lev][0],0,0,1,0);
		MultiFab::LinComb(rhoh0mac_edgey,0.5,rh0mac_old[lev][1],0,0.5,rh0mac_new[lev][1],0,0,1,0);
		MultiFab::LinComb(rhoh0mac_edgez,0.5,rh0mac_old[lev][2],0,0.5,rh0mac_new[lev][2],0,0,1,0);
	    } else {
		MultiFab::LinComb(rho0mac_edgex,0.5,r0mac_old[lev][0],0,0.5,r0mac_new[lev][0],0,0,1,0);
		MultiFab::LinComb(rho0mac_edgey,0.5,r0mac_old[lev][1],0,0.5,r0mac_new[lev][1],0,0,1,0);
		MultiFab::LinComb(rho0mac_edgez,0.5,r0mac_old[lev][2],0,0.5,r0mac_new[lev][2],0,0,1,0);
		MultiFab::LinComb(h0mac_edgex,0.5,h0mac_old[lev][0],0,0.5,h0mac_new[lev][0],0,0,1,0);
		MultiFab::LinComb(h0mac_edgey,0.5,h0mac_old[lev][1],0,0.5,h0mac_new[lev][1],0,0,1,0);
		MultiFab::LinComb(h0mac_edgez,0.5,h0mac_old[lev][2],0,0.5,h0mac_new[lev][2],0,0,1,0);
	    }
	}
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
	        // NOTE: not sure if this should be tiled (see make_rhoX_flux above)
        for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

	    // call fortran subroutine
	    // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
	    // lo/hi coordinates (including ghost cells), and/or the # of components
	    // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 0) {
#if (AMREX_SPACEDIM == 2)
		make_rhoh_flux_2d(
#elif (AMREX_SPACEDIM == 3)
                make_rhoh_flux_3d(
#endif
				  &lev, tileBox.loVect(), tileBox.hiVect(),
				  BL_TO_FORTRAN_FAB(sfluxx_mf[mfi]),
				  BL_TO_FORTRAN_FAB(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				  BL_TO_FORTRAN_FAB(sfluxz_mf[mfi]),
#endif
				  BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
				  BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				  BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
#endif
				  BL_TO_FORTRAN_3D(umac_mf[mfi]),
				  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
				  r0_old.dataPtr(), r0_edge_old.dataPtr(),
				  r0_new.dataPtr(), r0_edge_new.dataPtr(),
				  rh0_old.dataPtr(), rh0_edge_old.dataPtr(),
				  rh0_new.dataPtr(), rh0_edge_new.dataPtr(),
				  w0.dataPtr());
	    } else {

#if (AMREX_SPACEDIM == 3)
	        if (use_exact_base_state)
		{
		    make_rhoh_flux_3d_sphr_irreg(tileBox.loVect(), tileBox.hiVect(),
						 BL_TO_FORTRAN_FAB(sfluxx_mf[mfi]),
						 BL_TO_FORTRAN_FAB(sfluxy_mf[mfi]),
						 BL_TO_FORTRAN_FAB(sfluxz_mf[mfi]),
						 BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
						 BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
						 BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
						 BL_TO_FORTRAN_3D(umac_mf[mfi]),
						 BL_TO_FORTRAN_3D(vmac_mf[mfi]),
						 BL_TO_FORTRAN_3D(wmac_mf[mfi]),
						 BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
						 BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
						 BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
						 BL_TO_FORTRAN_3D(rhoh0mac_edgex[mfi]),
						 BL_TO_FORTRAN_3D(rhoh0mac_edgey[mfi]),
						 BL_TO_FORTRAN_3D(rhoh0mac_edgez[mfi]));
		}
		else
		{
		    make_rhoh_flux_3d_sphr(tileBox.loVect(), tileBox.hiVect(),
			               BL_TO_FORTRAN_FAB(sfluxx_mf[mfi]),
			               BL_TO_FORTRAN_FAB(sfluxy_mf[mfi]),
			               BL_TO_FORTRAN_FAB(sfluxz_mf[mfi]),
				       BL_TO_FORTRAN_FAB(sedgex_mf[mfi]),
				       BL_TO_FORTRAN_FAB(sedgey_mf[mfi]),
				       BL_TO_FORTRAN_FAB(sedgez_mf[mfi]),
				       BL_TO_FORTRAN_3D(umac_mf[mfi]),
				       BL_TO_FORTRAN_3D(vmac_mf[mfi]),
				       BL_TO_FORTRAN_3D(wmac_mf[mfi]),
				       BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
				       BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
				       BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
				       BL_TO_FORTRAN_3D(rho0mac_edgex[mfi]),
				       BL_TO_FORTRAN_3D(rho0mac_edgey[mfi]),
				       BL_TO_FORTRAN_3D(rho0mac_edgez[mfi]),
				       BL_TO_FORTRAN_3D(h0mac_edgex[mfi]),
				       BL_TO_FORTRAN_3D(h0mac_edgey[mfi]),
				       BL_TO_FORTRAN_3D(h0mac_edgez[mfi]));
		}
#else
	        Abort("MakeRhoHFlux: Spherical is not valid for DIM < 3");
#endif
	    }
	} // end MFIter loop


        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {

	    // Get the grid size
	    const Real* dx = geom[lev].CellSize();
            	    // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
	    const Real area[3] = {dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1]};
#else
	    const Real area[2] = {dx[1], dx[0]};
#endif

	    if (flux_reg_s[lev+1])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,RhoH,RhoH,1, -1.0*dt*area[i]);
                }
            }
	    if (flux_reg_s[lev])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,RhoH,RhoH,1, 1.0*dt*area[i]);
                }
            }
        }
    } // end loop over levels

    if (reflux_type == 1) {
	AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()

}

void
Maestro::UpdateScal(const Vector<MultiFab>& stateold,
                    Vector<MultiFab>& statenew,
                    const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                    const Vector<MultiFab>& force,
                    int start_comp, int num_comp,
                    const Real* p0, const Vector<MultiFab>& p0_cart)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateScal()",UpdateScal);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    Gpu::setLaunchRegion(true);
#endif

    // Make sure to pass in comp+1 for fortran indexing
    const int startcomp = start_comp + 1;
    const int endcomp = startcomp + num_comp;

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scalold_mf = stateold[lev];
	      MultiFab& scalnew_mf = statenew[lev];
        const MultiFab& sfluxx_mf = sflux[lev][0];
        const MultiFab& sfluxy_mf = sflux[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& sfluxz_mf = sflux[lev][2];

    	MultiFab p0_cart_dummy;
    	p0_cart_dummy.define(grids[lev], dmap[lev], 1, 1);
    	if (start_comp == RhoH && spherical == 1) {
    	    MultiFab::Copy(p0_cart_dummy,p0_cart[lev],0,0,1,1);
    	}
    	const MultiFab& p0cart_mf = p0_cart_dummy;
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(scalold_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

	    if (start_comp == RhoH)
	    {   // Enthalpy update

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
		if (spherical == 0) {
#pragma gpu box(tileBox)
            update_rhoh(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
                        lev,
    				    BL_TO_FORTRAN_ANYD(scalold_mf[mfi]),
    				    BL_TO_FORTRAN_ANYD(scalnew_mf[mfi]),
    				    BL_TO_FORTRAN_ANYD(sfluxx_mf[mfi]),
    				    BL_TO_FORTRAN_ANYD(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
    				    BL_TO_FORTRAN_ANYD(sfluxz_mf[mfi]),
#endif
    				    BL_TO_FORTRAN_ANYD(force_mf[mfi]),
    				    p0,
    				    AMREX_REAL_ANYD(dx), dt,
    				    NumSpec);
		} else {
#if (AMREX_SPACEDIM == 3)
#pragma gpu box(tileBox)
		    update_rhoh_3d_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
            					BL_TO_FORTRAN_ANYD(scalold_mf[mfi]),
            					BL_TO_FORTRAN_ANYD(scalnew_mf[mfi]),
            					BL_TO_FORTRAN_ANYD(sfluxx_mf[mfi]),
            					BL_TO_FORTRAN_ANYD(sfluxy_mf[mfi]),
            					BL_TO_FORTRAN_ANYD(sfluxz_mf[mfi]),
            					BL_TO_FORTRAN_ANYD(force_mf[mfi]),
            					BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]),
            					AMREX_REAL_ANYD(dx), dt,
            					NumSpec);
#else
		    Abort("UpdateScal: Spherical is not valid for DIM < 3");
#endif
		}

	    }
	    else if (start_comp == FirstSpec)
	    {   // RhoX update

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
                update_rhoX(AMREX_INT_ANYD(tileBox.loVect()),
                    AMREX_INT_ANYD(tileBox.hiVect()),
                    BL_TO_FORTRAN_ANYD(scalold_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(scalnew_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(sfluxx_mf[mfi]),
                    BL_TO_FORTRAN_ANYD(sfluxy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_ANYD(sfluxz_mf[mfi]),
#endif
                    BL_TO_FORTRAN_ANYD(force_mf[mfi]),
                    AMREX_REAL_ANYD(dx), dt,
                    startcomp, endcomp);
            }
	    else {
		Abort("Invalid scalar in UpdateScal().");
	    } // end if
	} // end MFIter loop
    } // end loop over levels

    // synchronize by refluxing and averaging down, starting from the finest_level-1/finest_level pair
    if (reflux_type == 2) {
        for (int lev=finest_level-1; lev>=0; --lev) {
            // update lev based on coarse-fine flux mismatch
            flux_reg_s[lev+1]->Reflux(statenew[lev], 1.0, start_comp, start_comp, num_comp, geom[lev]);
            if (start_comp == FirstSpec) {
                // do the same for density if we updated the species
                flux_reg_s[lev+1]->Reflux(statenew[lev], 1.0, Rho, Rho, 1, geom[lev]);
            }
        }
    }

    // average fine data onto coarser cells
    // fill ghost cells
    AverageDown(statenew,start_comp,num_comp);
    FillPatch(t_old, statenew, statenew, statenew, start_comp, start_comp, num_comp, start_comp, bcs_s);

    // do the same for density if we updated the species
    if (start_comp == FirstSpec) {
	AverageDown(statenew,Rho,1);
	FillPatch(t_old, statenew, statenew, statenew, Rho, Rho, 1, Rho, bcs_s);
    }

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::UpdateVel (const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                    const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                    const Vector<MultiFab>& force,
                    const Vector<MultiFab>& sponge,
                    const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateVel()",UpdateVel);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& uold_mf = uold[lev];
	      MultiFab& unew_mf = unew[lev];
	const MultiFab& umac_mf   = umac[lev][0];
        const MultiFab& uedgex_mf = uedge[lev][0];
	const MultiFab& vmac_mf   = umac[lev][1];
        const MultiFab& uedgey_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
	const MultiFab& wmac_mf   = umac[lev][2];
        const MultiFab& uedgez_mf = uedge[lev][2];

	// if spherical == 1
	const MultiFab& w0macx_mf = w0mac[lev][0];
	const MultiFab& w0macy_mf = w0mac[lev][1];
	const MultiFab& w0macz_mf = w0mac[lev][2];
#endif
        const MultiFab& force_mf = force[lev];
	const MultiFab& sponge_mf = sponge[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(force_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

	    if (spherical == 0) {
#pragma gpu box(tileBox)
		update_velocity(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()), lev,
				 BL_TO_FORTRAN_ANYD(uold_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(unew_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				 BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
#endif
				 BL_TO_FORTRAN_ANYD(uedgex_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(uedgey_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				 BL_TO_FORTRAN_ANYD(uedgez_mf[mfi]),
#endif
				 BL_TO_FORTRAN_ANYD(force_mf[mfi]),
				 BL_TO_FORTRAN_ANYD(sponge_mf[mfi]),
				 w0.dataPtr(),
				 AMREX_REAL_ANYD(dx), dt);
	    } else {
#if (AMREX_SPACEDIM == 3)
#pragma gpu box(tileBox)
		update_velocity_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
				      BL_TO_FORTRAN_ANYD(uold_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(unew_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(umac_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(vmac_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(wmac_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(uedgex_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(uedgey_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(uedgez_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(force_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(sponge_mf[mfi]),
				      w0.dataPtr(),
				      BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
				      BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
				      AMREX_REAL_ANYD(dx), dt);
#else
		Abort("UpdateVel: Spherical is not valid for DIM < 3");
#endif
	    }
        } // end MFIter loop
    } // end loop over levels

    // average fine data onto coarser cells
    AverageDown(unew,0,AMREX_SPACEDIM);

    // fill ghost cells
    FillPatch(t_old, unew, unew, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u, 1);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif

}
