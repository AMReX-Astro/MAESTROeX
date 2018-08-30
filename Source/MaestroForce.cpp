
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeVelForce (Vector<MultiFab>& vel_force,
                       int is_final_update,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                       const Vector<MultiFab>& rho,
                       const Vector<Real>& rho0,
                       const Vector<Real>& grav_cell,
                       const Vector<Real>& w0_force,
                       const Vector<MultiFab>& w0_force_cart,
#ifdef ROTATION
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
#endif
                       int do_add_utilde_force)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeVelForce()",MakeVelForce);

	// TODO: how do I properly do the w0_cart thing?

	// For spherical case
	Vector<MultiFab> w0_cart(finest_level+1);
	Vector<MultiFab> gradw0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		w0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
		w0_cart[lev].setVal(0.);
		gradw0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
		gradw0_cart[lev].setVal(0.);
	}

	if (spherical == 1) {
		Vector<Real> gradw0( (max_radial_level+1)*nr_fine );
		gradw0.shrink_to_fit();

		compute_grad_phi_rad(w0.dataPtr(), gradw0.dataPtr());

		Put1dArrayOnCart(gradw0,gradw0_cart,0,0,bcs_f,0);

#ifdef ROTATION
		Put1dArrayOnCart(w0,w0_cart,1,1,bcs_f,0);
#endif
	}

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		MultiFab& vel_force_mf = vel_force[lev];
		const MultiFab& gpi_mf = gpi[lev];
		const MultiFab& rho_mf = rho[lev];
		const MultiFab& uedge_mf = uedge[lev][0];
		const MultiFab& vedge_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
		const MultiFab& wedge_mf = uedge[lev][2];
		const MultiFab& gradw0_mf = gradw0_cart[lev];
		const MultiFab& normal_mf = normal[lev];
		const MultiFab& w0force_mf = w0_force_cart[lev];
#ifdef ROTATION
		const MultiFab& w0_mf = w0_cart[lev];
		const MultiFab& w0macx_mf = w0mac[lev][0];
		const MultiFab& w0macy_mf = w0mac[lev][1];
		const MultiFab& uold_mf = uold[lev];
#endif
#endif
		const MultiFab& cc_to_r = cell_cc_to_r[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		for ( MFIter mfi(vel_force_mf); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& validBox = mfi.validbox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			if (spherical == 0) {
				make_vel_force(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				               &is_final_update,
				               BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
				               BL_TO_FORTRAN_FAB(gpi_mf[mfi]),
				               BL_TO_FORTRAN_N_3D(rho_mf[mfi],Rho),
				               BL_TO_FORTRAN_3D(uedge_mf[mfi]),
				               BL_TO_FORTRAN_3D(vedge_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				               BL_TO_FORTRAN_3D(wedge_mf[mfi]),
#ifdef ROTATION
				               BL_TO_FORTRAN_FAB(uold_mf[mfi]),
#endif
#endif
				               w0.dataPtr(),
				               w0_force.dataPtr(),
				               rho0.dataPtr(),
				               grav_cell.dataPtr(),
				               &do_add_utilde_force);
			} else {

#if (AMREX_SPACEDIM == 3)
				make_vel_force_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				                    &is_final_update,
				                    BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
				                    BL_TO_FORTRAN_FAB(gpi_mf[mfi]),
				                    BL_TO_FORTRAN_N_3D(rho_mf[mfi],Rho),
				                    BL_TO_FORTRAN_3D(uedge_mf[mfi]),
				                    BL_TO_FORTRAN_3D(vedge_mf[mfi]),
				                    BL_TO_FORTRAN_3D(wedge_mf[mfi]),
				                    BL_TO_FORTRAN_FAB(normal_mf[mfi]),
				                    BL_TO_FORTRAN_3D(gradw0_mf[mfi]),
				                    BL_TO_FORTRAN_FAB(w0force_mf[mfi]),
#ifdef ROTATION
				                    BL_TO_FORTRAN_FAB(w0_mf[mfi]),
				                    BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
				                    BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
				                    BL_TO_FORTRAN_FAB(uold_mf[mfi]),
#endif
				                    rho0.dataPtr(),
				                    grav_cell.dataPtr(),
				                    dx,
				                    r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
				                    BL_TO_FORTRAN_3D(cc_to_r[mfi]),
				                    &do_add_utilde_force);
#else
				Abort("MakeVelForce: Spherical is not valid for DIM < 3");
#endif
			}
		}
	}

	// average fine data onto coarser cells
	AverageDown(vel_force,0,AMREX_SPACEDIM);

	// note - we need to reconsider the bcs type here
	// it matches fortran MAESTRO but is that correct?
	FillPatch(t_old, vel_force, vel_force, vel_force, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

}


void
Maestro::ModifyScalForce(Vector<MultiFab>& scal_force,
                         const Vector<MultiFab>& state,
                         const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                         const Vector<Real>& s0,
                         const Vector<Real>& s0_edge,
                         const Vector<MultiFab>& s0_cart,
                         int comp,
                         const Vector<BCRec>& bcs,
                         int fullform)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::ModifyScalForce()",ModifyScalForce);

	for (int lev=0; lev<=finest_level; ++lev) {

		// Get the index space and grid spacing of the domain
		const Box& domainBox = geom[lev].Domain();
		const Real* dx = geom[lev].CellSize();

		// get references to the MultiFabs at level lev
		MultiFab& scal_force_mf = scal_force[lev];
		const MultiFab& state_mf = state[lev];
		const MultiFab& umac_mf = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
		const MultiFab& vmac_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
		const MultiFab& wmac_mf = umac[lev][2];
		const MultiFab& s0cart_mf = s0_cart[lev];
		const MultiFab& cc_to_r = cell_cc_to_r[lev];
#endif
#endif

		// loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		for ( MFIter mfi(scal_force_mf); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& validBox = mfi.validbox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			if (spherical == 1) {
#if (AMREX_SPACEDIM == 3)
				modify_scal_force_sphr(domainBox.loVect(), domainBox.hiVect(),
				                       validBox.loVect(), validBox.hiVect(),
				                       scal_force_mf[mfi].dataPtr(comp),
				                       scal_force_mf[mfi].loVect(), scal_force_mf[mfi].hiVect(),
				                       state_mf[mfi].dataPtr(comp),
				                       state_mf[mfi].loVect(), state_mf[mfi].hiVect(),
				                       BL_TO_FORTRAN_3D(umac_mf[mfi]),
				                       BL_TO_FORTRAN_3D(vmac_mf[mfi]),
				                       BL_TO_FORTRAN_3D(wmac_mf[mfi]),
				                       BL_TO_FORTRAN_3D(s0cart_mf[mfi]),
				                       w0.dataPtr(), dx, &fullform,
				                       r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
				                       BL_TO_FORTRAN_3D(cc_to_r[mfi]));
#else
				Abort("ModifyScalForce: Spherical is not valid for DIM < 3");
#endif
			} else {
				modify_scal_force(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				                  scal_force_mf[mfi].dataPtr(comp),
				                  ARLIM_3D(scal_force_mf[mfi].loVect()), ARLIM_3D(scal_force_mf[mfi].hiVect()),
				                  state_mf[mfi].dataPtr(comp),
				                  ARLIM_3D(state_mf[mfi].loVect()), ARLIM_3D(state_mf[mfi].hiVect()),
				                  BL_TO_FORTRAN_3D(umac_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
				                  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
				                  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
#endif
				                  s0.dataPtr(), s0_edge.dataPtr(), w0.dataPtr(),
				                  dx, &fullform);
			}
		}
	}

	// average fine data onto coarser cells
	AverageDown(scal_force,comp,1);

	// fill ghost cells
	FillPatch(t_old, scal_force, scal_force, scal_force, comp, comp, 1, 0, bcs_f);

}

void
Maestro::MakeRhoHForce(Vector<MultiFab>& scal_force,
                       int is_prediction,
                       const Vector<MultiFab>& thermal,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       int add_thermal,
                       const int &which_step)

{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeRhoHForce()",MakeRhoHForce);

	// if we are doing the prediction, then it only makes sense to be in
	// this routine if the quantity we are predicting is rhoh', h, or rhoh
	if (is_prediction == 1 && !(enthalpy_pred_type == predict_rhohprime ||
	                            enthalpy_pred_type == predict_h ||
	                            enthalpy_pred_type == predict_rhoh) ) {
		Abort("ERROR: should only call mkrhohforce when predicting rhoh', h, or rhoh");
	}

	Vector<Real> rho0( (max_radial_level+1)*nr_fine );
	Vector<Real>   p0( (max_radial_level+1)*nr_fine );
	Vector<Real> grav( (max_radial_level+1)*nr_fine );
	Vector<Real> dpdt( (max_radial_level+1)*nr_fine );
	rho0.shrink_to_fit();
	p0.shrink_to_fit();
	grav.shrink_to_fit();
	dpdt.shrink_to_fit();

	if (which_step == 1) {
		rho0 = rho0_old;
		p0 =   p0_old;

		if (use_exact_base_state) {
			for (int i=0; i<p0_old.size(); ++i) {
				dpdt[i] = (p0_old[i] - p0_nm1[i])/dtold;
			}
		}
	}
	else {
		for(int i=0; i<rho0.size(); ++i) {
			rho0[i] = 0.5*(rho0_old[i]+rho0_new[i]);
			p0[i] = 0.5*(  p0_old[i]+  p0_new[i]);
		}

		if (use_exact_base_state) {
			for (int i=0; i<p0_old.size(); ++i) {
				dpdt[i] = (p0_new[i] - p0_old[i])/dt;
			}
		}
	}

	Vector<MultiFab> p0_cart(finest_level+1);
	Vector< std::array< MultiFab, AMREX_SPACEDIM > > p0mac(finest_level+1);
	if (spherical == 1) {
		for (int lev=0; lev<=finest_level; ++lev) {
			p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
			AMREX_D_TERM(p0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
			             p0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
			             p0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
		}

		Put1dArrayOnCart(p0, p0_cart, 0, 0, bcs_f, 0);
		MakeS0mac(p0, p0mac);
	}

	make_grav_cell(grav.dataPtr(),
	               rho0.dataPtr(),
	               r_cc_loc.dataPtr(),
	               r_edge_loc.dataPtr());


	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		MultiFab& scal_force_mf = scal_force[lev];
#if (AMREX_SPACEDIM >= 1)
		const MultiFab& umac_mf = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
		const MultiFab& vmac_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
		const MultiFab& wmac_mf = umac[lev][2];
		const MultiFab& p0cart_mf = p0_cart[lev];
		const MultiFab& p0macx_mf = p0mac[lev][0];
		const MultiFab& p0macy_mf = p0mac[lev][1];
		const MultiFab& p0macz_mf = p0mac[lev][2];
#endif
#endif
#endif
		const MultiFab& thermal_mf = thermal[lev];
		const MultiFab& cc_to_r = cell_cc_to_r[lev];

		// loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		for ( MFIter mfi(scal_force_mf); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& validBox = mfi.validbox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			if (spherical == 1) {
				const Real* dx = geom[lev].CellSize();
#if (AMREX_SPACEDIM == 3)
				if (use_exact_base_state) {
					// Use input parameter dpdt in place of psi
					mkrhohforce_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
					                 scal_force_mf[mfi].dataPtr(RhoH),
					                 ARLIM_3D(scal_force_mf[mfi].loVect()), ARLIM_3D(scal_force_mf[mfi].hiVect()),
					                 BL_TO_FORTRAN_3D(umac_mf[mfi]),
					                 BL_TO_FORTRAN_3D(vmac_mf[mfi]),
					                 BL_TO_FORTRAN_3D(wmac_mf[mfi]),
					                 BL_TO_FORTRAN_3D(thermal_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0cart_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0macx_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0macy_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0macz_mf[mfi]),
					                 dx, dpdt.dataPtr(),
					                 &is_prediction, &add_thermal,
					                 r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
					                 BL_TO_FORTRAN_3D(cc_to_r[mfi]));
				} else {
					mkrhohforce_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
					                 scal_force_mf[mfi].dataPtr(RhoH),
					                 ARLIM_3D(scal_force_mf[mfi].loVect()), ARLIM_3D(scal_force_mf[mfi].hiVect()),
					                 BL_TO_FORTRAN_3D(umac_mf[mfi]),
					                 BL_TO_FORTRAN_3D(vmac_mf[mfi]),
					                 BL_TO_FORTRAN_3D(wmac_mf[mfi]),
					                 BL_TO_FORTRAN_3D(thermal_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0cart_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0macx_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0macy_mf[mfi]),
					                 BL_TO_FORTRAN_3D(p0macz_mf[mfi]),
					                 dx, psi.dataPtr(),
					                 &is_prediction, &add_thermal,
					                 r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
					                 BL_TO_FORTRAN_3D(cc_to_r[mfi]));
				}
#else
				Abort("MakeRhoHForce: Spherical is not valid for DIM < 3");
#endif
			} else {
				mkrhohforce(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				            scal_force_mf[mfi].dataPtr(RhoH),
				            ARLIM_3D(scal_force_mf[mfi].loVect()), ARLIM_3D(scal_force_mf[mfi].hiVect()),
#if (AMREX_SPACEDIM == 1)
				            BL_TO_FORTRAN_3D(umac_mf[mfi]),
#elif (AMREX_SPACEDIM == 2)
				            BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#elif (AMREX_SPACEDIM == 3)
				            BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
				            BL_TO_FORTRAN_3D(thermal_mf[mfi]),
				            p0.dataPtr(), rho0.dataPtr(), grav.dataPtr(), psi.dataPtr(),
				            &is_prediction, &add_thermal);
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(scal_force,RhoH,1);
	FillPatch(t_old,scal_force,scal_force,scal_force,RhoH,RhoH,1,0,bcs_f);
}
