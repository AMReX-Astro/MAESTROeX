
#include <Maestro.H>

using namespace amrex;

// compute S at cell-centers
void
Maestro::Make_S_cc (Vector<MultiFab>& S_cc,
                    const Vector<MultiFab>& scal,
                    const Vector<MultiFab>& rho_omegadot,
                    const Vector<MultiFab>& rho_Hnuc,
                    const Vector<MultiFab>& rho_Hext,
                    const Vector<MultiFab>& thermal)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Make_S_cc()",Make_S_cc);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& S_cc_mf = S_cc[lev];
        const MultiFab& scal_mf = scal[lev];
        const MultiFab& rho_odot_mf = rho_omegadot[lev];
        const MultiFab& rho_Hnuc_mf = rho_Hnuc[lev];
        const MultiFab& rho_Hext_mf = rho_Hext[lev];
        const MultiFab& thermal_mf = thermal[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf,true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_S_cc(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                      BL_TO_FORTRAN_3D(S_cc_mf[mfi]),
                      BL_TO_FORTRAN_FAB(scal_mf[mfi]),
                      BL_TO_FORTRAN_FAB(rho_odot_mf[mfi]),
                      BL_TO_FORTRAN_3D(rho_Hnuc_mf[mfi]),
                      BL_TO_FORTRAN_3D(rho_Hext_mf[mfi]),
                      BL_TO_FORTRAN_3D(thermal_mf[mfi]));
        }

    }

    // average fine data onto coarser cells
    AverageDown(S_cc,0,1);

}

// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void
Maestro::MakeRHCCforNodalProj (Vector<MultiFab>& rhcc,
                               const Vector<MultiFab>& S_cc,
                               const Vector<Real>& Sbar,
                               const Vector<Real>& beta0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRHCCforNodalProj()",MakeRHCCforNodalProj);

    Vector<MultiFab> Sbar_cart(finest_level+1);
    Vector<MultiFab> beta0_cart(finest_level+1);

    if (spherical == 1) {
	for (int lev=0; lev<=finest_level; ++lev) {
	    Sbar_cart [lev].define(grids[lev], dmap[lev], 1, 0);
	    beta0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	    Sbar_cart [lev].setVal(0.);
	    beta0_cart[lev].setVal(0.);
	}

	Put1dArrayOnCart(Sbar,Sbar_cart,0,0,bcs_f,0);
	Put1dArrayOnCart(beta0,beta0_cart,0,0,bcs_f,0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // fill rhcc
        // get references to the MultiFabs at level lev
              MultiFab& rhcc_mf = rhcc[lev];
        const MultiFab& S_cc_mf = S_cc[lev];
	const MultiFab& Sbar_mf = Sbar_cart[lev];
	const MultiFab& beta0_mf = beta0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 1) {
		make_rhcc_for_nodalproj_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
					     BL_TO_FORTRAN_3D(rhcc_mf[mfi]),
					     BL_TO_FORTRAN_3D(S_cc_mf[mfi]),
					     BL_TO_FORTRAN_3D(Sbar_mf[mfi]),
					     BL_TO_FORTRAN_3D(beta0_mf[mfi]));
	    } else {
		make_rhcc_for_nodalproj(&lev, ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
					BL_TO_FORTRAN_3D(rhcc_mf[mfi]),
					BL_TO_FORTRAN_3D(S_cc_mf[mfi]),
					Sbar.dataPtr(), beta0.dataPtr());
	    }
        }
    }

    // averge down and fill ghost cells using first-order extrapolation
    AverageDown(rhcc,0,1);
    FillPatch(t_old, rhcc, rhcc, rhcc, 0, 0, 1, 0, bcs_f);
}

void
Maestro::CorrectRHCCforNodalProj(Vector<MultiFab>& rhcc,
				 const Vector<Real>& rho0,
				 const Vector<Real>& beta0,
				 const Vector<Real>& gamma1bar,
				 const Vector<Real>& p0,
				 const Vector<MultiFab>& delta_p_term)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CorrectRHCCforNodalProj()",CorrectRHCCforNodalProj);

    // Local variables
    Vector<MultiFab> correction_cc(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
	correction_cc[lev].define(grids[lev], dmap[lev], 1, 1);
	correction_cc[lev].setVal(0.);
    }

    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    Vector<MultiFab>        p0_cart(finest_level+1);
    Vector<MultiFab>     beta0_cart(finest_level+1);
    Vector<MultiFab>      rho0_cart(finest_level+1);

    if (spherical == 1) {
	for (int lev=0; lev<=finest_level; ++lev) {
	    gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	    p0_cart[lev]       .define(grids[lev], dmap[lev], 1, 0);
	    beta0_cart[lev]    .define(grids[lev], dmap[lev], 1, 0);
	    rho0_cart[lev]     .define(grids[lev], dmap[lev], 1, 0);
	    gamma1bar_cart[lev].setVal(0.);
	    p0_cart[lev]       .setVal(0.);
	    beta0_cart[lev]    .setVal(0.);
	    rho0_cart[lev]     .setVal(0.);
	}

	Put1dArrayOnCart(gamma1bar,gamma1bar_cart,0,0,bcs_f,0);
	Put1dArrayOnCart(p0,p0_cart,0,0,bcs_f,0);
	Put1dArrayOnCart(beta0,beta0_cart,0,0,bcs_f,0);
	Put1dArrayOnCart(rho0,rho0_cart,0,0,bcs_s,Rho);
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        const MultiFab& delta_p_mf = delta_p_term[lev];
	MultiFab& correction_cc_mf = correction_cc[lev];
	const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];
	const MultiFab& p0_mf = p0_cart[lev];
	const MultiFab& beta0_mf = beta0_cart[lev];
	const MultiFab& rho0_mf = rho0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(correction_cc_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 1) {
		create_correction_cc_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
					  BL_TO_FORTRAN_3D(correction_cc_mf[mfi]),
					  BL_TO_FORTRAN_3D(delta_p_mf[mfi]),
					  BL_TO_FORTRAN_3D(beta0_mf[mfi]),
					  BL_TO_FORTRAN_3D(gamma1bar_mf[mfi]),
					  BL_TO_FORTRAN_3D(p0_mf[mfi]),
					  BL_TO_FORTRAN_3D(rho0_mf[mfi]), &dt);
	    } else {
		create_correction_cc(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				     BL_TO_FORTRAN_3D(correction_cc_mf[mfi]),
				     BL_TO_FORTRAN_3D(delta_p_mf[mfi]),
				     beta0.dataPtr(), gamma1bar.dataPtr(),
				     p0.dataPtr(), &dt);
	    }
        }
    }

    // average down and fill ghost cells using first-order extrapolation
    AverageDown(correction_cc,0,1);
    FillPatch(t_old, correction_cc, correction_cc, correction_cc, 0, 0, 1, 0, bcs_f);

    // add correction term
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::Add(rhcc[lev],correction_cc[lev],0,0,1,1);
    }
}

// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void
Maestro::MakeRHCCforMacProj (Vector<MultiFab>& rhcc,
			     const Vector<Real>& rho0,
                             const Vector<MultiFab>& S_cc,
                             const Vector<Real>& Sbar,
                             const Vector<Real>& beta0,
                             const Vector<Real>& gamma1bar,
			     const Vector<Real>& p0,
                             const Vector<MultiFab>& delta_p_term,
			     Vector<MultiFab>& delta_chi,
			     int is_predictor)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRHCCforMacProj()",MakeRHCCforMacProj);

    for (int lev=0; lev<=finest_level; ++lev) {

        // fill rhcc
        // get references to the MultiFabs at level lev
              MultiFab& rhcc_mf = rhcc[lev];
        const MultiFab& S_cc_mf = S_cc[lev];
	const MultiFab& delta_p_mf = delta_p_term[lev];
	      MultiFab& delta_chi_mf = delta_chi[lev];
	const MultiFab& cc_to_r = cell_cc_to_r[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 1) {
		const Real* dx = geom[lev].CellSize();

		make_rhcc_for_macproj_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
					   BL_TO_FORTRAN_3D(rhcc_mf[mfi]),
					   BL_TO_FORTRAN_3D(S_cc_mf[mfi]),
					   Sbar.dataPtr(), beta0.dataPtr(),
					   rho0.dataPtr(), dx,
					   gamma1bar.dataPtr(), p0.dataPtr(),
					   BL_TO_FORTRAN_3D(delta_p_mf[mfi]),
					   BL_TO_FORTRAN_3D(delta_chi_mf[mfi]),
					   &dt, &is_predictor,
					   r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
					   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
	    } else {
		make_rhcc_for_macproj(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				      BL_TO_FORTRAN_3D(rhcc_mf[mfi]),
				      BL_TO_FORTRAN_3D(S_cc_mf[mfi]),
				      Sbar.dataPtr(), beta0.dataPtr(),
				      gamma1bar.dataPtr(), p0.dataPtr(),
				      BL_TO_FORTRAN_3D(delta_p_mf[mfi]),
				      BL_TO_FORTRAN_3D(delta_chi_mf[mfi]),
				      &dt, &is_predictor);
	    }
        }
    }
}
