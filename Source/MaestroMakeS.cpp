
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
        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_S_cc(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
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

    for (int lev=0; lev<=finest_level; ++lev) {

        // fill rhcc
        // get references to the MultiFabs at level lev
              MultiFab& rhcc_mf = rhcc[lev];
        const MultiFab& S_cc_mf = S_cc[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_rhcc_for_nodalproj(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                                    BL_TO_FORTRAN_3D(rhcc_mf[mfi]),
                                    BL_TO_FORTRAN_3D(S_cc_mf[mfi]),
                                    Sbar.dataPtr(), beta0.dataPtr());
        }
    }

    // averge down and fill ghost cells using first-order extrapolation
    AverageDown(rhcc,0,1);
    FillPatch(t_old, rhcc, rhcc, rhcc, 0, 0, 1, 0, bcs_f);
}


// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
void
Maestro::MakeRHCCforMacProj (Vector<MultiFab>& rhcc,
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

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
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
    Vector<MultiFab>    correction_cc(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
	correction_cc[lev]   .define(grids[lev], dmap[lev], 1, 1);
	correction_cc[lev]   .setVal(0.);
    }

    if (spherical == 1) {
	// build multifabs of 1d arrays
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        const MultiFab& delta_p_mf = delta_p_term[lev];
	MultiFab& correction_cc_mf = correction_cc[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(correction_cc_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            create_correction_cc(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				 BL_TO_FORTRAN_3D(correction_cc_mf[mfi]),
				 BL_TO_FORTRAN_3D(delta_p_mf[mfi]),
				 beta0.dataPtr(), gamma1bar.dataPtr(), 
				 p0.dataPtr(), &dt);
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
