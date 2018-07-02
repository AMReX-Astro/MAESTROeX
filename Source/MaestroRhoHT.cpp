
#include <Maestro.H>

using namespace amrex;

void
Maestro::TfromRhoH (Vector<MultiFab>& scal,
                    const Vector<Real>& p0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoH()",TfromRhoH);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];
	const MultiFab& cc_to_r = cell_cc_to_r[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
	    const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 1) {
		makeTfromRhoH_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				   BL_TO_FORTRAN_FAB(scal_mf[mfi]),
				   p0.dataPtr(), dx, 
				   r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
				   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
	    } else {
		makeTfromRhoH(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			      BL_TO_FORTRAN_FAB(scal_mf[mfi]),
			      p0.dataPtr());
	    }
        }

    }

    // average down and fill ghost cells
    AverageDown(scal,Temp,1);
    FillPatch(t_old,scal,scal,scal,Temp,Temp,1,Temp,bcs_s);
}

void
Maestro::TfromRhoP (Vector<MultiFab>& scal,
                    const Vector<Real>& p0,
                    int updateRhoH)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoP()",TfromRhoP);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];
	const MultiFab& cc_to_r = cell_cc_to_r[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
	    const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
	    if (spherical == 1) {
		makeTfromRhoP_sphr(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				   BL_TO_FORTRAN_FAB(scal_mf[mfi]),
				   p0.dataPtr(), dx, &updateRhoH, 
				   r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
				   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
	    } else {
		makeTfromRhoP(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			      BL_TO_FORTRAN_FAB(scal_mf[mfi]),
			      p0.dataPtr(),&updateRhoH);
	    }
        }

    }


    // average down and fill ghost cells (Temperature)
    AverageDown(scal,Temp,1);
    FillPatch(t_old,scal,scal,scal,Temp,Temp,1,Temp,bcs_s);

    // average down and fill ghost cells (Enthalpy)
    if (updateRhoH == 1) {
        AverageDown(scal,RhoH,1);
        FillPatch(t_old,scal,scal,scal,RhoH,RhoH,1,RhoH,bcs_s);
    }

}

void
Maestro::PfromRhoH (const Vector<MultiFab>& state,
                    const Vector<MultiFab>& s_old, 
		    Vector<MultiFab>& peos)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PfromRhoH()",PfromRhoH);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& state_mf = state[lev];
	const MultiFab& sold_mf = s_old[lev];
	MultiFab& peos_mf = peos[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(state_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            makePfromRhoH(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                          BL_TO_FORTRAN_FAB(state_mf[mfi]),
			  sold_mf[mfi].dataPtr(Temp), 
			  ARLIM_3D(sold_mf[mfi].loVect()), ARLIM_3D(sold_mf[mfi].hiVect()),
                          BL_TO_FORTRAN_3D(peos_mf[mfi]));
        }

    }

    // average down and fill ghost cells
    AverageDown(peos,0,1);
    FillPatch(t_old,peos,peos,peos,0,0,1,0,bcs_f);
}
