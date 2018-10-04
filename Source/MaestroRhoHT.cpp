
#include <Maestro.H>

using namespace amrex;

void
Maestro::TfromRhoH (Vector<MultiFab>& scal,
                    const RealVector& p0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoH()",TfromRhoH);

    // base state on cartesian grid
    Vector<MultiFab> p0_cart(finest_level+1);

    if (spherical) {
	for (int lev=0; lev<=finest_level; ++lev) {
	    p0_cart[lev].define(grids[lev], dmap[lev],1,0);
	}
	if (use_eos_e_instead_of_h) {
	    Put1dArrayOnCart(p0, p0_cart, 0, 0, bcs_f, 0);
	}
    }
    
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];
	const MultiFab& p0cart_mf = p0_cart[lev];

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

#pragma gpu
		makeTfromRhoH_sphr(AMREX_INT_ANYD(validBox.loVect()), AMREX_INT_ANYD(validBox.hiVect()),
				   BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf[mfi].nCompPtr(), 
				   BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]));
	    } else {

#pragma gpu
		makeTfromRhoH(AMREX_INT_ANYD(validBox.loVect()), AMREX_INT_ANYD(validBox.hiVect()), lev,
			      BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf[mfi].nCompPtr(),
			      p0.dataPtr());
	    }
        }

    }
    
#ifdef AMREX_USE_CUDA
    Device::synchronize();
#endif
    
    // average down and fill ghost cells
    AverageDown(scal,Temp,1);
    FillPatch(t_old,scal,scal,scal,Temp,Temp,1,Temp,bcs_s);
}

void
Maestro::TfromRhoP (Vector<MultiFab>& scal,
                    const RealVector& p0,
                    int updateRhoH)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TfromRhoP()",TfromRhoP);

    // base state on cartesian grid
    Vector<MultiFab> p0_cart(finest_level+1);

    if (spherical) {
	for (int lev=0; lev<=finest_level; ++lev) {
	    p0_cart[lev].define(grids[lev], dmap[lev],1,0);
	}

	Put1dArrayOnCart(p0, p0_cart, 0, 0, bcs_f, 0);
    }

    
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];
	const MultiFab& p0cart_mf = p0_cart[lev];

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

#pragma gpu
		makeTfromRhoP_sphr(AMREX_INT_ANYD(validBox.loVect()), AMREX_INT_ANYD(validBox.hiVect()),
				   BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf[mfi].nCompPtr(),
				   BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]), updateRhoH);
	    } else {

#pragma gpu
		makeTfromRhoP(AMREX_INT_ANYD(validBox.loVect()), AMREX_INT_ANYD(validBox.hiVect()), lev,
			      BL_TO_FORTRAN_ANYD(scal_mf[mfi]), scal_mf[mfi].nCompPtr(), 
			      p0.dataPtr(), updateRhoH);
	    }
        }

    }

#ifdef AMREX_USE_CUDA
    Device::synchronize();
#endif
    
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
#pragma gpu
	    makePfromRhoH(AMREX_INT_ANYD(validBox.loVect()), AMREX_INT_ANYD(validBox.hiVect()),
                          BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf[mfi].nCompPtr(), 
			  sold_mf[mfi].dataPtr(Temp), 
			  AMREX_INT_ANYD(sold_mf[mfi].loVect()), AMREX_INT_ANYD(sold_mf[mfi].hiVect()),
                          BL_TO_FORTRAN_ANYD(peos_mf[mfi]));
        }

    }

#ifdef AMREX_USE_CUDA
    Device::synchronize();
#endif
    
    // average down and fill ghost cells
    AverageDown(peos,0,1);
    FillPatch(t_old,peos,peos,peos,0,0,1,0,bcs_f);
}
