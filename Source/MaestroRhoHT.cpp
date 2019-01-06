
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
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			if (spherical == 1) {
				makeTfromRhoH_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				                   BL_TO_FORTRAN_FAB(scal_mf[mfi]),
				                   p0.dataPtr(), dx,
				                   r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
				                   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
			} else {
				makeTfromRhoH(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
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
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			if (spherical == 1) {
				makeTfromRhoP_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				                   BL_TO_FORTRAN_FAB(scal_mf[mfi]),
				                   p0.dataPtr(), dx, &updateRhoH,
				                   r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
				                   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
			} else {
				makeTfromRhoP(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
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
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			makePfromRhoH(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
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

void
Maestro::MachfromRhoH (const Vector<MultiFab>& scal,
                       const Vector<MultiFab>& vel,
                       const Vector<Real>& p0,
                       Vector<MultiFab>& mach)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MachfromRhoH()",MachfromRhoH);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& scal_mf = scal[lev];
		const MultiFab& vel_mf = vel[lev];
		MultiFab& mach_mf = mach[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			makeMachfromRhoH(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
			                 BL_TO_FORTRAN_FAB(scal_mf[mfi]),
			                 BL_TO_FORTRAN_3D(vel_mf[mfi]),
			                 p0.dataPtr(),w0.dataPtr(),
			                 BL_TO_FORTRAN_3D(mach_mf[mfi]));
		}

	}

	// average down and fill ghost cells
	AverageDown(mach,0,1);
	FillPatch(t_old,mach,mach,mach,0,0,1,0,bcs_f);
}

void
Maestro::MachfromRhoHSphr (const Vector<MultiFab>& scal,
                           const Vector<MultiFab>& vel,
                           const Vector<Real>& p0,
                           const Vector<MultiFab>& w0cart,
                           Vector<MultiFab>& mach)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MachfromRhoHSphr()",MachfromRhoHSphr);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& scal_mf = scal[lev];
		const MultiFab& vel_mf = vel[lev];
		const MultiFab& w0cart_mf = w0cart[lev];
		MultiFab& mach_mf = mach[lev];
		const MultiFab& cc_to_r = cell_cc_to_r[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
			makeMachfromRhoH_sphr(&lev,ARLIM_3D(tileBox.loVect()),
			                      ARLIM_3D(tileBox.hiVect()),
			                      BL_TO_FORTRAN_FAB(scal_mf[mfi]),
			                      BL_TO_FORTRAN_3D(vel_mf[mfi]),
			                      p0.dataPtr(),BL_TO_FORTRAN_3D(w0cart_mf[mfi]),
			                      dx,
			                      BL_TO_FORTRAN_3D(mach_mf[mfi]),
			                      r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
			                      BL_TO_FORTRAN_3D(cc_to_r[mfi]));
		}

	}

	// average down and fill ghost cells
	AverageDown(mach,0,1);
	FillPatch(t_old,mach,mach,mach,0,0,1,0,bcs_f);
}

void
Maestro::CsfromRhoH (const Vector<MultiFab>& scal,
                     const Vector<Real>& p0,
                     const Vector<MultiFab>& p0cart,
                     Vector<MultiFab>& cs)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::CsfromRhoH()",CsfromRhoH);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& scal_mf = scal[lev];
		MultiFab& cs_mf = cs[lev];

		if (spherical == 0) {

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				makeCsfromRhoH(&lev,ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				               BL_TO_FORTRAN_FAB(scal_mf[mfi]),
				               p0.dataPtr(),
				               BL_TO_FORTRAN_3D(cs_mf[mfi]));
			}
		} else {

			const MultiFab& p0cart_mf = p0cart[lev];

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(scal_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
				makeCsfromRhoH_sphr(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
				                    BL_TO_FORTRAN_FAB(scal_mf[mfi]),
				                    BL_TO_FORTRAN_3D(p0cart_mf[mfi]),
				                    BL_TO_FORTRAN_3D(cs_mf[mfi]));


			}
		}

	}

	// average down and fill ghost cells
	AverageDown(cs,0,1);
	FillPatch(t_old,cs,cs,cs,0,0,1,0,bcs_f);
}
