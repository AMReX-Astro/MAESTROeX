
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeAngularMomentum (const Vector<MultiFab>& state,
                              const Vector<MultiFab>& vel,
                              Vector<MultiFab>& angular_momentum)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeAngularMomentum()",MakeAngularMomentum);

	for (int lev=0; lev<=finest_level; ++lev) {

		// Get the index space and grid spacing of the domain
		const Box& domainBox = geom[lev].Domain();
		const Real* dx = geom[lev].CellSize();

		// get references to the MultiFabs at level lev
		const MultiFab& state_mf = state[lev];
		const MultiFab& vel_mf = vel[lev];
		MultiFab& angular_momentum_mf = angular_momentum[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(angular_momentum_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const RealBox temp    (tileBox,geom[lev].CellSize(),geom[lev].ProbLo());
			const Real* xlo     = temp.lo();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "tileBox", which specifies the "valid" region.
			derangmom(BL_TO_FORTRAN_3D(angular_momentum_mf[mfi]),
			           BL_TO_FORTRAN_FAB(state_mf[mfi]),
			           BL_TO_FORTRAN_3D(vel_mf[mfi]),
			           ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
			           domainBox.loVect(), domainBox.hiVect(),
			           dx, xlo);
		}

	}

	// average fine data onto coarser cells
	AverageDown(angular_momentum,0,AMREX_SPACEDIM);
	FillPatch(t_old,angular_momentum,angular_momentum,angular_momentum,0,0,AMREX_SPACEDIM,0,bcs_u);

}
