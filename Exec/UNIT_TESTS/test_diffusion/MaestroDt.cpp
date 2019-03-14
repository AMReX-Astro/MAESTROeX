
#include <Maestro.H>

using namespace amrex;

void
Maestro::FirstDt ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::FirstDt()",FirstDt);

	dt = 1.e20;

    if (fixed_dt != -1.0) {
		// fixed dt
		dt = fixed_dt;
		if (maestro_verbose > 0) {
			Print() << "Setting fixed dt = " << dt << std::endl;
		}
        return;
	}

	// build and compute vel_force
	Vector<MultiFab> vel_force(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev)
		vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

	for (int lev = 0; lev <= finest_level; ++lev) {
		Real dt_lev = 1.e99;
		Real umax_lev = 0.;

		// get references to the MultiFabs at level lev
		MultiFab& uold_mf = uold[lev];
		MultiFab& sold_mf = sold[lev];
		MultiFab& vel_force_mf = vel_force[lev];
		MultiFab& S_cc_old_mf = S_cc_old[lev];
		const MultiFab& cc_to_r = cell_cc_to_r[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
		for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

			Real dt_grid = 1.e99;
			Real umax_grid = 0.;

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.

			firstdt(&lev,&dt_grid,&umax_grid,
			        ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
			        ZFILL(dx),
			        BL_TO_FORTRAN_FAB(sold_mf[mfi]),
			        BL_TO_FORTRAN_FAB(uold_mf[mfi]),
			        BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
			        BL_TO_FORTRAN_3D(S_cc_old_mf[mfi]),
			        p0_old.dataPtr(),
			        gamma1bar_old.dataPtr());

			dt_lev = std::min(dt_lev,dt_grid);
		}

		// find the smallest dt over all processors
		ParallelDescriptor::ReduceRealMin(dt_lev);

		if (maestro_verbose > 0) {
			Print() << "Call to firstdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
		}

		// multiply by init_shrink
		dt_lev *= init_shrink;

		if (maestro_verbose > 0) {
			Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev << std::endl;
		}

		// update dt over all levels
		dt = std::min(dt,dt_lev);

	} // end loop over levels

	if (maestro_verbose > 0) {
		Print() << "Minimum firstdt over all levels = " << dt << std::endl;
	}

	if (dt < small_dt) {
		Abort("FirstDt: dt < small_dt");
	}

	if (dt > max_dt) {
		if (maestro_verbose > 0) {
			Print() << "max_dt limits the new dt = " << max_dt << std::endl;
		}
		dt = max_dt;
	}
}
