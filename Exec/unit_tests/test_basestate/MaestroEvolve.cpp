
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::Evolve()",Evolve);

	Print() << "Calling Evolve()" << std::endl;

	for (istep = start_step; istep <= max_step && t_old < stop_time; ++istep)
	{

		dtold = dt;

		// compute time step
		// if this is the first time step we already have a dt from either FirstDt()
		// or EstDt called during the divu_iters
		if (istep > 1) {

			EstDt();

			if (verbose > 0) {
				Print() << "Call to estdt at beginning of step " << istep
				        << " gives dt =" << dt << std::endl;
			}

			// fixme - add nuclear_dt_scalefac timestep limiter

			if (dt > max_dt_growth*dtold) {
				dt = max_dt_growth*dtold;
				if (verbose > 0) {
					Print() << "dt_growth factor limits the new dt = " << dt << std::endl;
				}
			}

			if (dt > max_dt) {
				if (verbose > 0) {
					Print() << "max_dt limits the new dt = " << max_dt << std::endl;
				}
				dt = max_dt;
			}

			if (stop_time >= 0. && t_old+dt > stop_time) {
				dt = std::min(dt,stop_time-t_old);
				Print() << "Stop time limits dt = " << dt << std::endl;
			}

			t_new = t_old + dt;
		}

		// advance the solution by dt

		AdvanceTimeStep(false);

		t_old = t_new;

		// move new state into old state by swapping pointers
		for (int lev=0; lev<=finest_level; ++lev) {
			std::swap(S_cc_old[lev], S_cc_new[lev]);
		}

		std::swap( rho0_old, rho0_new);
		rhoh0_old.swap(rhoh0_new);
		p0_nm1.swap(p0_old);
		p0_old.swap(p0_new);

		gamma1bar_old.swap(gamma1bar_new);
		grav_cell_old.swap(grav_cell_new);
	}

	// Now need to check the HSE-ness
	Real max_hse_error;
	check_hseness(rho0_old.dataPtr(), p0_old.dataPtr(), &max_hse_error);
	Print() << "Maximum HSE error = " << max_hse_error << std::endl;
}
