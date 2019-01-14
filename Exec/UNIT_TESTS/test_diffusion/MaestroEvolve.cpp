
#include <Maestro.H>

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

			if (fixed_dt != -1.) {
				dt = fixed_dt;
				if (maestro_verbose > 0) {
					Print() << "Setting fixed dt = " << dt;
				}
			}

			if (stop_time >= 0. && t_old+dt > stop_time) {
				dt = std::min(dt,stop_time-t_old);
				Print() << "Stop time limits dt = " << dt << std::endl;
			}

			t_new = t_old + dt;
		}

		AdvanceTimeStep(false);

		t_old = t_new;

		// move new state into old state by swapping pointers
		for (int lev=0; lev<=finest_level; ++lev)
			std::swap(    sold[lev],     snew[lev]);

	}
}
