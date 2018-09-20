
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Evolve()",Evolve);

    Print() << "Calling Evolve()" << std::endl;

    // check to make sure spherical is only used for 3d
    if (spherical == 1 && AMREX_SPACEDIM != 3) {
	Abort ("spherical = 1 and dm != 3");
    }

    // index for diag array buffer
    int diag_index=0;
    
    for (istep = start_step; istep <= max_step && t_old < stop_time; ++istep)
    {

        // check to see if we need to regrid, then regrid
        if (max_level > 0 && regrid_int > 0 && (istep-1) % regrid_int == 0) {
            Regrid();
        }

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

        // advance the solution by dt
	if (use_exact_base_state) {
	    AdvanceTimeStepIrreg(false);
	} else {
	    AdvanceTimeStep(false);
	}
	    
        t_old = t_new;

	// save diag output into buffer
	DiagFile(istep,t_new,rho0_new,p0_new,unew,snew,diag_index);

        // write a plotfile
        if (plot_int > 0 && ( (istep % plot_int == 0) || 
                              (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
                              (istep == max_step) ) )
        {
            Print() << "\nWriting plotfile " << istep << std::endl;
            WritePlotFile(istep,t_new,rho0_new,p0_new,unew,snew);
        }

        if (chk_int > 0 && (istep % chk_int == 0 || t_new >= stop_time || istep == max_step) )
        {
	    // write a checkpoint file
            Print() << "\nWriting checkpoint" << istep << std::endl;
            WriteCheckPoint(istep);
        }

	if (diag_index == diag_buf_size || istep == max_step) {
	    // write out any buffered diagnostic information
            WriteDiagFile(diag_index);
	}

        // move new state into old state by swapping pointers
        for (int lev=0; lev<=finest_level; ++lev) {
            std::swap(    sold[lev],     snew[lev]);
            std::swap(    uold[lev],     unew[lev]);
            std::swap(S_cc_old[lev], S_cc_new[lev]);
	}
	
	std::swap( rho0_old, rho0_new);
	std::swap(rhoh0_old,rhoh0_new);
	
	// Note std::swap() is incompatible with CudaManagedAllocator
	p0_nm1.swap(p0_old);
        p0_old.swap(p0_new);
        
	std::swap(beta0_old,beta0_new);
	std::swap(grav_cell_old,grav_cell_new);

    }
}
