
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    Print() << "Calling Evolve()" << endl;

    for (istep = 1; istep <= max_step && t_new < stop_time; ++istep)
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
        }

        // reset t_old and t_new
        t_old = t_new;
        t_new += dt;

        // advance the solution by dt
        AdvanceTimeStep(false);

        // write a plotfile
        if (plot_int > 0 && (istep % plot_int == 0 || istep == max_step) )
        {
            Print() << "\nWriting plotfile " << istep << endl;
            WritePlotFile(istep,t_new,rho0_new,p0_new,unew,snew);
        }

        // move new state into old state by swapping pointers
        for (int lev=0; lev<=finest_level; ++lev) {
            std::swap(    sold[lev],     snew[lev]);
            std::swap(    uold[lev],     unew[lev]);
            std::swap(S_cc_old[lev], S_cc_new[lev]);

            std::swap( rho0_old, rho0_new);
            std::swap(rhoh0_old,rhoh0_new);
            std::swap(   p0_old,   p0_new);

            std::swap(beta0_old,beta0_new);
            std::swap(grav_cell_old,grav_cell_new);
        }

    }
}
