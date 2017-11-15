
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    Print() << "Calling Evolve()" << endl;

    int last_plot_file_step = 0;

    for (istep = 1; istep <= max_step && t_new < stop_time; ++istep)
    {

        // check to see if we need to regrid, then regrid
        if (max_level > 0 && regrid_int > 0 && (istep-1) % regrid_int == 0) {
            Regrid();
        }

        // move new state into old state by swapping pointers
        for (int lev=0; lev<=finest_level; ++lev) {
            std::swap(    sold[lev],     snew[lev]);
            std::swap(    uold[lev],     unew[lev]);
            std::swap(S_cc_old[lev], S_cc_new[lev]);
        }

        // compute time step
        EstDt();

        // reset t_old and t_new
        t_old = t_new;
        t_new += dt;

        // advance the solution by dt
        AdvanceTimeStep(false);

        // write a plotfile
        if (plot_int > 0 && istep % plot_int == 0)
        {
            Print() << "\nWriting plotfile " << istep << std::endl;
            last_plot_file_step = istep;
            WritePlotFile(istep);
        }

    }

    // write a final plotfile if we haven't already
    if (plot_int > 0 && istep > last_plot_file_step)
    {
        Print() << "\nWriting plotfile " << istep-1 << std::endl;
        WritePlotFile(istep-1);
    }
}
