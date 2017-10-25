
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    Real cur_time = t_new;
    int last_plot_file_step = 0;

    for (istep = 1; istep <= max_step && cur_time < stop_time; ++istep)
    {

        if (regrid_int > 0)  // We may need to regrid
        {
            if ( (istep-1) % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                regrid(0, cur_time);
            }
        }
    
        // wallclock time
        const Real strt_total = ParallelDescriptor::second();

        // compute time step
        ComputeDt();

        Print() << "\nTimestep " << istep << " starts with TIME = " << cur_time 
                       << " DT = " << dt << std::endl << std::endl;

        AdvanceTimeStep(cur_time,false);

        cur_time += dt;

        Print() << "\nTimestep " << istep << " ends with TIME = " << cur_time 
                       << " DT = " << dt << std::endl;

        // wallclock time
        Real end_total = ParallelDescriptor::second() - strt_total;
	
        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        if (Verbose()) {
            Print() << "Time to advance time step: " << end_total << '\n';
        }

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
