
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    Real cur_time = t_new;
    int last_plot_file_step = 0;
    
    // this will change when we have restart capability
    int start_step = 1;

    ///////////////
    // initial projection
    ///////////////
    if (do_initial_projection) {
        InitProj();
    }

    ///////////////
    // divu iters
    ///////////////
    for (int iter=0; iter<init_divu_iter; iter++) {

        // compute time step
        //
        //

        DivuIter();
    }

    ///////////////
    // init iters (for gradpi)
    ///////////////
    for (int iter=0; iter<init_iter; iter++) {

        // compute time step
        //
        //

        InitIter();
        t_new = 0.;
    }

    ///////////////
    // regular time steps
    ///////////////
    for (istep = start_step; istep <= max_step && t_new < stop_time; ++istep)
    {

        if (regrid_int > 0)  // We may need to regrid
        {
            if ( (istep-1) % regrid_int == 0)
            {

                // wallclock time
                const Real strt_total = ParallelDescriptor::second();

                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                regrid(0, t_new);

                // wallclock time
                Real end_total = ParallelDescriptor::second() - strt_total;
                ParallelDescriptor::ReduceRealMax(end_total,
                                                  ParallelDescriptor::IOProcessorNumber());
	
                // print wallclock time
                if (Verbose()) {
                    Print() << "Time to regrid: " << end_total << '\n';
                }
            }
        }
    
        // wallclock time
        const Real strt_total = ParallelDescriptor::second();

        // compute time step
        ComputeDt();

        Print() << "\nTimestep " << istep << " starts with TIME = " << t_new
                << " DT = " << dt << std::endl << std::endl;

        AdvanceTimeStep(cur_time, false);

        Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
                << " DT = " << dt << std::endl;

        // wallclock time
        Real end_total = ParallelDescriptor::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
        
        // print wallclock time
        if (Verbose()) {
            Print() << "Time to advance time step: " << end_total << '\n';
        }

        // checkpoint file
        //
        //
        //
        
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

void
Maestro::InitProj()
{

}

void
Maestro::DivuIter()
{

}

void
Maestro::InitIter()
{

}
