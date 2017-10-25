
#include <Maestro.H>
#include <maestro_defaults.H>

using namespace amrex;

int Maestro::Rho       = -1;
int Maestro::RhoH      = -1;
int Maestro::FirstSpec = -1;
int Maestro::NumSpec   = -1;
int Maestro::Temp      = -1;
int Maestro::Pi        = -1;
int Maestro::NSCAL     = -1;

// constructor - reads in parameters from inputs file
//             - sizes multilevel vectors and data structures
Maestro::Maestro ()
{

    ///////
    // Geometry on all levels has been defined already.
    ///////

    // read in C++ parameters in maestro_queries.H using ParmParse pp("maestro");
    ReadParameters();

    // read in F90 parameters in meth_params.F90
    ca_set_maestro_method_params();

    // define (Rho, RhoH, etc.)
    // calls ca_network_init
    VariableSetup();

    // define additional module variables in meth_params.F90
    ca_set_method_params(Rho,RhoH,FirstSpec,Temp,Pi);

    // set up BCRec definitions for BC types
    BCSetup();

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep = 0;

    t_new = 0.0;
    t_old = -1.e100;

    // set this to a large number so change_max doesn't affect the first time step
    dt = 1.e100;

    snew.resize(nlevs_max);
    sold.resize(nlevs_max);
    unew.resize(nlevs_max);
    uold.resize(nlevs_max);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg_s.resize(nlevs_max);
    flux_reg_u.resize(nlevs_max);
}

Maestro::~Maestro ()
{

}

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

        AdvanceTimeStep(cur_time);

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
