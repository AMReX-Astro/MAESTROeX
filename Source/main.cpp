
#include <Maestro.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    // in AMReX.cpp
    Initialize(argc,argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {
        // declare an Maestro object to manage multilevel data
        Maestro maestro;
	
        // read in C++/F90 parameters
        // define global C++/F90 variables and initialize network
        // set up boundary conditions
        // initialize base state geometry parameters
        // set istep, t_new, t_old
        // allocate MultiFabs and base state arrays
        maestro.Setup();

        // initialize multifab and base state data
        // perform initial projection
        // perform divu iters
        // perform initial (pressure) iterations
        maestro.Init();

        // advance solution to final time
        maestro.Evolve();
	
        // wallclock time
        Real end_total = ParallelDescriptor::second() - strt_total;
	
        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        if (maestro.Verbose()) {
            Print() << "\nTotal Time: " << end_total << '\n';
        }
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    // in AMReX.cpp
    Finalize();
}
