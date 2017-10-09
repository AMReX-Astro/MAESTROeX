
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
	
        // initialize AMR data
        maestro.InitData();

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
