
#include <Maestro.H>
#include <MaestroPlot.H>

using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[]) {
    // in AMReX.cpp
    Initialize(argc, argv);

    // Refuse to continue if we did not provide an inputs file.

    if (argc <= 1) {
        Abort("Error: no inputs file provided on command line.");
    }

    // Save the inputs file name for later.

    if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
    }

    // timer for profiling
    BL_PROFILE_VAR("main()", main);

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

        // advance solution to final time
        maestro.Evolve();

        // wallclock time
        Real end_total = ParallelDescriptor::second() - strt_total;

        // print wallclock time
        ParallelDescriptor::ReduceRealMax(
            end_total, ParallelDescriptor::IOProcessorNumber());
        Print() << "\nTotal Time: " << end_total << '\n';
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(main);

    // in AMReX.cpp
    Finalize();
}
