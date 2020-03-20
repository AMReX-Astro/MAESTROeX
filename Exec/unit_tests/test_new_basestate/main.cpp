
#include <BaseState.H>
using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[])
{

    // in AMReX.cpp
    Initialize(argc, argv);

    // // Refuse to continue if we did not provide an inputs file.

    // if (argc <= 1) {
    //     Abort("Error: no inputs file provided on command line.");
    // }

    // // Save the inputs file name for later.

    // if (!strchr(argv[1], '=')) {
    //     inputs_name = argv[1];
    // }

    // timer for profiling
    BL_PROFILE_VAR("main()", main);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {
        BaseState<Real> base_state;

        const int nlevs = 3;
        const int len = 20;
        const int ncomp = 5;

        Print() << "Defining base state with " << nlevs << " levels, length " << len << " and " << ncomp << " components." << std::endl;

        base_state.define(nlevs, len, ncomp);

        Print() << "setting value to 1.0" << std:: endl;
        base_state.setVal(1.0);

        Print() << "print an element: " << base_state(0,0) << std::endl;

        Print() << "defining another base state" << std::endl;
        BaseState<Real> other_base_state(nlevs, len, ncomp);

        for (auto l = 0; l < nlevs; ++l) {
            std::fill(other_base_state[l].begin(), other_base_state[l].end(), 1.0);
        }

        Print() << "are the two base states equal? " << (base_state == other_base_state) << std::endl;
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(main);

    // in AMReX.cpp
    Finalize();
}
