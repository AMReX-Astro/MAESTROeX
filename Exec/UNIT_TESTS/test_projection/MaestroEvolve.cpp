
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Evolve()",Evolve);

    Print() << "Calling Evolve()" << std::endl;

}
