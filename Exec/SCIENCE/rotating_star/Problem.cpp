#include "Maestro.H"

#ifdef DO_PROBLEM_POST_TIMESTEP

using namespace amrex;

void
Maestro::problem_post_timestep() {
    Print() << "hiiiiiii" << std::endl;
}

#endif
