
#include <Maestro.H>
#include <AMReX_VisMF.H>

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

// compute S at cell-centers
void
Maestro::WriteCheckPoint () {

    std::string checkpointname = "test";
    bool call_barrier = true;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level+1, call_barrier);


    WriteCheckPointHeader(checkpointname);

}

void
Maestro::WriteCheckPointHeader(const std::string& name) {

}
