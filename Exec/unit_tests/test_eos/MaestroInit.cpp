
#include <AMReX_VisMF.H>
#include <Maestro.H>
#include <Maestro_F.H>
using namespace amrex;

void Maestro::Init() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Init()", Init);

    Print() << "Calling Init()" << std::endl;

    // fill in multifab and base state data
    InitData();
}

// fill in multifab and base state data
void Maestro::InitData() {
    // calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
    // that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
    InitFromScratch(t_old);

    // reset tagging array to include buffer zones
    TagArray();

    // compute numdisjointchunks, r_start_coord, r_end_coord
    BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                               base_geom.nr_fine);
    base_geom.InitMultiLevel(finest_level, tag_array_b.array());

    // average down data and fill ghost cells
    AverageDown(sold, 0, Nscal);
    FillPatch(t_old, sold, sold, sold, 0, 0, Nscal, 0, bcs_s);
}

// During initialization of a simulation, Maestro::InitData() calls
// AmrCore::InitFromScratch(), which calls
// a MakeNewGrids() function that repeatedly calls this function to build
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch(int lev, Real time, const BoxArray& ba,
                                      const DistributionMapping& dm) {
    sold[lev].define(ba, dm, Nscal, ng_s);
    sold[lev].setVal(0.);
}
