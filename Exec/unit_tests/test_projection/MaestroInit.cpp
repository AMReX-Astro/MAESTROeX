
#include <AMReX_VisMF.H>
#include <Maestro.H>
#include <Maestro_F.H>
#include <Problem_F.H>
using namespace amrex;

// initialize AMR data
// perform initial projection
// perform divu iters
// perform initial (pressure) iterations
void Maestro::Init() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Init()", Init);

    Print() << "Calling Init()" << std::endl;

    start_step = 1;

    // fill in multifab and base state data
    InitData();
}

// fill in multifab and base state data
void Maestro::InitData() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    Print() << "Calling InitData()" << std::endl;

    // calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
    // that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
    InitFromScratch(t_old);

    // reset tagging array to include buffer zones
    TagArray();

    // compute numdisjointchunks, r_start_coord, r_end_coord
    BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                               base_geom.nr_fine);
    base_geom.InitMultiLevel(finest_level, tag_array_b.array());
}

// During initialization of a simulation, Maestro::InitData() calls
// AmrCore::InitFromScratch(), which calls
// a MakeNewGrids() function that repeatedly calls this function to build
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch(int lev, Real time, const BoxArray& ba,
                                      const DistributionMapping& dm) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNewLevelFromScratch()",
                   MakeNewLevelFromScratch);

    sold[lev].define(ba, dm, Nscal, ng_s);
    snew[lev].define(ba, dm, Nscal, ng_s);
    uold[lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
    unew[lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
    gpi[lev].define(ba, dm, AMREX_SPACEDIM, 1);
    rhcc_for_nodalproj[lev].define(ba, dm, 1, 1);

    pi[lev].define(convert(ba, nodal_flag), dm, 1, 1);  // nodal

    sold[lev].setVal(0.);
    snew[lev].setVal(0.);
    uold[lev].setVal(0.);
    unew[lev].setVal(0.);
    gpi[lev].setVal(0.);
    rhcc_for_nodalproj[lev].setVal(0.);
    pi[lev].setVal(0.);
}
