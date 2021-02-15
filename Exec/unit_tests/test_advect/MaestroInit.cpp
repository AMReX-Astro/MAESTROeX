
#include <AMReX_VisMF.H>
#include <Maestro.H>
#include <Maestro_F.H>
using namespace amrex;

// initialize AMR data
void Maestro::Init() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Init()", Init);

    Print() << "Calling Init()" << std::endl;

    start_step = 1;

    // fill in multifab and base state data
    InitData();

    // compute numdisjointchunks, r_start_coord, r_end_coord
    BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                               base_geom.nr_fine);
    base_geom.InitMultiLevel(finest_level, tag_array_b.array());

    // compute initial time step
    FirstDt();

    if (stop_time >= 0. && t_old + dt > stop_time) {
        dt = amrex::min(dt, stop_time - t_old);
        Print() << "Stop time limits dt = " << dt << std::endl;
    }

    dtold = dt;
    t_new = t_old + dt;
}

// fill in multifab and base state data
void Maestro::InitData() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    Print() << "Calling InitData()" << std::endl;

    // read in model file and fill in s0_init and p0_init for all levels

    for (auto lev = 0; lev <= base_geom.max_radial_level; ++lev) {
        InitBaseState(rho0_old, rhoh0_old, p0_old, lev);
    }

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

    // first compute cutoff coordinates using initial density profile
    ComputeCutoffCoords(rho0_old);
    base_geom.ComputeCutoffCoords(rho0_old.array());

    // set rho0 to be the average
    Average(sold, rho0_old, Rho);
    ComputeCutoffCoords(rho0_old);
    base_geom.ComputeCutoffCoords(rho0_old.array());

    // call eos with r,p as input to recompute T,h
    TfromRhoP(sold, p0_old, 1);

    // set rhoh0 to be the average
    Average(sold, rhoh0_old, RhoH);
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
    S_cc_old[lev].define(ba, dm, 1, 0);
    S_cc_new[lev].define(ba, dm, 1, 0);
    gpi[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    w0_cart[lev].define(ba, dm, AMREX_SPACEDIM, 2);
    rhcc_for_nodalproj[lev].define(ba, dm, 1, 1);
    pi[lev].define(convert(ba, nodal_flag), dm, 1, 0);  // nodal

    sold[lev].setVal(0.);
    snew[lev].setVal(0.);
    uold[lev].setVal(0.);
    S_cc_old[lev].setVal(0.);
    S_cc_new[lev].setVal(0.);
    gpi[lev].setVal(0.);
    w0_cart[lev].setVal(0.);
    rhcc_for_nodalproj[lev].setVal(0.);
    pi[lev].setVal(0.);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        InitLevelData(lev, t_old, mfi, sold[lev].array(mfi),
                      uold[lev].array(mfi));
    }
}
