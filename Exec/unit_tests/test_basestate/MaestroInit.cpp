
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

    if (plot_int > 0) {
        // Need to fill normal vector to compute velrc in plotfile
        if (spherical) {
            MakeNormal();
        }

        Print() << "\nWriting plotfile plt_InitData after InitData"
                << std::endl;
        WritePlotFile(9999999, t_old, 0, rho0_old, rhoh0_old, p0_old,
                      gamma1bar_old, uold, sold, S_cc_old);
    }

    // compute numdisjointchunks, r_start_coord, r_end_coord
    BaseState<int> tag_array_b(tag_array, base_geom.max_radial_level + 1,
                               base_geom.nr_fine);
    base_geom.InitMultiLevel(finest_level, tag_array_b.array());

#if (AMREX_SPACEDIM == 3)
    if (spherical == 1) {
        MakeNormal();
        MakeCCtoRadii();
    }
#endif

    // make gravity
    MakeGravCell(grav_cell_old, rho0_old);

    // compute initial time step
    // FirstDt();
    dt = initial_dt;

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

    // first compute cutoff coordinates using initial density profile
    ComputeCutoffCoords(rho0_old);
    base_geom.ComputeCutoffCoords(rho0_old.array());

    // compute gravity
    MakeGravCell(grav_cell_old, rho0_old);

    // compute p0 with HSE
    EnforceHSE(rho0_old, p0_old, grav_cell_old);

    // set p0^{-1} = p0_old
    p0_nm1.copy(p0_old);

    // set some stuff to zero
    etarho_ec.setVal(0.0);
    etarho_cc.setVal(0.0);
    psi.setVal(0.0);
    w0.setVal(0.);
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
    S_cc_old[lev].define(ba, dm, 1, 0);
    S_cc_new[lev].define(ba, dm, 1, 0);
    gpi[lev].define(ba, dm, AMREX_SPACEDIM, 0);
    dSdt[lev].define(ba, dm, 1, 0);
    rhcc_for_nodalproj[lev].define(ba, dm, 1, 1);

    pi[lev].define(convert(ba, nodal_flag), dm, 1, 0);  // nodal

    sold[lev].setVal(0.);
    snew[lev].setVal(0.);
    uold[lev].setVal(0.);
    unew[lev].setVal(0.);
    S_cc_old[lev].setVal(0.);
    S_cc_new[lev].setVal(0.);
    gpi[lev].setVal(0.);
    dSdt[lev].setVal(0.);
    rhcc_for_nodalproj[lev].setVal(0.);
    pi[lev].setVal(0.);

    if (spherical == 1) {
        normal[lev].define(ba, dm, 3, 1);
        cell_cc_to_r[lev].define(ba, dm, 1, 0);
    }

    const Real* dx = geom[lev].CellSize();
    const Real* dx_fine = geom[max_level].CellSize();

    MultiFab& scal = sold[lev];
    MultiFab& vel = uold[lev];
    iMultiFab& cc_to_r = cell_cc_to_r[lev];

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tilebox = mfi.tilebox();
        const int* lo = tilebox.loVect();
        const int* hi = tilebox.hiVect();

        if (!spherical) {
            const Array4<Real> scal_arr = scal.array(mfi);
            const Array4<Real> vel_arr = vel.array(mfi);

            InitLevelData(lev, t_old, mfi, scal_arr, vel_arr);

        } else {
#if (AMREX_SPACEDIM == 3)
            const auto dx_fine_vec = geom[max_level].CellSizeArray();
            const auto dx_lev = geom[lev].CellSizeArray();

            InitBaseStateMapSphr(lev, mfi, dx_fine_vec, dx_lev);
#endif
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (spherical) {
        InitLevelDataSphr(lev, t_old, scal, vel);
    }
#endif
}

void Maestro::InitIter() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitIter()", InitIter);

    // advance the solution by dt

    AdvanceTimeStep(true);
}
