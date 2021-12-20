
#include <AMReX_VisMF.H>
#include <Maestro.H>
#include <extern_parameters.H>

using namespace amrex;

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
    // FirstDt();
    dt = min_time_step;

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

    // set tempbar to be the average
    Average(sold, tempbar, Temp);
    tempbar_init.copy(tempbar);

    for (int lev = 0; lev <= finest_level; ++lev)
        MultiFab::Copy(snew[lev], sold[lev], 0, 0, Nscal, ng_s);
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

    sold[lev].setVal(0.);
    snew[lev].setVal(0.);

    const Real* dx = geom[lev].CellSize();
    const Real* dx_fine = geom[max_level].CellSize();

    const Box& domainBox = geom[lev].Domain();
    const auto dom_lo = domainBox.loVect3d();
    const auto dom_hi = domainBox.hiVect3d();

    const auto xn_hi = dom_hi[2] - dom_lo[2];

    const auto dlogrho =
        (std::log10(dens_max / dens_min)) / (dom_hi[0] - dom_lo[0]);
    const auto dlogT =
        (std::log10(temp_max / temp_min)) / (dom_hi[1] - dom_lo[1]);

    const auto temp_min_l = temp_min;
    const auto dens_min_l = dens_min;

    GpuArray<Real, NumSpec> xn_zone;
    // FIXME: need to allocate the xns here

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(sold[lev], true); mfi.isValid(); ++mfi) {
        const Box& tilebox = mfi.tilebox();

        const Array4<Real> scal = sold[lev].array(mfi);

        ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // set the temperature
            const auto temp_zone =
                std::pow(10.0, std::log10(temp_min_l) + Real(j) * dlogT);

            // set the density
            const auto dens_zone =
                std::pow(10.0, std::log10(dens_min_l) + Real(i) * dlogrho);

            // call the eos with rho, temp & X as inputs
            eos_t eos_state;
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = xn_zone[comp];
            }

            eos(eos_input_rt, eos_state);

            // initialize this element of the state
            scal(i, j, k, Rho) = dens_zone;
            scal(i, j, k, RhoH) = dens_zone * eos_state.h;
            scal(i, j, k, Temp) = temp_zone;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) = dens_zone * xn_zone[comp];
            }
        });
    }
}
