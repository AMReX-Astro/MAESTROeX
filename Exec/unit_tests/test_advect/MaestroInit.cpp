
#include <Maestro.H>
#include <Maestro_F.H>
#include <AMReX_VisMF.H>
using namespace amrex;


// initialize AMR data
void
Maestro::Init ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::Init()",Init);

	Print() << "Calling Init()" << std::endl;

	start_step = 1;

	// fill in multifab and base state data
	InitData();

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	// compute initial time step
	FirstDt();

	if (stop_time >= 0. && t_old+dt > stop_time) {
		dt = std::min(dt,stop_time-t_old);
		Print() << "Stop time limits dt = " << dt << std::endl;
	}

	dtold = dt;
	t_new = t_old + dt;
}

// fill in multifab and base state data
void
Maestro::InitData ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::InitData()",InitData);

	Print() << "Calling InitData()" << std::endl;

	// read in model file and fill in s0_init and p0_init for all levels

	for (auto lev = 0; lev <= base_geom.max_radial_level; ++lev) {
	    InitBaseState(rho0_old, rhoh0_old,
			  p0_old, lev);
	}
	    
	// calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
	// that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
	InitFromScratch(t_old);

	// reset tagging array to include buffer zones
	TagArray();

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	// average down data and fill ghost cells
	AverageDown(sold,0,Nscal);
	FillPatch(t_old,sold,sold,sold,0,0,Nscal,0,bcs_s);

	// free memory in s0_init and p0_init by swapping it
	// with an empty vector that will go out of scope
	Vector<Real> s0_swap, p0_swap;
	std::swap(s0_swap,s0_init);
	std::swap(p0_swap,p0_init);

	// first compute cutoff coordinates using initial density profile
	compute_cutoff_coords(rho0_old.dataPtr());

	// set rho0 to be the average
	Average(sold,rho0_old,Rho);
	compute_cutoff_coords(rho0_old.dataPtr());

	// call eos with r,p as input to recompute T,h
	TfromRhoP(sold,p0_old,1);

	// set rhoh0 to be the average
	Average(sold,rhoh0_old,RhoH);

}

// During initialization of a simulation, Maestro::InitData() calls
// AmrCore::InitFromScratch(), which calls
// a MakeNewGrids() function that repeatedly calls this function to build
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                       const DistributionMapping& dm)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeNewLevelFromScratch()",MakeNewLevelFromScratch);

	sold              [lev].define(ba, dm,          Nscal, ng_s);
	snew              [lev].define(ba, dm,          Nscal, ng_s);
	uold              [lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
	S_cc_old          [lev].define(ba, dm,              1,    0);
	S_cc_new          [lev].define(ba, dm,              1,    0);
	gpi               [lev].define(ba, dm, AMREX_SPACEDIM,    0);
    w0_cart           [lev].define(ba, dm, AMREX_SPACEDIM,    2);
	rhcc_for_nodalproj[lev].define(ba, dm,              1,    1);
	pi[lev].define(convert(ba,nodal_flag), dm, 1, 0); // nodal

	sold              [lev].setVal(0.);
	snew              [lev].setVal(0.);
	uold              [lev].setVal(0.);
	S_cc_old          [lev].setVal(0.);
	S_cc_new          [lev].setVal(0.);
	gpi               [lev].setVal(0.);
    w0_cart           [lev].setVal(0.);
	rhcc_for_nodalproj[lev].setVal(0.);
	pi                [lev].setVal(0.);

	// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		InitLevelData(lev, t_old, mfi, sold[lev].array(mfi), uold[lev].array(mfi));
	}
}
