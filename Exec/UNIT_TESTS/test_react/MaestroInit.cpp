
#include <Maestro.H>
#include <AMReX_VisMF.H>
#include <Problem_F.H>
using namespace amrex;


// initialize AMR data
// perform initial projection
// perform divu iters
// perform initial (pressure) iterations
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
	// FirstDt();
    get_min_timestep(&dt);

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

	init_base_state(s0_init.dataPtr(),p0_init.dataPtr(),rho0_old.dataPtr(),
	                rhoh0_old.dataPtr(),p0_old.dataPtr(),tempbar.dataPtr(),
	                tempbar_init.dataPtr());

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

	// set tempbar to be the average
	Average(sold,tempbar,Temp);
	for (int i=0; i<tempbar.size(); ++i)
		tempbar_init[i] = tempbar[i];

    for (int lev=0; lev<=finest_level; ++lev)
		MultiFab::Copy(snew[lev],sold[lev],0,0,Nscal,ng_s);

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

	sold              [lev].setVal(0.);
	snew              [lev].setVal(0.);

	const Real* dx = geom[lev].CellSize();
	const Real* dx_fine = geom[max_level].CellSize();

	MultiFab& scal = sold[lev];

    const Box& domainBox = geom[lev].Domain();

	// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(scal, true); mfi.isValid(); ++mfi)
	{
		const Box& tilebox = mfi.tilebox();
		const int* lo  = tilebox.loVect();
		const int* hi  = tilebox.hiVect();

		initdata_thermal(ARLIM_3D(lo), ARLIM_3D(hi),
		         BL_TO_FORTRAN_FAB(scal[mfi]),
		         ARLIM_3D(domainBox.loVect()), ARLIM_3D(domainBox.hiVect()),
		         s0_init.dataPtr(), p0_init.dataPtr(),
		         ZFILL(dx));
	}
}
