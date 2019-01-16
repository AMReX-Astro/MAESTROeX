
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

}

// fill in multifab and base state data
void
Maestro::InitData ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::InitData()",InitData);

	Print() << "Calling InitData()" << std::endl;


	// init_base_state(s0_init.dataPtr(),p0_init.dataPtr(),rho0_old.dataPtr(),
	//                 rhoh0_old.dataPtr(),p0_old.dataPtr(),tempbar.dataPtr(),
	//                 tempbar_init.dataPtr());

	// calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
	// that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
	InitFromScratch(t_old);

	// reset tagging array to include buffer zones
	TagArray();

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	AverageDown(uold,0,AMREX_SPACEDIM);
	FillPatch(t_old,uold,uold,uold,0,0,AMREX_SPACEDIM,0,bcs_u);

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
	unew              [lev].define(ba, dm, AMREX_SPACEDIM, ng_s);
	gpi               [lev].define(ba, dm, AMREX_SPACEDIM,    1);
	rhcc_for_nodalproj[lev].define(ba, dm,              1,    1);

	pi[lev].define(convert(ba,nodal_flag), dm, 1, 1); // nodal

	sold              [lev].setVal(0.);
	snew              [lev].setVal(0.);
	uold              [lev].setVal(0.);
	unew              [lev].setVal(0.);
	gpi               [lev].setVal(0.);
	rhcc_for_nodalproj[lev].setVal(0.);
	pi                [lev].setVal(0.);

	const Real* dx = geom[lev].CellSize();
	const Real* dx_fine = geom[max_level].CellSize();

	int project_type;
	get_project_type(&project_type);

	if (project_type == 1) {

		MultiFab& vel = uold[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(vel, true); mfi.isValid(); ++mfi)
		{
			const Box& tilebox = mfi.tilebox();
			const int* lo  = tilebox.loVect();
			const int* hi  = tilebox.hiVect();

			init_vel(ARLIM_3D(lo), ARLIM_3D(hi),
			         BL_TO_FORTRAN_3D(vel[mfi]),
			         ZFILL(dx));
		}
	}
}
