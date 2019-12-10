
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

	if (plot_int > 0) {

		// Need to fill normal vector to compute velrc in plotfile
		if (spherical) { MakeNormal(); }

		Print() << "\nWriting plotfile plt_InitData after InitData" << std::endl;
		WritePlotFile(9999999,t_old,0,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);
	}

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	if (spherical == 1) {
		MakeNormal();
		MakeCCtoRadii();
	}

	// make gravity
	make_grav_cell(grav_cell_old.dataPtr(),
	               rho0_old.dataPtr(),
	               r_cc_loc.dataPtr(),
	               r_edge_loc.dataPtr());

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

	// first compute cutoff coordinates using initial density profile
	compute_cutoff_coords(rho0_old.dataPtr());

	// compute gravity
	make_grav_cell(grav_cell_old.dataPtr(),
	               rho0_old.dataPtr(),
	               r_cc_loc.dataPtr(),
	               r_edge_loc.dataPtr());

	// compute p0 with HSE
	enforce_HSE(rho0_old.dataPtr(),
	            p0_old.dataPtr(),
	            grav_cell_old.dataPtr(),
	            r_cc_loc.dataPtr(),
	            r_edge_loc.dataPtr());

	// set p0^{-1} = p0_old
	for (int i=0; i<p0_old.size(); ++i) {
		p0_nm1[i] = p0_old[i];
	}

    rhoX0_old.resize( (max_radial_level+1)*nr_fine*NumSpec);
    rhoX0_new.resize( (max_radial_level+1)*nr_fine*NumSpec);
	rhoX0_old.shrink_to_fit();
	rhoX0_new.shrink_to_fit();

    make_rhoX0(s0_init.dataPtr(), rhoX0_old.dataPtr());


	// set some stuff to zero
	std::fill(etarho_ec.begin(), etarho_ec.end(), 0.);
	std::fill(etarho_cc.begin(), etarho_cc.end(), 0.);
	std::fill(psi.begin(), psi.end(), 0.);
	std::fill(w0.begin(), w0.end(), 0.);


	// free memory in s0_init and p0_init by swapping it
	// with an empty vector that will go out of scope
	Vector<Real> s0_swap, p0_swap;
	std::swap(s0_swap,s0_init);
	std::swap(p0_swap,p0_init);

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
	S_cc_old          [lev].define(ba, dm,              1,    0);
	S_cc_new          [lev].define(ba, dm,              1,    0);
	gpi               [lev].define(ba, dm, AMREX_SPACEDIM,    0);
	dSdt              [lev].define(ba, dm,              1,    0);
	rhcc_for_nodalproj[lev].define(ba, dm,              1,    1);

	pi[lev].define(convert(ba,nodal_flag), dm, 1, 0); // nodal

	sold              [lev].setVal(0.);
	snew              [lev].setVal(0.);
	uold              [lev].setVal(0.);
	unew              [lev].setVal(0.);
	S_cc_old          [lev].setVal(0.);
	S_cc_new          [lev].setVal(0.);
	gpi               [lev].setVal(0.);
	dSdt              [lev].setVal(0.);
	rhcc_for_nodalproj[lev].setVal(0.);
	pi                [lev].setVal(0.);

	if (spherical == 1) {
		normal      [lev].define(ba, dm, 3, 1);
		cell_cc_to_r[lev].define(ba, dm, 1, 0);
	}

	const Real* dx = geom[lev].CellSize();
	const Real* dx_fine = geom[max_level].CellSize();

	MultiFab& scal = sold[lev];
	MultiFab& vel = uold[lev];
	MultiFab& cc_to_r = cell_cc_to_r[lev];

	// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(scal, true); mfi.isValid(); ++mfi)
	{
		const Box& tilebox = mfi.tilebox();
		const int* lo  = tilebox.loVect();
		const int* hi  = tilebox.hiVect();

		if (spherical == 0) {
			initdata(&lev, &t_old, ARLIM_3D(lo), ARLIM_3D(hi),
			         BL_TO_FORTRAN_FAB(scal[mfi]),
			         BL_TO_FORTRAN_FAB(vel[mfi]),
			         s0_init.dataPtr(), p0_init.dataPtr(),
			         ZFILL(dx));
		} else {
			init_base_state_map_sphr(BL_TO_FORTRAN_3D(cc_to_r[mfi]),
			                         ZFILL(dx_fine),
			                         ZFILL(dx));

			initdata_sphr(&t_old, ARLIM_3D(lo), ARLIM_3D(hi),
			              BL_TO_FORTRAN_FAB(scal[mfi]),
			              BL_TO_FORTRAN_FAB(vel[mfi]),
			              s0_init.dataPtr(), p0_init.dataPtr(),
			              ZFILL(dx),
			              r_cc_loc.dataPtr(), r_edge_loc.dataPtr(),
			              BL_TO_FORTRAN_3D(cc_to_r[mfi]));
		}
	}
}



void Maestro::InitIter ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::InitIter()",InitIter);

	// advance the solution by dt

		AdvanceTimeStep(true);


}
