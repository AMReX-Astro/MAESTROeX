
#include <Maestro.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStep (bool is_initIter) {

	// timer for profiling
	BL_PROFILE_VAR("Maestro::AdvanceTimeStep()",AdvanceTimeStep);

	// cell-centered MultiFabs needed within the AdvanceTimeStep routine
	Vector<MultiFab>      Tcoeff1(finest_level+1);
	Vector<MultiFab>      Tcoeff2(finest_level+1);
	Vector<MultiFab>      hcoeff1(finest_level+1);
	Vector<MultiFab>     Xkcoeff1(finest_level+1);
	Vector<MultiFab>      pcoeff1(finest_level+1);
	Vector<MultiFab>      hcoeff2(finest_level+1);
	Vector<MultiFab>     Xkcoeff2(finest_level+1);
	Vector<MultiFab>      pcoeff2(finest_level+1);

	// wallclock time
	const Real strt_total = ParallelDescriptor::second();

	Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
	        << " DT = " << dt << std::endl << std::endl;

	for (int lev=0; lev<=finest_level; ++lev) {
		// cell-centered MultiFabs
		Tcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
		Tcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
		hcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
		Xkcoeff1    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
		pcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
		hcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
		Xkcoeff2    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
		pcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
	}

	// thermal is the forcing for rhoh or temperature
	MakeThermalCoeffs(sold,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1);

    // on the first step, just copy coeffs for the time centering
	for (int lev=0; lev<=finest_level; ++lev) {
		MultiFab::Copy(Tcoeff2[lev],Tcoeff1[lev],1,1,1,1);
		MultiFab::Copy(hcoeff2[lev],hcoeff1[lev],1,1,1,1);
		MultiFab::Copy(Xkcoeff2[lev],Xkcoeff1[lev],1,1,NumSpec,1);
		MultiFab::Copy(pcoeff2[lev],pcoeff1[lev],1,1,1,1);
	}

    // diffuse the enthalpy
    Print() << '... conducting' << std::endl;
	ThermalConduct(sold,snew,hcoeff1,Xkcoeff1,pcoeff1,hcoeff2,Xkcoeff2,pcoeff2);

	// now update temperature
	TfromRhoH(snew,p0_old);

	Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
	        << " DT = " << dt << std::endl;

	// wallclock time
	Real end_total = ParallelDescriptor::second() - strt_total;

	// print wallclock time
	ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
	if (maestro_verbose > 0) {
		Print() << "Time to advance time step: " << end_total << '\n';
	}

}
