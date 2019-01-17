
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{

	Print() << "Running tests..." << std::endl;

	Vector<MultiFab> error(finest_level+1);

	// number of EoS tests that do_tests performs
	const int n_tests = 5;

	// setup error MultiFab
	for (int lev=0; lev<=finest_level; ++lev) {
		error[lev].define(grids[lev], dmap[lev], n_tests, ng_s);
		error[lev].setVal(0.);
	}


	// do the tests
	for (int lev=0; lev<=finest_level; ++lev) {

		MultiFab& scal = sold[lev];
		MultiFab& error_mf = error[lev];
		const Box& domainBox = geom[lev].Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(scal, true); mfi.isValid(); ++mfi)
		{
			const Box& tilebox = mfi.tilebox();
			const int* lo  = tilebox.loVect();
			const int* hi  = tilebox.hiVect();

			do_tests(ARLIM_3D(lo), ARLIM_3D(hi),
			         BL_TO_FORTRAN_3D(scal[mfi]),
			         BL_TO_FORTRAN_FAB(error_mf[mfi]),
			         domainBox.loVect(), domainBox.hiVect());
		}
	}

    // print the errors
    {
        int lev = finest_level;

        Print() << "\nError in T   from rho, h = " << error[lev].norm2(0) << std::endl;
        Print() << "Error in rho from T  , p = " << error[lev].norm2(1) << std::endl;
        Print() << "Error in T   from rho, p = " << error[lev].norm2(2) << std::endl;
        Print() << "Error in T   from rho, e = " << error[lev].norm2(3) << std::endl;
        Print() << "Error in T   from p  , s = " << error[lev].norm2(4) << std::endl;
    }

}
