
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::Evolve()",Evolve);

	Print() << "Calling Evolve()" << std::endl;

	// coefficients for diffusion
	Vector<MultiFab>      Tcoeff1(finest_level+1);
	Vector<MultiFab>      hcoeff1(finest_level+1);
	Vector<MultiFab>     Xkcoeff1(finest_level+1);
	Vector<MultiFab>      pcoeff1(finest_level+1);
	Vector<MultiFab>      Tcoeff2(finest_level+1);
	Vector<MultiFab>      hcoeff2(finest_level+1);
	Vector<MultiFab>     Xkcoeff2(finest_level+1);
	Vector<MultiFab>      pcoeff2(finest_level+1);
	Vector<MultiFab>     analytic(finest_level+1);
	Vector<MultiFab>        error(finest_level+1);

	for (int lev=0; lev<=finest_level; ++lev) {
		// coefficients for diffusion
		Tcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
		Tcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
		hcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
		hcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
		Xkcoeff1    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
		Xkcoeff2    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
		pcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
		pcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
		analytic    [lev].define(grids[lev], dmap[lev],       1,    0);
		error       [lev].define(grids[lev], dmap[lev],       1,    0);
		Tcoeff1[lev].setVal(0.);
		Tcoeff2[lev].setVal(0.);
		hcoeff1[lev].setVal(0.);
		hcoeff2[lev].setVal(0.);
		Xkcoeff1[lev].setVal(0.);
		Xkcoeff2[lev].setVal(0.);
		pcoeff1[lev].setVal(0.);
		pcoeff2[lev].setVal(0.);
		analytic[lev].setVal(0.);
		error[lev].setVal(0.);
	}

	// copy the state data to snew
	for (int lev=0; lev<=finest_level; ++lev)
		MultiFab::Copy(snew[lev],sold[lev],0,0,Nscal,ng_s);

	// calculate analytic solution
	for (int lev = 0; lev <= finest_level; ++lev) {

		MultiFab& analytic_mf = analytic[lev];
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(analytic_mf, true); mfi.isValid(); ++mfi ) {

			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			make_analytic_solution(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()), BL_TO_FORTRAN_3D(analytic_mf[mfi]), ZFILL(dx), t_new);
		}
	}

	// This problem uses a custom WritePlotFile, but as it's a member function of
	// the Maestro class, it has to use the same prototype as the original.
	// We shall therefore create a dummy variable to fill up all the variables
	// passed into the function that won't be used.
	Vector<Real> dummy;

	// dump initial state
	WritePlotFile(0,t_new,dt,dummy,dummy,dummy,dummy,sold,analytic,error);

	for (istep = start_step; istep <= max_step && t_old < stop_time; ++istep)
	{

		Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
		        << " DT = " << dt << std::endl << std::endl;

		t_new = t_old + dt;

		// AdvanceTimeStep
		{
			// thermal is the forcing for rhoh or temperature
			MakeThermalCoeffs(sold,Tcoeff1,hcoeff1,Xkcoeff1,pcoeff1);

			// on the first step, just copy coeffs for the time centering
			if (istep == 1) {
				for (int lev=0; lev<=finest_level; ++lev) {
					MultiFab::Copy(Tcoeff2[lev],Tcoeff1[lev],0,0,1,1);
					MultiFab::Copy(hcoeff2[lev],hcoeff1[lev],0,0,1,1);
					MultiFab::Copy(Xkcoeff2[lev],Xkcoeff1[lev],0,0,NumSpec,1);
					MultiFab::Copy(pcoeff2[lev],pcoeff1[lev],0,0,1,1);
				}
			}

			// diffuse the enthalpy
			Print() << "... conducting" << std::endl;
			ThermalConduct(sold,snew,hcoeff1,Xkcoeff1,pcoeff1,hcoeff2,Xkcoeff2,pcoeff2);

			// now update temperature
			TfromRhoH(snew,p0_new);
		}

		// calculate analytic solution
		for (int lev = 0; lev <= finest_level; ++lev) {

			MultiFab& analytic_mf = analytic[lev];
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(analytic_mf, true); mfi.isValid(); ++mfi ) {

				const Box& tileBox = mfi.tilebox();
				const Real* dx = geom[lev].CellSize();

				make_analytic_solution(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()), BL_TO_FORTRAN_3D(analytic_mf[mfi]), ZFILL(dx), t_new);
			}
		}

		// calculate the error
		MultiFab::Copy(error[finest_level],snew[finest_level],RhoH,0,1,0);
		MultiFab::Divide(error[finest_level],snew[finest_level],Rho,0,1,0);
		MultiFab::Subtract(error[finest_level],analytic[finest_level],0,0,1,0);

		Real L1norm = error[finest_level].norm1() / analytic[finest_level].norm1();
		Real L2norm = error[finest_level].norm2() / analytic[finest_level].norm2();

		Print() << "\nTime = " << t_new << ", L1norm = " << L1norm << ", L2norm = " << L2norm << '\n' << std::endl;

		// write a plotfile
        if (plot_int > 0 && ( (istep % plot_int == 0) ||
                              (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
                              (istep == max_step) ) || (t_old >= stop_time) )
        {
            Print() << "\nWriting plotfile " << istep << std::endl;
            WritePlotFile(istep,t_new,dt,dummy,dummy,dummy,dummy,sold,analytic,error);
        }

		t_old = t_new;
		p0_old = p0_new;

		// fill the mfs for the next timestep by switching pointers
		for (int lev=0; lev<=finest_level; ++lev) {
            std::swap(sold[lev], snew[lev]);
			std::swap(Tcoeff2[lev], Tcoeff1[lev]);
			std::swap(hcoeff2[lev], hcoeff1[lev]);
			std::swap(Xkcoeff2[lev], Xkcoeff1[lev]);
			std::swap(pcoeff2[lev], pcoeff1[lev]);
		}
	}
}
