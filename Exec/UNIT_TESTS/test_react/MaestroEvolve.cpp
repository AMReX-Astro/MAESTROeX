
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{

	Vector<MultiFab> rho_omegadot(finest_level+1);
	Vector<MultiFab>     rho_Hnuc(finest_level+1);
	Vector<MultiFab>     rho_Hext(finest_level+1);

	for (int lev=0; lev<=finest_level; ++lev) {
		// cell-centered MultiFabs
		rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec,    0);
		rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1,    0);
		rho_Hext    [lev].define(grids[lev], dmap[lev],       1,    0);
		rho_Hext[lev].setVal(0.);
	}

	auto dbo = do_burning;
	auto dho = do_heating;

	// Mode 1: No burning, no heating
	do_burning = false;
	do_heating = false;
	React(sold,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,dt);
	// write

	// Mode 2: Burning without heating
	do_burning = true;
	do_heating = false;
	React(sold,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,dt);
	// write

	// Mode 3: Heating without burning
	do_burning = false;
	do_heating = true;
	React(sold,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,dt);
	// write

	// Mode 4: Burning and heating
	do_burning = true;
	do_heating = true;
	React(sold,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,dt);
	// write

	// Explore ten orders of magnitude of the time domain using user inputs.
	do_burning = dbo;
	do_heating = dho;
	for (auto i=0; i < react_its; ++i) {
		React(sold,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,dt);
		// write
	}

}
