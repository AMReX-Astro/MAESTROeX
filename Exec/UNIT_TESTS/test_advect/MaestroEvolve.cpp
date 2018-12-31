
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::Evolve()",Evolve);

	Print() << "Calling Evolve()" << std::endl;

	// check to make sure spherical is only used for 3d
	if (spherical == 1 && AMREX_SPACEDIM != 3) {
		Abort ("spherical = 1 and dm != 3");
	}

	// there's only one level
	const int lev = 0;

	Vector<MultiFab>           s1(finest_level+1);
	Vector<MultiFab>   scal_force(finest_level+1);
	Vector<MultiFab>   etarhoflux(finest_level+1);
	Vector<std::array< MultiFab, AMREX_SPACEDIM > >  umac(finest_level+1);
	Vector<std::array< MultiFab, AMREX_SPACEDIM > > sedge(finest_level+1);
	Vector<std::array< MultiFab, AMREX_SPACEDIM > > sflux(finest_level+1);
	Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
	Vector<Real> rho0_predicted_edge( (max_radial_level+1)*(nr_fine+1) );

	rho0_predicted_edge.shrink_to_fit();

	s1          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
	s1[lev].setVal(0.);
	scal_force  [lev].define(grids[lev], dmap[lev],   Nscal,    1);

	AMREX_D_TERM(etarhoflux[lev].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
	             etarhoflux[lev].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
	             etarhoflux[lev].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

	// face-centered arrays of MultiFabs
	AMREX_D_TERM(umac [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1,     1); ,
	             umac [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1,     1); ,
	             umac [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1,     1); );
	AMREX_D_TERM(sedge[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0); ,
	             sedge[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0); ,
	             sedge[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0); );
	AMREX_D_TERM(sflux[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0); ,
	             sflux[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0); ,
	             sflux[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0); );

	// initialize umac
	for (int d=0; d < AMREX_SPACEDIM; ++d)
		umac[lev][d].setVal(0.);

#if (AMREX_SPACEDIM == 3)
	for (int lev=0; lev<=finest_level; ++lev) {
		w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
		w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
		w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
	}
	if (spherical == 1) {
		for (int lev=0; lev<=finest_level; ++lev) {
			w0_force_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
		}
	}
#endif

	// Store the initial density here
	for (int lev=0; lev<=finest_level; ++lev) {
		s1[lev].setVal(0.);
		MultiFab::Copy(s1[lev],sold[lev],Rho,Rho,1,ng_s);
	}

	const Real* dx = geom[lev].CellSize();

	dt = cfl * dx[0];

	for (istep = start_step; istep <= max_step && t_old < stop_time; ++istep)
	{

		t_old = t_new;

		DensityAdvance(1,sold,snew,sedge,sflux,scal_force,etarhoflux,umac,w0mac,rho0_predicted_edge);

		rho0_old = rho0_new;

		t_new = t_old + dt;

		// move new state into old state by swapping pointers
		for (int lev=0; lev<=finest_level; ++lev) {
			std::swap(    sold[lev],     snew[lev]);
			std::swap(    uold[lev],     unew[lev]);
			std::swap(S_cc_old[lev], S_cc_new[lev]);

			std::swap( rho0_old, rho0_new);
			std::swap(rhoh0_old,rhoh0_new);
			std::swap(   p0_nm1,   p0_old);
			std::swap(   p0_old,   p0_new);

			std::swap(beta0_old,beta0_new);
			std::swap(gamma1bar_old,gamma1bar_new);
			std::swap(grav_cell_old,grav_cell_new);
		}

	}
}
