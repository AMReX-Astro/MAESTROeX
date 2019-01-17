
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{

	// -------------------------------------------------------------------------
	//  loop over all possible directions and all ppm_types
	// -------------------------------------------------------------------------

	for (auto i = 0; i < AMREX_SPACEDIM; i++) {
		for (auto j = -1; j <= 1; j+=2) {

			t_old = 0.0;
			t_new = 0.0;
			istep = 0;

			// reset the density
			InitFromScratch(t_old);

			// average down data and fill ghost cells
			AverageDown(sold,0,Nscal);
			FillPatch(t_old,sold,sold,sold,0,0,Nscal,0,bcs_s);
			AverageDown(uold,0,AMREX_SPACEDIM);
			FillPatch(t_old,uold,uold,uold,0,0,AMREX_SPACEDIM,0,bcs_u);

			compute_cutoff_coords(rho0_old.dataPtr());

			// the base state will not carry any information in this test problem
			std::fill(rho0_old.begin(),  rho0_old.end(),  0.);
			std::fill(rhoh0_old.begin(), rhoh0_old.end(), 0.);

			// set tempbar to be the average
			Average(sold,tempbar,Temp);
			for (int i=0; i<tempbar.size(); ++i) {
				tempbar_init[i] = tempbar[i];
			}

			// set p0^{-1} = p0_old
			for (int i=0; i<p0_old.size(); ++i) {
				p0_nm1[i] = p0_old[i];
			}

			Print() << "\nInitial velocity = ";
			if (i == 0) {
				Print() << "(" << j << ", 0)" << std::endl;
			} else {
				Print() << "(0, " << j << ")" << std::endl;
			}

			// check to make sure spherical is only used for 3d
			if (spherical == 1 && AMREX_SPACEDIM != 3) {
				Abort ("spherical = 1 and dm != 3");
			}

			Vector<MultiFab>       s_orig(finest_level+1);
			Vector<MultiFab>      s_final(finest_level+1);
			Vector<MultiFab>        error(finest_level+1);
			Vector<MultiFab>   scal_force(finest_level+1);
			Vector<MultiFab>   etarhoflux(finest_level+1);
			Vector<std::array< MultiFab, AMREX_SPACEDIM > >  umac(finest_level+1);
			Vector<std::array< MultiFab, AMREX_SPACEDIM > > sedge(finest_level+1);
			Vector<std::array< MultiFab, AMREX_SPACEDIM > > sflux(finest_level+1);
			Vector<MultiFab> w0_force_cart(finest_level+1);
			Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
			Vector<Real> w0_force        ( (max_radial_level+1)*nr_fine );
			Vector<Real> rho0_predicted_edge( (max_radial_level+1)*(nr_fine+1) );

			Vector<Real> abs_norm(finest_level+1);
			Vector<Real> rel_norm(finest_level+1);

			rho0_predicted_edge.shrink_to_fit();
			w0_force.shrink_to_fit();

			for (int lev=0; lev<=finest_level; ++lev) {
				s_orig[lev].define(grids[lev], dmap[lev],   1, ng_s);
				s_orig[lev].setVal(0.);
				s_final[lev].define(grids[lev], dmap[lev],   1, ng_s);
				s_final[lev].setVal(0.);
				error[lev].define(grids[lev], dmap[lev],   1, ng_s);
				error[lev].setVal(0.);
				scal_force[lev].define(grids[lev], dmap[lev],   Nscal,    1);
				scal_force[lev].setVal(0.);

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


				etarhoflux[lev].setVal(0.);

				// initialize umac
				for (int d=0; d < AMREX_SPACEDIM; ++d) {
					umac[lev][d].setVal(0.);
					sedge[lev][d].setVal(0.);
					sflux[lev][d].setVal(0.);
				}
			}

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

			// initialize MultiFabs and Vectors to ZERO
			for (int lev=0; lev<=finest_level; ++lev) {
				for (int d=0; d<AMREX_SPACEDIM; ++d) {
					w0mac[lev][d].setVal(0.);
				}
			}
			if (spherical == 1) {
				for (int lev=0; lev<=finest_level; ++lev) {
					w0_force_cart[lev].setVal(0.);
				}
			}
#endif

			std::fill(w0_force.begin(), w0_force.end(), 0.);
			std::fill(rho0_predicted_edge.begin(), rho0_predicted_edge.end(), 0.);

			// Store the initial density here
			for (int lev=0; lev<=finest_level; ++lev)
				MultiFab::Copy(s_orig[lev],sold[lev],Rho,0,1,0);

			const Real* dx = geom[finest_level].CellSize();

			dt = cfl * dx[0];

			// initialize the velocity field -- it is unity in the
			// direction of propagation a negative itest_dir indicates
			// negative velocity
			for (int lev=0; lev<=finest_level; ++lev) {
				for (int d=0; d < AMREX_SPACEDIM; ++d)
					umac[lev][d].setVal(0.0);

				umac[lev][i].setVal(double(j));
			}

			AverageDownFaces(umac);

			Print() << "original density = " << s_orig[0].norm2() << std::endl;

			// advance the density using the constant velocity field
			for (istep = start_step; istep <= max_step && t_old < stop_time; ++istep)
			{

				t_old = t_new;

				// set etarhoflux to zero
				for (int lev=0; lev<=finest_level; ++lev)
					etarhoflux[lev].setVal(0.);

				// set sedge and sflux to zero
				for (int lev=0; lev<=finest_level; ++lev) {
					for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
						sedge[lev][idim].setVal(0.);
						sflux[lev][idim].setVal(0.);
					}
				}

				DensityAdvance(1,sold,snew,sedge,sflux,scal_force,etarhoflux,umac,w0mac,rho0_predicted_edge);

				p0_new = p0_old;

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

				t_new = t_old + dt;

				if (t_new + dt > stop_time)
					dt = stop_time - t_new;

			}

			// Store the final density here
			for (int lev=0; lev<=finest_level; ++lev)
				MultiFab::Copy(s_final[lev],snew[lev],Rho,0,1,0);

			Print() << "final density = " << s_final[0].norm2() << std::endl;

			// compare the initial and final density
			// compute dens_final - dens_orig
			for (int lev=0; lev<=finest_level; ++lev) {
				MultiFab::Copy(error[lev],s_final[lev],0,0,1,0);
				MultiFab::Subtract(error[lev],s_orig[lev],0,0,1,0);

				abs_norm[lev] = error[lev].norm2();

				MultiFab::Divide(error[lev],s_orig[lev],0,0,1,0);

				rel_norm[lev] = error[lev].norm2();

				Print() << "\tAbs norm = " << abs_norm[lev] << "  Rel norm = " << rel_norm[lev] << std::endl;

			}

		}
	}
}
