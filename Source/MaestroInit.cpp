
#include <Maestro.H>
#include <AMReX_VisMF.H>
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

	if (restart_file == "") {

		start_step = 1;

		// fill in multifab and base state data
		InitData();

		if (plot_int > 0 || plot_deltat > 0) {

			// Need to fill normal vector to compute velrc in plotfile
			if (spherical) { MakeNormal(); }

			Print() << "\nWriting plotfile "<< plot_base_name << "InitData after InitData" << std::endl;
			WritePlotFile(plotInitData,t_old,0,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);

		} else if (small_plot_int > 0 || small_plot_deltat > 0) {

			// Need to fill normal vector to compute velrc in plotfile
			if (spherical) { MakeNormal(); }

			Print() << "\nWriting small plotfile "<< small_plot_base_name << "InitData after InitData" << std::endl;
			WriteSmallPlotFile(plotInitData,t_old,0,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);

		}
	}
	else {
		Print() << "Initializing from checkpoint " << restart_file << std::endl;

		// read in checkpoint file
		// this builds (defines) and fills the following MultiFabs:
		//
		// snew, unew, gpi, dSdt, S_cc_new
		//
		// and also fills in the 1D arrays:
		//
		// rho0_new, p0_new, gamma1bar_new, rhoh0_new, beta0_new, psi, tempbar, etarho_cc, tempbar_init
		ReadCheckPoint();

		// build (define) the following MultiFabs (that weren't read in from checkpoint):
		// snew, unew, S_cc_new, rhcc_for_nodalproj, normal, pi
		for (int lev=0; lev<=finest_level; ++lev) {
			snew              [lev].define(grids[lev], dmap[lev],          Nscal, ng_s);
			unew              [lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
			S_cc_new          [lev].define(grids[lev], dmap[lev],              1,    0);
			rhcc_for_nodalproj[lev].define(grids[lev], dmap[lev],              1,    1);
			if (spherical == 1) {
				normal[lev].define(grids[lev], dmap[lev], 3, 1);
				cell_cc_to_r[lev].define(grids[lev], dmap[lev], 1, 0);
			}
			pi[lev].define(convert(grids[lev],nodal_flag), dmap[lev], 1, 0); // nodal

		}
	}

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	compute_cutoff_coords(rho0_old.dataPtr());

	if (spherical == 1) {
		MakeNormal();
		MakeCCtoRadii();
	}

	if (do_sponge) {
		if (use_exact_base_state) {
			init_sponge_irreg(rho0_old.dataPtr(),r_cc_loc.dataPtr(),r_edge_loc.dataPtr());
		} else {
			init_sponge(rho0_old.dataPtr());
		}
	}

	// make gravity
	make_grav_cell(grav_cell_old.dataPtr(),
	               rho0_old.dataPtr(),
	               r_cc_loc.dataPtr(),
	               r_edge_loc.dataPtr());

	if (restart_file == "") {

		// compute gamma1bar
		MakeGamma1bar(sold,gamma1bar_old,p0_old);

		// compute beta0
		if (use_exact_base_state) {
			make_beta0_irreg(beta0_old.dataPtr(),
			                 rho0_old.dataPtr(),
			                 p0_old.dataPtr(),
			                 gamma1bar_old.dataPtr(),
			                 grav_cell_old.dataPtr(),
			                 r_cc_loc.dataPtr(),
			                 r_edge_loc.dataPtr());
		} else {
			make_beta0(beta0_old.dataPtr(),
			           rho0_old.dataPtr(),
			           p0_old.dataPtr(),
			           gamma1bar_old.dataPtr(),
			           grav_cell_old.dataPtr());
		}

		// set beta0^{-1} = beta0_old
		for (int i=0; i<beta0_old.size(); ++i) {
			beta0_nm1[i] = beta0_old[i];
		}

		// initial projection
		if (do_initial_projection) {
			Print() << "Doing initial projection" << std::endl;
			InitProj();

			if (plot_int > 0 || plot_deltat > 0) {
				Print() << "\nWriting plotfile " << plot_base_name << "after_InitProj after InitProj" << std::endl;

				WritePlotFile(plotInitProj,t_old,0,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);

			} else if (small_plot_int > 0 || small_plot_deltat > 0) {
				Print() << "\nWriting small plotfile " << small_plot_base_name << "after_InitProj after InitProj" << std::endl;

				WriteSmallPlotFile(plotInitProj,t_old,0,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);
			}
		}

		// compute initial time step
		FirstDt();

		// divu iters - also update dt at end of each divu_iter
		if (init_divu_iter > 0) {
			for (int i=1; i<=init_divu_iter; ++i) {
				Print() << "Doing initial divu iteration #" << i << std::endl;
				DivuIter(i);
			}

			if (plot_int > 0 || plot_deltat > 0) {
				Print() << "\nWriting plotfile " << plot_base_name << "after_DivuIter after final DivuIter" << std::endl;
				WritePlotFile(plotDivuIter,t_old,dt,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);
			} else if (small_plot_int > 0 || small_plot_deltat > 0) {
				Print() << "\nWriting small plotfile " << small_plot_base_name << "after_DivuIter after final DivuIter" << std::endl;
				WriteSmallPlotFile(plotDivuIter,t_old,dt,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);
			}
		}

		if (stop_time >= 0. && t_old+dt > stop_time) {
			dt = std::min(dt,stop_time-t_old);
			Print() << "Stop time limits dt = " << dt << std::endl;
		}

		dtold = dt;
		t_new = t_old + dt;

		// copy S_cc_old into S_cc_new for the pressure iterations
		for (int lev=0; lev<=finest_level; ++lev) {
			MultiFab::Copy(S_cc_new[lev],S_cc_old[lev],0,0,1,0);
		}

		// initial (pressure) iters
		for (int i=1; i<= init_iter; ++i) {
			Print() << "Doing initial pressure iteration #" << i << std::endl;
			InitIter();
		}

		if (plot_int > 0 || plot_deltat > 0) {
			Print() << "\nWriting plotfile 0 after all initialization" << std::endl;
			WritePlotFile(0,t_old,dt,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);
		} else if (small_plot_int > 0 || small_plot_deltat > 0) {
			Print() << "\nWriting small plotfile 0 after all initialization" << std::endl;
			WriteSmallPlotFile(0,t_old,dt,rho0_old,rhoh0_old,p0_old,gamma1bar_old,uold,sold,S_cc_old);
		}

		if (chk_int > 0) {
			Print() << "\nWriting checkpoint 0 after all initialization" << std::endl;
			WriteCheckPoint(0);

			int index_dummy=0;
			Print() << "\nWriting diagnosis file after all initialization" << std::endl;
			DiagFile(0,t_old,rho0_old,p0_old,uold,sold,index_dummy);
		}
	}
}

// fill in multifab and base state data
void
Maestro::InitData ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::InitData()",InitData);

	Print() << "Calling InitData()" << std::endl;

	// read in model file and fill in s0_init and p0_init for all levels
	if (use_exact_base_state) {
		init_base_state_irreg(s0_init.dataPtr(),p0_init.dataPtr(),rho0_old.dataPtr(),
		                      rhoh0_old.dataPtr(),p0_old.dataPtr(),tempbar.dataPtr(),
		                      tempbar_init.dataPtr(),
		                      r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
		std::fill(psi.begin(), psi.end(), 0.);
	} else {
		init_base_state(s0_init.dataPtr(),p0_init.dataPtr(),rho0_old.dataPtr(),
		                rhoh0_old.dataPtr(),p0_old.dataPtr(),tempbar.dataPtr(),
		                tempbar_init.dataPtr());
	}

	// calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function
	// that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
	InitFromScratch(t_old);

	if (spherical == 0) {
		// reset tagging array to include buffer zones
		TagArray();
	}

	// set finest_radial_level in fortran
	// compute numdisjointchunks, r_start_coord, r_end_coord
	init_multilevel(tag_array.dataPtr(),&finest_level);

	// average down data and fill ghost cells
	AverageDown(sold,0,Nscal);
	FillPatch(t_old,sold,sold,sold,0,0,Nscal,0,bcs_s);
	AverageDown(uold,0,AMREX_SPACEDIM);
	FillPatch(t_old,uold,uold,uold,0,0,AMREX_SPACEDIM,0,bcs_u,1);

	// free memory in s0_init and p0_init by swapping it
	// with an empty vector that will go out of scope
	RealVector s0_swap, p0_swap;
	std::swap(s0_swap,s0_init);
	std::swap(p0_swap,p0_init);

	if (fix_base_state) {
		// compute cutoff coordinates
		compute_cutoff_coords(rho0_old.dataPtr());
		make_grav_cell(grav_cell_old.dataPtr(),
		               rho0_old.dataPtr(),
		               r_cc_loc.dataPtr(),
		               r_edge_loc.dataPtr());
	}
	else {

		// first compute cutoff coordinates using initial density profile
		compute_cutoff_coords(rho0_old.dataPtr());

		if (do_smallscale) {
			// set rho0_old = rhoh0_old = 0.
			std::fill(rho0_old.begin(),  rho0_old.end(),  0.);
			std::fill(rhoh0_old.begin(), rhoh0_old.end(), 0.);
		}
		else {
			// set rho0 to be the average
			Average(sold,rho0_old,Rho);
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

			// call eos with r,p as input to recompute T,h
			TfromRhoP(sold,p0_old,1);

			// set rhoh0 to be the average
			Average(sold,rhoh0_old,RhoH);
		}

		// set tempbar to be the average
		Average(sold,tempbar,Temp);
		for (int i=0; i<tempbar.size(); ++i) {
			tempbar_init[i] = tempbar[i];
		}
	}

	// set p0^{-1} = p0_old
	for (int i=0; i<p0_old.size(); ++i) {
		p0_nm1[i] = p0_old[i];
	}
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

	if (lev > 0 && reflux_type == 2) {
		flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
	}
}


void Maestro::InitProj ()
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::InitProj()",InitProj);

	Vector<MultiFab>       rho_omegadot(finest_level+1);
	Vector<MultiFab>            thermal(finest_level+1);
	Vector<MultiFab>           rho_Hnuc(finest_level+1);
	Vector<MultiFab>           rho_Hext(finest_level+1);
	Vector<MultiFab>            rhohalf(finest_level+1);
	Vector<MultiFab>             Tcoeff(finest_level+1);
	Vector<MultiFab>             hcoeff(finest_level+1);
	Vector<MultiFab>            Xkcoeff(finest_level+1);
	Vector<MultiFab>             pcoeff(finest_level+1);
	Vector<MultiFab>       delta_gamma1(finest_level+1);
	Vector<MultiFab>  delta_gamma1_term(finest_level+1);

	RealVector Sbar( (max_radial_level+1)*nr_fine );
	RealVector delta_gamma1_termbar( (max_radial_level+1)*nr_fine );
	Sbar.shrink_to_fit();
	delta_gamma1_termbar.shrink_to_fit();

	for (int lev=0; lev<=finest_level; ++lev) {
		rho_omegadot      [lev].define(grids[lev], dmap[lev], NumSpec, 0);
		thermal           [lev].define(grids[lev], dmap[lev],       1, 0);
		rho_Hnuc          [lev].define(grids[lev], dmap[lev],       1, 0);
		rho_Hext          [lev].define(grids[lev], dmap[lev],       1, 0);
		rhohalf           [lev].define(grids[lev], dmap[lev],       1, 1);
		Tcoeff            [lev].define(grids[lev], dmap[lev],       1,    1);
		hcoeff            [lev].define(grids[lev], dmap[lev],       1,    1);
		Xkcoeff           [lev].define(grids[lev], dmap[lev], NumSpec,    1);
		pcoeff            [lev].define(grids[lev], dmap[lev],       1,    1);
		delta_gamma1      [lev].define(grids[lev], dmap[lev],       1,    1);
		delta_gamma1_term [lev].define(grids[lev], dmap[lev],       1,    1);

		// we don't have a legit timestep yet, so we set rho_omegadot,
		// rho_Hnuc, and rho_Hext to 0
		rho_omegadot[lev].setVal(0.);
		rho_Hnuc[lev].setVal(0.);
		rho_Hext[lev].setVal(0.);
		delta_gamma1[lev].setVal(0.);
		delta_gamma1_term[lev].setVal(0.);

		// initial projection does not use density weighting
		rhohalf[lev].setVal(1.);
	}

	// compute thermal diffusion
	if (use_thermal_diffusion) {
		MakeThermalCoeffs(sold,Tcoeff,hcoeff,Xkcoeff,pcoeff);

		MakeExplicitThermal(thermal,sold,Tcoeff,hcoeff,Xkcoeff,pcoeff,p0_old,
		                    temp_diffusion_formulation);
	}
	else {
		for (int lev=0; lev<=finest_level; ++lev) {
			thermal[lev].setVal(0.);
		}
	}

	// compute S at cell-centers
	Make_S_cc(S_cc_old,delta_gamma1_term,delta_gamma1,sold,uold,rho_omegadot,rho_Hnuc,
	          rho_Hext,thermal,p0_old,gamma1bar_old,delta_gamma1_termbar,psi);

	// NOTE: not sure if valid for use_exact_base_state
	if (evolve_base_state && (use_exact_base_state == 0 && average_base_state == 0)) {
		// average S into Sbar
		Average(S_cc_old,Sbar,0);
	}

	// make the nodal rhs for projection beta0*(S_cc-Sbar) + beta0*delta_chi
	MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_old,Sbar,beta0_old,delta_gamma1_term);

	// perform a nodal projection
	NodalProj(initial_projection_comp,rhcc_for_nodalproj);
}


void Maestro::DivuIter (int istep_divu_iter)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::DivuIter()",DivuIter);

	Vector<MultiFab> stemp             (finest_level+1);
	Vector<MultiFab> rho_Hext          (finest_level+1);
	Vector<MultiFab> rho_omegadot      (finest_level+1);
	Vector<MultiFab> rho_Hnuc          (finest_level+1);
	Vector<MultiFab> thermal           (finest_level+1);
	Vector<MultiFab> rhohalf           (finest_level+1);
	Vector<MultiFab> Tcoeff            (finest_level+1);
	Vector<MultiFab> hcoeff            (finest_level+1);
	Vector<MultiFab> Xkcoeff           (finest_level+1);
	Vector<MultiFab> pcoeff            (finest_level+1);
	Vector<MultiFab> delta_gamma1      (finest_level+1);
	Vector<MultiFab> delta_gamma1_term (finest_level+1);

	RealVector Sbar                  ( (max_radial_level+1)*nr_fine );
	RealVector w0_force              ( (max_radial_level+1)*nr_fine );
	RealVector p0_minus_peosbar      ( (max_radial_level+1)*nr_fine );
	RealVector delta_chi_w0          ( (max_radial_level+1)*nr_fine );
	RealVector delta_gamma1_termbar  ( (max_radial_level+1)*nr_fine );

	Sbar.shrink_to_fit();
	w0_force.shrink_to_fit();
	p0_minus_peosbar.shrink_to_fit();
	delta_chi_w0.shrink_to_fit();
	delta_gamma1_termbar.shrink_to_fit();

	std::fill(Sbar.begin(),                 Sbar.end(),                 0.);
	std::fill(etarho_ec.begin(),            etarho_ec.end(),            0.);
	std::fill(w0_force.begin(),             w0_force.end(),             0.);
	std::fill(psi.begin(),                  psi.end(),                  0.);
	std::fill(etarho_cc.begin(),            etarho_cc.end(),            0.);
	std::fill(p0_minus_peosbar.begin(),     p0_minus_peosbar.end(),     0.);
	std::fill(delta_gamma1_termbar.begin(), delta_gamma1_termbar.end(), 0.);

	for (int lev=0; lev<=finest_level; ++lev) {
		stemp             [lev].define(grids[lev], dmap[lev],   Nscal, 0);
		rho_Hext          [lev].define(grids[lev], dmap[lev],       1, 0);
		rho_omegadot      [lev].define(grids[lev], dmap[lev], NumSpec, 0);
		rho_Hnuc          [lev].define(grids[lev], dmap[lev],       1, 0);
		thermal           [lev].define(grids[lev], dmap[lev],       1, 0);
		rhohalf           [lev].define(grids[lev], dmap[lev],       1, 1);
		Tcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
		hcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
		Xkcoeff           [lev].define(grids[lev], dmap[lev], NumSpec, 1);
		pcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
		delta_gamma1      [lev].define(grids[lev], dmap[lev],       1, 1);
		delta_gamma1_term [lev].define(grids[lev], dmap[lev],       1, 1);

		// divu_iters do not use density weighting
		rhohalf[lev].setVal(1.);
	}

	React(sold,stemp,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,0.5*dt);

	if (use_thermal_diffusion) {
		MakeThermalCoeffs(sold,Tcoeff,hcoeff,Xkcoeff,pcoeff);

		MakeExplicitThermal(thermal,sold,Tcoeff,hcoeff,Xkcoeff,pcoeff,p0_old,
		                    temp_diffusion_formulation);
	}
	else {
		for (int lev=0; lev<=finest_level; ++lev) {
			thermal[lev].setVal(0.);
		}
	}

	// compute S at cell-centers
	Make_S_cc(S_cc_old,delta_gamma1_term,delta_gamma1,sold,uold,rho_omegadot,rho_Hnuc,
	          rho_Hext,thermal,p0_old,gamma1bar_old,delta_gamma1_termbar,psi);

	// NOTE: not sure if valid for use_exact_base_state
	if (evolve_base_state) {
		if ((use_exact_base_state || average_base_state) && use_delta_gamma1_term) {
			for(int i=0; i<Sbar.size(); ++i) {
				Sbar[i] += delta_gamma1_termbar[i];
			}
		} else {
			Average(S_cc_old,Sbar,0);

			// compute Sbar = Sbar + delta_gamma1_termbar
			if (use_delta_gamma1_term) {
				for(int i=0; i<Sbar.size(); ++i) {
					Sbar[i] += delta_gamma1_termbar[i];
				}
			}

			int is_predictor = 1;
			make_w0(w0.dataPtr(), w0.dataPtr(), w0_force.dataPtr(),Sbar.dataPtr(),
			        rho0_old.dataPtr(), rho0_old.dataPtr(), p0_old.dataPtr(),
			        p0_old.dataPtr(), gamma1bar_old.dataPtr(), gamma1bar_old.dataPtr(),
			        p0_minus_peosbar.dataPtr(), etarho_ec.dataPtr(),
			        etarho_cc.dataPtr(), delta_chi_w0.dataPtr(), r_cc_loc.dataPtr(),
			        r_edge_loc.dataPtr(), &dt, &dt, &is_predictor);
		}
	}

	// make the nodal rhs for projection beta0*(S_cc-Sbar) + beta0*delta_chi
	MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_old,Sbar,beta0_old,delta_gamma1_term);

	// perform a nodal projection
	NodalProj(divu_iters_comp,rhcc_for_nodalproj,istep_divu_iter);

	Real dt_hold = dt;

	// compute new time step
	EstDt();

	if (maestro_verbose > 0) {
		Print() << "Call to estdt at end of istep_divu_iter = " << istep_divu_iter
		        << " gives dt = " << dt << std::endl;
	}

	dt *= init_shrink;
	if (maestro_verbose > 0) {
		Print() << "Multiplying dt by init_shrink; dt = " << dt << std::endl;
	}

	if (dt > dt_hold) {
		if (maestro_verbose > 0) {
			Print() << "Ignoring this new dt since it's larger than the previous dt = "
			        << dt_hold << std::endl;
		}
		dt = std::min(dt_hold,dt);
	}

	if (fixed_dt != -1.0) {
		// fixed dt
		dt = fixed_dt;
		if (maestro_verbose > 0) {
			Print() << "Setting fixed dt = " << dt << std::endl;
		}
	}
}

void Maestro::InitIter ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitIter()",InitIter);

    // wallclock time
    Real start_total = ParallelDescriptor::second();

    // advance the solution by dt
    if (use_exact_base_state) {
        AdvanceTimeStepIrreg(true);
    } else if (average_base_state) {
        AdvanceTimeStepAverage(true);
    } else {
        AdvanceTimeStep(true);
    }

    // wallclock time
    Real end_total = ParallelDescriptor::second() - start_total;
    ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());

    Print() << "Time to advance time step: " << end_total << '\n';

    // copy pi from snew to sold
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Copy(sold[lev],snew[lev],Pi,Pi,1,ng_s);
    }
}
