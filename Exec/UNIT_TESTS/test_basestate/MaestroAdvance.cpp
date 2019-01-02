
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStep (bool is_initIter) {

	// timer for profiling
	BL_PROFILE_VAR("Maestro::AdvanceTimeStep()",AdvanceTimeStep);

	// cell-centered MultiFabs needed within the AdvanceTimeStep routine
	Vector<MultiFab>           s1(finest_level+1);
	Vector<MultiFab>           s2(finest_level+1);
	Vector<MultiFab>    delta_chi(finest_level+1);

	// vectors store the multilevel 1D states as one very long array
	// these are cell-centered
	Vector<Real> p0_nph          ( (max_radial_level+1)*nr_fine );
	Vector<Real> p0_minus_peosbar( (max_radial_level+1)*nr_fine );
	Vector<Real> peosbar         ( (max_radial_level+1)*nr_fine );
	Vector<Real> w0_force        ( (max_radial_level+1)*nr_fine );
	Vector<Real> Sbar_old        ( (max_radial_level+1)*nr_fine );
	Vector<Real> Sbar_new        ( (max_radial_level+1)*nr_fine );
	Vector<Real> Sbar_nph        ( (max_radial_level+1)*nr_fine );
	Vector<Real> delta_chi_w0    ( (max_radial_level+1)*nr_fine );
	Vector<Real> Hext_bar        ( (max_radial_level+1)*nr_fine );
	Vector<Real> rhoX0_old       ( (max_radial_level+1)*nr_fine );
	Vector<Real> rhoX0_new       ( (max_radial_level+1)*nr_fine );
	Vector<Real> tempbar_new     ( (max_radial_level+1)*nr_fine );

	// vectors store the multilevel 1D states as one very long array
	// these are edge-centered
	Vector<Real> w0_old             ( (max_radial_level+1)*(nr_fine+1) );
	Vector<Real> rho0_predicted_edge( (max_radial_level+1)*(nr_fine+1) );

	// make sure C++ is as efficient as possible with memory usage
	p0_nph.shrink_to_fit();
	p0_minus_peosbar.shrink_to_fit();
	peosbar.shrink_to_fit();
	w0_force.shrink_to_fit();
	Sbar_old.shrink_to_fit();
	Sbar_new.shrink_to_fit();
	Sbar_nph.shrink_to_fit();
	delta_chi_w0.shrink_to_fit();
	w0_old.shrink_to_fit();
	rho0_predicted_edge.shrink_to_fit();
	Hext_bar.shrink_to_fit();
	rhoX0_old.shrink_to_fit();
	rhoX0_new.shrink_to_fit();

	int is_predictor;

	// wallclock time
	const Real strt_total = ParallelDescriptor::second();


	for (int lev=0; lev<=finest_level; ++lev) {
		// cell-centered MultiFabs
		s1          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
		s1[lev].setVal(0.);
		s2          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
		delta_chi   [lev].define(grids[lev], dmap[lev],       1,    0);
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute initial gamma1bar_old
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_gamma1bar(gamma1bar_old.dataPtr(), rho0_old.dataPtr(), tempbar.dataPtr(),
	                  rhoX0_old.dataPtr(), p0_old.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute the heating term and Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	get_heating(Hext_bar.dataPtr(), rho0_old.dataPtr(), tempbar.dataPtr(),
	            rhoX0_old.dataPtr(), t_old, dt, r_cc_loc.dataPtr());


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! make Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	make_Sbar(Sbar_old.dataPtr(), rho0_old.dataPtr(), tempbar.dataPtr(),
	          rhoX0_old.dataPtr(), Hext_bar.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute w_0
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	w0_old = w0;

	// compute w0, w0_force, and delta_chi_w0
	is_predictor = 1;
	make_w0(w0.dataPtr(),w0_old.dataPtr(),w0_force.dataPtr(),Sbar_old.dataPtr(),
	        rho0_old.dataPtr(),rho0_old.dataPtr(),p0_old.dataPtr(),p0_old.dataPtr(),
	        gamma1bar_old.dataPtr(),gamma1bar_old.dataPtr(),p0_minus_peosbar.dataPtr(),
	        psi.dataPtr(),etarho_ec.dataPtr(),etarho_cc.dataPtr(),delta_chi_w0.dataPtr(),
	        r_cc_loc.dataPtr(),r_edge_loc.dataPtr(),&dt,&dtold,&is_predictor);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update density and compute rho0_predicted_edge
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	advect_base_dens(w0.dataPtr(), rho0_old.dataPtr(), rho0_new.dataPtr(),
	                 rho0_predicted_edge.dataPtr(), &dt,
	                 r_cc_loc.dataPtr(), r_edge_loc.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! recompute cutoff coordinates now that rho0 has changed
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_new.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gravity
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	make_grav_cell(grav_cell_new.dataPtr(),
	               rho0_new.dataPtr(),
	               r_cc_loc.dataPtr(),
	               r_edge_loc.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update species
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	update_species(rho0_old.dataPtr(), rho0_predicted_edge.dataPtr(),
	               rhoX0_old.dataPtr(), rhoX0_new.dataPtr(), w0.dataPtr(),
                   r_edge_loc.dataPtr(), r_cc_loc.dataPtr(), dt);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update pressure
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// set new p0 through HSE
	p0_new = p0_old;

	enforce_HSE(rho0_new.dataPtr(),
	            p0_new.dataPtr(),
	            grav_cell_new.dataPtr(),
	            r_cc_loc.dataPtr(),
	            r_edge_loc.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gamma1bar_new
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_gamma1bar(gamma1bar_new.dataPtr(), rho0_new.dataPtr(), tempbar.dataPtr(),
	                  rhoX0_new.dataPtr(), p0_new.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update temperature
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::copy(tempbar.begin(), tempbar.end(), tempbar_new.begin());

	update_temp(rho0_new.dataPtr(), tempbar_new.dataPtr(), rhoX0_new.dataPtr(),
	            p0_new.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! reset cutoff coordinates
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_old.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! make Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	make_Sbar(Sbar_new.dataPtr(), rho0_new.dataPtr(), tempbar_new.dataPtr(),
	          rhoX0_new.dataPtr(), Hext_bar.dataPtr());

	for (int i=0; i<Sbar_nph.size(); ++i)
		Sbar_nph[i] = 0.5*(Sbar_old[i] + Sbar_new[i]);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute w_0
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w0_old = w0;

	is_predictor = 0;
	make_w0(w0.dataPtr(),w0_old.dataPtr(),w0_force.dataPtr(),Sbar_nph.dataPtr(),
	        rho0_old.dataPtr(),rho0_new.dataPtr(),p0_old.dataPtr(),p0_new.dataPtr(),
	        gamma1bar_old.dataPtr(),gamma1bar_new.dataPtr(),p0_minus_peosbar.dataPtr(),
	        psi.dataPtr(),etarho_ec.dataPtr(),etarho_cc.dataPtr(),delta_chi_w0.dataPtr(),
	        r_cc_loc.dataPtr(),r_edge_loc.dataPtr(),&dt,&dtold,&is_predictor);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update density and compute rho0_predicted_edge
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	advect_base_dens(w0.dataPtr(), rho0_old.dataPtr(), rho0_new.dataPtr(),
	                 rho0_predicted_edge.dataPtr(), &dt,
	                 r_cc_loc.dataPtr(), r_edge_loc.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! recompute cutoff coordinates now that rho0 has changed
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_new.dataPtr());


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gravity
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	make_grav_cell(grav_cell_new.dataPtr(),
	               rho0_new.dataPtr(),
	               r_cc_loc.dataPtr(),
	               r_edge_loc.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update species
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    update_species(rho0_old.dataPtr(), rho0_predicted_edge.dataPtr(),
	               rhoX0_old.dataPtr(), rhoX0_new.dataPtr(), w0.dataPtr(),
                   r_edge_loc.dataPtr(), r_cc_loc.dataPtr(), dt);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update pressure
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	p0_new = p0_old;

	enforce_HSE(rho0_new.dataPtr(),
	            p0_new.dataPtr(),
	            grav_cell_new.dataPtr(),
	            r_cc_loc.dataPtr(),
	            r_edge_loc.dataPtr());

	for (int i=0; i<p0_nph.size(); ++i)
		p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gamma1bar_new
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_gamma1bar(gamma1bar_new.dataPtr(), rho0_new.dataPtr(), tempbar.dataPtr(),
	                  rhoX0_new.dataPtr(), p0_new.dataPtr());


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update temperature
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::copy(tempbar.begin(), tempbar.end(), tempbar_new.begin());

    update_temp(rho0_new.dataPtr(), tempbar.dataPtr(), rhoX0_new.dataPtr(),
	            p0_new.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute the new timestep
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
	        << " DT = " << dt << std::endl;

	// wallclock time
	Real end_total = ParallelDescriptor::second() - strt_total;

	// print wallclock time
	ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
	if (maestro_verbose > 0)
		Print() << "Time to advance time step: " << end_total << '\n';

    std::swap(rhoX0_old, rhoX0_new);
    std::swap(tempbar, tempbar_new);

}
