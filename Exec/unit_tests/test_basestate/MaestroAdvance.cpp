
#include <Maestro.H>
#include <Maestro_F.H>
#include <Problem_F.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStep (bool is_initIter) {
	
	// // timer for profiling
	// BL_PROFILE_VAR("Maestro::AdvanceTimeStep()",AdvanceTimeStep);

	Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << std::endl << std::endl;

        const auto nr_fine = base_geom.nr_fine;
        const auto max_radial_level = base_geom.max_radial_level;

	// vectors store the multilevel 1D states as one very long array
	// these are cell-centered
        BaseState<Real> p0_minus_peosbar( max_radial_level+1, nr_fine );
	Vector<Real> w0_force        ( (max_radial_level+1)*nr_fine );
	Vector<Real> Sbar_old        ( (max_radial_level+1)*nr_fine );
	Vector<Real> Sbar_new        ( (max_radial_level+1)*nr_fine );
	Vector<Real> Sbar_nph        ( (max_radial_level+1)*nr_fine );
	Vector<Real> delta_chi_w0    ( (max_radial_level+1)*nr_fine );
	Vector<Real> Hext_bar        ( (max_radial_level+1)*nr_fine );
	Vector<Real> tempbar_new     ( (max_radial_level+1)*nr_fine );

	// vectors store the multilevel 1D states as one very long array
	// these are edge-centered
	Vector<Real> w0_old             ( (max_radial_level+1)*(nr_fine+1) );
	BaseState<Real> rho0_predicted_edge( max_radial_level+1, nr_fine+1 );

	// make sure C++ is as efficient as possible with memory usage
	w0_force.shrink_to_fit();
	Sbar_old.shrink_to_fit();
	Sbar_new.shrink_to_fit();
	Sbar_nph.shrink_to_fit();
	delta_chi_w0.shrink_to_fit();
	w0_old.shrink_to_fit();
	Hext_bar.shrink_to_fit();
	tempbar_new.shrink_to_fit();

	int is_predictor;

	p0_minus_peosbar.setVal(0.);
	std::fill(delta_chi_w0.begin(), delta_chi_w0.end(), 0.);

	// make Fortran-friendly RealVectors
	RealVector p0_old_vec( (max_radial_level+1)*nr_fine );
	RealVector gamma1bar_old_vec( (max_radial_level+1)*nr_fine );
	p0_old.toVector(p0_old_vec);
	gamma1bar_old.toVector(gamma1bar_old_vec);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute initial gamma1bar_old
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_gamma1bar(gamma1bar_old_vec.dataPtr(), rho0_old.dataPtr(), tempbar.dataPtr(),
	                  rhoX0_old.dataPtr(), p0_old_vec.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute the heating term and Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	get_heating(Hext_bar.dataPtr(), rho0_old.dataPtr(), tempbar.dataPtr(),
	            rhoX0_old.dataPtr(), t_old, dt, base_geom.r_cc_loc.dataPtr());


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
	Makew0(w0_old, w0_force, Sbar_old, rho0_old, rho0_old,
	       p0_old, p0_old, gamma1bar_old, gamma1bar_old,
	       p0_minus_peosbar, delta_chi_w0, dt, dtold, is_predictor);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update density and compute rho0_predicted_edge
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	AdvectBaseDens(rho0_predicted_edge);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! recompute cutoff coordinates now that rho0 has changed
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_new.dataPtr());
	ComputeCutoffCoords(rho0_new);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gravity
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	MakeGravCell(grav_cell_new, rho0_new);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update species
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// make a Fortran-friendly RealVector of rho0_predicted_edge
	RealVector rho0_predicted_edge_vec( (max_radial_level+1)*(nr_fine+1) );
	rho0_predicted_edge.toVector(rho0_predicted_edge_vec);
	
	update_species(rho0_old.dataPtr(), rho0_predicted_edge_vec.dataPtr(),
	               rhoX0_old.dataPtr(), rhoX0_new.dataPtr(), w0.dataPtr(),
	               base_geom.r_edge_loc.dataPtr(), base_geom.r_cc_loc.dataPtr(), dt);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update pressure
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// set new p0 through HSE
	p0_new.copy(p0_old);

	EnforceHSE(rho0_new, p0_new, grav_cell_new);

	// make Fortran-friendly RealVectors
	RealVector p0_new_vec( (max_radial_level+1)*nr_fine );
	RealVector gamma1bar_new_vec( (max_radial_level+1)*nr_fine );
	p0_new.toVector(p0_new_vec);
	gamma1bar_new.toVector(gamma1bar_new_vec);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gamma1bar_new
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_gamma1bar(gamma1bar_new_vec.dataPtr(), rho0_new.dataPtr(), tempbar.dataPtr(),
	                  rhoX0_new.dataPtr(), p0_new_vec.dataPtr());
	gamma1bar_new = BaseState<Real>(gamma1bar_new_vec, max_radial_level+1, nr_fine, 1);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update temperature
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	std::copy(tempbar.begin(), tempbar.end(), tempbar_new.begin());

	update_temp(rho0_new.dataPtr(), tempbar_new.dataPtr(), rhoX0_new.dataPtr(),
	            p0_new_vec.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! reset cutoff coordinates
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_old.dataPtr());
	ComputeCutoffCoords(rho0_old);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! make Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	make_Sbar(Sbar_new.dataPtr(), rho0_new.dataPtr(), tempbar_new.dataPtr(),
	          rhoX0_new.dataPtr(), Hext_bar.dataPtr());

	for (int i=0; i<Sbar_nph.size(); ++i)
		Sbar_nph[i] = 0.5*(Sbar_old[i] + Sbar_new[i]);

	p0_minus_peosbar.setVal(0.);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute w_0
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	w0_old = w0;

	is_predictor = 0;
	Makew0(w0_old, w0_force, Sbar_nph, rho0_old, rho0_new,
	       p0_old, p0_new, gamma1bar_old, gamma1bar_new,
	       p0_minus_peosbar, delta_chi_w0, dt, dtold, is_predictor);
	

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update density and compute rho0_predicted_edge
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	AdvectBaseDens(rho0_predicted_edge);
	rho0_predicted_edge.toVector(rho0_predicted_edge_vec);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! recompute cutoff coordinates now that rho0 has changed
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_new.dataPtr());
	ComputeCutoffCoords(rho0_new);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gravity
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	MakeGravCell(grav_cell_new, rho0_new);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update species
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	update_species(rho0_old.dataPtr(), rho0_predicted_edge_vec.dataPtr(),
	               rhoX0_old.dataPtr(), rhoX0_new.dataPtr(), w0.dataPtr(),
	               base_geom.r_edge_loc.dataPtr(), base_geom.r_cc_loc.dataPtr(), dt);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update pressure
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	p0_new.copy(p0_old);

	EnforceHSE(rho0_new, p0_new, grav_cell_new);
	p0_new.toVector(p0_new_vec);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute gamma1bar_new
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_gamma1bar(gamma1bar_new_vec.dataPtr(), rho0_new.dataPtr(), tempbar.dataPtr(),
	                  rhoX0_new.dataPtr(), p0_new_vec.dataPtr());
	gamma1bar_new = BaseState<Real>(gamma1bar_new_vec, max_radial_level+1, nr_fine, 1);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update temperature
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	std::copy(tempbar.begin(), tempbar.end(), tempbar_new.begin());

	update_temp(rho0_new.dataPtr(), tempbar_new.dataPtr(), rhoX0_new.dataPtr(),
	            p0_new.dataPtr());

	std::swap(rhoX0_old, rhoX0_new);
	std::swap(tempbar, tempbar_new);

}
