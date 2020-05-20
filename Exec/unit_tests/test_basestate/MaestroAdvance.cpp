
#include <Maestro.H>
#include <Maestro_F.H>
#include <Problem_F.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStep (bool is_initIter) {
	
	// // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStep()",AdvanceTimeStep);

	Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << std::endl << std::endl;

    const auto nr_fine = base_geom.nr_fine;
    const auto max_radial_level = base_geom.max_radial_level;

	// vectors store the multilevel 1D states as one very long array
	// these are cell-centered
	BaseState<Real> p0_minus_peosbar(max_radial_level+1, nr_fine);
	BaseState<Real> w0_force(max_radial_level+1, nr_fine);
	BaseState<Real> Sbar_old(max_radial_level+1, nr_fine);
	BaseState<Real> Sbar_new(max_radial_level+1, nr_fine);
	BaseState<Real> Sbar_nph(max_radial_level+1, nr_fine);
	BaseState<Real> delta_chi_w0(max_radial_level+1, nr_fine);
	BaseState<Real> Hext_bar(max_radial_level+1, nr_fine);
	BaseState<Real> tempbar_new(max_radial_level+1, nr_fine);
    BaseState<Real> rhoX0_old(max_radial_level+1, nr_fine, NumSpec);
    BaseState<Real> rhoX0_new(max_radial_level+1, nr_fine, NumSpec);

	// vectors store the multilevel 1D states as one very long array
	// these are edge-centered
	BaseState<Real> w0_old(max_radial_level+1, nr_fine+1);
	BaseState<Real> rho0_predicted_edge(max_radial_level+1, nr_fine+1);

	int is_predictor;

	p0_minus_peosbar.setVal(0.);
	delta_chi_w0.setVal(0.);

	// make Fortran-friendly RealVectors
	RealVector p0_old_vec( (max_radial_level+1)*nr_fine );
	RealVector gamma1bar_old_vec( (max_radial_level+1)*nr_fine );
	p0_old.toVector(p0_old_vec);
	gamma1bar_old.toVector(gamma1bar_old_vec);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute initial gamma1bar_old
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    auto rho0_old_arr = rho0_old.array();
	auto p0_old_arr = p0_old.array();
	auto rhoX0_old_arr = rhoX0_old.array();
	auto tempbar_arr = tempbar.array();
	auto gamma1bar_old_arr = gamma1bar_old.array();

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_old_arr(l,r);
            eos_state.p = p0_old_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_old_arr(l,r,n) / rho0_old_arr(l,r);
            } 
            eos_state.T = tempbar_arr(l,r);

            eos(eos_input_rp, eos_state);

            gamma1bar_old_arr(l,r) = eos_state.gam1;
        }
	}

    // copy rhoX0_old 
	auto s0_init_arr = s0_init.array();

    for (auto n = 0; n < base_geom.max_radial_level; ++n) {
        for (auto r = 0; r < base_geom.nr(n); ++r) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                rhoX0_old_arr(n,r,comp) = s0_init_arr(n,r,FirstSpec+comp);
            }
        }
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute the heating term and Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	get_heating(Hext_bar.dataPtr(), rho0_old.dataPtr(), tempbar.dataPtr(),
	            rhoX0_old.dataPtr(), t_old, dt, base_geom.r_cc_loc.dataPtr());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! make Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	auto Hext_bar_arr = Hext_bar.array();
	auto Sbar_old_arr = Sbar_old.array();

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_old_arr(l,r);
            eos_state.T = tempbar_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_old_arr(l,r,n) / rho0_old_arr(l,r);
            } 

            eos(eos_input_rt, eos_state);

            Sbar_old_arr(l,r) = Hext_bar_arr(l,r) * eos_state.dpdT / (eos_state.rho * eos_state.cp * eos_state.dpdr);
        }
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute w_0
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	w0_old.copy(w0);

	// compute w0, w0_force, and delta_chi_w0
	is_predictor = 1;
	Makew0(w0_old, w0_force, Sbar_old, rho0_old, rho0_old,
	       p0_old, p0_old, gamma1bar_old, gamma1bar_old,
	       p0_minus_peosbar, dt, dtold, is_predictor);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update density and compute rho0_predicted_edge
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	AdvectBaseDens(rho0_predicted_edge);
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! recompute cutoff coordinates now that rho0 has changed
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	auto rho0_new_arr = rho0_new.array();
	auto p0_new_arr = p0_new.array();
	auto rhoX0_new_arr = rhoX0_new.array();
	auto gamma1bar_new_arr = gamma1bar_new.array();
	auto tempbar_new_arr = tempbar_new.array();

	compute_cutoff_coords(rho0_new.dataPtr());
	ComputeCutoffCoords(rho0_new);
	base_geom.ComputeCutoffCoords(rho0_new.array());

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

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l,r);
            eos_state.p = p0_new_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l,r,n) / rho0_new_arr(l,r);
            } 
            eos_state.T = tempbar_arr(l,r);

            eos(eos_input_rp, eos_state);

            gamma1bar_new_arr(l,r) = eos_state.gam1;
        }
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update temperature
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l,r);
            eos_state.p = p0_new_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l,r,n) / rho0_new_arr(l,r);
            } 
            eos_state.T = tempbar_arr(l,r);

            eos(eos_input_rp, eos_state);

            tempbar_new_arr(l,r) = eos_state.T;
        }
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! reset cutoff coordinates
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	compute_cutoff_coords(rho0_old.dataPtr());
	ComputeCutoffCoords(rho0_old);
	base_geom.ComputeCutoffCoords(rho0_old.array());

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! make Sbar
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	auto Sbar_new_arr = Sbar_new.array();

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l,r);
            eos_state.T = tempbar_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l,r,n) / rho0_new_arr(l,r);
            } 
   
            eos(eos_input_rt, eos_state);

            Sbar_new_arr(l,r) = Hext_bar_arr(l,r) * eos_state.dpdT / (eos_state.rho * eos_state.cp * eos_state.dpdr);
        }
    }

	Sbar_nph.copy(0.5*(Sbar_old + Sbar_new));

	p0_minus_peosbar.setVal(0.);


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! compute w_0
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	w0_old.copy(w0);

	is_predictor = 0;
	Makew0(w0_old, w0_force, Sbar_nph, rho0_old, rho0_new,
	       p0_old, p0_new, gamma1bar_old, gamma1bar_new,
	       p0_minus_peosbar, dt, dtold, is_predictor);
	

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
	base_geom.ComputeCutoffCoords(rho0_new.array());

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

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l,r);
            eos_state.p = p0_new_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l,r,n) / rho0_new_arr(l,r);
            } 
            eos_state.T = tempbar_arr(l,r);

            eos(eos_input_rp, eos_state);

            gamma1bar_new_arr(l,r) = eos_state.gam1;
        }
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ! update temperature
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr(l); ++r) {
            eos_t eos_state;

            eos_state.rho = rho0_new_arr(l,r);
            eos_state.p = p0_new_arr(l,r);
            for (auto n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = rhoX0_new_arr(l,r,n) / rho0_new_arr(l,r);
            } 
            eos_state.T = tempbar_arr(l,r);

            eos(eos_input_rp, eos_state);

            tempbar_new_arr(l,r) = eos_state.T;
        }
	}

	// rhoX0_old.swap(rhoX0_new);
	tempbar.swap(tempbar_new);

    // copy rhoX0_new into s0_init_arr
    for (auto n = 0; n < base_geom.max_radial_level; ++n) {
        for (auto r = 0; r < base_geom.nr(n); ++r) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init_arr(n,r,FirstSpec+comp) = rhoX0_new_arr(n,r,comp);
            }
        }
	}
}
