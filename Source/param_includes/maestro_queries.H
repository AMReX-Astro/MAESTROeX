
// This file is automatically created by parse_maestro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

maestro::reflux_type = 1;
pp.query("reflux_type", maestro::reflux_type);

maestro::maestro_verbose = 1;
pp.query("maestro_verbose", maestro::maestro_verbose);

maestro::model_file = "";
pp.query("model_file", maestro::model_file);

maestro::perturb_model = false;
pp.query("perturb_model", maestro::perturb_model);

maestro::print_init_hse_diag = false;
pp.query("print_init_hse_diag", maestro::print_init_hse_diag);

maestro::stop_time = -1.0;
pp.query("stop_time", maestro::stop_time);

maestro::max_step = 1;
pp.query("max_step", maestro::max_step);

maestro::cfl = 0.5;
pp.query("cfl", maestro::cfl);

maestro::init_shrink = 1.0;
pp.query("init_shrink", maestro::init_shrink);

maestro::small_dt = 1.e-10;
pp.query("small_dt", maestro::small_dt);

maestro::max_dt_growth = 1.1;
pp.query("max_dt_growth", maestro::max_dt_growth);

maestro::max_dt = 1.e33;
pp.query("max_dt", maestro::max_dt);

maestro::fixed_dt = -1.0;
pp.query("fixed_dt", maestro::fixed_dt);

maestro::nuclear_dt_fac = -1.0;
pp.query("nuclear_dt_fac", maestro::nuclear_dt_fac);

maestro::use_soundspeed_firstdt = false;
pp.query("use_soundspeed_firstdt", maestro::use_soundspeed_firstdt);

maestro::use_divu_firstdt = false;
pp.query("use_divu_firstdt", maestro::use_divu_firstdt);

maestro::spherical = 0;
pp.query("spherical", maestro::spherical);

maestro::octant = false;
pp.query("octant", maestro::octant);

maestro::do_2d_planar_octant = 0;
pp.query("do_2d_planar_octant", maestro::do_2d_planar_octant);

maestro::regrid_int = -1;
pp.query("regrid_int", maestro::regrid_int);

maestro::amr_buf_width = -1;
pp.query("amr_buf_width", maestro::amr_buf_width);

maestro::drdxfac = 1;
pp.query("drdxfac", maestro::drdxfac);

maestro::minwidth = 8;
pp.query("minwidth", maestro::minwidth);

maestro::min_eff = 0.9;
pp.query("min_eff", maestro::min_eff);

maestro::use_tpert_in_tagging = false;
pp.query("use_tpert_in_tagging", maestro::use_tpert_in_tagging);

maestro::plot_int = 0;
pp.query("plot_int", maestro::plot_int);

maestro::small_plot_int = 0;
pp.query("small_plot_int", maestro::small_plot_int);

maestro::plot_deltat = -1.0;
pp.query("plot_deltat", maestro::plot_deltat);

maestro::small_plot_deltat = -1.0;
pp.query("small_plot_deltat", maestro::small_plot_deltat);

maestro::chk_int = 0;
pp.query("chk_int", maestro::chk_int);

maestro::chk_deltat = -1.0;
pp.query("chk_deltat", maestro::chk_deltat);

maestro::plot_h_with_use_tfromp = true;
pp.query("plot_h_with_use_tfromp", maestro::plot_h_with_use_tfromp);

maestro::plot_spec = true;
pp.query("plot_spec", maestro::plot_spec);

maestro::plot_omegadot = true;
pp.query("plot_omegadot", maestro::plot_omegadot);

maestro::plot_Hext = false;
pp.query("plot_Hext", maestro::plot_Hext);

maestro::plot_Hnuc = true;
pp.query("plot_Hnuc", maestro::plot_Hnuc);

maestro::plot_eta = false;
pp.query("plot_eta", maestro::plot_eta);

maestro::plot_trac = false;
pp.query("plot_trac", maestro::plot_trac);

maestro::plot_base_state = true;
pp.query("plot_base_state", maestro::plot_base_state);

maestro::plot_gpi = true;
pp.query("plot_gpi", maestro::plot_gpi);

maestro::plot_cs = false;
pp.query("plot_cs", maestro::plot_cs);

maestro::plot_grav = false;
pp.query("plot_grav", maestro::plot_grav);

maestro::plot_base_name = "plt";
pp.query("plot_base_name", maestro::plot_base_name);

maestro::small_plot_base_name = "smallplt";
pp.query("small_plot_base_name", maestro::small_plot_base_name);

maestro::check_base_name = "chk";
pp.query("check_base_name", maestro::check_base_name);

maestro::diag_buf_size = 10;
pp.query("diag_buf_size", maestro::diag_buf_size);

maestro::plot_ad_excess = false;
pp.query("plot_ad_excess", maestro::plot_ad_excess);

maestro::plot_processors = false;
pp.query("plot_processors", maestro::plot_processors);

maestro::plot_pidivu = false;
pp.query("plot_pidivu", maestro::plot_pidivu);

maestro::small_plot_vars = "rho p0 magvel";
pp.query("small_plot_vars", maestro::small_plot_vars);

maestro::init_iter = 4;
pp.query("init_iter", maestro::init_iter);

maestro::init_divu_iter = 4;
pp.query("init_divu_iter", maestro::init_divu_iter);

maestro::restart_file = "";
pp.query("restart_file", maestro::restart_file);

maestro::restart_into_finer = false;
pp.query("restart_into_finer", maestro::restart_into_finer);

maestro::do_initial_projection = true;
pp.query("do_initial_projection", maestro::do_initial_projection);

maestro::mg_verbose = 1;
pp.query("mg_verbose", maestro::mg_verbose);

maestro::cg_verbose = 0;
pp.query("cg_verbose", maestro::cg_verbose);

maestro::mg_cycle_type = 3;
pp.query("mg_cycle_type", maestro::mg_cycle_type);

maestro::hg_cycle_type = 3;
pp.query("hg_cycle_type", maestro::hg_cycle_type);

maestro::hg_bottom_solver = -1;
pp.query("hg_bottom_solver", maestro::hg_bottom_solver);

maestro::mg_bottom_solver = -1;
pp.query("mg_bottom_solver", maestro::mg_bottom_solver);

maestro::max_mg_bottom_nlevels = 1000;
pp.query("max_mg_bottom_nlevels", maestro::max_mg_bottom_nlevels);

maestro::mg_bottom_nu = 10;
pp.query("mg_bottom_nu", maestro::mg_bottom_nu);

maestro::mg_nu_1 = 2;
pp.query("mg_nu_1", maestro::mg_nu_1);

maestro::mg_nu_2 = 2;
pp.query("mg_nu_2", maestro::mg_nu_2);

maestro::hg_dense_stencil = true;
pp.query("hg_dense_stencil", maestro::hg_dense_stencil);

maestro::do_sponge = false;
pp.query("do_sponge", maestro::do_sponge);

maestro::sponge_kappa = 10.e0;
pp.query("sponge_kappa", maestro::sponge_kappa);

maestro::sponge_center_density = 3.e6;
pp.query("sponge_center_density", maestro::sponge_center_density);

maestro::sponge_start_factor = 3.333e0;
pp.query("sponge_start_factor", maestro::sponge_start_factor);

maestro::plot_sponge_fdamp = false;
pp.query("plot_sponge_fdamp", maestro::plot_sponge_fdamp);

maestro::anelastic_cutoff_density = -1.0;
pp.query("anelastic_cutoff_density", maestro::anelastic_cutoff_density);

maestro::base_cutoff_density = -1.0;
pp.query("base_cutoff_density", maestro::base_cutoff_density);

maestro::burning_cutoff_density_lo = -1.0;
pp.query("burning_cutoff_density_lo", maestro::burning_cutoff_density_lo);

maestro::burning_cutoff_density_hi = 1.e100;
pp.query("burning_cutoff_density_hi", maestro::burning_cutoff_density_hi);

maestro::heating_cutoff_density_lo = -1.0;
pp.query("heating_cutoff_density_lo", maestro::heating_cutoff_density_lo);

maestro::heating_cutoff_density_hi = 1.e100;
pp.query("heating_cutoff_density_hi", maestro::heating_cutoff_density_hi);

maestro::buoyancy_cutoff_factor = 5.0;
pp.query("buoyancy_cutoff_factor", maestro::buoyancy_cutoff_factor);

maestro::dpdt_factor = 0.0;
pp.query("dpdt_factor", maestro::dpdt_factor);

maestro::do_planar_invsq_grav = false;
pp.query("do_planar_invsq_grav", maestro::do_planar_invsq_grav);

maestro::planar_invsq_mass = 0.0;
pp.query("planar_invsq_mass", maestro::planar_invsq_mass);

maestro::evolve_base_state = true;
pp.query("evolve_base_state", maestro::evolve_base_state);

maestro::use_exact_base_state = false;
pp.query("use_exact_base_state", maestro::use_exact_base_state);

maestro::fix_base_state = false;
pp.query("fix_base_state", maestro::fix_base_state);

maestro::average_base_state = false;
pp.query("average_base_state", maestro::average_base_state);

maestro::do_smallscale = false;
pp.query("do_smallscale", maestro::do_smallscale);

maestro::do_eos_h_above_cutoff = true;
pp.query("do_eos_h_above_cutoff", maestro::do_eos_h_above_cutoff);

maestro::enthalpy_pred_type = 1;
pp.query("enthalpy_pred_type", maestro::enthalpy_pred_type);

maestro::species_pred_type = 1;
pp.query("species_pred_type", maestro::species_pred_type);

maestro::use_delta_gamma1_term = true;
pp.query("use_delta_gamma1_term", maestro::use_delta_gamma1_term);

maestro::use_etarho = true;
pp.query("use_etarho", maestro::use_etarho);

maestro::slope_order = 4;
pp.query("slope_order", maestro::slope_order);

maestro::grav_const = -1.5e10;
pp.query("grav_const", maestro::grav_const);

maestro::ppm_type = 1;
pp.query("ppm_type", maestro::ppm_type);

maestro::bds_type = 0;
pp.query("bds_type", maestro::bds_type);

maestro::ppm_trace_forces = 0;
pp.query("ppm_trace_forces", maestro::ppm_trace_forces);

maestro::beta0_type = 1;
pp.query("beta0_type", maestro::beta0_type);

maestro::use_linear_grav_in_beta0 = false;
pp.query("use_linear_grav_in_beta0", maestro::use_linear_grav_in_beta0);

maestro::rotational_frequency = 0.0;
pp.query("rotational_frequency", maestro::rotational_frequency);

maestro::co_latitude = 0.0;
pp.query("co_latitude", maestro::co_latitude);

maestro::rotation_radius = 1.0e6;
pp.query("rotation_radius", maestro::rotation_radius);

maestro::mach_max_abort = -1.0;
pp.query("mach_max_abort", maestro::mach_max_abort);

maestro::drive_initial_convection = false;
pp.query("drive_initial_convection", maestro::drive_initial_convection);

maestro::stop_initial_convection = -1;
pp.query("stop_initial_convection", maestro::stop_initial_convection);

maestro::restart_with_vel_field = false;
pp.query("restart_with_vel_field", maestro::restart_with_vel_field);

maestro::use_alt_energy_fix = true;
pp.query("use_alt_energy_fix", maestro::use_alt_energy_fix);

maestro::use_omegadot_terms_in_S = true;
pp.query("use_omegadot_terms_in_S", maestro::use_omegadot_terms_in_S);

maestro::use_thermal_diffusion = false;
pp.query("use_thermal_diffusion", maestro::use_thermal_diffusion);

maestro::temp_diffusion_formulation = 2;
pp.query("temp_diffusion_formulation", maestro::temp_diffusion_formulation);

maestro::thermal_diffusion_type = 1;
pp.query("thermal_diffusion_type", maestro::thermal_diffusion_type);

maestro::limit_conductivity = false;
pp.query("limit_conductivity", maestro::limit_conductivity);

maestro::do_burning = true;
pp.query("do_burning", maestro::do_burning);

maestro::burner_threshold_species = "";
pp.query("burner_threshold_species", maestro::burner_threshold_species);

maestro::burner_threshold_cutoff = 1.e-10;
pp.query("burner_threshold_cutoff", maestro::burner_threshold_cutoff);

maestro::do_subgrid_burning = false;
pp.query("do_subgrid_burning", maestro::do_subgrid_burning);

maestro::reaction_sum_tol = 1.e-10;
pp.query("reaction_sum_tol", maestro::reaction_sum_tol);

maestro::small_temp = 5.e6;
pp.query("small_temp", maestro::small_temp);

maestro::small_dens = 1.e-5;
pp.query("small_dens", maestro::small_dens);

maestro::use_tfromp = false;
pp.query("use_tfromp", maestro::use_tfromp);

maestro::use_eos_e_instead_of_h = false;
pp.query("use_eos_e_instead_of_h", maestro::use_eos_e_instead_of_h);

maestro::use_pprime_in_tfromp = false;
pp.query("use_pprime_in_tfromp", maestro::use_pprime_in_tfromp);

maestro::s0_interp_type = 3;
pp.query("s0_interp_type", maestro::s0_interp_type);

maestro::w0_interp_type = 2;
pp.query("w0_interp_type", maestro::w0_interp_type);

maestro::s0mac_interp_type = 1;
pp.query("s0mac_interp_type", maestro::s0mac_interp_type);

maestro::w0mac_interp_type = 1;
pp.query("w0mac_interp_type", maestro::w0mac_interp_type);

#ifdef AMREX_DEBUG
maestro::print_fortran_warnings = 1;
#else
maestro::print_fortran_warnings = 0;
#endif
pp.query("print_fortran_warnings", maestro::print_fortran_warnings);

#ifdef AMREX_DEBUG
maestro::print_update_diagnostics = 1;
#else
maestro::print_update_diagnostics = 0;
#endif
pp.query("print_update_diagnostics", maestro::print_update_diagnostics);

maestro::track_grid_losses = 0;
pp.query("track_grid_losses", maestro::track_grid_losses);

maestro::sum_interval = -1;
pp.query("sum_interval", maestro::sum_interval);

maestro::sum_per = -1.0e0;
pp.query("sum_per", maestro::sum_per);

maestro::show_center_of_mass = 0;
pp.query("show_center_of_mass", maestro::show_center_of_mass);

maestro::hard_cfl_limit = 1;
pp.query("hard_cfl_limit", maestro::hard_cfl_limit);

maestro::job_name = "";
pp.query("job_name", maestro::job_name);

maestro::output_at_completion = 1;
pp.query("output_at_completion", maestro::output_at_completion);

maestro::reset_checkpoint_time = -1.e200;
pp.query("reset_checkpoint_time", maestro::reset_checkpoint_time);

maestro::reset_checkpoint_step = -1;
pp.query("reset_checkpoint_step", maestro::reset_checkpoint_step);

maestro::use_particles = false;
pp.query("use_particles", maestro::use_particles);

maestro::store_particle_vels = false;
pp.query("store_particle_vels", maestro::store_particle_vels);

maestro::do_heating = false;
pp.query("do_heating", maestro::do_heating);

maestro::sdc_iters = 1;
pp.query("sdc_iters", maestro::sdc_iters);

maestro::sdc_couple_mac_velocity = false;
pp.query("sdc_couple_mac_velocity", maestro::sdc_couple_mac_velocity);

maestro::deterministic_nodal_solve = false;
pp.query("deterministic_nodal_solve", maestro::deterministic_nodal_solve);

maestro::eps_init_proj_cart = 1.e-12;
pp.query("eps_init_proj_cart", maestro::eps_init_proj_cart);

maestro::eps_init_proj_sph = 1.e-10;
pp.query("eps_init_proj_sph", maestro::eps_init_proj_sph);

maestro::eps_divu_cart = 1.e-12;
pp.query("eps_divu_cart", maestro::eps_divu_cart);

maestro::eps_divu_sph = 1.e-10;
pp.query("eps_divu_sph", maestro::eps_divu_sph);

maestro::divu_iter_factor = 100.;
pp.query("divu_iter_factor", maestro::divu_iter_factor);

maestro::divu_level_factor = 10.;
pp.query("divu_level_factor", maestro::divu_level_factor);

maestro::eps_mac = 1.e-10;
pp.query("eps_mac", maestro::eps_mac);

maestro::eps_mac_max = 1.e-8;
pp.query("eps_mac_max", maestro::eps_mac_max);

maestro::mac_level_factor = 10.;
pp.query("mac_level_factor", maestro::mac_level_factor);

maestro::eps_mac_bottom = 1.e-3;
pp.query("eps_mac_bottom", maestro::eps_mac_bottom);

maestro::eps_hg = 1.e-12;
pp.query("eps_hg", maestro::eps_hg);

maestro::eps_hg_max = 1.e-10;
pp.query("eps_hg_max", maestro::eps_hg_max);

maestro::hg_level_factor = 10.;
pp.query("hg_level_factor", maestro::hg_level_factor);

maestro::eps_hg_bottom = 1.e-4;
pp.query("eps_hg_bottom", maestro::eps_hg_bottom);
