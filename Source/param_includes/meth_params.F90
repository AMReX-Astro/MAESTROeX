
! This file is automatically created by parse_maestro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in read_method_params().

module meth_params_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! variables in the module

  integer, save :: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp
  integer, save :: nscal
  double precision, save :: prob_lo(3), prob_hi(3)
  double precision, save :: rel_eps


  ! Begin the declarations of the ParmParse parameters

  integer                       , save :: maestro_verbose
  character (len=:), allocatable, save :: model_file
  logical                       , save :: perturb_model
  logical                       , save :: print_init_hse_diag
  double precision              , save :: cfl
  logical                       , save :: use_soundspeed_firstdt
  logical                       , save :: use_divu_firstdt
  integer                       , save :: spherical
  logical                       , save :: octant
  integer                       , save :: do_2d_planar_octant
  integer                       , save :: drdxfac
  logical                       , save :: use_tpert_in_tagging
  logical                       , save :: do_sponge
  double precision              , save :: sponge_kappa
  double precision              , save :: sponge_center_density
  double precision              , save :: sponge_start_factor
  double precision              , save :: anelastic_cutoff
  double precision              , save :: base_cutoff_density
  double precision              , save :: burning_cutoff_density
  double precision              , save :: buoyancy_cutoff_factor
  double precision              , save :: dpdt_factor
  logical                       , save :: do_planar_invsq_grav
  double precision              , save :: planar_invsq_mass
  logical                       , save :: evolve_base_state
  logical                       , save :: use_exact_base_state
  logical                       , save :: average_base_state
  logical                       , save :: do_smallscale
  logical                       , save :: do_eos_h_above_cutoff
  integer                       , save :: enthalpy_pred_type
  integer                       , save :: species_pred_type
  logical                       , save :: use_delta_gamma1_term
  integer                       , save :: slope_order
  double precision              , save :: grav_const
  integer                       , save :: ppm_type
  integer                       , save :: ppm_trace_forces
  integer                       , save :: beta0_type
  logical                       , save :: use_linear_grav_in_beta0
  double precision              , save :: rotational_frequency
  double precision              , save :: co_latitude
  logical                       , save :: drive_initial_convection
  logical                       , save :: use_alt_energy_fix
  integer                       , save :: temp_diffusion_formulation
  integer                       , save :: thermal_diffusion_type
  logical                       , save :: limit_conductivity
  character (len=:), allocatable, save :: burner_threshold_species
  double precision              , save :: burner_threshold_cutoff
  double precision              , save :: reaction_sum_tol
  double precision              , save :: small_temp
  double precision              , save :: small_dens
  logical                       , save :: use_eos_e_instead_of_h
  logical                       , save :: use_pprime_in_tfromp
  integer                       , save :: s0_interp_type
  integer                       , save :: w0_interp_type
  integer                       , save :: s0mac_interp_type
  integer                       , save :: w0mac_interp_type
  integer                       , save :: track_grid_losses

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_method_params() bind(C, name="read_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp

    maestro_verbose = 1;
    allocate(character(len=1)::model_file)
    model_file = "";
    perturb_model = .false.;
    print_init_hse_diag = .false.;
    cfl = 0.5d0;
    use_soundspeed_firstdt = .false.;
    use_divu_firstdt = .false.;
    spherical = 0;
    octant = .false.;
    do_2d_planar_octant = 0;
    drdxfac = 1;
    use_tpert_in_tagging = .false.;
    do_sponge = .false.;
    sponge_kappa = 10.d0;
    sponge_center_density = 3.d6;
    sponge_start_factor = 3.333d0;
    anelastic_cutoff = 3.d6;
    base_cutoff_density = 3.d6;
    burning_cutoff_density = 3.d6;
    buoyancy_cutoff_factor = 5.0d0;
    dpdt_factor = 0.0d0;
    do_planar_invsq_grav = .false.;
    planar_invsq_mass = 0.0d0;
    evolve_base_state = .true.;
    use_exact_base_state = .false.;
    average_base_state = .false.;
    do_smallscale = .false.;
    do_eos_h_above_cutoff = .true.;
    enthalpy_pred_type = 1;
    species_pred_type = 1;
    use_delta_gamma1_term = .false.;
    slope_order = 4;
    grav_const = -1.5d10;
    ppm_type = 1;
    ppm_trace_forces = 0;
    beta0_type = 1;
    use_linear_grav_in_beta0 = .false.;
    rotational_frequency = 0.0d0;
    co_latitude = 0.0d0;
    drive_initial_convection = .false.;
    use_alt_energy_fix = .true.;
    temp_diffusion_formulation = 2;
    thermal_diffusion_type = 1;
    limit_conductivity = .false.;
    allocate(character(len=1)::burner_threshold_species)
    burner_threshold_species = "";
    burner_threshold_cutoff = 1.d-10;
    reaction_sum_tol = 1.d-10;
    small_temp = 5.d6;
    small_dens = 1.d-5;
    use_eos_e_instead_of_h = .false.;
    use_pprime_in_tfromp = .false.;
    s0_interp_type = 3;
    w0_interp_type = 2;
    s0mac_interp_type = 1;
    w0mac_interp_type = 1;
    track_grid_losses = 0;

    call amrex_parmparse_build(pp, "maestro")
    call pp%query("maestro_verbose", maestro_verbose)
    call pp%query("model_file", model_file)
    call pp%query("perturb_model", perturb_model)
    call pp%query("print_init_hse_diag", print_init_hse_diag)
    call pp%query("cfl", cfl)
    call pp%query("use_soundspeed_firstdt", use_soundspeed_firstdt)
    call pp%query("use_divu_firstdt", use_divu_firstdt)
    call pp%query("spherical", spherical)
    call pp%query("octant", octant)
    call pp%query("do_2d_planar_octant", do_2d_planar_octant)
    call pp%query("drdxfac", drdxfac)
    call pp%query("use_tpert_in_tagging", use_tpert_in_tagging)
    call pp%query("do_sponge", do_sponge)
    call pp%query("sponge_kappa", sponge_kappa)
    call pp%query("sponge_center_density", sponge_center_density)
    call pp%query("sponge_start_factor", sponge_start_factor)
    call pp%query("anelastic_cutoff", anelastic_cutoff)
    call pp%query("base_cutoff_density", base_cutoff_density)
    call pp%query("burning_cutoff_density", burning_cutoff_density)
    call pp%query("buoyancy_cutoff_factor", buoyancy_cutoff_factor)
    call pp%query("dpdt_factor", dpdt_factor)
    call pp%query("do_planar_invsq_grav", do_planar_invsq_grav)
    call pp%query("planar_invsq_mass", planar_invsq_mass)
    call pp%query("evolve_base_state", evolve_base_state)
    call pp%query("use_exact_base_state", use_exact_base_state)
    call pp%query("average_base_state", average_base_state)
    call pp%query("do_smallscale", do_smallscale)
    call pp%query("do_eos_h_above_cutoff", do_eos_h_above_cutoff)
    call pp%query("enthalpy_pred_type", enthalpy_pred_type)
    call pp%query("species_pred_type", species_pred_type)
    call pp%query("use_delta_gamma1_term", use_delta_gamma1_term)
    call pp%query("slope_order", slope_order)
    call pp%query("grav_const", grav_const)
    call pp%query("ppm_type", ppm_type)
    call pp%query("ppm_trace_forces", ppm_trace_forces)
    call pp%query("beta0_type", beta0_type)
    call pp%query("use_linear_grav_in_beta0", use_linear_grav_in_beta0)
    call pp%query("rotational_frequency", rotational_frequency)
    call pp%query("co_latitude", co_latitude)
    call pp%query("drive_initial_convection", drive_initial_convection)
    call pp%query("use_alt_energy_fix", use_alt_energy_fix)
    call pp%query("temp_diffusion_formulation", temp_diffusion_formulation)
    call pp%query("thermal_diffusion_type", thermal_diffusion_type)
    call pp%query("limit_conductivity", limit_conductivity)
    call pp%query("burner_threshold_species", burner_threshold_species)
    call pp%query("burner_threshold_cutoff", burner_threshold_cutoff)
    call pp%query("reaction_sum_tol", reaction_sum_tol)
    call pp%query("small_temp", small_temp)
    call pp%query("small_dens", small_dens)
    call pp%query("use_eos_e_instead_of_h", use_eos_e_instead_of_h)
    call pp%query("use_pprime_in_tfromp", use_pprime_in_tfromp)
    call pp%query("s0_interp_type", s0_interp_type)
    call pp%query("w0_interp_type", w0_interp_type)
    call pp%query("s0mac_interp_type", s0mac_interp_type)
    call pp%query("w0mac_interp_type", w0mac_interp_type)
    call pp%query("track_grid_losses", track_grid_losses)
    call amrex_parmparse_destroy(pp)



  end subroutine read_method_params


  subroutine finalize_meth_params() bind(C, name="finalize_meth_params")
    implicit none

    if (allocated(model_file)) then
        deallocate(model_file)
    end if
    if (allocated(burner_threshold_species)) then
        deallocate(burner_threshold_species)
    end if


    
  end subroutine finalize_meth_params

end module meth_params_module
