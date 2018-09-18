
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

  integer, allocatable, save :: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp
  integer, allocatable, save :: nscal
  double precision, allocatable, save :: prob_lo(:), prob_hi(:)
  double precision, allocatable, save :: rel_eps

#ifdef AMREX_USE_CUDA
  attributes(managed) :: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp
  attributes(managed) :: nscal
  attributes(managed) :: prob_lo, prob_hi
  attributes(managed) :: rel_eps
#endif


  ! Begin the declarations of the ParmParse parameters

  integer,           allocatable, save :: maestro_verbose
  character (len=:), allocatable, save :: model_file
  logical,           allocatable, save :: perturb_model
  logical,           allocatable, save :: print_init_hse_diag
  double precision,  allocatable, save :: cfl
  logical,           allocatable, save :: use_soundspeed_firstdt
  logical,           allocatable, save :: use_divu_firstdt
  integer,           allocatable, save :: spherical
  logical,           allocatable, save :: octant
  integer,           allocatable, save :: do_2d_planar_octant
  integer,           allocatable, save :: drdxfac
  logical,           allocatable, save :: do_sponge
  double precision,  allocatable, save :: sponge_kappa
  double precision,  allocatable, save :: sponge_center_density
  double precision,  allocatable, save :: sponge_start_factor
  double precision,  allocatable, save :: anelastic_cutoff
  double precision,  allocatable, save :: base_cutoff_density
  double precision,  allocatable, save :: burning_cutoff_density
  double precision,  allocatable, save :: buoyancy_cutoff_factor
  double precision,  allocatable, save :: dpdt_factor
  logical,           allocatable, save :: do_planar_invsq_grav
  double precision,  allocatable, save :: planar_invsq_mass
  logical,           allocatable, save :: evolve_base_state
  logical,           allocatable, save :: use_exact_base_state
  logical,           allocatable, save :: do_eos_h_above_cutoff
  integer,           allocatable, save :: enthalpy_pred_type
  integer,           allocatable, save :: species_pred_type
  logical,           allocatable, save :: use_delta_gamma1_term
  integer,           allocatable, save :: slope_order
  double precision,  allocatable, save :: grav_const
  integer,           allocatable, save :: ppm_type
  integer,           allocatable, save :: beta0_type
  logical,           allocatable, save :: use_linear_grav_in_beta0
  double precision,  allocatable, save :: rotational_frequency
  double precision,  allocatable, save :: co_latitude
  logical,           allocatable, save :: drive_initial_convection
  logical,           allocatable, save :: use_alt_energy_fix
  integer,           allocatable, save :: temp_diffusion_formulation
  integer,           allocatable, save :: thermal_diffusion_type
  logical,           allocatable, save :: limit_conductivity
  character (len=:), allocatable, save :: burner_threshold_species
  double precision,  allocatable, save :: burner_threshold_cutoff
  double precision,  allocatable, save :: reaction_sum_tol
  double precision,  allocatable, save :: small_temp
  double precision,  allocatable, save :: small_dens
  logical,           allocatable, save :: use_eos_e_instead_of_h
  logical,           allocatable, save :: use_pprime_in_tfromp
  integer,           allocatable, save :: s0_interp_type
  integer,           allocatable, save :: w0_interp_type
  integer,           allocatable, save :: s0mac_interp_type
  integer,           allocatable, save :: w0mac_interp_type
  integer,           allocatable, save :: track_grid_losses

#ifdef AMREX_USE_CUDA
  attributes(managed) :: maestro_verbose
  
  attributes(managed) :: perturb_model
  attributes(managed) :: print_init_hse_diag
  attributes(managed) :: cfl
  attributes(managed) :: use_soundspeed_firstdt
  attributes(managed) :: use_divu_firstdt
  attributes(managed) :: spherical
  attributes(managed) :: octant
  attributes(managed) :: do_2d_planar_octant
  attributes(managed) :: drdxfac
  attributes(managed) :: do_sponge
  attributes(managed) :: sponge_kappa
  attributes(managed) :: sponge_center_density
  attributes(managed) :: sponge_start_factor
  attributes(managed) :: anelastic_cutoff
  attributes(managed) :: base_cutoff_density
  attributes(managed) :: burning_cutoff_density
  attributes(managed) :: buoyancy_cutoff_factor
  attributes(managed) :: dpdt_factor
  attributes(managed) :: do_planar_invsq_grav
  attributes(managed) :: planar_invsq_mass
  attributes(managed) :: evolve_base_state
  attributes(managed) :: use_exact_base_state
  attributes(managed) :: do_eos_h_above_cutoff
  attributes(managed) :: enthalpy_pred_type
  attributes(managed) :: species_pred_type
  attributes(managed) :: use_delta_gamma1_term
  attributes(managed) :: slope_order
  attributes(managed) :: grav_const
  attributes(managed) :: ppm_type
  attributes(managed) :: beta0_type
  attributes(managed) :: use_linear_grav_in_beta0
  attributes(managed) :: rotational_frequency
  attributes(managed) :: co_latitude
  attributes(managed) :: drive_initial_convection
  attributes(managed) :: use_alt_energy_fix
  attributes(managed) :: temp_diffusion_formulation
  attributes(managed) :: thermal_diffusion_type
  attributes(managed) :: limit_conductivity
  
  attributes(managed) :: burner_threshold_cutoff
  attributes(managed) :: reaction_sum_tol
  attributes(managed) :: small_temp
  attributes(managed) :: small_dens
  attributes(managed) :: use_eos_e_instead_of_h
  attributes(managed) :: use_pprime_in_tfromp
  attributes(managed) :: s0_interp_type
  attributes(managed) :: w0_interp_type
  attributes(managed) :: s0mac_interp_type
  attributes(managed) :: w0mac_interp_type
  attributes(managed) :: track_grid_losses
#endif

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
    do_eos_h_above_cutoff = .true.;
    enthalpy_pred_type = 1;
    species_pred_type = 1;
    use_delta_gamma1_term = .false.;
    slope_order = 4;
    grav_const = -1.5d10;
    ppm_type = 1;
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
    call pp%query("do_eos_h_above_cutoff", do_eos_h_above_cutoff)
    call pp%query("enthalpy_pred_type", enthalpy_pred_type)
    call pp%query("species_pred_type", species_pred_type)
    call pp%query("use_delta_gamma1_term", use_delta_gamma1_term)
    call pp%query("slope_order", slope_order)
    call pp%query("grav_const", grav_const)
    call pp%query("ppm_type", ppm_type)
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
