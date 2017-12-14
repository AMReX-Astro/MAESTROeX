
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
  double precision              , save :: pert_temp_factor
  double precision              , save :: pert_rad_factor
  logical                       , save :: do_small_domain
  double precision              , save :: cfl
  logical                       , save :: use_soundspeed_firstdt
  integer                       , save :: spherical
  logical                       , save :: octant
  integer                       , save :: do_2d_planar_octant
  integer                       , save :: drdxfac
  double precision              , save :: sponge_kappa
  double precision              , save :: sponge_center_density
  double precision              , save :: sponge_start_factor
  double precision              , save :: anelastic_cutoff
  double precision              , save :: base_cutoff_density
  double precision              , save :: burning_cutoff_density
  logical                       , save :: do_planar_invsq_grav
  double precision              , save :: planar_invsq_mass
  double precision              , save :: grav_const
  integer                       , save :: beta0_type
  logical                       , save :: use_linear_grav_in_beta0
  double precision              , save :: rotational_frequency
  double precision              , save :: co_latitude
  logical                       , save :: drive_initial_convection
  character (len=:), allocatable, save :: burner_threshold_species
  double precision              , save :: burner_threshold_cutoff
  double precision              , save :: reaction_sum_tol
  double precision              , save :: small_temp
  double precision              , save :: small_dens
  logical                       , save :: use_eos_e_instead_of_h
  logical                       , save :: use_pprime_in_tfromp

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_method_params() bind(C, name="read_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp

    maestro_verbose = 1;
    allocate(character(len=1)::model_file)
    model_file = "model.hse";
    perturb_model = .false.;
    print_init_hse_diag = .false.;
    pert_temp_factor = 1.0d0;
    pert_rad_factor = 1.0d0;
    do_small_domain = .false.;
    cfl = 0.5d0;
    use_soundspeed_firstdt = .false.;
    spherical = 0;
    octant = .false.;
    do_2d_planar_octant = 0;
    drdxfac = 1;
    sponge_kappa = 10.d0;
    sponge_center_density = 3.d6;
    sponge_start_factor = 3.333d0;
    anelastic_cutoff = 3.d6;
    base_cutoff_density = 3.d6;
    burning_cutoff_density = 3.d6;
    do_planar_invsq_grav = .false.;
    planar_invsq_mass = 0.0d0;
    grav_const = -1.5d10;
    beta0_type = 1;
    use_linear_grav_in_beta0 = .false.;
    rotational_frequency = 0.0d0;
    co_latitude = 0.0d0;
    drive_initial_convection = .false.;
    allocate(character(len=1)::burner_threshold_species)
    burner_threshold_species = "";
    burner_threshold_cutoff = 1.d-10;
    reaction_sum_tol = 1.d-10;
    small_temp = 5.d6;
    small_dens = 1.d-5;
    use_eos_e_instead_of_h = .false.;
    use_pprime_in_tfromp = .false.;

    call amrex_parmparse_build(pp, "maestro")
    call pp%query("maestro_verbose", maestro_verbose)
    call pp%query("model_file", model_file)
    call pp%query("perturb_model", perturb_model)
    call pp%query("print_init_hse_diag", print_init_hse_diag)
    call pp%query("pert_temp_factor", pert_temp_factor)
    call pp%query("pert_rad_factor", pert_rad_factor)
    call pp%query("do_small_domain", do_small_domain)
    call pp%query("cfl", cfl)
    call pp%query("use_soundspeed_firstdt", use_soundspeed_firstdt)
    call pp%query("spherical", spherical)
    call pp%query("octant", octant)
    call pp%query("do_2d_planar_octant", do_2d_planar_octant)
    call pp%query("drdxfac", drdxfac)
    call pp%query("sponge_kappa", sponge_kappa)
    call pp%query("sponge_center_density", sponge_center_density)
    call pp%query("sponge_start_factor", sponge_start_factor)
    call pp%query("anelastic_cutoff", anelastic_cutoff)
    call pp%query("base_cutoff_density", base_cutoff_density)
    call pp%query("burning_cutoff_density", burning_cutoff_density)
    call pp%query("do_planar_invsq_grav", do_planar_invsq_grav)
    call pp%query("planar_invsq_mass", planar_invsq_mass)
    call pp%query("grav_const", grav_const)
    call pp%query("beta0_type", beta0_type)
    call pp%query("use_linear_grav_in_beta0", use_linear_grav_in_beta0)
    call pp%query("rotational_frequency", rotational_frequency)
    call pp%query("co_latitude", co_latitude)
    call pp%query("drive_initial_convection", drive_initial_convection)
    call pp%query("burner_threshold_species", burner_threshold_species)
    call pp%query("burner_threshold_cutoff", burner_threshold_cutoff)
    call pp%query("reaction_sum_tol", reaction_sum_tol)
    call pp%query("small_temp", small_temp)
    call pp%query("small_dens", small_dens)
    call pp%query("use_eos_e_instead_of_h", use_eos_e_instead_of_h)
    call pp%query("use_pprime_in_tfromp", use_pprime_in_tfromp)
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
