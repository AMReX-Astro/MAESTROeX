
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

  double precision, allocatable, save, public :: prob_lo(:)
  double precision, allocatable, save, public :: prob_hi(:)

  ! Begin the declarations of the ParmParse parameters

  character (len=:), allocatable, save :: model_file
  integer                       , save :: spherical
  logical                       , save :: octant
  real(rt)                      , save :: anelastic_cutoff
  real(rt)                      , save :: base_cutoff_density
  real(rt)                      , save :: burning_cutoff_density
  real(rt)                      , save :: rotational_frequency
  real(rt)                      , save :: co_latitude
  real(rt)                      , save :: small_temp
  real(rt)                      , save :: small_dens

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_method_params() bind(C, name="read_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp

    allocate(character(len=1)::model_file)
    model_file = "model.hse";
    spherical = 0;
    octant = .false.;
    anelastic_cutoff = 3.d6;
    base_cutoff_density = 3.d6;
    burning_cutoff_density = 3.d6;
    rotational_frequency = 0.0d0;
    co_latitude = 0.0d0;
    small_temp = 5.d6;
    small_dens = 1.d-5;

    call amrex_parmparse_build(pp, "maestro")
    call pp%query("model_file", model_file)
    call pp%query("spherical", spherical)
    call pp%query("octant", octant)
    call pp%query("anelastic_cutoff", anelastic_cutoff)
    call pp%query("base_cutoff_density", base_cutoff_density)
    call pp%query("burning_cutoff_density", burning_cutoff_density)
    call pp%query("rotational_frequency", rotational_frequency)
    call pp%query("co_latitude", co_latitude)
    call pp%query("small_temp", small_temp)
    call pp%query("small_dens", small_dens)
    call amrex_parmparse_destroy(pp)



  end subroutine read_method_params


  subroutine finalize_meth_params() bind(C, name="finalize_meth_params")
    implicit none

    if (allocated(model_file)) then
        deallocate(model_file)
    end if


    
  end subroutine finalize_meth_params

end module meth_params_module
