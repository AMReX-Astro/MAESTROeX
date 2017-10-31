
! This file is automatically created by parse_maestro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in ca_set_maestro_method_params().

module meth_params_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! variables in the module
  integer, save :: RHO, RHOH, RHOX, PRES, TEMP, NSTATE
  integer, save :: VELX, VELY, VELZ

  ! Begin the declarations of the ParmParse parameters

  real(rt), save :: cflfac
  real(rt), save :: small_temp
  real(rt), save :: small_dens
  integer         , save :: use_tfromp
  integer         , save :: use_eos_e_instead_of_h
  integer         , save :: use_pprime_in_tfromp

  ! End the declarations of the ParmParse parameters

contains

  subroutine ca_set_maestro_method_params() bind(C, name="ca_set_maestro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp

    cflfac = 0.8d0;
    small_temp = 5.d6;
    small_dens = 1.d-5;
    use_tfromp = 0;
    use_eos_e_instead_of_h = 0;
    use_pprime_in_tfromp = 0;

    call amrex_parmparse_build(pp, "maestro")
    call pp%query("cflfac", cflfac)
    call pp%query("small_temp", small_temp)
    call pp%query("small_dens", small_dens)
    call pp%query("use_tfromp", use_tfromp)
    call pp%query("use_eos_e_instead_of_h", use_eos_e_instead_of_h)
    call pp%query("use_pprime_in_tfromp", use_pprime_in_tfromp)
    call amrex_parmparse_destroy(pp)



  end subroutine ca_set_maestro_method_params


  subroutine ca_finalize_meth_params() bind(C, name="ca_finalize_meth_params")
    implicit none



    
  end subroutine ca_finalize_meth_params

end module meth_params_module
