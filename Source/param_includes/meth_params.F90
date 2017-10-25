
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
  integer, save :: RHO, RHOH, RHOX, TEMP, PRES, NUM_STATE
  integer, save :: VELX, VELY, VELZ

  ! Begin the declarations of the ParmParse parameters

  real(rt), save :: cfl

  ! End the declarations of the ParmParse parameters

contains

  subroutine ca_set_maestro_method_params() bind(C, name="ca_set_maestro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp

    cfl = 0.8d0;

    call amrex_parmparse_build(pp, "maestro")
    call pp%query("cfl", cfl)
    call amrex_parmparse_destroy(pp)



  end subroutine ca_set_maestro_method_params


  subroutine ca_finalize_meth_params() bind(C, name="ca_finalize_meth_params")
    implicit none



    
  end subroutine ca_finalize_meth_params

end module meth_params_module
