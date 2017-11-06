
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
  integer, save :: RHO, RHOH, RHOX, PRES, TEMP, NSTATE
  integer, save :: VELX, VELY, VELZ

  ! Begin the declarations of the ParmParse parameters

  character (len=:), allocatable, save :: model_file
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
    small_temp = 5.d6;
    small_dens = 1.d-5;

    call amrex_parmparse_build(pp, "maestro")
    call pp%query("model_file", model_file)
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
