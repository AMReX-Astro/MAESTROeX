! inlet_bc_module is a simple container module that holds the parameters
! that are used by physbc to implement the inlet boundary conditions.
! As these are problem-specific, any problem needing inlet boundary
! conditions should create its own version of this module, using this
! outline.

module inlet_bc_module

  use amrex_constants_module
  use network
  use eos_module, only: eos, eos_input_rp
  use eos_type_module

  implicit none

  ! parameters that would be used by physbc in the EXT_DIR sections
  ! would be stored here with the 'save' attribute
  double precision, save :: INLET_RHO
  double precision, save :: INLET_RHOH
  double precision, save :: INLET_TEMP
  double precision, save :: INLET_CS

  logical, save :: inlet_bc_initialized = .false.

contains

  ! here we would initialize the parameters that are module variables.
  ! this routine is called when the base state is defined initially,
  ! and upon restart, just after the base state is read in.
  subroutine set_inlet_bcs()

    type (eos_t) :: eos_state

    eos_state%T     = 10.d0
    eos_state%rho   = 1.d-3
    eos_state%p     = 1.d6
    eos_state%xn(:) = 1.d0

    call eos(eos_input_rp, eos_state)

    INLET_CS   = eos_state%cs

    eos_state%T     = 10.d0
    eos_state%rho   = 5.d-4
    eos_state%p     = 1.d6
    eos_state%xn(:) = 1.d0

    call eos(eos_input_rp, eos_state)

    INLET_RHO  = eos_state%rho
    INLET_RHOH = eos_state%rho * eos_state%h
    INLET_TEMP = eos_state%T

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
