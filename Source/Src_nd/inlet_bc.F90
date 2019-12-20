module inlet_bc_module
  ! inlet_bc_module is a simple container module that holds the parameters
  ! that are used by physbc to implement the inlet boundary conditions.
  ! As these are problem-specific, any problem needing inlet boundary
  ! conditions should create its own version of this module, using this
  ! outline.

  implicit none

  ! parameters that would be used by physbc in the EXT_DIR sections
  ! would be stored here with the 'save' attribute

  logical, save :: inlet_bc_initialized = .false.

contains

  subroutine set_inlet_bcs()
    ! here we would initialize the parameters that are module variables.
    ! this routine is called when the base state is defined initially,
    ! and upon restart, just after the base state is read in.
    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
