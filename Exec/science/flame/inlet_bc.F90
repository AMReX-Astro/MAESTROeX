! inlet_bc_module is a simple container module that holds the parameters
! that are used by physbc to implement the inlet boundary conditions.
! As these are problem-specific, any problem needing inlet boundary
! conditions should create its own version of this module, using this
! outline.

module inlet_bc_module

  use amrex_constants_module
  use amrex_error_module
  use network
  use eos_type_module
  use eos_module
  use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel

  implicit none

  ! parameters that would be used by physbc in the EXT_DIR sections
  ! would be stored here with the 'save' attribute
  double precision, save :: INLET_RHO
  double precision, save :: INLET_RHOH
  double precision, save :: INLET_TEMP
  double precision, save :: INLET_RHOX(nspec)
  double precision, save :: INLET_VEL

  logical, save :: inlet_bc_initialized = .false.

contains

  ! here we would initialize the parameters that are module variables.
  ! this routine is called when the base state is defined initially,
  ! and upon restart, just after the base state is read in.
  subroutine set_inlet_bcs() bind(C, name="set_inlet_bcs")

    integer :: ic12, io16

    type (eos_t) :: eos_state

    ! figure out the indices for different species
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")

    if (ic12 < 0 .or. io16 < 0) then
       call amrex_error("ERROR: species indices undefined in inlet_bc")
    endif

    eos_state%rho = dens_fuel
    eos_state%T   = temp_fuel

    eos_state%xn(:)    = ZERO
    eos_state%xn(ic12) = xc12_fuel
    eos_state%xn(io16) = 1.d0 - xc12_fuel

    call eos(eos_input_rt, eos_state)

    INLET_RHO     = dens_fuel
    INLET_RHOH    = dens_fuel*eos_state%h
    INLET_TEMP    = temp_fuel
    INLET_RHOX(:) = dens_fuel*eos_state%xn(:)
    INLET_VEL     = vel_fuel

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

  subroutine get_inlet_bcs(params) bind(C, name="get_inlet_bcs")
    
    double precision, intent(inout) :: params(4+nspec)
    
    params(1) = INLET_RHO
    params(2) = INLET_RHOH
    params(3) = INLET_TEMP
    params(4:nspec+3) = INLET_RHOX(:)
    params(nspec+4) = INLET_VEL
    
  end subroutine get_inlet_bcs
  
end module inlet_bc_module
