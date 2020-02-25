module sdc_type_module

  use amrex_fort_module, only : rt => amrex_real
  use actual_network, only: nspec

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn
  ! in the SDC formalism.

#if defined(SDC_EVOLVE_ENERGY)

  ! these indicies represent the order that the conserved state comes
  ! into the ODE integration from the hydro code.
  !
  ! they also represent the order of the advective sources
  !
  ! integrate rho*X, internal energy, total energy
  ! carry momentum as an unevolved variable

  integer, parameter :: SEDEN = 1
  integer, parameter :: SEINT = 2
  integer, parameter :: SFS   = 3
  integer, parameter :: SRHO  = SFS + nspec
  integer, parameter :: SMX   = SRHO + 1
  integer, parameter :: SMY   = SRHO + 2
  integer, parameter :: SMZ   = SRHO + 3

  integer, parameter :: SVAR  = SMZ
  integer, parameter :: SVAR_EVOLVE = SRHO - 1

#elif defined(SDC_EVOLVE_ENTHALPY)

  ! integrate rho*X (species masses) and rho*h (enthalpy)
  ! carry pressure for EOS calls in RHS

  integer, parameter :: SFS = 1
  integer, parameter :: SENTH = SFS + nspec
  integer, parameter :: SVAR  = SENTH
  integer, parameter :: SVAR_EVOLVE = SVAR

#endif
  
  type :: sdc_t

     real(rt) :: y(SVAR)
     real(rt) :: ydot_a(SVAR)

#if defined(SDC_EVOLVE_ENERGY)
     logical :: T_from_eden
#elif defined(SDC_EVOLVE_ENTHALPY)
     ! Pressure in case we wish to use it for EOS calls
     real(rt) :: p0
     ! Density is defined by sum(rho*X) = rho in this method
     real(rt) :: rho
#endif

     integer :: i
     integer :: j
     integer :: k

     integer :: n_rhs
     integer :: n_jac

     integer :: sdc_iter

     logical :: success
  end type sdc_t

end module sdc_type_module
