module burn_type_module

  use bl_types, only: dp_t
#ifdef REACT_SPARSE_JACOBIAN
  use actual_network, only: nspec, nspec_evolve, naux, NETWORK_SPARSE_JAC_NNZ
#else
  use actual_network, only: nspec, nspec_evolve, naux
#endif
  
  implicit none

  ! A generic structure holding data necessary to do a nuclear burn.

  ! Set the number of independent variables -- this should be
  ! temperature, enuc + the number of species which participate
  ! in the evolution equations.

  integer, parameter :: neqs = 2 + nspec_evolve

  ! Indices of the temperature and energy variables in the work arrays.

  integer, parameter :: net_itemp = nspec_evolve + 1
  integer, parameter :: net_ienuc = nspec_evolve + 2

  type :: burn_t

    real(dp_t) :: rho
    real(dp_t) :: T
    real(dp_t) :: e
    real(dp_t) :: xn(nspec)
#if naux > 0
    real(dp_t) :: aux(naux)
#endif

    real(dp_t) :: cv
    real(dp_t) :: cp
    real(dp_t) :: y_e
    real(dp_t) :: eta
    real(dp_t) :: cs
    real(dp_t) :: dx
    real(dp_t) :: abar
    real(dp_t) :: zbar

    ! Last temperature we evaluated the EOS at
    real(dp_t) :: T_old

    ! Temperature derivatives of specific heat
    real(dp_t) :: dcvdT
    real(dp_t) :: dcpdT

    ! The following are the actual integration data.
    ! To avoid potential incompatibilities we won't
    ! include the integration vector y itself here.
    ! It can be reconstructed from all of the above
    ! data, particularly xn, e, and T.

    real(dp_t) :: ydot(neqs)
    real(dp_t) :: jac(neqs, neqs)

    ! Whether we are self-heating or not.

    logical          :: self_heat

    ! Zone index information.

    integer          :: i
    integer          :: j
    integer          :: k

    ! diagnostics
    integer :: n_rhs
    integer :: n_jac

    ! Integration time.

    real(dp_t) :: time

    ! Was the burn successful?

    logical :: success

  end type burn_t

contains

  ! Implement a manual copy routine since CUDA Fortran doesn't
  ! (yet) support derived type copying on the device.
  subroutine copy_burn_t(to_state, from_state)

    implicit none

    type (burn_t), intent(in   ) :: from_state
    type (burn_t), intent(  out) :: to_state

    to_state % rho = from_state % rho
    to_state % T   = from_state % T
    to_state % e   = from_state % e
    to_state % xn(1:nspec) = from_state % xn(1:nspec)

#if naux > 0
    to_state % aux(1:naux) = from_state % aux(1:naux)
#endif

    to_state % cv  = from_state % cv
    to_state % cp  = from_state % cp
    to_state % y_e = from_state % y_e
    to_state % eta = from_state % eta
    to_state % cs  = from_state % cs
    to_state % dx  = from_state % dx

    to_state % abar = from_state % abar
    to_state % zbar = from_state % zbar

    to_state % T_old = from_state % T_old

    to_state % dcvdT = from_state % dcvdT
    to_state % dcpdT = from_state % dcpdT

    to_state % ydot(1:neqs) = from_state % ydot(1:neqs)
    to_state % jac(1:neqs, 1:neqs) = from_state % jac(1:neqs, 1:neqs)

    to_state % self_heat = from_state % self_heat

    to_state % i = from_state % i
    to_state % j = from_state % j
    to_state % k = from_state % k

    to_state % n_rhs = from_state % n_rhs
    to_state % n_jac = from_state % n_jac

    to_state % time = from_state % time

    to_state % success = from_state % success

  end subroutine copy_burn_t


  ! Given an eos type, copy the data relevant to the burn type.

  subroutine eos_to_burn(eos_state, burn_state)

    !$acc routine seq

    use eos_type_module, only: eos_t

    implicit none

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    burn_state % rho  = eos_state % rho
    burn_state % T    = eos_state % T
    burn_state % e    = eos_state % e
    burn_state % xn   = eos_state % xn
#if naux > 0
    burn_state % aux  = eos_state % aux
#endif
    burn_state % cv   = eos_state % cv
    burn_state % cp   = eos_state % cp
    burn_state % y_e  = eos_state % y_e
    burn_state % eta  = eos_state % eta
    burn_state % cs   = eos_state % cs
    burn_state % abar = eos_state % abar
    burn_state % zbar = eos_state % zbar

  end subroutine eos_to_burn



  ! Given a burn type, copy the data relevant to the eos type.

  subroutine burn_to_eos(burn_state, eos_state)

    !$acc routine seq

    use eos_type_module, only: eos_t

    implicit none

    type (burn_t) :: burn_state
    type (eos_t)  :: eos_state

    eos_state % rho  = burn_state % rho
    eos_state % T    = burn_state % T
    eos_state % e    = burn_state % e
    eos_state % xn   = burn_state % xn
#if naux > 0
    eos_state % aux  = burn_state % aux
#endif
    eos_state % cv   = burn_state % cv
    eos_state % cp   = burn_state % cp
    eos_state % y_e  = burn_state % y_e
    eos_state % eta  = burn_state % eta
    eos_state % cs   = burn_state % cs
    eos_state % abar = burn_state % abar
    eos_state % zbar = burn_state % zbar

  end subroutine burn_to_eos


  subroutine normalize_abundances_burn(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use extern_probin_module, only: small_x

    implicit none

    type (burn_t), intent(inout) :: state

    state % xn(:) = max(small_x, min(ONE, state % xn(:)))
    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances_burn

  
#ifdef REACT_SPARSE_JACOBIAN
  subroutine lookup_csr_jac_loc(row, col, csr_loc)

    !$acc routine seq

    use actual_network, only: csr_jac_col_index, csr_jac_row_count

    implicit none

    integer, intent(in   ) :: row, col
    integer, intent(out  ) :: csr_loc

    integer :: num_in_row, row_start_loc, row_end_loc, i

    !$gpu

    ! Looks up the index into a CSR-formatted Jacobian
    ! matrix given row and col indices into the
    ! equivalent dense matrix.
    !
    ! Assumes the base in first element of CSR row count array is 1

    num_in_row = csr_jac_row_count(row+1) - csr_jac_row_count(row)
    row_start_loc = csr_jac_row_count(row)
    row_end_loc   = row_start_loc + num_in_row - 1

    csr_loc = -1
    do i = row_start_loc, row_end_loc
       if (csr_jac_col_index(i) == col) then
          csr_loc = i
          exit
       endif
    enddo
  end subroutine lookup_csr_jac_loc


  subroutine set_csr_jac_entry(csr_jac, row, col, val)

    !$acc routine seq

    implicit none

    real(dp_t), intent(inout) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)
    integer   , intent(in   ) :: row, col
    real(dp_t), intent(in   ) :: val

    integer :: csr_loc

    !$gpu

    ! Get index into the CSR Jacobian
    call lookup_csr_jac_loc(row, col, csr_loc)

    ! Set value in CSR Jacobian only if row, col entry exists
    if (csr_loc /= -1) then
       csr_jac(csr_loc) = val
    endif

  end subroutine set_csr_jac_entry


  subroutine scale_csr_jac_entry(csr_jac, row, col, val)

    !$acc routine seq

    implicit none

    real(dp_t), intent(inout) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)
    integer   , intent(in   ) :: row, col
    real(dp_t), intent(in   ) :: val

    integer :: csr_loc

    !$gpu

    ! Get index into the CSR Jacobian
    call lookup_csr_jac_loc(row, col, csr_loc)

    ! Scale value in CSR Jacobian only if row, col entry exists
    if (csr_loc /= -1) then
       csr_jac(csr_loc) = csr_jac(csr_loc) * val
    endif

  end subroutine scale_csr_jac_entry


  subroutine get_csr_jac_entry(csr_jac, row, col, val)

    !$acc routine seq

    implicit none

    real(dp_t), intent(in   ) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)
    integer   , intent(in   ) :: row, col
    real(dp_t), intent(out  ) :: val

    integer :: csr_loc

    !$gpu

    ! Get index into the CSR Jacobian
    call lookup_csr_jac_loc(row, col, csr_loc)

    ! Get value from CSR Jacobian only if row, col entry exists
    if (csr_loc /= -1) then
       val = csr_jac(csr_loc)
    endif

  end subroutine get_csr_jac_entry
#endif


  subroutine set_jac_entry(state, row, col, val)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state
    integer      , intent(in   ) :: row, col
    real(dp_t)   , intent(in   ) :: val

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    call set_csr_jac_entry(state % sparse_jac, row, col, val)
#else
    state % jac(row, col) = val
#endif

  end subroutine set_jac_entry


  subroutine scale_jac_entry(state, row, col, val)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state
    integer      , intent(in   ) :: row, col
    real(dp_t)   , intent(in   ) :: val

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    call scale_csr_jac_entry(state % sparse_jac, row, col, val)
#else
    state % jac(row, col) = state % jac(row, col) * val
#endif

  end subroutine scale_jac_entry  


  subroutine get_jac_entry(state, row, col, val)

    !$acc routine seq

    implicit none

    type (burn_t), intent(in   ) :: state
    integer      , intent(in   ) :: row, col
    real(dp_t)   , intent(out  ) :: val

    integer :: csr_loc

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    call get_csr_jac_entry(state % sparse_jac, row, col, val)
#else
    val = state % jac(row, col)
#endif

  end subroutine get_jac_entry


  subroutine set_jac_zero(state)

    !$acc routine seq

    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t), intent(inout) :: state

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    state % sparse_jac(:) = ZERO
#else
    state % jac(:,:) = ZERO
#endif

  end subroutine set_jac_zero
  
end module burn_type_module
