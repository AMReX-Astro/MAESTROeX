module actual_rhs_module

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state, ydot)

    !$acc routine seq

    use burn_type_module, only: burn_t, neqs
    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state
    double precision :: ydot(neqs)
    !$gpu

    ! Do nothing in this RHS.

    ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state, jac)

    !$acc routine seq

    use burn_type_module, only: burn_t, njrows, njcols
    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state
    double precision :: jac(njrows, njcols)

    !$gpu

    ! Do nothing in this RHS.

    jac(:,:) = ZERO

  end subroutine actual_jac



  subroutine update_unevolved_species(state)

    !$acc routine seq

    use burn_type_module, only: burn_t

    implicit none

    type (burn_t)    :: state

    !$gpu

  end subroutine update_unevolved_species

end module actual_rhs_module
