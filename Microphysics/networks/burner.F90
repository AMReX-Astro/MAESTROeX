module burner_module

  use amrex_constants_module
  use network
  use eos_module
  use burn_type_module
  use actual_burner_module

  logical :: burner_initialized = .false.

contains

  subroutine burner_init() bind(C, name="burner_init")

    implicit none

    call actual_burner_init()

    burner_initialized = .true.

  end subroutine burner_init



  subroutine burner(state_in, state_out, dt, time_in)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    double precision, intent(in) :: dt
    double precision, intent(in) :: time_in

    !$gpu

    ! Make sure the network and burner have been initialized.

#if !defined(ACC) && !defined(AMREX_USE_CUDA)
    if (.NOT. network_initialized) then
       call amrex_error("ERROR in burner: must initialize network first.")
    endif

    if (.NOT. burner_initialized) then
       call amrex_error("ERROR in burner: must initialize burner first.")
    endif
#endif

    ! Initialize the final state by assuming it does not change.

    call copy_burn_t(state_out, state_in)

    ! Do the burning.

    call actual_burner(state_in, state_out, dt, time_in)

  end subroutine burner

end module burner_module
