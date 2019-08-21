module burner_module

  use amrex_constants_module
  use network
  use eos_module
#ifndef SDC
  use burn_type_module, only: burn_t
#else
#if (SDC_METHOD == 1)
  use integrator_module
#elif (SDC_METHOD == 2)
  use sdc_type_module, only: sdc_t
#endif
  use actual_burner_module

  logical :: burner_initialized = .false.

contains

  subroutine burner_init() bind(C, name="burner_init")

    implicit none

#ifdef SDC
#if (SDC_METHOD == 1)
    call integrator_init()
#elif (SDC_METHOD == 2)
    call actual_burner_init()
#endif
#else
    call actual_burner_init()
#endif

    burner_initialized = .true.

  end subroutine burner_init



#ifndef SDC
  subroutine burner(state_in, state_out, dt)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    double precision, intent(in) :: dt

    double precision :: time

    !$gpu

    time = 0.d0

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

    call actual_burner(state_in, state_out, dt, time)

  end subroutine burner

#else
    subroutine burner(state_in, state_out, dt)

    !$acc routine seq

    implicit none

    type (sdc_t), intent(inout) :: state_in
    type (sdc_t), intent(inout) :: state_out
    double precision, intent(in) :: dt

    double precision :: time

    !$gpu

    time = 0.d0

    ! Make sure the network and burner have been initialized.

#if !defined(ACC) && !defined(AMREX_USE_CUDA)
    if (.NOT. network_initialized) then
       call amrex_error("ERROR in burner: must initialize network first.")
    endif

    if (.NOT. burner_initialized) then
       call amrex_error("ERROR in burner: must initialize burner first.")
    endif
#endif

    ! Do the burning.

    call actual_burner(state_in, state_out, dt, time)

  end subroutine burner
#endif

end module burner_module
