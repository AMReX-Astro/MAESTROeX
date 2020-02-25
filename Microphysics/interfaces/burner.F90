module burner_module

  use amrex_fort_module, only: rt => amrex_real
  use burn_type_module, only: burn_t

  logical, save :: burner_initialized = .false.

contains

  subroutine burner_init() bind(C, name="burner_init")

#ifdef SIMPLIFIED_SDC
    use integrator_module, only: integrator_init
#else
    use actual_burner_module, only: actual_burner_init
#endif

    implicit none

#ifdef SIMPLIFIED_SDC
    call integrator_init()
#else
    call actual_burner_init()
#endif

    burner_initialized = .true.

  end subroutine burner_init



#ifndef SIMPLIFIED_SDC
  subroutine burner(state_in, state_out, dt, time)

    use actual_burner_module, only: actual_burner

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt), intent(in) :: dt, time

    !$gpu

    call actual_burner(state_in, state_out, dt, time)

  end subroutine burner
#endif

end module burner_module
