module microphysics_module

  use network
  use eos_module, only : eos_init
#ifdef REACTIONS
  use actual_rhs_module, only : actual_rhs_init
#ifndef SIMPLIFIED_SDC
  use actual_burner_module, only : actual_burner_init
#endif
#endif

#ifdef CONDUCTIVITY
  use actual_conductivity_module, only: actual_conductivity_init
#endif

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine microphysics_init(small_temp, small_dens)

    real(rt)        , optional :: small_temp
    real(rt)        , optional :: small_dens

    call network_init()

    if (present(small_temp) .and. present(small_dens)) then
       call eos_init(small_temp=small_temp, small_dens=small_dens)
    else if (present(small_temp)) then
       call eos_init(small_temp=small_temp)
    else if (present(small_dens)) then
       call eos_init(small_dens=small_dens)
    else
       call eos_init()
    endif

#ifdef REACTIONS
    call actual_rhs_init()
#ifndef SIMPLIFIED_SDC
    call actual_burner_init()
#endif
#endif

#ifdef CONDUCTIVITY
    call actual_conductivity_init()
#endif

  end subroutine microphysics_init

  subroutine microphysics_finalize()

    use eos_module, only: eos_finalize
#ifdef USE_SCREENING
    use screening_module, only: screening_finalize
    call screening_finalize()
#endif
    call eos_finalize()
    call network_finalize()

  end subroutine microphysics_finalize

end module microphysics_module

