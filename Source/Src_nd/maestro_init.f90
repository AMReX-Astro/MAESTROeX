
module maestro_init_module

  use amrex_error_module
  use network, only: network_init, nspec, short_spec_names
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use amrex_fort_module, only: amrex_spacedim
  use model_parser_module
  use meth_params_module, only: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp, &
       nscal, small_dens, small_temp, prob_lo, prob_hi, rel_eps
  use eos_module, only: eos_init
  use runtime_init_module
  implicit none

  private

contains

  subroutine maestro_network_init() bind(C, name="maestro_network_init")

    use actual_rhs_module, only: actual_rhs_init
    use burner_loop_module, only: burner_loop_init

    call network_init()
    call actual_rhs_init()
    call burner_loop_init()

  end subroutine maestro_network_init

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine maestro_extern_init() bind(C, name="maestro_extern_init")

    ! initialize the external runtime parameters in
    ! extern_probin_module

    use amrex_fort_module, only: rt => amrex_real
    !
    call runtime_init()

  end subroutine maestro_extern_init

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine maestro_conductivity_init() bind(C, name="maestro_conductivity_init")

    ! initialize the external runtime parameters in
    ! extern_probin_module

    use conductivity_module, only: conductivity_init

    call conductivity_init()

  end subroutine maestro_conductivity_init

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_num_spec(nspec_out) bind(C, name="get_num_spec")

    integer, intent(out) :: nspec_out

    nspec_out = nspec

  end subroutine get_num_spec

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_spec_names(spec_names,ispec,len) bind(C, name="get_spec_names")

    integer, intent(in   ) :: ispec
    integer, intent(inout) :: len
    integer, intent(inout) :: spec_names(len)

    ! Local variables
    integer   :: i

    len = len_trim(short_spec_names(ispec+1))

    do i = 1,len
       spec_names(i) = ichar(short_spec_names(ispec+1)(i:i))
    end do

  end subroutine get_spec_names

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_spec_az(ispec,A,Z) bind(C, name="get_spec_az")

    use network, only: nspec, aion, zion
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: ispec
    real(rt), intent(inout) :: A, Z

    ! C++ is 0-based indexing, so increment
    A = aion(ispec+1)
    Z = zion(ispec+1)

  end subroutine get_spec_az

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine set_method_params(Density,Enthalpy,FirstSpec,Temperature, &
       Pressure,Nscalars,prob_lo_in,prob_hi_in) &
       bind(C, name="set_method_params")

    integer         , intent(in) :: Density, Enthalpy, FirstSpec, Temperature
    integer         , intent(in) :: Pressure, Nscalars
    double precision, intent(in) :: prob_lo_in(3), prob_hi_in(3)

    integer :: ioproc

    !---------------------------------------------------------------------
    ! conserved state components
    !---------------------------------------------------------------------

    if (parallel_IOProcessor()) then
       print*,'Calling set_method_params()'
    end if

    rho_comp  = Density+1
    rhoh_comp = Enthalpy+1
    spec_comp = FirstSpec+1
    temp_comp = Temperature+1
    pi_comp   = Pressure+1

    ! nscal = Nscalars

    prob_lo(1:3) = prob_lo_in(1:3)
    prob_hi(1:3) = prob_hi_in(1:3)

    !---------------------------------------------------------------------
    ! other initializations
    !---------------------------------------------------------------------

    ! This is a routine which links to the C++ ParallelDescriptor class

    call bl_pd_is_ioproc(ioproc)

    !---------------------------------------------------------------------
    ! safety checks
    !---------------------------------------------------------------------

    if (small_dens <= 0.d0) then
       if (ioproc == 1) then
          call amrex_warning("Warning:: small_dens has not been set, defaulting to 1.d-200.")
       endif
       small_dens = 1.d-200
    endif

    if (small_temp <= 0.d0) then
       if (ioproc == 1) then
          call amrex_warning("Warning:: small_temp has not been set, defaulting to 1.d-200.")
       endif
       small_temp = 1.d-200
    endif

    ! Note that the EOS may modify our choices because of its
    ! internal limitations, so the small_dens and small_temp
    ! may be modified coming back out of this routine.

    call eos_init(small_temp, small_dens)

  end subroutine set_method_params

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine set_rel_eps(rel_eps_in) bind(C,name="set_rel_eps")

    double precision, intent(in) :: rel_eps_in

    rel_eps = rel_eps_in

  end subroutine set_rel_eps

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_rel_eps(rel_eps_in) bind(C,name="get_rel_eps")

    double precision, intent(inout) :: rel_eps_in

    rel_eps_in = rel_eps

  end subroutine get_rel_eps

end module maestro_init_module
