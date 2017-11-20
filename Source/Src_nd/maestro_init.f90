
module maestro_init_module

  use network, only: network_init, nspec, short_spec_names
  use parallel, only: parallel_IOProcessor
  use amrex_fort_module, only: amrex_spacedim
  use model_parser_module
  use meth_params_module, only: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp, &
                                nscal, small_dens, small_temp
  use eos_module, only: eos_init

  implicit none

  private

contains

  subroutine maestro_network_init() bind(C, name="maestro_network_init")


    call network_init()

  end subroutine maestro_network_init


  subroutine get_num_spec(nspec_out) bind(C, name="get_num_spec")

    integer, intent(out) :: nspec_out

    nspec_out = nspec

  end subroutine get_num_spec


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

  subroutine set_method_params(Density,Enthalpy,FirstSpec,Temperature,Pressure,Nscalars) &
       bind(C, name="set_method_params")

    integer, intent(in) :: Density, Enthalpy, FirstSpec, Temperature, Pressure, Nscalars

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

    nscal = Nscalars

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
          call bl_warning("Warning:: small_dens has not been set, defaulting to 1.d-200.")
       endif
       small_dens = 1.d-200
    endif

    if (small_temp <= 0.d0) then
       if (ioproc == 1) then
          call bl_warning("Warning:: small_temp has not been set, defaulting to 1.d-200.")
       endif
       small_temp = 1.d-200
    endif

    ! Note that the EOS may modify our choices because of its
    ! internal limitations, so the small_dens and small_temp
    ! may be modified coming back out of this routine.

    call eos_init(small_temp, small_dens)

  end subroutine set_method_params

end module maestro_init_module
