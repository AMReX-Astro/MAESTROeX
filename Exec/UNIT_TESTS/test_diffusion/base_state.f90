! init_base_state is used to initialize the base state arrays from the
! model file.  The actual reading of the model file is handled by the
! model_parser_module in Util/
!
! Note: The initial base state quantities returned from this routine
! are only a temporary base state.  These quantities are mapped onto
! the full 2- or 3-d state in initscaldata.f90 and a new base state is
! created after initialization by averaging the density and calling
! enforce_HSE in initialize.f90.

module base_state_module

  use eos_type_module
  use eos_module
  use amrex_constants_module
  use simple_log_module
  use inlet_bc_module
  use network, only: nspec, network_species_index, network_init
  use meth_params_module, only: nscal, rho_comp, &
       rhoh_comp, spec_comp, temp_comp
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use probin_module, only: ambient_h, ambient_dens, &
       ambient_he4, ambient_c12, ambient_fe56
  use extern_probin_module, only: const_conductivity

  implicit none

  private

contains

  subroutine init_base_state(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init) &
       bind(C, name="init_base_state")

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)

    ! local variables
    integer :: ihe4, ic12, ife56
    integer :: n,r
    double precision :: diffusion_coefficient
    type (eos_t) :: eos_state

    call network_init()

    ihe4  = network_species_index("helium-4")
    ic12  = network_species_index("carbon-12")
    ife56 = network_species_index("iron-56")

    if (ihe4 < 0 .or. ic12 < 0 .or. ife56 < 0) then
       print *, ihe4, ic12, ife56
       call bl_error("Invalid species in init_base_state.")
    endif

    eos_state%h         = ambient_h
    eos_state%rho       = ambient_dens

    eos_state%xn(:) = ZERO
    eos_state%xn(ihe4)  = ambient_he4
    eos_state%xn(ic12)  = ambient_c12
    eos_state%xn(ife56) = ambient_fe56

    call eos(eos_input_rh, eos_state)

    diffusion_coefficient = const_conductivity / (eos_state%cp * ambient_dens)
    tempbar(0:max_radial_level,0:nr_fine-1) = diffusion_coefficient

    do n=0,max_radial_level

       do r=0,nr(n)-1

          s0_init(n,r,rho_comp) = eos_state%rho
          s0_init(n,r,rhoh_comp) = eos_state%rho * eos_state%h
          s0_init(n,r,spec_comp:spec_comp+nspec-1) = eos_state%rho * eos_state%xn(1:nspec)
          p0_init(n,r) = eos_state%p
          s0_init(n,r,temp_comp) = eos_state%T

       end do


       ! initialize any inlet BC parameters
       call set_inlet_bcs()

    end do ! end loop over levels

  end subroutine init_base_state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_base_state_irreg(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init, &
       r_cc_loc, r_edge_loc) &
       bind(C, name="init_base_state_irreg")

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )


  end subroutine init_base_state_irreg
end module base_state_module
