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

  use model_parser_module
  use amrex_constants_module
  use simple_log_module
  use meth_params_module, only: nscal, spherical
  use base_state_geometry_module, only: nr_fine,  max_radial_level

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

    if (spherical .eq. 1) then
       call amrex_error("ERROR: test_advect base_state is not valid for spherical")
    endif

    s0_init(0:max_radial_level,0:nr_fine-1,1:nscal) = 0.0d0
    p0_init(0:max_radial_level,0:nr_fine-1) = 0.0d0
    rho0(0:max_radial_level,0:nr_fine-1) = 0.0d0
    rhoh0(0:max_radial_level,0:nr_fine-1) = 0.0d0
    p0(0:max_radial_level,0:nr_fine-1) = 0.0d0
    tempbar(0:max_radial_level,0:nr_fine-1) = 0.0d0
    tempbar_init(0:max_radial_level,0:nr_fine-1) = 0.0d0

  end subroutine init_base_state


end module base_state_module
