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
  use amrex_fort_module, only: amrex_spacedim
  use network, only: nspec
  use meth_params_module, only: nscal, rho_comp
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level

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
    integer         :: i,j,n, r,comp
    double precision :: dist


    do n=0,max_radial_level

       do r=0,nr(n)-1

          ! height above the bottom of the domain
          dist = (dble(r) + HALF)*dr(n)

          s0_init(n, r, rho_comp) = exp(-dist**2/0.1d0)

       end do

    end do ! end loop over levels

    p0(:,:) = s0_init(:,:,rho_comp)

  end subroutine init_base_state

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

  end subroutine  init_base_state_irreg

end module base_state_module
