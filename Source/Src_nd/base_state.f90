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
  use network, only: nspec
  use meth_params_module, only: nscal, model_file
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(s0_init,p0_init,nlevs,nr_fine) bind(C, name="init_base_state")

    integer         , intent(in   ) :: nlevs, nr_fine
    double precision, intent(in   ) :: s0_init(1:nlevs,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(1:nlevs,0:nr_fine-1)

    ! local variables
    double precision :: base_cutoff_density_loc
    double precision :: rho_above_cutoff, rhoh_above_cutoff
    double precision :: spec_above_cutoff(nspec), p_above_cutoff
    double precision :: temp_above_cutoff

    call read_model_file(model_file)

  end subroutine init_base_state

end module base_state_module
