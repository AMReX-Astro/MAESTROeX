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
  use eos_type_module
  use eos_module
  use amrex_constants_module
  use simple_log_module
  use inlet_bc_module
  use fundamental_constants_module, only: Gconst
  use amrex_fort_module, only: amrex_spacedim
  use network, only: nspec
  use meth_params_module, only: nscal, model_file, spherical, base_cutoff_density, &
       do_2d_planar_octant, do_planar_invsq_grav, rho_comp, &
       rhoh_comp, spec_comp, temp_comp, grav_const, &
       planar_invsq_mass, print_init_hse_diag, prob_lo, &
       prob_hi, small_dens, small_temp, &
       anelastic_cutoff_density, buoyancy_cutoff_factor
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use probin_module, only: rho_base, p_base

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
    integer         :: n,i,j,r,comp
    integer         :: ia, ib
    double precision :: rloc,rshr1,rshr2
    double precision :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    double precision :: t_guess
    double precision :: xn_jet(nspec), xn_still(nspec)
    double precision :: min_dens, max_dens, min_temp, max_temp
    double precision :: dpdr, rhog
    double precision :: max_hse_error
    double precision,     parameter :: SMALL   = 1.d-12
    double precision,     parameter :: TINY    = 1.0d-10
    character(len=*),    parameter :: FMT_SEP = "(78('-'))"
    character(len=*),    parameter :: FMT_MSG = "(a60,g18.10)"

    type (eos_t) :: eos_state

    if (spherical .eq. 1) then
       call amrex_error("ERROR: Incompressible shear jet base_state is not valid for spherical")
    endif

    if ( parallel_IOProcessor()) then
       ! Output block for cutoff densities
       write (*,FMT_SEP)
       write (*,*)   'cutoff densities:'
       write (*,FMT_MSG) '    low density cutoff (for mapping the model) =      ', &
            base_cutoff_density
       write (*,FMT_MSG) '    buoyancy cutoff density                           '
       write (*,FMT_MSG) '        (for zeroing rho - rho_0, centrifugal term) = ', &
            buoyancy_cutoff_factor*base_cutoff_density
       write (*,FMT_MSG) '    anelastic cutoff =                                ', &
            anelastic_cutoff_density
       write (*,FMT_MSG) ' '
    end if

    !Check min density
    min_dens = rho_base
    if (min_dens < small_dens) then
       if ( parallel_IOProcessor()) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than the EOS cutoff'
          print *, '         density, small_dens'
       endif
    endif

    if ( parallel_IOProcessor()) then
       ! Close the cutoff density output block
       write (*,FMT_SEP)
       write (*,*)   ' '
    end if

    !-- Initialize s0_init and p0_init --
    !Location of the two thin shear layers
    rshr1 = 0.25d0*(prob_lo(amrex_spacedim) + prob_hi(amrex_spacedim))
    rshr2 = 0.75d0*(prob_lo(amrex_spacedim) + prob_hi(amrex_spacedim))

    ! A and B act as tags.
    ! A tags the part of the fluid initially in the "jet zone," i.e. the part of
    ! the fluid with positive initial velocity.
    ! B tags the part of the fluid with negative initial velocity
    ia = network_species_index("A")
    ib = network_species_index("B")

    ! xn is an array where the ith element is the mass fraction (percent) of the
    ! ith species.  Initially, the middle half of the domain is composed entirely of
    ! "A" particles and negligible "B" particles, while the reverse is true of
    ! the outer regions of the domain.

    ! As we evolve in time these arrays track how the two fluids mix.
    xn_jet(:)      = SMALL
    xn_jet(ia)     = ONE - (nspec-1)*SMALL

    xn_still(:)  = SMALL
    xn_still(ib) = ONE - (nspec-1)*SMALL

    ! set a guess for the temperature of the EOS calls
    t_guess = 1.e-8

    !Iterate through each base state cell and initialize
    !   -The components of the fluid state 's':
    !       density, enthalpy, species mass fractions, and temperature
    !   -The pressure (note, pressure is NOT a component of the 's' multifab)
    do n=0,max_radial_level

       do r=0,nr(n)-1

          ! Current height above the bottom of the domain
          rloc = (dble(r) + HALF)*dr(n)

          !Init density, pressure, and temp
          d_ambient = rho_base
          p_ambient = p_base
          t_ambient = t_guess

          ! Depending on location, initialize the mass fraction
          if (rloc > rshr1 .and. rloc < rshr2) then
             ! Middle -- jet
             xn_ambient(:) = xn_jet(:)
          else
             ! Outer regions -- still fluid
             xn_ambient(:) = xn_still(:)
          endif

          ! use the EOS to make the state consistent
          ! We set density and pressure, and from this the EoS yields many
          ! thermodynamic quantities (temperature and enthalpy being the two we
          ! care about in this problem).
          eos_state%T     = t_ambient
          eos_state%rho   = d_ambient
          eos_state%p     = p_ambient
          eos_state%xn(:) = xn_ambient(:)

          ! (rho,p) --> T, h
          call eos(eos_input_rp, eos_state)

          !Now that we've calculated all of the ambient values and churned them
          !through the EoS we can finally initialize the fluid state.
          s0_init(n, r, rho_comp)                    = d_ambient
          s0_init(n, r, rhoh_comp)                   = d_ambient * eos_state%h
          s0_init(n, r, spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
          s0_init(n, r, temp_comp)                   = eos_state%T
          p0_init(n, r)                              = p_base

       end do

       !-- Post State Initialization --
       !Check that the temperature is consistent with the EoS
       min_temp = minval(s0_init(n, :,temp_comp))
       if (min_temp < small_temp) then
          if ( parallel_IOProcessor() .and. n == 1) then
             print *, ' '
             print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
             print *, '         temperature, small_temp'
          endif
       endif

       ! initialize any inlet BC parameters
       call set_inlet_bcs()

    end do ! end loop over levels

    ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    rho0 = s0_init(:,:,rho_comp)
    rhoh0 = s0_init(:,:,rhoh_comp)
    tempbar = s0_init(:,:,temp_comp)
    tempbar_init = s0_init(:,:,temp_comp)
    p0 = p0_init

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
