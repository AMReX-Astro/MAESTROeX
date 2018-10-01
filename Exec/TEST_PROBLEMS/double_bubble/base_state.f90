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
       anelastic_cutoff, buoyancy_cutoff_factor
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use probin_module, only: dens_base, pres_base, do_isentropic

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
    integer         :: n,r,j
    real(kind=dp_t) :: H,z,z0, gamma_const
    real(kind=dp_t) :: dens_zone, temp_zone
    real(kind=dp_t) :: xn_zone(nspec)

    type (eos_t) :: eos_state

    if (spherical .eq. 1) then
       call bl_error("ERROR: rt base_state is not valid for spherical")
    endif

    ! only initialize the first species
    xn_zone(:) = ZERO
    xn_zone(1) = 1.d0

    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(grav_const)

    do n=0,max_radial_level

       ! for isentropic, we satisfy p ~ rho^gamma, but we'll need to get gamma
       eos_state % rho = dens_base
       eos_state % p = pres_base
       eos_state % xn(:) = xn_zone(:)

       ! initial guess
       eos_state % T = 1000.0d0

       call eos(eos_input_rp, eos_state)

       gamma_const = pres_base/(dens_base * eos_state % e) + 1.0d0

       p0_init(n,0) = pres_base
       s0_init(n,0, rho_comp) = dens_base

       s0_init(n,0,rhoh_comp) = dens_base * eos_state%h

       s0_init(n,0,spec_comp:spec_comp-1+nspec) = dens_base*xn_zone(:)

       s0_init(n,0,temp_comp) = eos_state%T

       z0 = 0.5d0*dr(n)

       ! set an initial guess for the temperature -- this will be reset
       ! by the EOS
       temp_zone = 1000.d0

       do r=1,nr(n)-1

          z = (dble(r)+HALF) * dr(n)

          if (do_isentropic) then

             ! we can integrate HSE with p = K rho^gamma analytically
             dens_zone = dens_base*(grav_const*dens_base*(gamma_const - 1.0)* &
                  (z-z0)/ &
                  (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))

          else

             ! the density of an isothermal gamma-law atm is exponential
             dens_zone = dens_base * exp(-z/H)

          end if

          s0_init(n,r, rho_comp) = dens_zone

          ! compute the pressure by discretizing HSE
          p0_init(n,r) = p0_init(n,r-1) - &
               dr(n) * HALF * (s0_init(n,r,rho_comp) + s0_init(n,r-1,rho_comp)) * &
               abs(grav_const)

          ! use the EOS to make the state consistent
          eos_state%rho   = dens_zone
          eos_state%p     = p0_init(n,r)
          eos_state%T     = temp_zone
          eos_state%xn(:) = xn_zone(:)

          ! (rho,p) --> T, h
          call eos(eos_input_rp, eos_state)

          s0_init(n,r, rho_comp) = dens_zone
          s0_init(n,r,rhoh_comp) = dens_zone * eos_state%h

          s0_init(n,r,spec_comp:spec_comp-1+nspec) = dens_zone*xn_zone(:)

          s0_init(n,r,temp_comp) = eos_state%T

       end do

       ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
       rho0 = s0_init(:,:,rho_comp)
       rhoh0 = s0_init(:,:,rhoh_comp)
       tempbar = s0_init(:,:,temp_comp)
       tempbar_init = s0_init(:,:,temp_comp)
       p0 = p0_init

       ! initialize any inlet BC parameters
       call set_inlet_bcs()

    end do ! end loop over levels

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
