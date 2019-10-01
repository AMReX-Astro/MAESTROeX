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
  use probin_module, only: dens_base, pres_base, do_isentropic, do_uniform, &
        use_p_dens_g, use_p_H_g, scale_height

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
    double precision :: H,z,z0, gamma_const
    double precision :: dens_zone
    double precision :: xn_zone(nspec)

    type (eos_t) :: eos_state

    if (spherical .eq. 1) then
       call amrex_error("ERROR: base_state is not valid for spherical")
    endif

    ! strictly single species
    xn_zone(:) = ZERO
    xn_zone(1) = 1.d0

    ! calculate H or dens_base, as necessary
    if (use_p_dens_g .and. use_p_H_g) then
      call amrex_error("ERROR: use_p_dens_g AND use_p_H_g cannot both be true")
    elseif (use_p_dens_g) then
      H = pres_base / dens_base / abs(grav_const)
    elseif (use_p_H_g) then
      H = scale_height
      dens_base = pres_base / H / abs(grav_const)
    else
      call amrex_error("ERROR: either use_p_dens_g or use_p_H_g must be true")
    endif


    do n=0,max_radial_level

       ! Set the state in the first cell in the vertical direction
       ! ("the bottom")
       z0 = 0.5d0*dr(n) ! generic "z", actually 2d/3d agnostic

       ! set p0 and rho0 based on analytical value at the cc coord
       p0_init(n,0) = pres_base * exp(-z0/H)
       s0_init(n,0, rho_comp) = dens_base * exp(-z0/H)

       ! use eos call to be consistent
       eos_state % rho = p0_init(n,0)
       eos_state % p = s0_init(n,0, rho_comp)
       eos_state % xn(:) = xn_zone(:)
       eos_state % T = 1000.0d0 ! just a guess... needed? no guess on e.g. h
       call eos(eos_input_rp, eos_state) ! (rho, p) --> T, h

       ! set other base vars
       s0_init(n,0,rhoh_comp) = dens_base * eos_state%h
       s0_init(n,0,spec_comp:spec_comp-1+nspec) = dens_base*xn_zone(:)
       s0_init(n,0,temp_comp) = eos_state%T

      
      do r=1,nr(n)-1

          z = (dble(r)+HALF) * dr(n)

          ! set rho analytically          
          dens_zone = dens_base * exp(-z/H)
          s0_init(n,r, rho_comp) = dens_zone ! needs to be set before pressure 

          ! compute the pressure by discretizing HSE
          p0_init(n,r) = p0_init(n,r-1) - &
               dr(n) * HALF * (s0_init(n,r,rho_comp) + s0_init(n,r-1,rho_comp)) * &
               abs(grav_const)

          ! use the EOS to make the state consistent
          eos_state%rho   = dens_zone
          eos_state%p     = p0_init(n,r)
          eos_state%T     = 1000.d0 ! just a  guess
          eos_state%xn(:) = xn_zone(:)
          call eos(eos_input_rp, eos_state) ! (rho,p) --> T, h

          ! save the state
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
