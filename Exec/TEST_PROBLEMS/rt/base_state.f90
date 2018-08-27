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
  use probdata_module, only: p0_base, rho_1, rho_2

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
    real(kind=dp_t) :: rloc,rmid
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)

    real(kind=dp_t) :: p0_light, p0_heavy, t_guess
    real(kind=dp_t) :: xn_light(nspec), xn_heavy(nspec)
    integer :: ia, ib

    real(kind=dp_t) :: min_dens, max_dens, min_temp, max_temp

    real(kind=dp_t), parameter :: TINY = 1.0d-10

    real(kind=dp_t) :: dpdr, rhog
    real(kind=dp_t) :: max_hse_error

    real(kind=dp_t), parameter :: SMALL = 1.d-12

    type (eos_t) :: eos_state

887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)

    if (spherical .eq. 1) then
       call bl_error("ERROR: rt base_state is not valid for spherical")
    endif

    if ( parallel_IOProcessor()) then
       ! output block for cutoff density information
       write (*,887)
       write (*,*)   'cutoff densities:'
       write (*,888) '    low density cutoff (for mapping the model) =      ', &
            base_cutoff_density
       write (*,888) '    buoyancy cutoff density                           '
       write (*,888) '        (for zeroing rho - rho_0, centrifugal term) = ', &
            buoyancy_cutoff_factor*base_cutoff_density
       write (*,888) '    anelastic cutoff =                                ', &
            anelastic_cutoff
       write (*,888) ' '
    end if

    min_dens = min(rho_1, rho_2)

    if (anelastic_cutoff > min_dens .or. base_cutoff_density > min_dens) then
       call bl_error("ERROR: for the RT problem, the anelastic and base cutoff densities > min(rho)")
    endif

    if (min_dens < small_dens) then
       if ( parallel_IOProcessor()) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than the EOS cutoff'
          print *, '         density, small_dens'
       endif
    endif

    if ( parallel_IOProcessor()) then
       ! close the cutoff density output block
       write (*,887)
       write (*,*)   ' '
    end if

    ! rmid is the middle of the domain
    rmid = HALF*(prob_lo(amrex_spacedim) + prob_hi(amrex_spacedim))

    ! p0_light is the pressure at the base of the light fluid (lower
    ! half of the domain)
    p0_light = p0_base

    ! p0heavy is the pressure at the base of the heavy fluid (top half
    ! of the domain).  To find this, we integrate dp/dr = rho g from the
    ! bottom of the domain (p = p0light) to the interface (y = rmid).  The
    ! density is constant in that region, rho = rho_1
    p0_heavy = p0_base + rho_1*grav_const*(rmid - prob_lo(amrex_spacedim))

    ! set the compositions
    ia = network_species_index("A")
    ib = network_species_index("B")

    xn_light(:) = SMALL
    xn_light(ia) = ONE - (nspec-1)*SMALL

    xn_heavy(:) = SMALL
    xn_heavy(ib) = ONE - (nspec-1)*SMALL

    ! set a guess for the temperature for the EOS calls
    t_guess = 1.e-8

    do n=0,max_radial_level

       do r=0,nr(n)-1

           ! height above the bottom of the domain
           rloc = (dble(r) + HALF)*dr(n)

           if (rloc > rmid) then

              ! top half -- heavy fluid
              d_ambient = rho_2
              p_ambient = p0_heavy + rho_2*grav_const*(rloc - rmid)
              t_ambient = t_guess
              xn_ambient(:) = xn_heavy(:)

           else

              ! lower half -- light fluid
              d_ambient = rho_1
              p_ambient = p0_light + rho_1*grav_const*(rloc - prob_lo(amrex_spacedim))
              t_ambient = t_guess
              xn_ambient(:) = xn_light(:)

           endif

           ! use the EOS to make the state consistent
           eos_state%T     = t_ambient
           eos_state%rho   = d_ambient
           eos_state%p     = p_ambient
           eos_state%xn(:) = xn_ambient(:)

           ! (rho,p) --> T, h
           call eos(eos_input_rp, eos_state)

           s0_init(n, r, rho_comp) = d_ambient
           s0_init(n, r,rhoh_comp) = d_ambient * eos_state%h
           s0_init(n, r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
           p0_init(n, r) = eos_state%p

           s0_init(n, r,temp_comp) = eos_state%T

       end do

       ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
       rho0 = s0_init(:,:,rho_comp)
       rhoh0 = s0_init(:,:,rhoh_comp)
       tempbar = s0_init(:,:,temp_comp)
       tempbar_init = s0_init(:,:,temp_comp)
       p0 = p0_init

       min_temp = minval(s0_init(n, :,temp_comp))

       if (min_temp < small_temp) then
          if ( parallel_IOProcessor() .and. n == 1) then
             print *, ' '
             print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
             print *, '         temperature, small_temp'
          endif
       endif

       max_hse_error = -1.d30

       do r=1,nr(n)-1

           rloc = prob_lo(amrex_spacedim) + (dble(r) + HALF)*dr(n)

           dpdr = (p0_init(n, r) - p0_init(n, r-1))/dr(n)
           rhog = HALF*(s0_init(n, r,rho_comp) + s0_init(n, r-1,rho_comp))*grav_const

           max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))

          !
          !    max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))
          !

       enddo

       if ( parallel_IOProcessor() ) then
          write (*,*) ' '
          write (*,*) 'Maximum HSE Error = ', max_hse_error
          write (*,*) '   (after putting initial model into base state arrays, and'
          write (*,*) '    for density < base_cutoff_density)'
          write (*,887)
          write (*,*) ' '
       endif

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
