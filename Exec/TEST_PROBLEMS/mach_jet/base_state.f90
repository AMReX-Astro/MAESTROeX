module base_state_module

  use model_parser_module
  use eos_type_module
  use eos_module, only: eos, eos_input_rp
  use amrex_constants_module
  use amrex_error_module
  use simple_log_module
  use inlet_bc_module
  use amrex_fort_module, only: amrex_spacedim
  use inlet_bc_module
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use network, only: nspec
  use meth_params_module, only: nscal, spherical, base_cutoff_density, &
       rho_comp, rhoh_comp, spec_comp, temp_comp, grav_const, &
       prob_lo, prob_hi, small_dens, small_temp, &
       anelastic_cutoff, buoyancy_cutoff_factor
  use probin_module, only:  do_stratified, do_isentropic

  implicit none

  private

contains

  subroutine init_base_state(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init) &
       bind(C, name="init_base_state")

    ! use bl_prof_module
    ! use parallel
    ! use bl_error_module
    ! use bl_constants_module
    ! use eos_module, only: eos, eos_input_rp
    ! use eos_type_module
    ! use network, only: spec_names, network_species_index
    ! use probin_module, only: grav_const, do_stratified, do_isentropic
    ! use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    ! use geometry, only: dr, spherical, nr
    ! use inlet_bc_module
    ! use extern_probin_module, only : eos_gamma

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)

    ! local variables
    integer         :: i,j,n,r,comp
    double precision :: z, H, eos_gamma
    double precision :: dens_zone, temp_zone, pres_zone
    double precision :: xn_zone(nspec)

    type (eos_t) :: eos_state

887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)

    if (spherical .eq. 1) then
       call amrex_error("ERROR: mach_jet base_state is not valid for spherical")
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

    if ( parallel_IOProcessor()) then
       ! close the cutoff density output block
       write (*,887)
       write (*,*)   ' '
    end if

    do n=0,max_radial_level

       if (do_stratified) then

          ! use the EOS to make the state consistent
          dens_zone = 1.d-3
          pres_zone = 1.d6

          ! only initialize the first species
          xn_zone(:) = ZERO
          xn_zone(1) = 1.d0

          p0_init(0) = pres_zone

          ! H = pres_base / dens_base / abs(grav_const)
          H = 1.d6 / 1.d-3 / abs(grav_const)

          ! set an initial guess for the temperature -- this will be reset
          ! by the EOS
          temp_zone = 10.d0

          eos_gamma = eos_state%gam1

          do j=0,nr(n)-1

             z = (dble(j)+HALF) * dr(n)

             if (do_isentropic) then
                dens_zone = 1.d-3*(grav_const*1.d-3*(eos_gamma - 1.0)*z/ &
                     (eos_gamma*1.d6) + 1.d0)**(1.d0/(eos_gamma - 1.d0))
             else
                dens_zone = 1.d-3*exp(-z/H)
             end if

             s0_init(j, rho_comp) = dens_zone

             if (j.eq.0) then
                p0_init(j) = p0_init(j) - &
                     dr(1) * HALF * s0_init(j,rho_comp) * &
                     abs(grav_const)
             else if (j.gt.0) then
                p0_init(j) = p0_init(j-1) - &
                     dr(1) * HALF * (s0_init(j,rho_comp) + s0_init(j-1,rho_comp)) * &
                     abs(grav_const)
             end if

             pres_zone = p0_init(j)

             ! use the EOS to make the state consistent
             eos_state%rho   = dens_zone
             eos_state%T     = temp_zone
             eos_state%p     = pres_zone
             eos_state%xn(:) = xn_zone(:)

             ! (rho,p) --> T, h
             call eos(eos_input_rp, eos_state)

             s0_init(j, rho_comp) = dens_zone
             s0_init(j,rhoh_comp) = dens_zone*eos_state%h

             s0_init(j,spec_comp:spec_comp-1+nspec) = ZERO
             s0_init(j,spec_comp) = dens_zone

             s0_init(j,temp_comp) = eos_state%T
             s0_init(j,trac_comp) = ZERO

          end do

       else

          ! use the EOS to make the state consistent
          eos_state%T    = 10.d0
          eos_state%rho   = 1.d-3
          eos_state%p     = 1.d6
          eos_state%xn(:) = 1.d0

          ! (rho,p) --> T, h
          call eos(eos_input_rp, eos_state)

          s0_init(0:nr(n)-1, rho_comp) = eos_state%rho
          s0_init(0:nr(n)-1,rhoh_comp) = eos_state%rho * eos_state%h
          s0_init(0:nr(n)-1,spec_comp) = eos_state%rho
          s0_init(0:nr(n)-1,temp_comp) = eos_state%T
          s0_init(0:nr(n)-1,trac_comp) = ZERO

          p0_init(0:nr(n)-1) = eos_state%p

       end if

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

  end subroutine init_base_state_irreg

end module base_state_module
