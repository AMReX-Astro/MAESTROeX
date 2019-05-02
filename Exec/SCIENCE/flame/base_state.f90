module base_state_module

  use model_parser_module
  use eos_type_module
  use eos_module
  use amrex_constants_module
  use amrex_error_module
  use simple_log_module
  use inlet_bc_module
  use amrex_fort_module, only: amrex_spacedim
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use network, only: nspec, network_species_index
  use meth_params_module, only: nscal, base_cutoff_density, &
       rho_comp, rhoh_comp, spec_comp, temp_comp, grav_const, &
       prob_lo, prob_hi, small_dens, small_temp, &
       anelastic_cutoff_density, buoyancy_cutoff_factor
  use probin_module, only: dens_fuel, temp_fuel, xc12_fuel, vel_fuel, &
       interface_pos_frac, smooth_len_frac, &
       temp_ash, temp_fuel

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
    integer         :: i,j,n,r, ic12, io16, img24
    double precision :: p_ambient, dens_ash, rhoh_fuel, rhoh_ash
    double precision :: xn_fuel(nspec), xn_ash(nspec), xn_smooth(nspec)
    double precision :: rlen, rloc

    type (eos_t) :: eos_state

887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)

    ! sanity check
    if (dens_fuel < base_cutoff_density .or.  dens_fuel < anelastic_cutoff_density) then
       call amrex_error('ERROR: fuel density < (base_cutoff_density or anelastic_cutoff_density)')
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
            anelastic_cutoff_density
       write (*,888) ' '
    end if

    if ( parallel_IOProcessor()) then
       ! close the cutoff density output block
       write (*,887)
       write (*,*)   ' '
    end if

    ! figure out the indices for different species
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")
    img24  = network_species_index("magnesium-24")

    if (ic12 < 0 .or. io16 < 0 .or. img24 < 0) then
       call amrex_error("ERROR: species indices not defined")
    end if

    ! length of the domain
    rlen = (prob_hi(amrex_spacedim) - prob_lo(amrex_spacedim))

    ! figure out the thermodynamics of the fuel and ash state

    ! fuel
    xn_fuel(:)    = ZERO
    xn_fuel(ic12) = xc12_fuel
    xn_fuel(io16) = 1.d0 - xc12_fuel

    eos_state%rho   = dens_fuel
    eos_state%T     = temp_fuel
    eos_state%xn(:) = xn_fuel(:)

    call eos(eos_input_rt, eos_state)

    ! note: p_ambient should be = p0_init
    p_ambient = eos_state%p
    rhoh_fuel = dens_fuel*eos_state%h

    ! ash
    xn_ash(:)     = ZERO
    xn_ash(ic12)  = ZERO
    xn_ash(io16)  = 1.d0 - xc12_fuel
    xn_ash(img24) = xc12_fuel

    eos_state%rho   = dens_fuel    ! initial guess
    eos_state%T     = temp_ash
    eos_state%xn(:) = xn_ash(:)
    eos_state%p     = p_ambient

    call eos(eos_input_tp, eos_state)

    dens_ash = eos_state%rho
    rhoh_ash = dens_ash*eos_state%h

    do n=0,max_radial_level

       do j=0,nr(n)-1

          rloc = prob_lo(amrex_spacedim) + (dble(j) + 0.5d0) * dr(n)

          ! the flame propagates in the -y direction.  The fuel/ash division
          ! is interface_pos_frac through the domain
          if (rloc < prob_lo(amrex_spacedim) + interface_pos_frac*rlen) then

             ! fuel
             s0_init(n,j,rho_comp)  = dens_fuel
             s0_init(n,j,rhoh_comp) = rhoh_fuel
             !s0_init(r,temp_comp) = temp_fuel
             s0_init(n,j,spec_comp:spec_comp+nspec-1) = dens_fuel*xn_fuel(:)

          else

             ! ash
             s0_init(n,j,rho_comp)  = dens_ash
             s0_init(n,j,rhoh_comp) = rhoh_ash
             !s0_init(r,temp_comp) = temp_ash
             s0_init(n,j,spec_comp:spec_comp+nspec-1) = dens_ash*xn_ash(:)

          endif

          ! give the temperature a smooth profile
          s0_init(n,j,temp_comp) = temp_fuel + (temp_ash - temp_fuel) * &
               HALF * (ONE + &
               tanh( (rloc - (prob_lo(amrex_spacedim) + interface_pos_frac*rlen)) / &
               (smooth_len_frac*rlen) ) )

          ! give the carbon mass fraction a smooth profile too
          xn_smooth(:) = ZERO
          xn_smooth(ic12) = xn_fuel(ic12) + (xn_ash(ic12) - xn_fuel(ic12)) * &
               HALF * (ONE + &
               tanh( (rloc - (prob_lo(amrex_spacedim) + interface_pos_frac*rlen)) / &
               (smooth_len_frac*rlen) ) )

          xn_smooth(io16) = xn_fuel(io16)
          xn_smooth(img24) = 1.d0 - xn_smooth(ic12) - xn_smooth(io16)

          ! get the new density and enthalpy
          eos_state%rho   = s0_init(n,j,rho_comp)
          eos_state%T     = s0_init(n,j,temp_comp)
          eos_state%xn(:) = xn_smooth(:)
          eos_state%p     = p_ambient

          call eos(eos_input_tp, eos_state)

          s0_init(n,j,rho_comp)  = eos_state%rho
          s0_init(n,j,spec_comp:spec_comp+nspec-1) = eos_state%rho*xn_smooth(:)
          s0_init(n,j,rhoh_comp) = eos_state%rho * eos_state%h
       enddo

       ! set the base state pressure to be the ambient pressure.
       p0_init(:,:) = p_ambient

    end do ! end loop over levels

    ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    rho0 = s0_init(:,:,rho_comp)
    rhoh0 = s0_init(:,:,rhoh_comp)
    tempbar = s0_init(:,:,temp_comp)
    tempbar_init = s0_init(:,:,temp_comp)
    p0 = p0_init

    call set_inlet_bcs()

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
