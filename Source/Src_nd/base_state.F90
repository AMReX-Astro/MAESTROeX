module base_state_module
  ! init_base_state is used to initialize the base state arrays from the
  ! model file.  The actual reading of the model file is handled by the
  ! model_parser_module in Util/
  !
  ! Note: The initial base state quantities returned from this routine
  ! are only a temporary base state.  These quantities are mapped onto
  ! the full 2- or 3-d state in initscaldata.f90 and a new base state is
  ! created after initialization by averaging the density and calling
  ! enforce_HSE in initialize.f90.

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
                                planar_invsq_mass, print_init_hse_diag, prob_lo
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level

  implicit none

  private

contains

  subroutine init_base_state(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init) &
       bind(C, name="init_base_state")
    ! Binds to C function ``init_base_state``

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)

    ! local variables
    double precision :: base_cutoff_density_loc
    double precision :: rho_above_cutoff, rhoh_above_cutoff
    double precision :: spec_above_cutoff(nspec), p_above_cutoff
    double precision :: temp_above_cutoff

    integer         :: n,r,comp
    double precision :: rloc,model_dr,rmax,starting_rad,mod_dr
    double precision :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    double precision :: sumX

    double precision, parameter :: TINY = 1.0d-10

    double precision :: mencl, g, r_l, r_r, dpdr, rhog
    double precision :: max_hse_error

    type (eos_t) :: eos_state

887 format(78('-'))

    base_cutoff_density_loc = 1.d99

    ! only need to read in the initial model once --
    ! model_parser_module stores the model data
    call read_model_file(model_file)

    model_dr = (model_r(npts_model) - model_r(1)) / dble(npts_model-1)
    rmax = model_r(npts_model)

    do n=0,max_radial_level

       if ( parallel_IOProcessor() ) then
          write (*,887)
          if (spherical .ne. 1) then
             call log('model file mapping, level: ', n)
          else
             call log('model file mapping (spherical base state)')
          endif

          call log('dr of MAESTRO base state =                            ', dr(n))
          call log('dr of input file data =                               ', model_dr)
          call log(' ')
          call log('maximum radius (cell-centered) of input model =       ', rmax)

          if (dr(n) .lt. model_dr) then
             mod_dr = mod(model_dr,dr(n))
          else
             mod_dr = mod(dr(n),model_dr)
          endif

          if (mod_dr .gt. TINY) then
             call log(' ')
             call log("WARNING: resolution of base state array is not an integer")
             call log("         multiple of the initial model's resolution.     ")
             call log("         make sure this is a desired property as this    ")
             call log("         could lead to aliasing when performing the      ")
             call log("         interpolation.                                  ")
             call log(" ")
             call log("modulus = ", mod_dr)
          endif
       end if

       if (spherical .eq. 0) then
          starting_rad = prob_lo(amrex_spacedim)
       else
          starting_rad = ZERO
       endif

       do r=0,nr(n)-1

          rloc = starting_rad + (dble(r) + HALF)*dr(n)

          ! here we account for r > rmax of the model.hse array, assuming
          ! that the state stays constant beyond rmax
          rloc = min(rloc, rmax)

          ! also, if we've falled below the cutoff density, just keep the
          ! model constant
          if (rloc .gt. base_cutoff_density_loc) then

             s0_init(n,r,rho_comp) = rho_above_cutoff
             s0_init(n,r,rhoh_comp) = rhoh_above_cutoff
             s0_init(n,r,spec_comp:spec_comp+nspec-1) = spec_above_cutoff(1:nspec)
             p0_init(n,r) = p_above_cutoff
             s0_init(n,r,temp_comp) = temp_above_cutoff

          else

             d_ambient = interpolate(rloc, idens_model)
             t_ambient = interpolate(rloc, itemp_model)
             p_ambient = interpolate(rloc, ipres_model)

             sumX = ZERO
             do comp = 1, nspec
                xn_ambient(comp) = max(ZERO, min(ONE, &
                     interpolate(rloc, ispec_model-1+comp)))
                sumX = sumX + xn_ambient(comp)
             enddo
             xn_ambient = xn_ambient/sumX

             ! use the EOS to make the state consistent
             eos_state%T     = t_ambient
             eos_state%rho   = d_ambient
             eos_state%p     = p_ambient
             eos_state%xn(:) = xn_ambient(:)

             ! (rho,T) --> p,h
             call eos(eos_input_rt, eos_state)

             s0_init(n,r, rho_comp ) = d_ambient
             s0_init(n,r,rhoh_comp ) = d_ambient * eos_state%h
             s0_init(n,r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
             p0_init(n,r) = eos_state%p ! p_ambient !
             s0_init(n,r,temp_comp) = t_ambient

             ! keep track of the height where we drop below the cutoff density
             if (s0_init(n,r,rho_comp) .le. base_cutoff_density .and. &
                  base_cutoff_density_loc .eq. 1.d99) then

                if ( parallel_IOProcessor() ) then
                   write (*,*) ' '
                   write (*,*) 'setting r_cutoff to ', r
                   write (*,*) 'radius at r_cutoff ', rloc
                end if

                base_cutoff_density_loc = rloc

                rho_above_cutoff           = s0_init(n,r,rho_comp)
                rhoh_above_cutoff          = s0_init(n,r,rhoh_comp)
                spec_above_cutoff(1:nspec) = s0_init(n,r,spec_comp:spec_comp+nspec-1)
                temp_above_cutoff          = s0_init(n,r,temp_comp)
                p_above_cutoff             = p0_init(n,r)

             end if

          end if

       end do

       ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
       rho0 = s0_init(:,:,rho_comp)
       rhoh0 = s0_init(:,:,rhoh_comp)
       tempbar = s0_init(:,:,temp_comp)
       tempbar_init = s0_init(:,:,temp_comp)
       p0 = p0_init

       ! check whether we are in HSE

       mencl = zero

       if (spherical .eq. 1 .OR. do_2d_planar_octant .eq. 1) then
          mencl = four3rd*m_pi*dr(n)**3*s0_init(n,0,rho_comp)
       endif

       max_hse_error = -1.d30

       do r=1,nr(n)-1

          rloc = starting_rad + (dble(r) + HALF)*dr(n)
          rloc = min(rloc, rmax)

          if (rloc .lt. base_cutoff_density_loc) then

             r_r = starting_rad + dble(r+1)*dr(n)
             r_l = starting_rad + dble(r)*dr(n)

             if (spherical .eq. 1 .OR. do_2d_planar_octant .eq. 1) then
                g = -Gconst*mencl/r_l**2
                mencl = mencl &
                     + four3rd*m_pi*dr(n)*(r_l**2+r_l*r_r+r_r**2)*s0_init(n,r,rho_comp)
             else
                if (.not. do_planar_invsq_grav) then
                   g = grav_const
                else
                   g = -Gconst*planar_invsq_mass / r_l**2
                endif
             endif

             dpdr = (p0_init(n,r) - p0_init(n,r-1))/dr(n)
             rhog = HALF*(s0_init(n,r,rho_comp) + s0_init(n,r-1,rho_comp))*g

             if (print_init_hse_diag) then
                if ( parallel_IOProcessor() ) then
                   print *, 'r, dpdr, rhog, err: ', rloc, dpdr, rhog, &
                        abs(dpdr - rhog)/abs(rhog)
                endif
             endif

             max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))

          end if

       enddo

       if ( parallel_IOProcessor() ) then
          call log(' ')
          call log('Maximum HSE Error = ', max_hse_error)
          call log('   (after putting initial model into base state arrays, and')
          call log('    for density < base_cutoff_density)')
          write (*,887)
          call log(' ')
       endif

       ! initialize any inlet BC parameters
       call set_inlet_bcs()

    end do ! end loop over levels

  end subroutine init_base_state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_base_state_irreg(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init, &
                                     r_cc_loc, r_edge_loc) &
       bind(C, name="init_base_state_irreg")
    ! Binds to C function ``init_base_state_irreg``

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

    ! local variables
    double precision :: base_cutoff_density_loc
    double precision :: rho_above_cutoff, rhoh_above_cutoff
    double precision :: spec_above_cutoff(nspec), p_above_cutoff
    double precision :: temp_above_cutoff

    integer         :: n,r,comp,i
    double precision :: rloc,rmax,starting_rad,mod_dr
    double precision :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    double precision :: sumX
    double precision :: dr_irreg, rfrac

    double precision, allocatable :: model_dr(:)

    double precision, parameter :: TINY = 1.0d-10

    double precision :: mencl, g, r_l, r_r, dpdr, rhog
    double precision :: max_hse_error

    type (eos_t) :: eos_state

887 format(78('-'))

    base_cutoff_density_loc = 1.d99

    ! only need to read in the initial model once --
    ! model_parser_module stores the model data
    call read_model_file(model_file)

    allocate (model_dr(0:npts_model-1))

    model_dr(0) = model_r(1)
    do i=1,npts_model-1
       model_dr(i) = (model_r(i+1) - model_r(i))
    end do
    rmax = model_r(npts_model)

    do n=0,max_radial_level

       if ( parallel_IOProcessor() ) then
          write (*,887)
          if (spherical .ne. 1) then
             call log('model file mapping, level: ', n)
          else
             call log('model file mapping (spherical base state)')
          endif

          call log('coarsest dr of MAESTRO base state =                   ', r_cc_loc(n,0))
          call log('coarsest dr of input file data =                      ', model_dr(0))
          call log(' ')
          call log('maximum radius (cell-centered) of input model =       ', rmax)

          if (r_cc_loc(n,0) .lt. model_dr(0)) then
             mod_dr = mod(model_dr(0),r_cc_loc(n,0))
          else
             mod_dr = mod(r_cc_loc(n,0),model_dr(0))
          endif

          if (mod_dr .gt. TINY) then
             call log(' ')
             call log("WARNING: resolution of base state array is not an integer")
             call log("         multiple of the initial model's resolution.     ")
             call log("         make sure this is a desired property as this    ")
             call log("         could lead to aliasing when performing the      ")
             call log("         interpolation.                                  ")
             call log(" ")
             call log("modulus = ", mod_dr)
          endif
       end if

       if (spherical .eq. 0) then
          call amrex_error('init_base_state_irreg() not valid for planar')
       else
          starting_rad = ZERO
       endif

       do r=0,nr(n)-1

          rloc = starting_rad + r_cc_loc(n,r)

          ! here we account for r > rmax of the model.hse array, assuming
          ! that the state stays constant beyond rmax
          rloc = min(rloc, rmax)

          ! also, if we've falled below the cutoff density, just keep the
          ! model constant
          if (rloc .gt. base_cutoff_density_loc) then

             s0_init(n,r,rho_comp) = rho_above_cutoff
             s0_init(n,r,rhoh_comp) = rhoh_above_cutoff
             s0_init(n,r,spec_comp:spec_comp+nspec-1) = spec_above_cutoff(1:nspec)
             p0_init(n,r) = p_above_cutoff
             s0_init(n,r,temp_comp) = temp_above_cutoff

          else

             d_ambient = interpolate(rloc, idens_model)
             t_ambient = interpolate(rloc, itemp_model)
             p_ambient = interpolate(rloc, ipres_model)

             sumX = ZERO
             do comp = 1, nspec
                xn_ambient(comp) = max(ZERO, min(ONE, &
                     interpolate(rloc, ispec_model-1+comp)))
                sumX = sumX + xn_ambient(comp)
             enddo
             xn_ambient = xn_ambient/sumX

             ! use the EOS to make the state consistent
             eos_state%T     = t_ambient
             eos_state%rho   = d_ambient
             eos_state%p     = p_ambient
             eos_state%xn(:) = xn_ambient(:)

             ! (rho,T) --> p,h
             call eos(eos_input_rt, eos_state)

             s0_init(n,r, rho_comp ) = d_ambient
             s0_init(n,r,rhoh_comp ) = d_ambient * eos_state%h
             s0_init(n,r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
             p0_init(n,r) = eos_state%p ! p_ambient !
             s0_init(n,r,temp_comp) = t_ambient

             ! keep track of the height where we drop below the cutoff density
             if (s0_init(n,r,rho_comp) .le. base_cutoff_density .and. &
                  base_cutoff_density_loc .eq. 1.d99) then

                if ( parallel_IOProcessor() ) then
                   write (*,*) ' '
                   write (*,*) 'setting r_cutoff to ', r
                   write (*,*) 'radius at r_cutoff ', rloc
                end if

                base_cutoff_density_loc = rloc

                rho_above_cutoff           = s0_init(n,r,rho_comp)
                rhoh_above_cutoff          = s0_init(n,r,rhoh_comp)
                spec_above_cutoff(1:nspec) = s0_init(n,r,spec_comp:spec_comp+nspec-1)
                temp_above_cutoff          = s0_init(n,r,temp_comp)
                p_above_cutoff             = p0_init(n,r)

             end if

          end if

       end do

       ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
       rho0 = s0_init(:,:,rho_comp)
       rhoh0 = s0_init(:,:,rhoh_comp)
       tempbar = s0_init(:,:,temp_comp)
       tempbar_init = s0_init(:,:,temp_comp)
       p0 = p0_init

       ! check whether we are in HSE

       mencl = zero

       dr_irreg = r_edge_loc(n,1) - r_edge_loc(n,0)  ! edge-to-edge

       if (spherical .eq. 1 .OR. do_2d_planar_octant .eq. 1) then
          mencl = four3rd*m_pi*dr_irreg**3*s0_init(n,0,rho_comp)
       endif

       max_hse_error = -1.d30

       do r=1,nr(n)-1

          rloc = starting_rad + r_cc_loc(n,r)
          rloc = min(rloc, rmax)

          if (rloc .lt. base_cutoff_density_loc) then

             r_r = starting_rad + r_edge_loc(n,r+1)
             r_l = starting_rad + r_edge_loc(n,r)

             if (spherical .eq. 1 .OR. do_2d_planar_octant .eq. 1) then
                g = -Gconst*mencl/r_l**2
                mencl = mencl &
                     + four3rd*m_pi*(r_r-r_l)*(r_l**2+r_l*r_r+r_r**2)*s0_init(n,r,rho_comp)
             else
                if (.not. do_planar_invsq_grav) then
                   g = grav_const
                else
                   g = -Gconst*planar_invsq_mass / r_l**2
                endif
             endif

             ! dpdr and rhog are computed consistently with model file generator
             dr_irreg = (r_cc_loc(n,r)-r_cc_loc(n,r-1))
             dpdr = (p0_init(n,r) - p0_init(n,r-1))/dr_irreg

             rfrac = (r_edge_loc(n,r) - r_cc_loc(n,r-1))/dr_irreg
             rhog = ((1.d0-rfrac)*s0_init(n,r,rho_comp) + rfrac*s0_init(n,r-1,rho_comp))*g

             if (print_init_hse_diag) then
                if ( parallel_IOProcessor() ) then
                   print *, 'r, dpdr, rhog, err: ', rloc, dpdr, rhog, &
                        abs(dpdr - rhog)/abs(rhog)
                endif
             endif

             max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))

          end if

       enddo

       if ( parallel_IOProcessor() ) then
          call log(' ')
          call log('Maximum HSE Error = ', max_hse_error)
          call log('   (after putting initial model into base state arrays, and')
          call log('    for density < base_cutoff_density)')
          write (*,887)
          call log(' ')
       endif

       ! initialize any inlet BC parameters
       call set_inlet_bcs()

    end do ! end loop over levels

  end subroutine init_base_state_irreg
end module base_state_module
