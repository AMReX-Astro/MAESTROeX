module test_basestate_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use probin_module, ONLY : heating_time, heating_rad, heating_peak, &
       heating_sigma, prob_type, &
       cooling_rad, cooling_peak, cooling_sigma
  use network
  use eos_module
  use eos_type_module
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use meth_params_module, only: spherical

  implicit none

  public get_heating

contains

  subroutine get_heating(Hbar,rho0,tempbar,rhoX0,time,dt,r_cc_loc) bind(C, name="get_heating")

    double precision, intent(inout) :: Hbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rhoX0(0:max_radial_level,0:nr_fine-1,1:nspec)
    double precision, value, intent(in   ) :: time,dt
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)

    double precision :: fac
    integer :: r, n

    integer         :: h1_comp
    integer         :: he4_comp
    integer         :: c12_comp
    integer         :: n14_comp
    integer         :: o16_comp
    double precision :: rho, T_6_third, X_CNO, X_1, g14
    double precision :: tmp1, tmp2, tmp3

    Hbar(:,:) = 0.d0

    if (prob_type .eq. 1) then

       if (time .le. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do n=0,max_radial_level
             do r = 0, nr(1)-1
                if (spherical .eq. 0) then
                   ! plane-parallel -- do the heating term in paper II (section 4)
                   Hbar(n,r) = fac * heating_peak * &
                        exp(-((r_cc_loc(n,r) - heating_rad)**2)/ heating_sigma)
                else
                   ! spherical -- lower amplitude heating term
                   Hbar(n,r) = fac * heating_peak * &
                        exp(-((r_cc_loc(n,r) - heating_rad)**2)/ heating_sigma)
                endif
             enddo
          enddo
       end if

    elseif (prob_type .eq. 2) then

       ! analytic heating modeling CNO cycle

       h1_comp = network_species_index("hydrogen-1")
       c12_comp = network_species_index("carbon-12")
       n14_comp = network_species_index("nitrogen-14")
       o16_comp = network_species_index("oxygen-16")

       do n=0,max_radial_level
          do r = 0, nr(1)-1
             rho = rho0(n,r)
             T_6_third = (tempbar(n,r) / 1.0d6) ** THIRD
             tmp1 = rhoX0(n,r,c12_comp)
             tmp2 = rhoX0(n,r,n14_comp)
             tmp3 = rhoX0(n,r,o16_comp)
             X_CNO = (tmp1 + tmp2 + tmp3) / rho
             X_1 = rhoX0(n,r,h1_comp) / rho
             tmp1 =   2.7d-3 * T_6_third
             tmp2 = -7.78d-3 * T_6_third**2
             tmp3 = -1.49d-4 * T_6_third**3
             g14 = 1.0_dp_t + tmp1 + tmp2 + tmp3
             tmp1 = 8.67d27 * g14 * X_CNO * X_1 * rho / T_6_third**2
             tmp2 = dexp(-1.5228d2 / T_6_third)
             Hbar(n,r) = tmp1 * tmp2
          enddo
       enddo

    elseif (prob_type .eq. 3) then

       he4_comp = network_species_index("helium-4")

       ! off-center heating for sub_chandra
       if (time .le. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do n=0,max_radial_level
             do r = 0, nr(1)-1
                if (spherical .eq. 0) then
                   call bl_error("ERROR: heating not supported")
                else
                   ! spherical -- lower amplitude heating term
                   Hbar(n,r) = fac * heating_peak * &
                        exp(-((r_cc_loc(n,r) - heating_rad)**2)/ heating_sigma**2)

                   ! only heat if there is He-4
                   Hbar(n,r) = (rhoX0(n,r,he4_comp)/rho0(n,r)) * Hbar(n,r)
                endif
             enddo
          enddo
       end if

    elseif (prob_type .eq. 4) then
       ! Apply both heating and cooling for an Urca process

       if (time .le. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do n=0,max_radial_level
             do r = 0, nr(1)-1
                if (spherical .eq. 0) then
                   ! plane-parallel -- do the heating term in paper II (section 4)
                   ! plus a similar cooling term for Urca
                   Hbar(n,r) = fac * (&
                        heating_peak * exp(-((r_cc_loc(n,r) - heating_rad)**2)/ heating_sigma) + &
                        cooling_peak * exp(-((r_cc_loc(n,r) - cooling_rad)**2)/ cooling_sigma))

                else
                   ! spherical -- lower amplitude heating/cooling term
                   Hbar(n,r) = fac * (&
                        heating_peak * exp(-((r_cc_loc(n,r) - heating_rad)**2)/ heating_sigma) + &
                        cooling_peak * exp(-((r_cc_loc(n,r) - cooling_rad)**2)/ cooling_sigma))
                endif
             enddo
          enddo
       end if

    else

       write(*,*) prob_type
       call bl_error("prob_type not yet supported.")

    endif

    return
  end subroutine get_heating

  subroutine make_Sbar(Sbar, rho0, tempbar, rhoX0, Hbar) bind(C, name="make_Sbar")

    double precision, intent(inout) :: Sbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rhoX0(0:max_radial_level,0:nr_fine-1, 1:nspec)
    double precision, intent(in   ) :: Hbar(0:max_radial_level,0:nr_fine-1)

    integer :: r, n
    type (eos_t) :: eos_state

    do n=0,max_radial_level
       do r=0,nr(1)-1

          ! (rho, T) --> p,h, etc
          eos_state%rho   = rho0(n,r)
          eos_state%T     = tempbar(n,r)
          eos_state%xn(:) = rhoX0(n,r,1:nspec)/rho0(n,r)

          call eos(eos_input_rt, eos_state)

          Sbar(n,r) = Hbar(n,r) * eos_state%dpdt / (eos_state%rho * eos_state%cp * eos_state%dpdr)

       enddo
    enddo

  end subroutine make_Sbar

  subroutine compute_gamma1bar(gamma1bar, rho0, tempbar, rhoX0, p0) bind(C, name="compute_gamma1bar")

    double precision, intent(inout) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rhoX0(0:max_radial_level,0:nr_fine-1,1:nspec)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    integer :: r, n
    type (eos_t) :: eos_state

    do n=0,max_radial_level
       do r=0,nr(1)-1

          ! (rho, p) --> gamma1bar
          eos_state%rho   = rho0(n,r)
          eos_state%p     = p0(n,r)
          eos_state%xn(:) = rhoX0(n,r,:)/rho0(n,r)

          eos_state%T     = tempbar(n,r)

          call eos(eos_input_rp, eos_state)

          gamma1bar(n,r) = eos_state%gam1

       enddo
    enddo

  end subroutine compute_gamma1bar

  subroutine update_temp(rho0, tempbar, rhoX0, p0) bind(C, name="update_temp")

    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rhoX0(0:max_radial_level,0:nr_fine-1,1:nspec)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    integer :: r, n
    type (eos_t) :: eos_state

    do n=0,max_radial_level
       do r=0,nr(1)-1

          ! (rho,p) --> T,h, etc
          eos_state%rho   = rho0(n,r)
          eos_state%p     = p0(1,r)
          eos_state%xn(:) = rhoX0(n,r,:)/rho0(n,r)

          eos_state%T     = tempbar(n,r)

          call eos(eos_input_rp, eos_state)

          tempbar(n,r) = eos_state%T

       enddo
    enddo

  end subroutine update_temp

  subroutine update_species(rho0, rho0_predicted_edge, rhoX0_old, rhoX0_new, &
       w0, r_cc_loc, r_edge_loc, dt) bind(C, name="update_species")

    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in) :: rhoX0_old(0:max_radial_level,0:nr_fine-1,1:nspec)
    double precision, intent(inout) :: rhoX0_new(0:max_radial_level,0:nr_fine-1,1:nspec)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine-1)
    double precision, value, intent(in   ) :: dt

    integer :: r, n, comp, comp2
    double precision :: X0(0:max_radial_level,0:nr_fine-1)
    double precision :: edge(0:max_radial_level,0:nr_fine-1)
    double precision :: force, delta, sumX, frac

    do comp=1,nspec

       do n=0,max_radial_level
          do r=0,nr(1)-1

             X0(n, r) = max(rhoX0_old(n, r, comp) / rho0(n, r), ZERO)

             force = 0

          enddo
       enddo

       do n=0,max_radial_level
          do r=0,nr(1)-1

             ! make edge state

             edge(n, r) = rho0_predicted_edge(n, r) * edge(n, r)

             ! update rhoX_0
             if (spherical .eq. 0) then
                rhoX0_new(n, r, comp) = rhoX0_old(n, r, comp) - dt / dr(n) *(edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r))
             else
                rhoX0_new(n, r, comp) = rhoX0_old(n, r, comp) - dt / dr(n) /r_cc_loc(n,r)**2* &
                     (r_edge_loc(n,r+1)**2 * edge(n,r+1) * w0(n,r+1) - &
                     r_edge_loc(n,r  )**2 * edge(n,r  ) * w0(n,r  ))
             endif

          enddo
       enddo
    enddo

    ! don't let the species leave here negative
    do comp=1,nspec

       if (minval(rhoX0_new(:,:,comp)) .lt. ZERO) then
          do n=0,max_radial_level
             do r=0,nr(1)-1
                if (rhoX0_new(n,r,comp) .lt. ZERO) then
                   delta = -rhoX0_new(n,r,comp)
                   sumX = ZERO
                   do comp2 = 1, nspec
                      if (comp2 .ne. comp .and. rhoX0_new(n,r,comp2) .ge. ZERO) then
                         sumX = sumX + rhoX0_new(n,r,comp2)
                      endif
                   enddo
                   do comp2 = 1, nspec
                      if (comp2 .ne. comp .and. rhoX0_new(n,r,comp2) .ge. ZERO) then
                         frac = rhoX0_new(n,r,comp2) / sumX
                         rhoX0_new(n,r,comp2) = rhoX0_new(n,r,comp2) - frac * delta
                      endif
                   enddo
                   rhoX0_new(n,r,comp) = ZERO
                endif
             enddo
          enddo
       endif

    enddo

  end subroutine update_species


end module test_basestate_module
