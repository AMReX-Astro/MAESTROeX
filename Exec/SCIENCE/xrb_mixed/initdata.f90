
module initdata_module

  use network, only: nspec
  use amrex_error_module
  use amrex_constants_module
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, &
       spec_comp, pi_comp, prob_lo, prob_hi, perturb_model
  use probin_module, only: xrb_pert_height, xrb_pert_size, xrb_pert_factor, &
       xrb_pert_type, num_vortices, apply_vel_field, velpert_scale, &
       velpert_amplitude, velpert_height_loc
  use eos_module
  use eos_type_module

  implicit none

  private

contains

  subroutine initdata(lev, time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, dx) bind(C, name="initdata")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), 1:nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), 1:nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k,r,vortex
    double precision :: x, y, z, xcen, ycen, zcen, dist, xdist, ydist
    double precision :: xloc(3), upert(3), xloc_vortices(num_vortices), offset
    double precision :: rhoX_pert(nspec), dens_pert, rhoh_pert, temp_pert

    offset = (prob_hi(1) - prob_lo(1)) / (num_vortices + 1)

    do i = 1, num_vortices
       xloc_vortices(i) = dble(i) * offset
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (amrex_spacedim .eq. 1) then
                r = i
             else if (amrex_spacedim .eq. 2) then
                r = j
             else if (amrex_spacedim .eq. 3) then
                r = k
             end if

             ! set scalars using s0
             scal(i,j,k,rho_comp)  = s0_init(lev,r,rho_comp)
             scal(i,j,k,rhoh_comp) = s0_init(lev,r,rhoh_comp)
             scal(i,j,k,temp_comp) = s0_init(lev,r,temp_comp)
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  s0_init(lev,r,spec_comp:spec_comp+nspec-1)

             ! initialize pi to zero for now
             scal(i,j,k,pi_comp) = 0.d0

          end do
       end do
    end do

    if (perturb_model) then

       xcen = center(1)
       if (amrex_spacedim .eq. 2) then
          ycen = xrb_pert_height
          zcen = ZERO
       else
          ycen = center(2)
          zcen = xrb_pert_height
       endif

       do k = lo(3), hi(3)
          z = prob_lo(3) + (dble(k)+HALF) * dx(3)
          do j = lo(2), hi(2)
             y = prob_lo(2) + (dble(j)+HALF) * dx(2)
             do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i)+HALF) * dx(1)

                if (amrex_spacedim .eq. 1) then
                   r = i
                else if (amrex_spacedim .eq. 2) then
                   r = j
                else if (amrex_spacedim .eq. 3) then
                   r = k
                end if

                dist = sqrt((x-xcen)**2 + (y-ycen)**2 + (z-zcen)**2)

                call perturb(dist, p0_init(lev,r), s0_init(lev,r,:), &
                     dens_pert, rhoh_pert, rhoX_pert, temp_pert)

                scal(i,j,k,rho_comp) = dens_pert
                scal(i,j,k,rhoh_comp) = rhoh_pert
                scal(i,j,k,temp_comp) = temp_pert
                scal(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
             enddo
          enddo
       enddo
    endif

    ! set velocity to zero
    vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3),1:nc_v) = 0.d0

    if (apply_vel_field) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             xloc(2) = prob_lo(2) + (dble(j)+HALF)*dx(2)
             ydist = xloc(2) - velpert_height_loc

             do i = lo(1), hi(1)

                upert(:) = ZERO

                xloc(1) = prob_lo(1) + (dble(i)+HALF)*dx(1)

                ! loop over each vortex
                do vortex = 1, num_vortices

                   xdist = xloc(1) - xloc_vortices(vortex)

                   r = sqrt(xdist**2 + ydist**2)

                   ! e.g. Calder et al. ApJSS 143, 201-229 (2002)
                   ! we set things up so that every other vortex has the same
                   ! orientation
                   upert(1) = upert(1) - ydist * &
                        velpert_amplitude * exp( -r**2/(TWO*velpert_scale)) &
                        * (-ONE)**vortex

                   upert(2) = upert(2) + xdist * &
                        velpert_amplitude * exp(-r**2/(TWO*velpert_scale)) &
                        * (-ONE)**vortex
                enddo

                vel(i,j,k,:) = vel(i,j,k,:) + upert(:)

             enddo
          enddo
       enddo
    endif


  end subroutine initdata

  subroutine initdata_sphr(time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, &
       dx, r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind(C, name="initdata_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc  (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: cc_to_r   (ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

  end subroutine initdata_sphr


  subroutine perturb(distance, p0_init, s0_init,  &
       dens_pert, rhoh_pert, rhoX_pert, temp_pert)
    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: distance
    double precision, intent(in ) :: p0_init, s0_init(:)
    double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
    double precision, intent(out) :: rhoX_pert(:)

    double precision :: temp,dens,t0,d0,rad_pert

    integer, parameter :: perturb_temp = 1, perturb_dens = 2
    integer :: eos_input_flag

    type (eos_t) :: eos_state

    rad_pert = -xrb_pert_size**2 / (FOUR*log(HALF))

    select case (xrb_pert_type)
    case(perturb_temp)

       t0 = s0_init(temp_comp)
       temp = t0 * (ONE + xrb_pert_factor * exp(-distance**2 / rad_pert) )
       dens = s0_init(rho_comp)
       eos_input_flag = eos_input_tp

    case(perturb_dens)

       d0 = s0_init(rho_comp)
       dens = d0 * (ONE + xrb_pert_factor * exp(-distance**2 / rad_pert) )
       temp = s0_init(temp_comp)
       eos_input_flag = eos_input_rp

    end select

    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = dens
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_flag, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert(:) = dens_pert*eos_state%xn(:)

    temp_pert = eos_state%T

  end subroutine perturb

end module initdata_module
