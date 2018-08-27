
module initdata_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
                                perturb_model, prob_lo, prob_hi
  use probdata_module, only: pert_factor, y_pert_center, pert_width, single, do_isentropic
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

  private

contains

  subroutine initdata(lev, time, lo, hi, &
                      scal, scal_lo, scal_hi, nc_s, &
                      vel, vel_lo, vel_hi, nc_v, &
                      s0_init, p0_init, &
                      dx) bind(C, name="initdata")

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

    integer          :: i,j,k,r
    double precision :: x,y,z

    double precision :: dens_pert, rhoh_pert, temp_pert
    double precision :: rhoX_pert(nspec)

    ! set velocity to zero
    vel = 0.d0

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)

          if (amrex_spacedim .eq. 2) then
             r = j
          else
             r = k
          end if

          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

             scal(i,j,k,rho_comp)  = s0_init(lev,r,rho_comp)
             scal(i,j,k,rhoh_comp) = s0_init(lev,r,rhoh_comp)
             scal(i,j,k,temp_comp) = s0_init(lev,r,temp_comp)
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  s0_init(lev,r,spec_comp:spec_comp+nspec-1)

             ! initialize pi_comp to zero for now
             scal(i,j,k,pi_comp) = 0.d0

          end do
       end do
    end do

    ! add an optional perturbation
    if (perturb_model) then

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
          do j = lo(2), hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
             do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                if (amrex_spacedim .eq. 2) then
                   call perturb_2d(x, y, p0_init(lev,j), s0_init(lev,j,:), &
                        dens_pert, rhoh_pert, rhoX_pert, temp_pert)
                else
                   call perturb_3d(x, y, z, p0_init(lev,k), s0_init(lev,k,:), &
                                   dens_pert, rhoh_pert, rhoX_pert, temp_pert)
                end if

                scal(i,j,k,rho_comp) = dens_pert
                scal(i,j,k,rhoh_comp) = rhoh_pert
                scal(i,j,k,temp_comp) = temp_pert
                scal(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:nspec)

             end do
          end do
       end do

    end if

  end subroutine initdata

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: x, y
    double precision, intent(in ) :: p0_init, s0_init(1:nscal)
    double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
    double precision, intent(out) :: rhoX_pert(1:nspec)

    double precision :: xn(nspec)
    double precision :: dens
    double precision :: x1, y1, r1, x2, y2, r2

    type (eos_t) :: eos_state

    if (.not. single) then
       x1 = prob_lo(1) + (prob_hi(1)-prob_lo(1))/3.d0
       x2 = prob_lo(1) + 2.d0*(prob_hi(1)-prob_lo(1))/3.d0

       y1 = y_pert_center
       y2 = y_pert_center

       r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width
       r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / pert_width

       if (r1 < 2.0d0) then
          dens = s0_init(rho_comp) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r1))))
          xn(:) = 0.d0
          xn(2) = 1.d0

       else if (r2 < 2.0d0) then
          dens = s0_init(rho_comp) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r2))))
          xn(:) = 0.d0
          xn(3) = 1.d0

       else
          dens = s0_init(rho_comp)
          xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)
       endif

    else

       x1 = prob_lo(1) + 0.5d0*(prob_hi(1)-prob_lo(1))
       y1 = y_pert_center
       r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width

       if (r1 < 2.0d0) then
          dens = s0_init(rho_comp) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r1))))
          xn(:) = 0.d0
          xn(2) = 1.d0

       else
          dens = s0_init(rho_comp)
          xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)
       endif

    endif

    ! Use the EOS to make this temperature perturbation occur at constant
    ! pressure
    eos_state%T     = 10000.d0   ! guess
    eos_state%p     = p0_init
    eos_state%rho   = dens
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rp, eos_state)

    dens_pert = dens
    rhoh_pert = dens_pert * eos_state%h
    rhoX_pert = dens_pert * xn(:)

    temp_pert = eos_state % T

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: x, y, z
    double precision, intent(in ) :: p0_init, s0_init(1:nscal)
    double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
    double precision, intent(out) :: rhoX_pert(1:nspec)

    call bl_error("ERROR: perturb_3d not implemented")

  end subroutine perturb_3d


  subroutine initdata_sphr(time, lo, hi, &
                           scal, scal_lo, scal_hi, nc_s, &
                           vel, vel_lo, vel_hi, nc_v, &
                           s0_init, p0_init, &
                           dx, r_cc_loc, r_edge_loc, &
                           cc_to_r, ccr_lo, ccr_hi) bind(C, name="initdata_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
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
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
                                               ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    call bl_error("ERROR: initdata_sphr not implemented")

  end subroutine initdata_sphr

  subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles
    double precision, intent(in ) :: x, y, z
    double precision, intent(in ) :: p0_init, s0_init(nscal)
    double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
    double precision, intent(out) :: rhoX_pert(nspec)

    call bl_error("ERROR: perturb_3d_sphr not implemented")

  end subroutine perturb_3d_sphr

end module initdata_module
