
module initdata_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use eos_type_module
  use eos_module
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       perturb_model, prob_lo, prob_hi
  use probin_module, only: rho_1, rho_2, vel_amplitude, vel_width, nmodes
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

    integer          :: i,j,k,r,n
    double precision :: x,y,z

    double precision :: dens_pert, rhoh_pert, temp_pert, rand
    double precision :: rhoX_pert(nspec), vel_pert(3)
    double precision :: alpha(nmodes), phi(nmodes)

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

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

       do n = 1, nmodes
          call random_number(rand)
          alpha(n) = 2.0d0*rand - 1.0d0
          call random_number(rand)
          phi(n) = 2.0d0*M_PI*rand
       enddo

       ! write(*,*) "alpha = ", alpha

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
          do j = lo(2), hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
             do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                if (amrex_spacedim .eq. 2) then
                   call perturb_2d(x, y, p0_init(lev,j), s0_init(lev,j,:), &
                        alpha, phi, &
                        dens_pert, rhoX_pert, vel_pert)
                else
                   call perturb_3d(x, y, z, p0_init(lev,k), s0_init(lev,k,:), &
                        alpha, phi, &
                        dens_pert, rhoX_pert, vel_pert)
                end if

                scal(i,j,k,rho_comp) = dens_pert
                scal(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:nspec)

                vel(i,j,k,1:nc_v) = vel_pert(1:nc_v)

             end do
          end do
       end do

    end if

  end subroutine initdata

  subroutine perturb_2d(x, y, p0_init, s0_init, alpha, phi, dens_pert, rhoX_pert, vel_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: x, y
    double precision, intent(in ) :: p0_init, s0_init(1:nscal), alpha(1:nmodes), phi(1:nmodes)
    double precision, intent(out) :: dens_pert
    double precision, intent(out) :: rhoX_pert(1:nspec)
    double precision, intent(out) :: vel_pert(1:3)

    double precision :: L_x, pertheight, y0, pert
    integer :: n

    type (eos_t) :: eos_state

    L_x = prob_hi(1) - prob_lo(1)

    pertheight = 0.01d0*0.5d0*(cos(2.0d0*M_PI*x/L_x)+cos(2.0d0*M_PI*(L_x-x)/L_x)) + 0.5d0

    dens_pert = rho_1 + ((rho_2-rho_1)/2.0d0)*(1.0d0+tanh((y-pertheight)/0.005d0))

    rhoX_pert(:) = dens_pert*s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    ! y0 = 0.5d0 * (prob_lo(2) + prob_hi(2))
    !
    ! pert = 0.0d0
    !
    ! if (nmodes == 1) then
    !     pert = pert + vel_amplitude*0.5d0*(cos(2.d0*M_PI*x/L_x) + &
    !     cos(2.d0*M_PI*(L_x- x)/L_x))
    ! else
    !     do n = 1, nmodes
    !         pert = pert + vel_amplitude*alpha(n)* &
    !              cos(2.d0*M_PI*x/L_x + phi(n))
    !     enddo
    ! endif

    vel_pert(:) = 0.0d0
    !vel_pert(2) = exp(-(y-y0)**2/vel_width**2)*pert

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, alpha, phi, dens_pert, rhoX_pert, vel_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: x, y, z
    double precision, intent(in ) :: p0_init, s0_init(1:nscal), alpha(1:nmodes), phi(1:nmodes)
    double precision, intent(out) :: dens_pert
    double precision, intent(out) :: rhoX_pert(1:nspec)
    double precision, intent(out) :: vel_pert(1:3)

    double precision :: L_x, pertheight, pert, y0
    integer :: n

    type (eos_t) :: eos_state

    L_x = prob_hi(1) - prob_lo(1)

    pertheight = 0.01d0*0.5d0*(cos(2.0d0*M_PI*x/L_x)+cos(2.0d0*M_PI*(L_x-x)/L_x)) + 0.5d0

    dens_pert = rho_1 + ((rho_2-rho_1)/2.0d0)*(1.0d0+tanh((y-pertheight)/0.005d0))

    rhoX_pert = dens_pert*s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    y0 = 0.5d0 * (prob_lo(2) + prob_hi(2))

    pert = 0.0d0

    if (nmodes == 1) then
       pert = pert + vel_amplitude*0.5d0*(cos(2.d0*M_PI*x/L_x) + &
            cos(2.d0*M_PI*(L_x- x)/L_x))
    else
       do n = 1, nmodes
          pert = pert + vel_amplitude*alpha(n)* &
               cos(2.d0*M_PI*x/L_x + phi(n))
       enddo
    endif

    vel_pert(:) = 0.0d0
    vel_pert(2) = exp(-(y-y0)**2/vel_width**2)*pert

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

    !     Local variables
    integer          :: i,j,k,comp
    double precision :: x,y,z
    double precision :: dens_pert, rhoh_pert, temp_pert
    double precision :: rhoX_pert(nspec)
    double precision, pointer :: p0_cart(:,:,:,:)

    type (eos_t) :: eos_state

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

    ! initialize the domain with the base state
    scal(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_s) = 0.d0

    ! if we are spherical, we want to make sure that p0 is good, since that is
    ! what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    ! then initialize h from rho, X, and p0.
    call bl_allocate(p0_cart,lo,hi,1)

    ! initialize temp
    call put_1d_array_on_cart_sphr(lo,hi,scal(:,:,:,temp_comp),scal_lo,scal_hi,1, &
         s0_init(0,:,temp_comp),dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)

    ! initialize p0_cart
    call put_1d_array_on_cart_sphr(lo,hi,p0_cart,lo,hi,1,p0_init,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)

    ! initialize species
    do comp = spec_comp, spec_comp+nspec-1
       call put_1d_array_on_cart_sphr(lo,hi,scal(:,:,:,comp),scal_lo,scal_hi,1, &
            s0_init(0,:,comp),dx,0,0,r_cc_loc,r_edge_loc, &
            cc_to_r,ccr_lo,ccr_hi)
    end do

    ! initialize rho as sum of partial densities rho*X_i
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             scal(i,j,k,rho_comp) = 0.d0
             do comp = spec_comp, spec_comp+nspec-1
                scal(i,j,k,rho_comp) = scal(i,j,k,rho_comp) + scal(i,j,k,comp)
             enddo
          enddo
       enddo
    enddo

    ! initialize (rho h) and T using the EOS
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%T     = scal(i,j,k,temp_comp)
             eos_state%p     = p0_cart(i,j,k,1)
             eos_state%rho   = scal(i,j,k,rho_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             call eos(eos_input_rp, eos_state)

             scal(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h
             scal(i,j,k,temp_comp) = eos_state%T

          enddo
       enddo
    enddo

    ! if (perturb_model) then
    !
    !    ! add an optional perturbation
    !    do k = lo(3), hi(3)
    !       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
    !
    !       do j = lo(2), hi(2)
    !          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
    !
    !          do i = lo(1), hi(1)
    !             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
    !
    !             call perturb_3d_sphr(x, y, z, p0_cart(i,j,k,1), scal(i,j,k,:), &
    !                                  dens_pert, rhoh_pert, rhoX_pert, temp_pert)
    !
    !             scal(i,j,k,rho_comp) = dens_pert
    !             scal(i,j,k,rhoh_comp) = rhoh_pert
    !             scal(i,j,k,temp_comp) = temp_pert
    !             scal(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
    !          enddo
    !       enddo
    !    enddo
    !
    ! end if

    call bl_deallocate(p0_cart)

  end subroutine initdata_sphr
  !
  ! subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
  !                            rhoX_pert, temp_pert)
  !
  !   ! apply an optional perturbation to the initial temperature field
  !   ! to see some bubbles
  !   double precision, intent(in ) :: x, y, z
  !   double precision, intent(in ) :: p0_init, s0_init(nscal)
  !   double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
  !   double precision, intent(out) :: rhoX_pert(nspec)
  !
  !   double precision :: temp, t0
  !   double precision :: x0, y0, z0, r0
  !
  !   type (eos_t) :: eos_state
  !
  !   t0 = s0_init(temp_comp)
  !
  !   ! center of star is at 2.5d8
  !   x0 = 2.5d8
  !   y0 = 2.5d8
  !   z0 = 3.0d8
  !
  !   ! Tanh bubbles
  !   r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 2.5d6
  !
  !   ! This case works
  !   temp = t0 * (1.0d0 + 2.0d0*(.150d0 * 0.5d0 * (1.0d0 + tanh((2.0d0-r0)))))
  !
  !   ! Use the EOS to make this temperature perturbation occur at constant
  !   ! pressure
  !   eos_state%T     = temp
  !   eos_state%p     = p0_init
  !   eos_state%rho   = s0_init(rho_comp)
  !   eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)
  !
  !   call eos(eos_input_tp, eos_state)
  !
  !   dens_pert = eos_state%rho
  !   rhoh_pert = eos_state%rho * eos_state%h
  !   rhoX_pert = dens_pert*eos_state%xn(:)
  !
  !   temp_pert = temp
  !
  ! end subroutine perturb_3d_sphr

end module initdata_module
