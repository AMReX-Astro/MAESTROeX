
module initdata_module

  use parallel, only: parallel_IOProcessor
  use amrex_fort_module, only : amrex_spacedim, amrex_random
  use amrex_constants_module
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       prob_lo, prob_hi
  use probin_module, only: velpert_amplitude, velpert_radius, velpert_steep, velpert_scale, velpert_min

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
    integer :: ip, jp, kp
    double precision :: xloc(3), upert(3)

    double precision :: rand

    ! random numbers between -1 and 1
    double precision :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    double precision :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k
    double precision :: normk(3,3,3)

    ! cos and sin of (2*pi*kx/L + phix), etc
    double precision :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    double precision :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    double precision :: rloc

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

    ! generate random numbers
    ! random numbers are not currently supported
    ! use functions that result in numbers between (-1,1)
    do k=1,3
       do j=1,3
          do i=1,3
             rand = amrex_random()
             rand = 2.0d0*rand - 1.0d0
             ! rand = (.5)**i * (.7)**j * (.3)**k * (-1.)**i
             alpha(i,j,k) = rand
             rand = amrex_random()
             rand = 2.0d0*rand - 1.0d0
             ! rand = (.5)**i * (.3)**j * (.7)**k * (-1.)**j
             beta(i,j,k) = rand
             rand = amrex_random()
             rand = 2.0d0*rand - 1.0d0
             ! rand = (.3)**i * (.5)**j * (.7)**k * (-1.)**k
             gamma(i,j,k) = rand
             rand = amrex_random()
             ! rand = (.3)**i * (.7)**j * (.5)**k
             rand = 2.0d0*M_PI*rand
             phix(i,j,k) = rand
             rand = amrex_random()
             ! rand = (.7)**i * (.3)**j * (.5)**k
             rand = 2.0d0*M_PI*rand
             phiy(i,j,k) = rand
             rand = amrex_random()
             ! rand = (.7)**i * (.5)**j * (.3)**k
             rand = 2.0d0*M_PI*rand
             phiz(i,j,k) = rand
          enddo
       enddo
    enddo

    ! compute the norm of k
    do k=1,3
       do j=1,3
          do i=1,3
             normk(i,j,k) = sqrt(dble(i)**2+dble(j)**2+dble(k)**2)
          enddo
       enddo
    enddo

    do k=lo(3), hi(3)
       xloc(3) = prob_lo(3) + (dble(k)+0.5d0)*dx(3)

       do j=lo(2),hi(2)
          xloc(2) = prob_lo(2) + (dble(j)+0.5d0)*dx(2)

          do i=lo(1),hi(1)
             xloc(1) = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

             upert = 0.d0

             rloc = xloc(2)

             ! loop over the 9 combinations of fourier components
             do kp = 1, 3
                do jp = 1, 3
                   do ip = 1, 3
                      ! compute cosines and sines
                      cx(ip,jp,kp) = cos(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(ip,jp,kp))
                      cy(ip,jp,kp) = cos(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(ip,jp,kp))
                      sx(ip,jp,kp) = sin(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(ip,jp,kp))
                      sy(ip,jp,kp) = sin(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(ip,jp,kp))
                   end do
                end do
             end do

             ! loop over the 9 combinations of fourier components
             do kp = 1, 3
                do jp = 1, 3
                   do ip = 1, 3
                      ! compute contribution from perturbation velocity from each mode
                      upert(1) = upert(1) + &
                           (-gamma(ip,jp,kp)*dble(j)*cx(ip,jp,kp)*sy(ip,jp,kp) &
                             +beta(ip,jp,kp)*dble(j)*cy(ip,jp,kp)*sx(ip,jp,kp)) &
                             / normk(ip,jp,kp)

                      upert(2) = upert(2) + &
                           (gamma(ip,jp,kp)*dble(i)*cy(ip,jp,kp)*sx(ip,jp,kp) &
                           -alpha(ip,jp,kp)*dble(j)*cx(ip,jp,kp)*sy(ip,jp,kp)) &
                             / normk(ip,jp,kp)

                   end do
                end do
             end do

             ! apply the cutoff function to the perturbational velocity
             do ip = 1, nc_v
                upert(ip) = velpert_amplitude*upert(ip) &
                     *(0.5d0+0.5d0*tanh((velpert_radius-rloc)/velpert_steep))
             end do

             ! add perturbational velocity to background velocity
             vel(i,j,k,1:nc_v) = vel(i,j,k,1:nc_v) + upert(1:nc_v)

          end do
       end do
    end do

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

  end subroutine initdata_sphr

end module initdata_module
