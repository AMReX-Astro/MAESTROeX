
module initdata_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use amrex_constants_module
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim, amrex_random
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       prob_lo, prob_hi, model_file
  use probin_module, only: velpert_amplitude, velpert_radius, velpert_steep, velpert_scale
  use eos_module
  use eos_type_module
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr
  use model_parser_module
  ! use rotation_module, only: omega

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

    integer          :: i,j,k,r

    if (parallel_IOProcessor()) then
       print*,"Create a local copy of initdata.f90 in your build directory"
       print*,"Here is a sample that initializes v=0 and the scalars using s0"
    end if

    ! abort program
    call bl_error("Planar initdata not written")

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

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

    !     Local variables
    integer          :: i,j,k,comp
    integer :: iloc, jloc, kloc
    double precision, pointer :: p0_cart(:,:,:,:)

    type (eos_t) :: eos_state
    integer :: pt_index(3)

    ! radius, or distance, to center of star
    double precision :: rloc

    ! the point we're at
    double precision :: x, y, z

    ! normal
    double precision :: normal(3)

    ! cos and sin of (2*pi*kx/L + phix), etc
    double precision :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    double precision :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    double precision :: rand

    ! random numbers between -1 and 1
    double precision :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    double precision :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! perturbational velocity to add
    double precision :: vpert(3)

    ! L2 norm of k
    double precision :: normk(3,3,3)

    ! density pertubation and amplitude
    ! double precision :: delta_rho
    ! double precision, parameter :: amplitude = 2.d-4

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
    !! add pertubation and renormalize partial densities
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             scal(i,j,k,rho_comp) = 0.d0
             do comp = spec_comp, spec_comp+nspec-1
                scal(i,j,k,rho_comp) = scal(i,j,k,rho_comp) + scal(i,j,k,comp)
             enddo

             ! delta_rho = amrex_random()
             ! delta_rho = amplitude * (2.0d0*delta_rho - 1.0d0) * scal(i,j,k,rho_comp)

             ! do comp = spec_comp, spec_comp+nspec-1
             !    scal(i,j,k,comp) = scal(i,j,k,comp)  * (scal(i,j,k,rho_comp) + delta_rho) / scal(i,j,k,rho_comp)
             ! enddo
             !
             ! scal(i,j,k,rho_comp) = scal(i,j,k,rho_comp) + delta_rho

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

             pt_index = (/ i, j, k /)

             call eos(eos_input_rp, eos_state, pt_index)

             scal(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h
             scal(i,j,k,temp_comp) = eos_state%T

          enddo
       enddo
    enddo

    call bl_deallocate(p0_cart)

    ! initialize the velocity to zero everywhere
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

    if (iconvel_model > 0) then

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

       ! now do the big loop over all points in the domain
       do kloc = lo(3),hi(3)
          z = prob_lo(3) + (dble(kloc)+HALF)*dx(3) 
          do jloc = lo(2),hi(2)
             y = prob_lo(2) + (dble(jloc)+HALF)*dx(2)
             do iloc = lo(1),hi(1)
                x = prob_lo(1) + (dble(iloc)+HALF)*dx(1) 

                ! set perturbational velocity to zero
                vpert = ZERO

                rloc = sqrt((x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2)

                ! loop over the 27 combinations of fourier components
                do i=1,3
                   do j=1,3
                      do k=1,3
                         ! compute cosines and sines
                         cx(i,j,k) = cos(2.0d0*M_PI*dble(i)*x/velpert_scale + phix(i,j,k))
                         cy(i,j,k) = cos(2.0d0*M_PI*dble(j)*y/velpert_scale + phiy(i,j,k))
                         cz(i,j,k) = cos(2.0d0*M_PI*dble(k)*z/velpert_scale + phiz(i,j,k))
                         sx(i,j,k) = sin(2.0d0*M_PI*dble(i)*x/velpert_scale + phix(i,j,k))
                         sy(i,j,k) = sin(2.0d0*M_PI*dble(j)*y/velpert_scale + phiy(i,j,k))
                         sz(i,j,k) = sin(2.0d0*M_PI*dble(k)*z/velpert_scale + phiz(i,j,k))
                      enddo
                   enddo
                enddo

                ! loop over the 27 combinations of fourier components
                do k=1,3
                   do j=1,3
                      do i=1,3
                         ! compute contribution from perturbation velocity from each mode
                         vpert(1) = vpert(1) + &
                              (-gamma(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k) &
                              +beta(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
                              / normk(i,j,k)

                         vpert(2) = vpert(2) + &
                              (gamma(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                              -alpha(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
                              / normk(i,j,k)

                         vpert(3) = vpert(3) + &
                              ( -beta(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                              +alpha(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k)) &
                              / normk(i,j,k)
                      enddo
                   enddo
                enddo

                ! apply the cutoff function to the perturbational velocity
                do i=1,3
                   vpert(i) = velpert_amplitude * vpert(i) &
                        *(0.5d0+0.5d0*tanh((velpert_radius-rloc)/velpert_steep))
                enddo

                if (rloc > 0.0d0) then
                   normal(1) = x * ONE / rloc
                   normal(2) = y * ONE / rloc
                   normal(3) = z * ONE / rloc

                   vel(i,j,k,1:3) = normal(1:3) * interpolate(rloc, iconvel_model)
                endif

                ! add perturbational velocity to background velocity
                do i=1,3
                   vel(iloc,jloc,kloc,i) = vel(iloc,jloc,kloc,i) * (1.0d0 + vpert(i))
                enddo

             enddo
          enddo
       enddo

       ! ! now do the big loop over all points in the domain
       ! do k = lo(3),hi(3)
       !    z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
       !    do j = lo(2),hi(2)
       !       y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
       !       do i = lo(1),hi(1)
       !          x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
       !
       !          ! compute perpendicular distance to the rotation axis
       !          rloc = sqrt(x**2 + y**2 + z**2)
       !
       !          if (rloc > 0.0d0) then
       !             normal(1) = x * ONE / rloc
       !             normal(2) = y * ONE / rloc
       !             normal(3) = z * ONE / rloc
       !
       !             vel(i,j,k,1:3) = normal(1:3) * interpolate(rloc, iconvel_model)
       !          endif
       !
       !       enddo
       !    enddo
       ! enddo

    endif

    !
    ! #ifdef ROTATION
    !
    !     ! initialize solid body rotation
    !
    !     ! define where center of star is
    !     ! this currently assumes the star is at the center of the domain
    !     xc(1) = 0.5d0*(prob_lo(1)+prob_hi(1))
    !     xc(2) = 0.5d0*(prob_lo(2)+prob_hi(2))
    !     xc(3) = 0.5d0*(prob_lo(3)+prob_hi(3))
    !
    !     ! now do the big loop over all points in the domain
    !     do iloc = lo(1),hi(1)
    !        do jloc = lo(2),hi(2)
    !           do kloc = lo(3),hi(3)
    !
    !              ! compute where we physically are
    !              xloc(1) = prob_lo(1) + (dble(iloc)+0.5d0)*dx(1)
    !              xloc(2) = prob_lo(2) + (dble(jloc)+0.5d0)*dx(2)
    !              xloc(3) = prob_lo(3) + (dble(kloc)+0.5d0)*dx(3)
    !
    !              ! compute perpendicular distance to the rotation axis
    !              rloc = sqrt(sum(xloc(1:2) - xc(1:2))**2)
    !
    !              ! write(*,*) "omega =", omega, "rloc =", rloc
    !
    !              ! vx = -2 pi r omega sin(theta)
    !              ! vx = -2 pi r omega cos(theta)
    !              ! sin(theta) = y/r, cos(theta) = x/r
    !              if (rloc > 0.0d0) then
    !                 vel(iloc,jloc,kloc,1) = -2 * M_PI * rloc * omega * xloc(2) / rloc
    !                 vel(iloc,jloc,kloc,2) = 2 * M_PI * rloc * omega * xloc(1) / rloc
    !              endif
    !
    !           enddo
    !        enddo
    !     enddo
    !
    !
    ! #endif

  end subroutine initdata_sphr

end module initdata_module
