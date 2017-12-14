
module sponge_module

  use bl_constants_module
  use parallel, only: parallel_IOProcessor
  use base_state_geometry_module, only: dr, r_end_coord, max_radial_level, nr_fine, center
  use meth_params_module, only: maestro_verbose, sponge_start_factor, &
                                sponge_center_density, sponge_kappa, spherical, &
                                drdxfac, prob_lo

  implicit none

  private

  double precision, save :: r_sp
  double precision, save :: r_md
  double precision, save :: r_tp
  double precision, save :: r_sp_outer ! outer sponge parameters used for spherical problems
  double precision, save :: r_tp_outer ! outer sponge parameters used for spherical problems
  
  ! the sponge_start_density should be the density below which the
  ! sponge first turns on.  Different problems may compute this in
  ! different ways (i.e. not using sponge_center_density and
  ! sponge_start_factor), so we provide this public module variable to
  ! ensure that the rest of the code always knows at what density the
  ! sponge begins.
  double precision, save, public :: sponge_start_density

contains

 subroutine init_sponge(rho0)

    ! The sponge has a HALF * ( 1 - cos( (r - r_sp)/L)) profile, where
    ! the width, L, is r_tp - r_sp.
    !
    ! The center of the sponge, r_md, is set to the radius where r =
    ! sponge_center_density
    !
    ! The start of the sponge, r_sp, (moving outward from the center)
    ! is the radius where r = sponge_start_factor * sponge_center_density
    ! 
    ! The top of the sponge is then 2 * r_md - r_tp

    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)

    double precision :: prob_lo_r,r_top
    integer :: r

    if (spherical .eq. 1) then
       prob_lo_r = 0.d0
    else
       prob_lo_r = prob_lo(AMREX_SPACEDIM)
    end if

    r_top = prob_lo_r + dble(r_end_coord(0,1)+1) * dr(0)
    r_sp = r_top

    sponge_start_density = sponge_start_factor*sponge_center_density

    ! set r_sp
    do r=0,r_end_coord(0,1)
       if (rho0(0,r) < sponge_start_density) then
          r_sp = prob_lo_r + (dble(r)+HALF) * dr(0)
          exit
       endif
    enddo

    ! set r_md
    r_md = r_top
    do r=0,r_end_coord(0,1)
       if (rho0(0,r) < sponge_center_density) then
          r_md = prob_lo_r + (dble(r)+HALF) * dr(0)
          exit
       endif
    enddo

    ! set r_tp
    r_tp = TWO * r_md - r_sp

    ! outer sponge parameters used for spherical problems
    if (spherical .eq. 1) then
       r_sp_outer = r_tp
       r_tp_outer = r_sp_outer + 4.d0 * drdxfac * dr(0)
    end if

    if (parallel_IOProcessor() .and. maestro_verbose .ge. 1) &
         write(6,1000) r_sp, r_tp
    if (parallel_IOProcessor() .and. maestro_verbose .ge. 1 .and. spherical .eq. 1) &
         write(6,1001) r_sp_outer, r_tp_outer
    if (parallel_IOProcessor() .and. maestro_verbose .ge. 1) &
         print*,""

1000 format('inner sponge: r_sp      , r_tp      : ',e20.12,2x,e20.12)
1001 format('outer sponge: r_sp_outer, r_tp_outer: ',e20.12,2x,e20.12)

  end subroutine init_sponge

  subroutine mk_sponge(lo,hi,sponge,s_lo,s_hi,dx,dt) bind(C, name="mk_sponge")

    integer        , intent(in   ) :: lo(3),hi(3),s_lo(3),s_hi(3)
    real(kind=dp_t), intent(inout) :: sponge(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    real(kind=dp_t), intent(in   ) :: dx(3),dt

    integer         :: i,j,k
    real(kind=dp_t) :: x,y,z,r,smdamp

    sponge = ONE

    if (spherical .eq. 0) then

#if (AMREX_SPACEDIM == 1)
       do i=lo(1),hi(1)
          r = prob_lo(1) + (dble(i)+HALF)*dx(1)
#elif (AMREX_SPACEDIM == 2)
       do j=lo(2),hi(2)
          r = prob_lo(2) + (dble(2)+HALF)*dx(2)
#else
       do k = lo(3),hi(3)
          r = prob_lo(3) + (dble(k)+HALF)*dx(3)
#endif
          if (r >= r_sp) then
             if (r < r_tp) then
                smdamp = HALF*(ONE - cos(M_PI*(r - r_sp)/(r_tp - r_sp)))
             else
                smdamp = ONE
             endif
#if (AMREX_SPACEDIM == 1)
             sponge(i,:,:) = ONE / (ONE + dt * smdamp* sponge_kappa)
#elif (AMREX_SPACEDIM == 2)
             sponge(:,j,:) = ONE / (ONE + dt * smdamp* sponge_kappa)
#else
             sponge(:,:,k) = ONE / (ONE + dt * smdamp* sponge_kappa)
#endif
          endif
       end do

    else

       do k = lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+HALF)*dx(3)

          do j = lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+HALF)*dx(2)
             
             do i = lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+HALF)*dx(1)

                r = sqrt( (x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2 )
                
                ! Inner sponge: damps velocities at edge of star
                if (r >= r_sp) then
                   if (r < r_tp) then
                      smdamp = HALF*(ONE - cos(M_PI*(r - r_sp)/(r_tp - r_sp)))
                   else
                      smdamp = ONE
                   endif
                   sponge(i,j,k) = ONE / (ONE + dt * smdamp * sponge_kappa)
                endif

                ! Outer sponge: damps velocities in the corners of the domain
                if (r >= r_sp_outer) then
                   if (r < r_tp_outer) then
                      smdamp = HALF * &
                           (ONE - cos(M_PI*(r - r_sp_outer)/(r_tp_outer - r_sp_outer)))
                   else
                      smdamp = ONE
                   endif
                   sponge(i,j,k) = sponge(i,j,k) / &
                        (ONE + dt * smdamp * 10.d0 * sponge_kappa)
                endif

             end do
          end do
       end do

    end if

  end subroutine mk_sponge


end module sponge_module
