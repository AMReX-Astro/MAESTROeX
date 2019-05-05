
module sponge_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : amrex_spacedim
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use base_state_geometry_module, only: dr, r_end_coord, max_radial_level, nr_fine, center
  use meth_params_module, only: maestro_verbose, sponge_start_factor, &
       sponge_center_density, prob_lo, anelastic_cutoff_density
  use probin_module, only: xrb_use_bottom_sponge, sponge_min

  implicit none

  private

  double precision, save :: topsponge_lo_r, topsponge_hi_r
  double precision, save :: botsponge_lo_r, botsponge_hi_r

  ! the sponge_start_density should be the density below which the
  ! sponge first turns on.  Different problems may compute this in
  ! different ways (i.e. not using sponge_center_density and
  ! sponge_start_factor), so we provide this public module variable to
  ! ensure that the rest of the code always knows at what density the
  ! sponge begins.
  double precision, save, public :: sponge_start_density

contains

  subroutine init_sponge(rho0) bind(C, name="init_sponge")

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

    prob_lo_r = prob_lo(AMREX_SPACEDIM)

    r_top = prob_lo_r + dble(r_end_coord(0,1)+1) * dr(0)
    topsponge_lo_r = r_top

    sponge_start_density = sponge_start_factor*sponge_center_density

    ! set topsponge_lo_r
    do r=0,r_end_coord(0,1)
       if (rho0(0,r) < sponge_start_density) then
          topsponge_lo_r = prob_lo_r + (dble(r)+HALF) * dr(0)
          exit
       endif
    enddo

    ! set topsponge_hi_r
    topsponge_hi_r = r_top
    do r=0,r_end_coord(0,1)
       if (rho0(0,r) < anelastic_cutoff_density) then
          topsponge_hi_r = prob_lo_r + (dble(r)+HALF) * dr(0)
          exit
       endif
    enddo

    ! set botsponge_lo_r
    do r=0,r_end_coord(0,1)
       if (rho0(0,r) < 6.d7) then
          botsponge_lo_r = prob_lo_r + (dble(r)+HALF) * dr(0)
          exit
       endif
    enddo

    ! set botsponge_hi_r
    do r=0,r_end_coord(0,1)
       if (rho0(0,r) < 5.d7) then
          botsponge_hi_r = prob_lo_r + (dble(r)+HALF) * dr(0)
          exit
       endif
    enddo

    if (parallel_IOProcessor() .and. maestro_verbose .ge. 1) &
         write(6,1000) topsponge_lo_r, topsponge_hi_r
    if (parallel_IOProcessor() .and. maestro_verbose .ge. 1) &
         print*,""

1000 format('inner sponge: r_sp      , r_tp      : ',e20.12,2x,e20.12)
1001 format('outer sponge: r_sp_outer, r_tp_outer: ',e20.12,2x,e20.12)

  end subroutine init_sponge

  subroutine init_sponge_irreg(rho0,r_cc_loc,r_edge_loc) &
       bind(C, name="init_sponge_irreg")

    double precision, intent(in   ) ::       rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

    call amrex_error("ERROR: irregular sponge not supported for xrb_mixed")

  end subroutine init_sponge_irreg


  subroutine mk_sponge(lo,hi,sponge,s_lo,s_hi,dx,dt) bind(C, name="mk_sponge")

    integer        , intent(in   ) :: lo(3),hi(3),s_lo(3),s_hi(3)
    double precision, intent(inout) :: sponge(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: dx(3),dt

    integer         :: j
    double precision :: y

    if (amrex_spacedim .ne. 2) call amrex_error("ERROR: sponge only supported for 2d in xrb_mixed")

    sponge(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ONE

    if (xrb_use_bottom_sponge) then
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF)*dx(2)

          if(y .le. botsponge_lo_r) then
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = sponge_min
          else if(y .le. botsponge_hi_r) then
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = -HALF*(ONE-sponge_min) &
                  * cos(M_PI*(y-botsponge_lo_r)/(botsponge_hi_r-botsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else if(y .le. topsponge_lo_r) then
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = ONE
          else if (y .le. topsponge_hi_r) then
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = HALF*(ONE-sponge_min) &
                  * cos(M_PI*(y-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = sponge_min
          end if

       enddo
    else

       do j = lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+HALF)*dx(2)
          if(y .le. topsponge_lo_r) then
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = ONE
          else if (y .le. topsponge_hi_r) then
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = HALF*(ONE-sponge_min) &
                  * cos(M_PI*(y-topsponge_lo_r)/(topsponge_hi_r-topsponge_lo_r)) &
                  + HALF*(ONE+sponge_min)
          else
             sponge(lo(1):hi(1),j,lo(3):hi(3)) = sponge_min
          end if
       end do

    endif

  end subroutine mk_sponge

end module sponge_module
