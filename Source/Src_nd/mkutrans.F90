#include "AMReX_BC_TYPES.H"

module mkutrans_module


#ifndef AMREX_USE_CUDA
  use amrex_error_module
#endif
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use slope_module
  use ppm_module
  use meth_params_module, only: rel_eps, ppm_type, spherical
  use base_state_geometry_module, only: nr_fine, max_radial_level

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)
  subroutine mkutrans_2d(lo, hi, idir, domlo, domhi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       Ip, ip_lo, ip_hi, &
       Im, im_lo, im_hi, &
       w0_cart, w_lo, w_hi, &
       dx,dt,adv_bc,phys_bc) bind(C,name="mkutrans_2d")

    implicit none

    integer, value, intent(in   ) :: idir
    integer, intent(in) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
    integer, value  , intent(in   ) :: nc_ut, ng_ut
    integer         , intent(in   ) :: uf_lo(3), uf_hi(3)
    integer, value  , intent(in   ) :: nc_uf, ng_uf
    integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
    integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
    integer         , intent(in   ) :: ip_lo(3), ip_hi(3)
    integer         , intent(in   ) :: im_lo(3), im_hi(3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
    double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),1:2)
    double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),1:2)
    double precision, intent(in   ) :: w0_cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),AMREX_SPACEDIM)
    double precision, intent(in   ) :: dx(3)
    double precision, value, intent(in   ) :: dt
    integer         , intent(in   ) :: adv_bc(2,2,2), phys_bc(2,2) ! dim, lohi, (comp)

    double precision hx,hy,dt2,uavg
    double precision :: ulx, urx, vly, vry

    integer :: i,j,k

    logical :: test

    !$gpu

    ! for ppm 0, slopex = Ip, slopey = Im

    dt2 = HALF*dt

    hx = dx(1)
    hy = dx(2)

    k = lo(3)

    if (idir == 1) then



       !******************************************************************
       ! create utrans
       !******************************************************************

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (ppm_type .eq. 0) then

                ulx = utilde(i-1,j,k,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,k,1)))*Ip(i-1,j,k,1)
                urx = utilde(i  ,j,k,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,k,1)))*Ip(i  ,j,k,1)

             else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                ! extrapolate to edges
                ulx = Ip(i-1,j,k,1)
                urx = Im(i  ,j,k,1)
             end if

             ! impose lo i side bc's
             if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                select case(phys_bc(1,1))
                case (Inflow)
                   ulx = utilde(i-1,j,k,1)
                   urx = utilde(i-1,j,k,1)
                case (SlipWall, NoSlipWall, Symmetry)
                   ulx = ZERO
                   urx = ZERO
                case (Outflow)
                   ulx = min(urx,ZERO)
                   urx = ulx
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_2d: invalid boundary type phys_bc(1,1)")
#endif
                end select
             end if

             ! impose hi i side bc's
             if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                select case(phys_bc(1,2))
                case (Inflow)
                   ulx = utilde(i,j,k,1)
                   urx = utilde(i,j,k,1)
                case (SlipWall, NoSlipWall, Symmetry)
                   ulx = ZERO
                   urx = ZERO
                case (Outflow)
                   ulx = max(ulx,ZERO)
                   urx = ulx
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_2d: invalid boundary type phys_bc(1,2)")
#endif
                end select
             end if

             ! solve Riemann problem using full velocity
             uavg = HALF*(ulx+urx)
             test = ((ulx .le. ZERO .and. urx .ge. ZERO) .or. &
                  (abs(ulx+urx) .lt. rel_eps))
             utrans(i,j,k) = merge(ulx,urx,uavg .gt. ZERO)
             utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
          end do
       end do

    else ! idir == 2

       !******************************************************************
       ! create vtrans
       !******************************************************************

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (ppm_type .eq. 0) then
                ! ! extrapolate to edges
                vly = utilde(i,j-1,k,2) + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,k,2)))*Im(i,j-1,k,1)
                vry = utilde(i,j  ,k,2) - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,k,2)))*Im(i,j  ,k,1)

             else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                ! extrapolate to edges
                vly = Ip(i,j-1,k,2)
                vry = Im(i,j  ,k,2)
             end if

             ! impose lo side bc's
             if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                select case(phys_bc(2,1))
                case (Inflow)
                   vly = utilde(i,j-1,k,2)
                   vry = utilde(i,j-1,k,2)
                case (SlipWall, NoSlipWall, Symmetry)
                   vry = ZERO
                   vry = ZERO
                case (Outflow)
                   vly = min(vry,ZERO)
                   vry = vly
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_2d: invalid boundary type phys_bc(2,1)")
#endif
                end select
             end if

             ! impose hi side bc's
             if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                select case(phys_bc(2,2))
                case (Inflow)
                   vly = utilde(i,j,k,2)
                   vry = utilde(i,j,k,2)
                case (SlipWall, NoSlipWall, Symmetry)
                   vly = ZERO
                   vry = ZERO
                case (Outflow)
                   vly = max(vly,ZERO)
                   vry = vly
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_2d: invalid boundary type phys_bc(2,2)")
#endif
                end select
             end if

             ! solve Riemann problem using full velocity
             uavg = HALF*(vly+vry)
             test = ((vly+w0_cart(i,j,k,AMREX_SPACEDIM) .le. ZERO .and. vry+w0_cart(i,j,k,AMREX_SPACEDIM) .ge. ZERO) .or. &
                  (abs(vly+vry+TWO*w0_cart(i,j,k,AMREX_SPACEDIM)) .lt. rel_eps))
             vtrans(i,j,k) = merge(vly,vry,uavg+w0_cart(i,j,k,AMREX_SPACEDIM) .gt. ZERO)
             vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
          enddo
       enddo

    endif

  end subroutine mkutrans_2d
#endif

#if (AMREX_SPACEDIM == 3)
  subroutine mkutrans_3d(lo, hi, idir, domlo, domhi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       wtrans, uw_lo, uw_hi, &
       Ip, ip_lo, ip_hi, &
       Im, im_lo, im_hi, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       w0_cart, w_lo, w_hi, &
       dx,dt,adv_bc,phys_bc) bind(C,name="mkutrans_3d")

    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
    integer, value  , intent(in   ) :: idir, ng_ut, nc_ut
    integer         , intent(in   ) :: uf_lo(3), uf_hi(3)
    integer, value  , intent(in   ) :: ng_uf, nc_uf
    integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
    integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
    integer         , intent(in   ) :: uw_lo(3), uw_hi(3)
    integer         , intent(in   ) :: ip_lo(3), ip_hi(3)
    integer         , intent(in   ) :: im_lo(3), im_hi(3)
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
    double precision, intent(inout) :: wtrans(uw_lo(1):uw_hi(1),uw_lo(2):uw_hi(2),uw_lo(3):uw_hi(3))
    double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),AMREX_SPACEDIM)
    double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),AMREX_SPACEDIM)
    double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    double precision, intent(in   ) :: w0_cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),AMREX_SPACEDIM)
    double precision, intent(in   ) :: dx(3)
    double precision, value, intent(in   ) :: dt
    integer         , intent(in   ) :: adv_bc(AMREX_SPACEDIM,2,AMREX_SPACEDIM), phys_bc(AMREX_SPACEDIM,2) ! dim, lohi, (comp)

    double precision hx,hy,hz,dt2,uavg

    logical :: test

    integer :: i,j,k

    double precision :: ulx,urx,vly,vry
    double precision :: wlz,wrz

    !$gpu

    dt2 = HALF*dt

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    if (idir == 1) then

       !******************************************************************
       ! create utrans
       !******************************************************************

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                if (ppm_type .eq. 0) then
                   ! extrapolate to edges
                   ulx = utilde(i-1,j,k,1) &
                        + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,k,1)))*Ip(i-1,j,k,1)
                   urx = utilde(i  ,j,k,1) &
                        - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,k,1)))*Ip(i  ,j,k,1)
                else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                   ! extrapolate to edges
                   ulx = Ip(i-1,j,k,1)
                   urx = Im(i  ,j,k,1)
                end if

                ! impose lo side bc's
                if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                   select case(phys_bc(1,1))
                   case (Inflow)
                      ulx = utilde(i-1,j,k,1)
                      urx = utilde(i-1,j,k,1)
                   case (SlipWall, NoSlipWall, Symmetry)
                      ulx = ZERO
                      urx = ZERO
                   case (Outflow)
                      ulx = min(urx,ZERO)
                      urx = ulx
                   case (Interior)
                   case  default
#ifndef AMREX_USE_CUDA
                      call amrex_error("mkutrans_3d: invalid boundary type phys_bc(1,1)")
#endif
                   end select
                end if

                ! impose hi side bc's
                if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                   select case(phys_bc(1,2))
                   case (Inflow)
                      ulx = utilde(i+1,j,k,1)
                      urx = utilde(i+1,j,k,1)
                   case (SlipWall, NoSlipWall, Symmetry)
                      ulx = ZERO
                      urx = ZERO
                   case (Outflow)
                      ulx = max(ulx,ZERO)
                      urx = ulx
                   case (Interior)
                   case  default
#ifndef AMREX_USE_CUDA
                      call amrex_error("mkutrans_3d: invalid boundary type phys_bc(1,2)")
#endif
                   end select
                end if

                if (spherical .eq. 1) then

                   ! solve Riemann problem using full velocity
                   uavg = HALF*(ulx+urx)
                   test = ((ulx+w0macx(i,j,k) .le. ZERO .and. &
                        urx+w0macx(i,j,k) .ge. ZERO) .or. &
                        (abs(ulx+urx+TWO*w0macx(i,j,k)) .lt. rel_eps))
                   utrans(i,j,k) = merge(ulx,urx,uavg+w0macx(i,j,k) .gt. ZERO)
                   utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)

                else

                   ! solve Riemann problem using full velocity
                   uavg = HALF*(ulx+urx)
                   test = ((ulx .le. ZERO .and. urx .ge. ZERO) .or. &
                        (abs(ulx+urx) .lt. rel_eps))
                   utrans(i,j,k) = merge(ulx,urx,uavg .gt. ZERO)
                   utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)

                end if
             end do
          end do
       end do

    elseif (idir == 2) then

       !******************************************************************
       ! create vtrans
       !******************************************************************

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                if (ppm_type .eq. 0) then
                   ! extrapolate to edges
                   vly = utilde(i,j-1,k,2) &
                        + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,k,2)))*Im(i,j-1,k,1)
                   vry = utilde(i,j  ,k,2) &
                        - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,k,2)))*Im(i,j  ,k,1)

                else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

                   ! extrapolate to edges
                   vly = Ip(i,j-1,k,2)
                   vry = Im(i,j  ,k,2)

                end if

                ! impose lo side bc's
                if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                   select case(phys_bc(2,1))
                   case (Inflow)
                      vly = utilde(i,j-1,k,2)
                      vry = utilde(i,j-1,k,2)
                   case (SlipWall, NoSlipWall, Symmetry)
                      vly = ZERO
                      vry = ZERO
                   case (Outflow)
                      vly = min(vry,ZERO)
                      vry = vly
                   case (Interior)
                   case  default
#ifndef AMREX_USE_CUDA
                      call amrex_error("mkutrans_3d: invalid boundary type phys_bc(2,1)")
#endif
                   end select
                end if

                ! impose hi side bc's
                if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                   select case(phys_bc(2,2))
                   case (Inflow)
                      vly = utilde(i,j+1,k,2)
                      vry = utilde(i,j+1,k,2)
                   case (SlipWall, NoSlipWall, Symmetry)
                      vly = ZERO
                      vry = ZERO
                   case (Outflow)
                      vly = max(vly,ZERO)
                      vry = vly
                   case (Interior)
                   case  default
#ifndef AMREX_USE_CUDA
                      call amrex_error("mkutrans_3d: invalid boundary type phys_bc(2,2)")
#endif
                   end select
                end if

                if (spherical .eq. 1) then
                   ! solve Riemann problem using full velocity
                   uavg = HALF*(vly+vry)
                   test = ((vly+w0macy(i,j,k) .le. ZERO .and. &
                        vry+w0macy(i,j,k) .ge. ZERO) .or. &
                        (abs(vly+vry+TWO*w0macy(i,j,k)) .lt. rel_eps))
                   vtrans(i,j,k) = merge(vly,vry,uavg+w0macy(i,j,k) .gt. ZERO)
                   vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
                else
                   ! solve Riemann problem using full velocity
                   uavg = HALF*(vly+vry)
                   test = ((vly .le. ZERO .and. vry .ge. ZERO) .or. &
                        (abs(vly+vry) .lt. rel_eps))
                   vtrans(i,j,k) = merge(vly,vry,uavg .gt. ZERO)
                   vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
                end if

             enddo
          enddo
       enddo

    else

       !******************************************************************
       ! create wtrans
       !******************************************************************

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                if (ppm_type .eq. 0) then
                   ! extrapolate to edges
                   wlz = utilde(i,j,k-1,3) &
                        + (HALF-(dt2/hz)*max(ZERO,ufull(i,j,k-1,3)))*Im(i,j,k-1,1)
                   wrz = utilde(i,j,k  ,3) &
                        - (HALF+(dt2/hz)*min(ZERO,ufull(i,j,k  ,3)))*Im(i,j,k  ,1)
                else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                   ! extrapolate to edges
                   wlz = Ip(i,j,k-1,3)
                   wrz = Im(i,j,k  ,3)
                end if

                ! impose lo side bc's
                if (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                   select case(phys_bc(3,1))
                   case (Inflow)
                      wlz = utilde(i,j,k-1,3)
                      wrz = utilde(i,j,k-1,3)
                   case (SlipWall, NoSlipWall, Symmetry)
                      wlz = ZERO
                      wrz = ZERO
                   case (Outflow)
                      wlz = min(wrz,ZERO)
                      wrz = wlz
                   case (Interior)
                   case  default
#ifndef AMREX_USE_CUDA
                      call amrex_error("mkutrans_3d: invalid boundary type phys_bc(3,1)")
#endif
                   end select
                end if

                ! impose hi side bc's
                if (k .eq. hi(3) .and. hi(3)-1 .eq. domhi(3)) then
                   select case(phys_bc(3,2))
                   case (Inflow)
                      wlz = utilde(i,j,k+1,3)
                      wrz = utilde(i,j,k+1,3)
                   case (SlipWall, NoSlipWall, Symmetry)
                      wlz = ZERO
                      wrz = ZERO
                   case (Outflow)
                      wlz = max(wlz,ZERO)
                      wrz = wlz
                   case (Interior)
                   case  default
#ifndef AMREX_USE_CUDA
                      call amrex_error("mkutrans_3d: invalid boundary type phys_bc(3,2)")
#endif
                   end select
                end if

                if (spherical .eq. 1) then
                   ! solve Riemann problem using full velocity
                   uavg = HALF*(wlz+wrz)
                   test = ((wlz+w0macz(i,j,k) .le. ZERO .and. &
                        wrz+w0macz(i,j,k) .ge. ZERO) .or. &
                        (abs(wlz+wrz+TWO*w0macz(i,j,k)) .lt. rel_eps))
                   wtrans(i,j,k) = merge(wlz,wrz,uavg+w0macz(i,j,k) .gt. ZERO)
                   wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
                else
                   ! solve Riemann problem using full velocity
                   uavg = HALF*(wlz+wrz)
                   test = ((wlz+w0_cart(i,j,k,AMREX_SPACEDIM).le.ZERO .and. wrz+w0_cart(i,j,k,AMREX_SPACEDIM).ge.ZERO) .or. &
                        (abs(wlz+wrz+TWO*w0_cart(i,j,k,AMREX_SPACEDIM)) .lt. rel_eps))
                   wtrans(i,j,k) = merge(wlz,wrz,uavg+w0_cart(i,j,k,AMREX_SPACEDIM) .gt. ZERO)
                   wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
                end if
             enddo
          enddo
       enddo

    endif

  end subroutine mkutrans_3d
#endif

end module mkutrans_module
