#include "AMReX_BC_TYPES.H"

module mkutrans_module


#ifndef AMREX_USE_CUDA
  use amrex_error_module
#endif
  use amrex_constants_module
  use slope_module
  use ppm_module
  use meth_params_module, only: rel_eps, ppm_type, spherical
  use base_state_geometry_module, only: nr_fine, max_radial_level

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 1)
  subroutine mkutrans_1d(lev, domlo, domhi, lo, hi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       w0,dx,dt,adv_bc,phys_bc) bind(C,name="mkutrans_1d")

    integer         , intent(in   ) :: lev, domlo(1), domhi(1), lo(1), hi(1)
    integer         , intent(in   ) :: ut_lo(1), ut_hi(1), nc_ut
    integer, value  , intent(in   ) :: ng_ut
    integer         , intent(in   ) :: uf_lo(1), uf_hi(1), nc_uf
    integer, value  , intent(in   ) :: ng_uf
    integer         , intent(in   ) :: uu_lo(1), uu_hi(1)
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1))
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(1), dt
    integer         , intent(in   ) :: adv_bc(1,2,1), phys_bc(1,2) ! dim, lohi, (comp)

    double precision, pointer :: slopex(:,:)

    double precision, pointer :: Ip(:)
    double precision, pointer :: Im(:)

    double precision, pointer :: ulx(:),urx(:)

    double precision hx,dt2,uavg

    integer :: i,is,ie

    logical :: test

    allocate(slopex(lo(1)-1:hi(1)+1,1,1))

    allocate(ulx(lo(1):hi(1)+1))
    allocate(urx(lo(1):hi(1)+1))

    allocate(Ip(lo(1)-1:hi(1)+1))
    allocate(Im(lo(1)-1:hi(1)+1))

    is = lo(1)
    ie = hi(1)

    dt2 = HALF*dt

    hx = dx(1)

    if (ppm_type .eq. 0) then
       call slopex_1d(utilde(:,1:1),slopex,domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,1:1))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_1d(utilde(:,1),ng_ut,ufull(:,1),ng_uf,Ip,Im, &
            domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
    end if


    !******************************************************************
    ! create utrans
    !******************************************************************

    if (ppm_type .eq. 0) then
       do i=lo(1),hi(1)+1
          ! extrapolate to edges
          ulx(i) = utilde(i-1,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,1)))*slopex(i-1,1)
          urx(i) = utilde(i  ,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,1)))*slopex(i  ,1)
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do i=lo(1),hi(1)+1
          ! extrapolate to edges
          ulx(i) = Ip(i-1)
          urx(i) = Im(i  )
       end do
    end if

    ! impose lo i side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          ulx(lo(1)) = utilde(lo(1)-1,1)
          urx(lo(1)) = utilde(lo(1)-1,1)
       case (SlipWall, NoSlipWall, Symmetry)
          ulx(lo(1)) = ZERO
          urx(lo(1)) = ZERO
       case (Outflow)
          ulx(lo(1)) = min(urx(lo(1)),ZERO)
          urx(lo(1)) = ulx(lo(1))
       case (Interior)
       case  default
#ifndef AMREX_USE_CUDA
          call amrex_error("mkutrans_1d: invalid boundary type phys_bc(1,1)")
#endif
       end select
    end if

    ! impose hi i side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          ulx(hi(1)+1) = utilde(hi(1)+1,1)
          urx(hi(1)+1) = utilde(hi(1)+1,1)
       case (SlipWall, NoSlipWall, Symmetry)
          ulx(hi(1)+1) = ZERO
          urx(hi(1)+1) = ZERO
       case (Outflow)
          ulx(hi(1)+1) = max(ulx(hi(1)+1),ZERO)
          urx(hi(1)+1) = ulx(hi(1)+1)
       case (Interior)
       case  default
#ifndef AMREX_USE_CUDA
          call amrex_error("mkutrans_1d: invalid boundary type phys_bc(1,2)")
#endif
       end select
    end if

    do i=lo(1),hi(1)+1
       ! solve Riemann problem using full velocity
       uavg = HALF*(ulx(i)+urx(i))
       test = ((ulx(i)+w0(lev,i) .le. ZERO .and. urx(i)+w0(lev,i) .ge. ZERO) .or. &
            (abs(ulx(i)+urx(i)+TWO*w0(lev,i)) .lt. rel_eps))
       utrans(i) = merge(ulx(i),urx(i),uavg+w0(lev,i) .gt. ZERO)
       utrans(i) = merge(ZERO,utrans(i),test)
    end do

    deallocate(slopex)
    deallocate(Ip)
    deallocate(Im)
    deallocate(ulx)
    deallocate(urx)

  end subroutine mkutrans_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine mkutrans_2d(lo, hi, lev, domlo, domhi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       Ip, ip_lo, ip_hi, &
       Im, im_lo, im_hi, &
       w0,dx,dt,adv_bc,phys_bc) bind(C,name="mkutrans_2d")

    implicit none

    integer, value, intent(in   ) :: lev
    integer, intent(in) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
    integer, value  , intent(in   ) :: nc_ut, ng_ut
    integer         , intent(in   ) :: uf_lo(3), uf_hi(3)
    integer, value  , intent(in   ) :: nc_uf, ng_uf
    integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
    integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
    integer         , intent(in   ) :: ip_lo(3), ip_hi(3)
    integer         , intent(in   ) :: im_lo(3), im_hi(3)
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
    double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),1:2)
    double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),1:2)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3)
    double precision, value, intent(in   ) :: dt
    integer         , intent(in   ) :: adv_bc(2,2,2), phys_bc(2,2) ! dim, lohi, (comp)

    double precision hx,hy,dt2,uavg
    double precision :: ul, ur, vl, vr

    integer :: i,j,k

    logical :: test

    !$gpu

    ! for ppm 0, slopex = Ip, slopey = Im

    dt2 = HALF*dt

    hx = dx(1)
    hy = dx(2)

    k = lo(3)

    if (ppm_type .eq. 0) then
       ! call slopex_2d(utilde(:,:,1:1),slopex,domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,1:1))
       ! call slopey_2d(utilde(:,:,2:2),slopey,domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,2:2))

       ! call slopex_2d(utilde(:,:,k,1:1),Ip(:,:,k,1:1),domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,1:1))
       ! call slopey_2d(utilde(:,:,k,2:2),Im(:,:,k,1:1),domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,2:2))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

       call ppm_2d(lo,hi,utilde,ut_lo,ut_hi,nc_ut, &
            ufull(:,:,:,1),uf_lo,uf_hi,ufull(:,:,:,2),uf_lo,uf_hi, &
            Ip,ip_lo,ip_hi,Im,im_lo,im_hi,domlo,domhi,adv_bc,dx,dt,.false.,1,1)

    end if

    !******************************************************************
    ! create utrans
    !******************************************************************

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          if (ppm_type .eq. 0) then
             ! extrapolate to edges
             ! ul = utilde(i-1,j,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,1)))*slopex(i-1,j,1)
             ! ur = utilde(i  ,j,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,1)))*slopex(i  ,j,1)

             ul = utilde(i-1,j,k,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,k,1)))*Ip(i-1,j,k,1)
             ur = utilde(i  ,j,k,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,k,1)))*Ip(i  ,j,k,1)
          else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
             ! extrapolate to edges
             ul = Ip(i-1,j,k,1)
             ur = Im(i  ,j,k,1)
          end if
          ! end do

          ! impose lo i side bc's
          if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
             select case(phys_bc(1,1))
             case (Inflow)
                ul = utilde(lo(1)-1,j,k,1)
                ur = utilde(lo(1)-1,j,k,1)
             case (SlipWall, NoSlipWall, Symmetry)
                ul = ZERO
                ur = ZERO
             case (Outflow)
                ul = min(ur,ZERO)
                ur = ul
             case (Interior)
             case  default
#ifndef AMREX_USE_CUDA
                call amrex_error("mkutrans_2d: invalid boundary type phys_bc(1,1)")
#endif
             end select
          end if

          ! impose hi i side bc's
          if (i .eq. hi(1)+1 .and. hi(1) .eq. domhi(1)) then
             select case(phys_bc(1,2))
             case (Inflow)
                ul = utilde(hi(1)+1,j,k,1)
                ur = utilde(hi(1)+1,j,k,1)
             case (SlipWall, NoSlipWall, Symmetry)
                ul = ZERO
                ur = ZERO
             case (Outflow)
                ul = max(ul,ZERO)
                ur = ul
             case (Interior)
             case  default
#ifndef AMREX_USE_CUDA
                call amrex_error("mkutrans_2d: invalid boundary type phys_bc(1,2)")
#endif
             end select
          end if

          ! solve Riemann problem using full velocity
          uavg = HALF*(ul+ur)
          test = ((ul .le. ZERO .and. ur .ge. ZERO) .or. &
               (abs(ul+ur) .lt. rel_eps))
          utrans(i,j,k) = merge(ul,ur,uavg .gt. ZERO)
          utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
       end do
    end do

    !******************************************************************
    ! create vtrans
    !******************************************************************
    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

       call ppm_2d(lo,hi,utilde,ut_lo,ut_hi,nc_ut, &
            ufull(:,:,:,1),uf_lo,uf_hi,&
            ufull(:,:,:,2),uf_lo,uf_hi, &
            Ip,ip_lo,ip_hi,Im,im_lo,im_hi,domlo,domhi,adv_bc,dx,dt,.false.,2,2)

    end if

    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          if (ppm_type .eq. 0) then
             ! ! extrapolate to edges
             ! vl = utilde(i,j-1,2) + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,2)))*slopey(i,j-1,1)
             ! vr = utilde(i,j  ,2) - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,2)))*slopey(i,j  ,1)

             ! extrapolate to edges
             vl = utilde(i,j-1,k,2) + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,k,2)))*Im(i,j-1,k,1)
             vr = utilde(i,j  ,k,2) - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,k,2)))*Im(i,j  ,k,1)

          else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
             ! extrapolate to edges
             vl = Ip(i,j-1,k,2)
             vr = Im(i,j  ,k,2)
          end if

          ! impose lo side bc's
          if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
             select case(phys_bc(2,1))
             case (Inflow)
                vl = utilde(i,lo(2)-1,k,2)
                vr = utilde(i,lo(2)-1,k,2)
             case (SlipWall, NoSlipWall, Symmetry)
                vr = ZERO
                vr = ZERO
             case (Outflow)
                vl = min(vr,ZERO)
                vr = vl
             case (Interior)
             case  default
#ifndef AMREX_USE_CUDA
                call amrex_error("mkutrans_2d: invalid boundary type phys_bc(2,1)")
#endif
             end select
          end if

          ! impose hi side bc's
          if (j .eq. hi(2)+1 .and. hi(2) .eq. domhi(2)) then
             select case(phys_bc(2,2))
             case (Inflow)
                vl = utilde(i,hi(2)+1,k,2)
                vr = utilde(i,hi(2)+1,k,2)
             case (SlipWall, NoSlipWall, Symmetry)
                vl = ZERO
                vr = ZERO
             case (Outflow)
                vl = max(vl,ZERO)
                vr = vl
             case (Interior)
             case  default
#ifndef AMREX_USE_CUDA
                call amrex_error("mkutrans_2d: invalid boundary type phys_bc(2,2)")
#endif
             end select
          end if

          ! solve Riemann problem using full velocity
          uavg = HALF*(vl+vr)
          test = ((vl+w0(lev,j) .le. ZERO .and. vr+w0(lev,j) .ge. ZERO) .or. &
               (abs(vl+vr+TWO*w0(lev,j)) .lt. rel_eps))
          vtrans(i,j,k) = merge(vl,vr,uavg+w0(lev,j) .gt. ZERO)
          vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
       enddo
    enddo

  end subroutine mkutrans_2d
#endif

#if (AMREX_SPACEDIM == 3)
  subroutine mkutrans_3d(lev, domlo, domhi, lo, hi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       wtrans, uw_lo, uw_hi, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       w0, &
       dx,dt,adv_bc,phys_bc) bind(C,name="mkutrans_3d")

    integer         , intent(in   ) :: lev, domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: ut_lo(3), ut_hi(3), nc_ut
    integer, value  , intent(in   ) :: ng_ut
    integer         , intent(in   ) :: uf_lo(3), uf_hi(3), nc_uf
    integer, value  , intent(in   ) :: ng_uf
    integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
    integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
    integer         , intent(in   ) :: uw_lo(3), uw_hi(3)
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
    double precision, intent(inout) :: wtrans(uw_lo(1):uw_hi(1),uw_lo(2):uw_hi(2),uw_lo(3):uw_hi(3))
    double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: adv_bc(3,2,3), phys_bc(3,2) ! dim, lohi, (comp)

    double precision, pointer :: slopex(:,:,:,:)
    double precision, pointer :: slopey(:,:,:,:)
    double precision, pointer :: slopez(:,:,:,:)

    double precision, pointer :: Ip(:,:,:,:)
    double precision, pointer :: Im(:,:,:,:)

    double precision hx,hy,hz,dt2,uavg

    logical :: test

    integer :: i,j,k,is,js,ks,ie,je,ke

    double precision, pointer:: ulx(:,:,:),urx(:,:,:)
    double precision, pointer:: vly(:,:,:),vry(:,:,:)
    double precision, pointer:: wlz(:,:,:),wrz(:,:,:)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))

    allocate(Ip(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Im(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)

    dt2 = HALF*dt

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    if (ppm_type .eq. 0) then
       do k = lo(3)-1,hi(3)+1
          call slopex_2d(utilde(:,:,k,1:1),slopex(:,:,k,:),domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,1:1))
          call slopey_2d(utilde(:,:,k,2:2),slopey(:,:,k,:),domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,2:2))
       end do
       call slopez_3d(utilde(:,:,:,3:3),slopez,domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,3:3))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(utilde(:,:,:,1),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
    end if

    !******************************************************************
    ! create utrans
    !******************************************************************

    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             if (ppm_type .eq. 0) then
                !$OMP PARALLEL DO PRIVATE(i,j,k)
                ! extrapolate to edges
                ulx(i,j,k) = utilde(i-1,j,k,1) &
                     + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,k,1)))*slopex(i-1,j,k,1)
                urx(i,j,k) = utilde(i  ,j,k,1) &
                     - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,k,1)))*slopex(i  ,j,k,1)
                !$OMP END PARALLEL DO
             else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                ! extrapolate to edges
                ulx(i,j,k) = Ip(i-1,j,k,1)
                urx(i,j,k) = Im(i  ,j,k,1)
             end if

             ! impose lo side bc's
             if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                select case(phys_bc(1,1))
                case (Inflow)
                   ulx(lo(1),j,k) = utilde(lo(1)-1,j,k,1)
                   urx(lo(1),j,k) = utilde(lo(1)-1,j,k,1)
                case (SlipWall, NoSlipWall, Symmetry)
                   ulx(lo(1),j,k) = ZERO
                   urx(lo(1),j,k) = ZERO
                case (Outflow)
                   ulx(lo(1),j,k) = min(urx(lo(1),j,k),ZERO)
                   urx(lo(1),j,k) = ulx(lo(1),j,k)
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_3d: invalid boundary type phys_bc(1,1)")
#endif
                end select
             end if

             ! impose hi side bc's
             if (i .eq. hi(1) .and. hi(1) .eq. domhi(1)) then
                select case(phys_bc(1,2))
                case (Inflow)
                   ulx(hi(1)+1,j,k) = utilde(hi(1)+1,j,k,1)
                   urx(hi(1)+1,j,k) = utilde(hi(1)+1,j,k,1)
                case (SlipWall, NoSlipWall, Symmetry)
                   ulx(hi(1)+1,j,k) = ZERO
                   urx(hi(1)+1,j,k) = ZERO
                case (Outflow)
                   ulx(hi(1)+1,j,k) = max(ulx(hi(1)+1,j,k),ZERO)
                   urx(hi(1)+1,j,k) = ulx(hi(1)+1,j,k)
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_3d: invalid boundary type phys_bc(1,2)")
#endif
                end select
             end if

             if (spherical .eq. 1) then

                ! solve Riemann problem using full velocity
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k)+w0macx(i,j,k) .le. ZERO .and. &
                     urx(i,j,k)+w0macx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)+TWO*w0macx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg+w0macx(i,j,k) .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)

             else

                ! solve Riemann problem using full velocity
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k) .le. ZERO .and. urx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)

             end if
          end do
       end do
    end do

    deallocate(slopex)
    deallocate(ulx)
    deallocate(urx)

    !******************************************************************
    ! create vtrans
    !******************************************************************

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(utilde(:,:,:,2),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
    end if

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             if (ppm_type .eq. 0) then
                !$OMP PARALLEL DO PRIVATE(i,j,k)
                ! extrapolate to edges
                vly(i,j,k) = utilde(i,j-1,k,2) &
                     + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,k,2)))*slopey(i,j-1,k,1)
                vry(i,j,k) = utilde(i,j  ,k,2) &
                     - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,k,2)))*slopey(i,j  ,k,1)

             else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

                ! extrapolate to edges
                vly(i,j,k) = Ip(i,j-1,k,2)
                vry(i,j,k) = Im(i,j  ,k,2)

             end if

             ! impose lo side bc's
             if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                select case(phys_bc(2,1))
                case (Inflow)
                   vly(i,lo(2),k) = utilde(i,lo(2)-1,k,2)
                   vry(i,lo(2),k) = utilde(i,lo(2)-1,k,2)
                case (SlipWall, NoSlipWall, Symmetry)
                   vly(i,lo(2),k) = ZERO
                   vry(i,lo(2),k) = ZERO
                case (Outflow)
                   vly(i,lo(2),k) = min(vry(i,lo(2),k),ZERO)
                   vry(i,lo(2),k) = vly(i,lo(2),k)
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_3d: invalid boundary type phys_bc(2,1)")
#endif
                end select
             end if

             ! impose hi side bc's
             if (j .eq. hi(2) .and. hi(2) .eq. domhi(2)) then
                select case(phys_bc(2,2))
                case (Inflow)
                   vly(i,hi(2)+1,k) = utilde(i,hi(2)+1,k,2)
                   vry(i,hi(2)+1,k) = utilde(i,hi(2)+1,k,2)
                case (SlipWall, NoSlipWall, Symmetry)
                   vly(i,hi(2)+1,k) = ZERO
                   vry(i,hi(2)+1,k) = ZERO
                case (Outflow)
                   vly(i,hi(2)+1,k) = max(vly(i,hi(2)+1,k),ZERO)
                   vry(i,hi(2)+1,k) = vly(i,hi(2)+1,k)
                case (Interior)
                case  default
#ifndef AMREX_USE_CUDA
                   call amrex_error("mkutrans_3d: invalid boundary type phys_bc(2,2)")
#endif
                end select
             end if

             if (spherical .eq. 1) then
                ! solve Riemann problem using full velocity
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k)+w0macy(i,j,k) .le. ZERO .and. &
                     vry(i,j,k)+w0macy(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)+TWO*w0macy(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg+w0macy(i,j,k) .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             else
                ! solve Riemann problem using full velocity
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k) .le. ZERO .and. vry(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             end if

          enddo
       enddo
    enddo

    deallocate(slopey)
    deallocate(vly)
    deallocate(vry)

    !******************************************************************
    ! create wtrans
    !******************************************************************

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(utilde(:,:,:,3),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,3),dx,dt,.false.)
    end if

    allocate(wlz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(wrz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                ! extrapolate to edges
                wlz(i,j,k) = utilde(i,j,k-1,3) &
                     + (HALF-(dt2/hz)*max(ZERO,ufull(i,j,k-1,3)))*slopez(i,j,k-1,1)
                wrz(i,j,k) = utilde(i,j,k  ,3) &
                     - (HALF+(dt2/hz)*min(ZERO,ufull(i,j,k  ,3)))*slopez(i,j,k  ,1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                ! extrapolate to edges
                wlz(i,j,k) = Ip(i,j,k-1,3)
                wrz(i,j,k) = Im(i,j,k  ,3)
             end do
          end do
       end do
    end if

    deallocate(slopez)
    deallocate(Ip)
    deallocate(Im)

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       select case(phys_bc(3,1))
       case (Inflow)
          wlz(lo(1):hi(1),lo(2):hi(2),lo(3)) = utilde(lo(1):hi(1),lo(2):hi(2),lo(3)-1,3)
          wrz(lo(1):hi(1),lo(2):hi(2),lo(3)) = utilde(lo(1):hi(1),lo(2):hi(2),lo(3)-1,3)
       case (SlipWall, NoSlipWall, Symmetry)
          wlz(lo(1):hi(1),lo(2):hi(2),lo(3)) = ZERO
          wrz(lo(1):hi(1),lo(2):hi(2),lo(3)) = ZERO
       case (Outflow)
          wlz(lo(1):hi(1),lo(2):hi(2),lo(3)) = min(wrz(lo(1):hi(1),lo(2):hi(2),lo(3)),ZERO)
          wrz(lo(1):hi(1),lo(2):hi(2),lo(3)) = wlz(lo(1):hi(1),lo(2):hi(2),lo(3))
       case (Interior)
       case  default
#ifndef AMREX_USE_CUDA
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(3,1)")
#endif
       end select
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       select case(phys_bc(3,2))
       case (Inflow)
          wlz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = utilde(lo(1):hi(1),lo(2):hi(2),hi(3)+1,3)
          wrz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = utilde(lo(1):hi(1),lo(2):hi(2),hi(3)+1,3)
       case (SlipWall, NoSlipWall, Symmetry)
          wlz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = ZERO
          wrz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = ZERO
       case (Outflow)
          wlz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = max(wlz(lo(1):hi(1),lo(2):hi(2),hi(3)+1),ZERO)
          wrz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = wlz(lo(1):hi(1),lo(2):hi(2),hi(3)+1)
       case (Interior)
       case  default
#ifndef AMREX_USE_CUDA
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(3,2)")
#endif
       end select
    end if

    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                ! solve Riemann problem using full velocity
                uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
                test = ((wlz(i,j,k)+w0macz(i,j,k) .le. ZERO .and. &
                     wrz(i,j,k)+w0macz(i,j,k) .ge. ZERO) .or. &
                     (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0macz(i,j,k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0macz(i,j,k) .gt. ZERO)
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                ! solve Riemann problem using full velocity
                uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
                test = ((wlz(i,j,k)+w0(lev,k).le.ZERO .and. wrz(i,j,k)+w0(lev,k).ge.ZERO) .or. &
                     (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0(lev,k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0(lev,k) .gt. ZERO)
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    deallocate(wlz)
    deallocate(wrz)

  end subroutine mkutrans_3d
#endif

end module mkutrans_module
