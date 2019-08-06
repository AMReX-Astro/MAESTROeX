! velpred is called by advance_premac -- it is used to predict the
! normal velocities to the interfaces.  We don't care about the
! transverse velocities here.  The prediction is done piecewise linear (for now)

#include "AMReX_BC_TYPES.H"

module velpred_module

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use slope_module
  use ppm_module
  use meth_params_module, only: ppm_type, rel_eps, spherical, ppm_trace_forces
  use base_state_geometry_module, only: nr_fine, max_radial_level

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)
  subroutine velpred_2d(lev, domlo, domhi, lo, hi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       umac  , mu_lo, mu_hi, &
       vmac  , mv_lo, mv_hi, &
       force,   f_lo,  f_hi, nc_f, ng_f, &
       w0,dx,dt,adv_bc,phys_bc) bind(C,name="velpred_2d")

    integer         , intent(in   ) :: lev, domlo(2), domhi(2), lo(2), hi(2)
    integer         , intent(in   ) :: ut_lo(2), ut_hi(2), nc_ut
    integer, value,   intent(in   ) :: ng_ut
    integer         , intent(in   ) :: uf_lo(2), uf_hi(2), nc_uf
    integer, value,   intent(in   ) :: ng_uf
    integer         , intent(in   ) :: uu_lo(2), uu_hi(2)
    integer         , intent(in   ) :: uv_lo(2), uv_hi(2)
    integer         , intent(in   ) :: mu_lo(2), mu_hi(2)
    integer         , intent(in   ) :: mv_lo(2), mv_hi(2)
    integer         , intent(in   ) ::  f_lo(2),  f_hi(2), nc_f
    integer, value,   intent(in   ) :: ng_f
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2))
    double precision, intent(inout) :: umac  (mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2))
    double precision, intent(inout) :: vmac  (mv_lo(1):mv_hi(1),mv_lo(2):mv_hi(2))
    double precision, intent(in   ) :: force ( f_lo(1): f_hi(1), f_lo(2): f_hi(2),nc_f)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(2), dt
    integer         , intent(in   ) :: adv_bc(2,2,2), phys_bc(2,2) ! dim, lohi, (comp)

    ! Local variables
    double precision, pointer :: slopex(:,:,:)
    double precision, pointer :: slopey(:,:,:)

    double precision, pointer :: Ipu(:,:,:), Ipfx(:,:,:)
    double precision, pointer :: Imu(:,:,:), Imfx(:,:,:)
    double precision, pointer :: Ipv(:,:,:), Ipfy(:,:,:)
    double precision, pointer :: Imv(:,:,:), Imfy(:,:,:)

    ! these correspond to u_L^x, etc.
    double precision, pointer :: ulx(:,:,:),urx(:,:,:),uimhx(:,:,:)
    double precision, pointer :: uly(:,:,:),ury(:,:,:),uimhy(:,:,:)

    ! these correspond to umac_L, etc.
    double precision, pointer :: umacl(:,:),umacr(:,:)
    double precision, pointer :: vmacl(:,:),vmacr(:,:)

    double precision :: hx, hy, dt2, dt4, uavg, maxu, minu
    double precision :: fl, fr

    integer :: i,j,is,js,ie,je

    logical :: test

    call bl_allocate(slopex,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(slopey,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    call bl_allocate(Ipu,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Imu,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Ipv,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Imv,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    call bl_allocate(Ipfx,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Imfx,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Ipfy,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Imfy,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    call bl_allocate(  ulx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(  urx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(uimhx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    call bl_allocate(  uly,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,1,2)
    call bl_allocate(  ury,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,1,2)
    call bl_allocate(uimhy,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,1,2)

    call bl_allocate(umacl,lo(1),hi(1)+1,lo(2),hi(2))
    call bl_allocate(umacr,lo(1),hi(1)+1,lo(2),hi(2))

    call bl_allocate(vmacl,lo(1),hi(1),lo(2),hi(2)+1)
    call bl_allocate(vmacr,lo(1),hi(1),lo(2),hi(2)+1)

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    if (ppm_type .eq. 0) then
       call slopex_2d(utilde,slopex,domlo,domhi,lo,hi,ng_ut,2,adv_bc)
       call slopey_2d(utilde,slopey,domlo,domhi,lo,hi,ng_ut,2,adv_bc)
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(utilde(:,:,1),ng_ut, &
            ufull(:,:,1),ufull(:,:,2),ng_uf, &
            Ipu,Imu,domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
       call ppm_2d(utilde(:,:,2),ng_ut, &
            ufull(:,:,1),ufull(:,:,2),ng_uf, &
            Ipv,Imv,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)

       ! trace forces, if necessary.  Note by default the ppm routines
       ! will trace each component to each interface in all coordinate
       ! directions, but we really only need the force traced along
       ! its respective dimension.  This should be simplified later.
       if (ppm_trace_forces .eq. 1) then
          call ppm_2d(force(:,:,1),ng_f, &
               ufull(:,:,1),ufull(:,:,2),ng_uf, &
               Ipfx,Imfx,domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
          call ppm_2d(force(:,:,2),ng_f, &
               ufull(:,:,1),ufull(:,:,2),ng_uf, &
               Ipfy,Imfy,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
       endif
    end if

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    if (ppm_type .eq. 0) then
       do j=js-1,je+1
          do i=is,ie+1
             maxu = max(ZERO,ufull(i-1,j,1))
             minu = min(ZERO,ufull(i  ,j,1))
             ! extrapolate both components of velocity to left face
             ulx(i,j,1) = utilde(i-1,j,1) + (HALF - (dt2/hx)*maxu)*slopex(i-1,j,1)
             ulx(i,j,2) = utilde(i-1,j,2) + (HALF - (dt2/hx)*maxu)*slopex(i-1,j,2)
             ! extrapolate both components of velocity to right face
             urx(i,j,1) = utilde(i  ,j,1) - (HALF + (dt2/hx)*minu)*slopex(i  ,j,1)
             urx(i,j,2) = utilde(i  ,j,2) - (HALF + (dt2/hx)*minu)*slopex(i  ,j,2)
          end do
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate both components of velocity to left face
             ulx(i,j,1) = Ipu(i-1,j,1)
             ulx(i,j,2) = Ipv(i-1,j,1)
             ! extrapolate both components of velocity to right face
             urx(i,j,1) = Imu(i,j,1)
             urx(i,j,2) = Imv(i,j,1)
          end do
       end do
    end if

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          ulx(is,js-1:je+1,1:2) = utilde(is-1,js-1:je+1,1:2)
          urx(is,js-1:je+1,1:2) = utilde(is-1,js-1:je+1,1:2)
       case (SlipWall, Symmetry)
          ulx(is,js-1:je+1,1) = ZERO
          urx(is,js-1:je+1,1) = ZERO
          ulx(is,js-1:je+1,2) = urx(is,js-1:je+1,2)
       case (NoSlipWall)
          ulx(is,js-1:je+1,1:2) = ZERO
          urx(is,js-1:je+1,1:2) = ZERO
       case (Outflow)
          urx(is,js-1:je+1,1) = min(urx(is,js-1:je+1,1),ZERO)
          urx(is,js-1:je+1,1:2) = ulx(is,js-1:je+1,1:2)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          ulx(ie+1,js-1:je+1,1:2) = utilde(ie+1,js-1:je+1,1:2)
          urx(ie+1,js-1:je+1,1:2) = utilde(ie+1,js-1:je+1,1:2)
       case (SlipWall, Symmetry)
          ulx(ie+1,js-1:je+1,1) = ZERO
          urx(ie+1,js-1:je+1,1) = ZERO
          urx(ie+1,js-1:je+1,2) = ulx(ie+1,js-1:je+1,2)
       case (NoSlipWall)
          ulx(ie+1,js-1:je+1,1:2) = ZERO
          urx(ie+1,js-1:je+1,1:2) = ZERO
       case (Outflow)
          ulx(ie+1,js-1:je+1,1) = max(ulx(ie+1,js-1:je+1,1),ZERO)
          urx(ie+1,js-1:je+1,1:2) = ulx(ie+1,js-1:je+1,1:2)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    do j=js-1,je+1
       do i=is,ie+1
          ! No need to compute uimhx(:,:,1) since it's equal to utrans-w0
          ! upwind using full velocity to get transverse component of uimhx
          ! Note: utrans already contains w0
          uimhx(i,j,2) = merge(ulx(i,j,2),urx(i,j,2),utrans(i,j).gt.ZERO)
          uavg = HALF*(ulx(i,j,2)+urx(i,j,2))
          uimhx(i,j,2) = merge(uavg,uimhx(i,j,2),abs(utrans(i,j)).lt.rel_eps)
       enddo
    enddo

    if (ppm_type .eq. 0) then
       do j=js,je+1
          do i=is-1,ie+1
             maxu = max(ZERO,ufull(i,j-1,2))
             minu = min(ZERO,ufull(i,j  ,2))
             ! extrapolate both components of velocity to left face
             uly(i,j,1) = utilde(i,j-1,1) + (HALF-(dt2/hy)*maxu)*slopey(i,j-1,1)
             uly(i,j,2) = utilde(i,j-1,2) + (HALF-(dt2/hy)*maxu)*slopey(i,j-1,2)
             ! extrapolate both components of velocity to right face
             ury(i,j,1) = utilde(i,j  ,1) - (HALF+(dt2/hy)*minu)*slopey(i,j  ,1)
             ury(i,j,2) = utilde(i,j  ,2) - (HALF+(dt2/hy)*minu)*slopey(i,j  ,2)
          end do
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate both components of velocity to left face
             uly(i,j,1) = Ipu(i,j-1,2)
             uly(i,j,2) = Ipv(i,j-1,2)
             ! extrapolate both components of velocity to right face
             ury(i,j,1) = Imu(i,j,2)
             ury(i,j,2) = Imv(i,j,2)
          end do
       end do
    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          uly(is-1:ie+1,js,1:2) = utilde(is-1:ie+1,js-1,1:2)
          ury(is-1:ie+1,js,1:2) = utilde(is-1:ie+1,js-1,1:2)
       case (SlipWall, Symmetry)
          uly(is-1:ie+1,js,1) = ury(is-1:ie+1,js,1)
          uly(is-1:ie+1,js,2) = ZERO
          ury(is-1:ie+1,js,2) = ZERO
       case (NoSlipWall)
          uly(is-1:ie+1,js,1:2) = ZERO
          ury(is-1:ie+1,js,1:2) = ZERO
       case (Outflow)
          ury(is-1:ie+1,js,2) = min(ury(is-1:ie+1,js,2),ZERO)
          uly(is-1:ie+1,js,1:2) = ury(is-1:ie+1,js,1:2)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          uly(is-1:ie+1,je+1,1:2) = utilde(is-1:ie+1,je+1,1:2)
          ury(is-1:ie+1,je+1,1:2) = utilde(is-1:ie+1,je+1,1:2)
       case (SlipWall, Symmetry)
          ury(is-1:ie+1,je+1,1) = uly(is-1:ie+1,je+1,1)
          uly(is-1:ie+1,je+1,2) = ZERO
          ury(is-1:ie+1,je+1,2) = ZERO
       case (NoSlipWall)
          uly(is-1:ie+1,je+1,1:2) = ZERO
          ury(is-1:ie+1,je+1,1:2) = ZERO
       case (Outflow)
          uly(is-1:ie+1,je+1,2)   = max(uly(is-1:ie+1,je+1,2),ZERO)
          ury(is-1:ie+1,je+1,1:2) = uly(is-1:ie+1,je+1,1:2)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    do j=js,je+1
       do i=is-1,ie+1
          ! No need to compute uimhy(:,:,2) since it's equal to vtrans-w0
          ! upwind using full velocity to get transverse component of uimhy
          ! Note: vtrans already contains w0
          uimhy(i,j,1) = merge(uly(i,j,1),ury(i,j,1),vtrans(i,j).gt.ZERO)
          uavg = HALF*(uly(i,j,1)+ury(i,j,1))
          uimhy(i,j,1) = merge(uavg,uimhy(i,j,1),abs(vtrans(i,j)).lt.rel_eps)
       enddo
    enddo

    !******************************************************************
    ! Create umac and vmac
    !******************************************************************

    do j=js,je
       do i=is,ie+1
          ! use the traced force if ppm_trace_forces = 1
          fl = merge(force(i-1,j,1), Ipfx(i-1,j,1), ppm_trace_forces == 0)
          fr = merge(force(i,j  ,1), Imfx(i,  j,1), ppm_trace_forces == 0)

          ! extrapolate to edges
          umacl(i,j) = ulx(i,j,1) &
               - (dt4/hy)*(vtrans(i-1,j+1)+vtrans(i-1,j)) &
               * (uimhy(i-1,j+1,1)-uimhy(i-1,j,1)) + dt2*fl
          umacr(i,j) = urx(i,j,1) &
               - (dt4/hy)*(vtrans(i  ,j+1)+vtrans(i  ,j)) &
               * (uimhy(i  ,j+1,1)-uimhy(i  ,j,1)) + dt2*fr

          ! solve Riemann problem using full velocity
          uavg = HALF*(umacl(i,j)+umacr(i,j))
          test = ((umacl(i,j) .le. ZERO .and. umacr(i,j) .ge. ZERO) .or. &
               (abs(umacl(i,j)+umacr(i,j)) .lt. rel_eps))
          umac(i,j) = merge(umacl(i,j),umacr(i,j),uavg .gt. ZERO)
          umac(i,j) = merge(ZERO,umac(i,j),test)
       enddo
    enddo

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          umac(is,js:je) = utilde(is-1,js:je,1)
       case (SlipWall, NoSlipWall, Symmetry)
          umac(is,js:je) = ZERO
       case (Outflow)
          umac(is,js:je) = min(umacr(is,js:je),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          umac(ie+1,js:je) = utilde(ie+1,js:je,1)
       case (SlipWall, NoSlipWall, Symmetry)
          umac(ie+1,js:je) = ZERO
       case (Outflow)
          umac(ie+1,js:je) = max(umacl(ie+1,js:je),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(1,2)")
       end select
    end if


    do j=js,je+1
       do i=is,ie
          ! use the traced force if ppm_trace_forces = 1
          fl = merge(force(i,j-1,2), Ipfy(i,j-1,2), ppm_trace_forces == 0)
          fr = merge(force(i,j  ,2), Imfy(i,j  ,2), ppm_trace_forces == 0)

          ! extrapolate to edges
          vmacl(i,j) = uly(i,j,2) &
               - (dt4/hx)*(utrans(i+1,j-1)+utrans(i,j-1)) &
               * (uimhx(i+1,j-1,2)-uimhx(i,j-1,2)) + dt2*fl
          vmacr(i,j) = ury(i,j,2) &
               - (dt4/hx)*(utrans(i+1,j  )+utrans(i,j  )) &
               * (uimhx(i+1,j  ,2)-uimhx(i,j  ,2)) + dt2*fr

          ! solve Riemann problem using full velocity
          uavg = HALF*(vmacl(i,j)+vmacr(i,j))
          test = ((vmacl(i,j)+w0(lev,j) .le. ZERO .and. vmacr(i,j)+w0(lev,j) .ge. ZERO) .or. &
               (abs(vmacl(i,j)+vmacr(i,j)+TWO*w0(lev,j)) .lt. rel_eps))
          vmac(i,j) = merge(vmacl(i,j),vmacr(i,j),uavg+w0(lev,j) .gt. ZERO)
          vmac(i,j) = merge(ZERO,vmac(i,j),test)

       enddo
    enddo


    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          vmac(is:ie,js) = utilde(is:ie,js-1,2)
       case (SlipWall, NoSlipWall, Symmetry)
          vmac(is:ie,js) = ZERO
       case (Outflow)
          vmac(is:ie,js) = min(vmacr(is:ie,js),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          vmac(is:ie,je+1) = utilde(is:ie,je+1,2)
       case (SlipWall, NoSlipWall, Symmetry)
          vmac(is:ie,je+1) = ZERO
       case (Outflow)
          vmac(is:ie,je+1) = max(vmacl(is:ie,je+1),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_2d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    call bl_deallocate(slopex)
    call bl_deallocate(slopey)

    call bl_deallocate(Ipu)
    call bl_deallocate(Imu)
    call bl_deallocate(Ipfx)
    call bl_deallocate(Imfx)
    call bl_deallocate(Ipv)
    call bl_deallocate(Imv)
    call bl_deallocate(Ipfy)
    call bl_deallocate(Imfy)

    call bl_deallocate(ulx)
    call bl_deallocate(urx)
    call bl_deallocate(uimhx)

    call bl_deallocate(uly)
    call bl_deallocate(ury)
    call bl_deallocate(uimhy)

    call bl_deallocate(umacl)
    call bl_deallocate(umacr)
    call bl_deallocate(vmacl)
    call bl_deallocate(vmacr)

  end subroutine velpred_2d
#endif

#if (AMREX_SPACEDIM == 3)
  subroutine velpred_3d(lev, domlo, domhi, lo, hi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       wtrans, uw_lo, uw_hi, &
       umac,   mu_lo, mu_hi, &
       vmac,   mv_lo, mv_hi, &
       wmac,   mw_lo, mw_hi, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       force,   f_lo,  f_hi, nc_f, ng_f, &
       w0,dx,dt,adv_bc,phys_bc) bind(C,name="velpred_3d")

    integer         , intent(in   ) :: lev, domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: ut_lo(3), ut_hi(3), nc_ut
    integer, value,   intent(in   ) :: ng_ut
    integer         , intent(in   ) :: uf_lo(3), uf_hi(3), nc_uf
    integer, value,   intent(in   ) :: ng_uf
    integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
    integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
    integer         , intent(in   ) :: uw_lo(3), uw_hi(3)
    integer         , intent(in   ) :: mu_lo(3), mu_hi(3)
    integer         , intent(in   ) :: mv_lo(3), mv_hi(3)
    integer         , intent(in   ) :: mw_lo(3), mw_hi(3)
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    integer         , intent(in   ) ::  f_lo(3),  f_hi(3), nc_f
    integer, value,   intent(in   ) :: ng_f
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
    double precision, intent(inout) :: wtrans(uw_lo(1):uw_hi(1),uw_lo(2):uw_hi(2),uw_lo(3):uw_hi(3))
    double precision, intent(inout) :: umac  (mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))
    double precision, intent(inout) :: vmac  (mv_lo(1):mv_hi(1),mv_lo(2):mv_hi(2),mv_lo(3):mv_hi(3))
    double precision, intent(inout) :: wmac  (mw_lo(1):mw_hi(1),mw_lo(2):mw_hi(2),mw_lo(3):mw_hi(3))
    double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    double precision, intent(in   ) :: force ( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3),nc_f)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: adv_bc(3,2,3), phys_bc(3,2) ! dim, lohi, (comp)

    ! local variables
    double precision, pointer :: slopex(:,:,:,:)
    double precision, pointer :: slopey(:,:,:,:)
    double precision, pointer :: slopez(:,:,:,:)

    double precision, pointer :: Ipu(:,:,:,:), Ipfx(:,:,:,:)
    double precision, pointer :: Imu(:,:,:,:), Imfx(:,:,:,:)
    double precision, pointer :: Ipv(:,:,:,:), Ipfy(:,:,:,:)
    double precision, pointer :: Imv(:,:,:,:), Imfy(:,:,:,:)
    double precision, pointer :: Ipw(:,:,:,:), Ipfz(:,:,:,:)
    double precision, pointer :: Imw(:,:,:,:), Imfz(:,:,:,:)

    ! these correspond to u_L^x, etc.
    double precision, pointer:: ulx(:,:,:,:),urx(:,:,:,:),uimhx(:,:,:,:)
    double precision, pointer:: uly(:,:,:,:),ury(:,:,:,:),uimhy(:,:,:,:)
    double precision, pointer:: ulz(:,:,:,:),urz(:,:,:,:),uimhz(:,:,:,:)

    ! these correspond to u_L^{y|z}, etc.
    double precision, pointer:: ulyz(:,:,:)
    double precision, pointer:: uryz(:,:,:)
    double precision, pointer:: uimhyz(:,:,:)

    double precision, pointer:: ulzy(:,:,:)
    double precision, pointer:: urzy(:,:,:)
    double precision, pointer:: uimhzy(:,:,:)

    double precision, pointer:: vlxz(:,:,:)
    double precision, pointer:: vrxz(:,:,:)
    double precision, pointer:: vimhxz(:,:,:)

    double precision, pointer:: vlzx(:,:,:)
    double precision, pointer:: vrzx(:,:,:)
    double precision, pointer:: vimhzx(:,:,:)

    double precision, pointer:: wlxy(:,:,:)
    double precision, pointer:: wrxy(:,:,:)
    double precision, pointer:: wimhxy(:,:,:)

    double precision, pointer:: wlyx(:,:,:)
    double precision, pointer:: wryx(:,:,:)
    double precision, pointer:: wimhyx(:,:,:)

    ! these correspond to umac_L, etc.
    double precision, pointer:: umacl(:,:,:),umacr(:,:,:)
    double precision, pointer:: vmacl(:,:,:),vmacr(:,:,:)
    double precision, pointer:: wmacl(:,:,:),wmacr(:,:,:)

    double precision :: hx, hy, hz, dt2, dt4, dt6, uavg, maxu, minu
    double precision :: fl, fr

    integer :: i,j,k,is,js,ie,je,ks,ke

    logical :: test

    call bl_allocate(slopex,lo-1,hi+1,3)
    call bl_allocate(slopey,lo-1,hi+1,3)
    call bl_allocate(slopez,lo-1,hi+1,3)

    call bl_allocate(Ipu,lo-1,hi+1,3)
    call bl_allocate(Imu,lo-1,hi+1,3)
    call bl_allocate(Ipv,lo-1,hi+1,3)
    call bl_allocate(Imv,lo-1,hi+1,3)
    call bl_allocate(Ipw,lo-1,hi+1,3)
    call bl_allocate(Imw,lo-1,hi+1,3)

    call bl_allocate(Ipfx,lo-1,hi+1,3)
    call bl_allocate(Imfx,lo-1,hi+1,3)
    call bl_allocate(Ipfy,lo-1,hi+1,3)
    call bl_allocate(Imfy,lo-1,hi+1,3)
    call bl_allocate(Ipfz,lo-1,hi+1,3)
    call bl_allocate(Imfz,lo-1,hi+1,3)

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)

    dt2 = HALF*dt
    dt4 = dt/4.0d0
    dt6 = dt/6.0d0

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(k)
       do k = lo(3)-1,hi(3)+1
          call slopex_2d(utilde(:,:,k,:),slopex(:,:,k,:),domlo,domhi,lo,hi,ng_ut,3,adv_bc)
          call slopey_2d(utilde(:,:,k,:),slopey(:,:,k,:),domlo,domhi,lo,hi,ng_ut,3,adv_bc)
       end do
       !$OMP END PARALLEL DO
       call slopez_3d(utilde,slopez,domlo,domhi,lo,hi,ng_ut,3,adv_bc)
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(utilde(:,:,:,1),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ipu,Imu,domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
       call ppm_3d(utilde(:,:,:,2),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ipv,Imv,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
       call ppm_3d(utilde(:,:,:,3),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ipw,Imw,domlo,domhi,lo,hi,adv_bc(:,:,3),dx,dt,.false.)

       ! trace forces, if necessary.  Note by default the ppm routines
       ! will trace each component to each interface in all coordinate
       ! directions, but we really only need the force traced along
       ! its respective dimension.  This should be simplified later.
       if (ppm_trace_forces .eq. 1) then
          call ppm_3d(force(:,:,:,1),ng_f, &
               ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
               Ipfx,Imfx,domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
          call ppm_3d(force(:,:,:,2),ng_f, &
               ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
               Ipfy,Imfy,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
          call ppm_3d(force(:,:,:,3),ng_f, &
               ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
               Ipfz,Imfz,domlo,domhi,lo,hi,adv_bc(:,:,3),dx,dt,.false.)
       endif
    end if

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    call bl_allocate(ulx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1,1,3)
    call bl_allocate(urx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1,1,3)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,maxu,minu)
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                maxu = (HALF - dt2*max(ZERO,ufull(i-1,j,k,1))/hx)
                minu = (HALF + dt2*min(ZERO,ufull(i  ,j,k,1))/hx)

                ! extrapolate all components of velocity to left face
                ulx(i,j,k,1) = utilde(i-1,j,k,1) + maxu * slopex(i-1,j,k,1)
                ulx(i,j,k,2) = utilde(i-1,j,k,2) + maxu * slopex(i-1,j,k,2)
                ulx(i,j,k,3) = utilde(i-1,j,k,3) + maxu * slopex(i-1,j,k,3)

                ! extrapolate all components of velocity to right face
                urx(i,j,k,1) = utilde(i,j,k,1) - minu * slopex(i,j,k,1)
                urx(i,j,k,2) = utilde(i,j,k,2) - minu * slopex(i,j,k,2)
                urx(i,j,k,3) = utilde(i,j,k,3) - minu * slopex(i,j,k,3)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! extrapolate all components of velocity to left face
                ulx(i,j,k,1) = Ipu(i-1,j,k,1)
                ulx(i,j,k,2) = Ipv(i-1,j,k,1)
                ulx(i,j,k,3) = Ipw(i-1,j,k,1)

                ! extrapolate all components of velocity to right face
                urx(i,j,k,1) = Imu(i,j,k,1)
                urx(i,j,k,2) = Imv(i,j,k,1)
                urx(i,j,k,3) = Imw(i,j,k,1)
             end do
          end do
       end do
    end if

    call bl_deallocate(slopex)

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          ulx(is,js-1:je+1,ks-1:ke+1,1:3) = utilde(is-1,js-1:je+1,ks-1:ke+1,1:3)
          urx(is,js-1:je+1,ks-1:ke+1,1:3) = utilde(is-1,js-1:je+1,ks-1:ke+1,1:3)
       case (SlipWall, Symmetry)
          ulx(is,js-1:je+1,ks-1:ke+1,1) = ZERO
          urx(is,js-1:je+1,ks-1:ke+1,1) = ZERO
          ulx(is,js-1:je+1,ks-1:ke+1,2) = urx(is,js-1:je+1,ks-1:ke+1,2)
          ulx(is,js-1:je+1,ks-1:ke+1,3) = urx(is,js-1:je+1,ks-1:ke+1,3)
       case (NoSlipWall)
          ulx(is,js-1:je+1,ks-1:ke+1,1:3) = ZERO
          urx(is,js-1:je+1,ks-1:ke+1,1:3) = ZERO
       case (Outflow)
          urx(is,js-1:je+1,ks-1:ke+1,1) = min(urx(is,js-1:je+1,ks-1:ke+1,1),ZERO)
          ulx(is,js-1:je+1,ks-1:ke+1,1:3) = ulx(is,js-1:je+1,ks-1:ke+1,1:3)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          ulx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = utilde(ie+1,js-1:je+1,ks-1:ke+1,1:)
          urx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = utilde(ie+1,js-1:je+1,ks-1:ke+1,1:3)
       case (SlipWall, Symmetry)
          ulx(ie+1,js-1:je+1,ks-1:ke+1,1) = ZERO
          urx(ie+1,js-1:je+1,ks-1:ke+1,1) = ZERO
          urx(ie+1,js-1:je+1,ks-1:ke+1,2) = ulx(ie+1,js-1:je+1,ks-1:ke+1,2)
          urx(ie+1,js-1:je+1,ks-1:ke+1,3) = ulx(ie+1,js-1:je+1,ks-1:ke+1,3)
       case (NoSlipWall)
          ulx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = ZERO
          urx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = ZERO
       case (Outflow)
          ulx(ie+1,js-1:je+1,ks-1:ke+1,1) = max(ulx(ie+1,js-1:je+1,ks-1:ke+1,1),ZERO)
          urx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = ulx(ie+1,js-1:je+1,ks-1:ke+1,1:3)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    call bl_allocate(uimhx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1,1,3)

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             ! No need to compute uimhx(:,:,:,1) since it's equal to utrans-w0
             ! upwind using full velocity to get transverse components of uimhx
             ! Note: utrans already contains w0
             uimhx(i,j,k,2) = merge(ulx(i,j,k,2),urx(i,j,k,2),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,2)+urx(i,j,k,2))
             uimhx(i,j,k,2) = merge(uavg,uimhx(i,j,k,2),abs(utrans(i,j,k)).lt.rel_eps)

             uimhx(i,j,k,3) = merge(ulx(i,j,k,3),urx(i,j,k,3),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,3)+urx(i,j,k,3))
             uimhx(i,j,k,3) = merge(uavg,uimhx(i,j,k,3),abs(utrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! normal predictor states

    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    call bl_allocate(uly,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1,1,3)
    call bl_allocate(ury,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1,1,3)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,minu,maxu)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                maxu = (HALF - dt2*max(ZERO,ufull(i,j-1,k,2))/hy)
                minu = (HALF + dt2*min(ZERO,ufull(i,j  ,k,2))/hy)

                ! extrapolate all components of velocity to left face
                uly(i,j,k,1) = utilde(i,j-1,k,1) + maxu * slopey(i,j-1,k,1)
                uly(i,j,k,2) = utilde(i,j-1,k,2) + maxu * slopey(i,j-1,k,2)
                uly(i,j,k,3) = utilde(i,j-1,k,3) + maxu * slopey(i,j-1,k,3)

                ! extrapolate all components of velocity to right face
                ury(i,j,k,1) = utilde(i,j,k,1) - minu * slopey(i,j,k,1)
                ury(i,j,k,2) = utilde(i,j,k,2) - minu * slopey(i,j,k,2)
                ury(i,j,k,3) = utilde(i,j,k,3) - minu * slopey(i,j,k,3)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                ! extrapolate all components of velocity to left face
                uly(i,j,k,1) = Ipu(i,j-1,k,2)
                uly(i,j,k,2) = Ipv(i,j-1,k,2)
                uly(i,j,k,3) = Ipw(i,j-1,k,2)

                ! extrapolate all components of velocity to right face
                ury(i,j,k,1) = Imu(i,j,k,2)
                ury(i,j,k,2) = Imv(i,j,k,2)
                ury(i,j,k,3) = Imw(i,j,k,2)
             enddo
          enddo
       enddo
    end if

    call bl_deallocate(slopey)

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          uly(is-1:ie+1,js,ks-1:ke+1,1:3) = utilde(is-1:ie+1,js-1,ks-1:ke+1,1:3)
          ury(is-1:ie+1,js,ks-1:ke+1,1:3) = utilde(is-1:ie+1,js-1,ks-1:ke+1,1:3)
       case (SlipWall, Symmetry)
          uly(is-1:ie+1,js,ks-1:ke+1,1) = ury(is-1:ie+1,js,ks-1:ke+1,1)
          uly(is-1:ie+1,js,ks-1:ke+1,2) = ZERO
          ury(is-1:ie+1,js,ks-1:ke+1,2) = ZERO
          uly(is-1:ie+1,js,ks-1:ke+1,3) = ury(is-1:ie+1,js,ks-1:ke+1,3)
       case (NoSlipWall)
          uly(is-1:ie+1,js,ks-1:ke+1,1:3) = ZERO
          ury(is-1:ie+1,js,ks-1:ke+1,1:3) = ZERO
       case (Outflow)
          ury(is-1:ie+1,js,ks-1:ke+1,2) = min(ury(is-1:ie+1,js,ks-1:ke+1,2),ZERO)
          uly(is-1:ie+1,js,ks-1:ke+1,1:3) = ury(is-1:ie+1,js,ks-1:ke+1,1:3)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          uly(is-1:ie+1,je+1,ks-1:ke+1,1:3) = utilde(is-1:ie+1,je+1,ks-1:ke+1,1:3)
          ury(is-1:ie+1,je+1,ks-1:ke+1,1:3) = utilde(is-1:ie+1,je+1,ks-1:ke+1,1:3)
       case (SlipWall, Symmetry)
          ury(is-1:ie+1,je+1,ks-1:ke+1,1) = uly(is-1:ie+1,je+1,ks-1:ke+1,1)
          uly(is-1:ie+1,je+1,ks-1:ke+1,2) = ZERO
          ury(is-1:ie+1,je+1,ks-1:ke+1,2) = ZERO
          ury(is-1:ie+1,je+1,ks-1:ke+1,3) = uly(is-1:ie+1,je+1,ks-1:ke+1,3)
       case (NoSlipWall)
          uly(is-1:ie+1,je+1,ks-1:ke+1,1:3) = ZERO
          ury(is-1:ie+1,je+1,ks-1:ke+1,1:3) = ZERO
       case (Outflow)
          uly(is-1:ie+1,je+1,ks-1:ke+1,2) = max(uly(is-1:ie+1,je+1,ks-1:ke+1,2),ZERO)
          ury(is-1:ie+1,je+1,ks-1:ke+1,1:3) = uly(is-1:ie+1,je+1,ks-1:ke+1,1:3)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    call bl_allocate(uimhy,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1,1,3)

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             ! No need to compute uimhy(:,:,:,2) since it's equal to vtrans-w0
             ! upwind using full velocity to get transverse components of uimhy
             ! Note: vtrans already contains w0
             uimhy(i,j,k,1) = merge(uly(i,j,k,1),ury(i,j,k,1),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(uly(i,j,k,1)+ury(i,j,k,1))
             uimhy(i,j,k,1) = merge(uavg,uimhy(i,j,k,1),abs(vtrans(i,j,k)).lt.rel_eps)

             uimhy(i,j,k,3) = merge(uly(i,j,k,3),ury(i,j,k,3),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(uly(i,j,k,3)+ury(i,j,k,3))
             uimhy(i,j,k,3) = merge(uavg,uimhy(i,j,k,3),abs(vtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    call bl_allocate(ulz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1,1,3)
    call bl_allocate(urz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1,1,3)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,minu,maxu)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                maxu = (HALF - dt2*max(ZERO,ufull(i,j,k-1,3))/hz)
                minu = (HALF + dt2*min(ZERO,ufull(i,j,k  ,3))/hz)

                ! extrapolate all components of velocity to left face
                ulz(i,j,k,1) = utilde(i,j,k-1,1) + maxu * slopez(i,j,k-1,1)
                ulz(i,j,k,2) = utilde(i,j,k-1,2) + maxu * slopez(i,j,k-1,2)
                ulz(i,j,k,3) = utilde(i,j,k-1,3) + maxu * slopez(i,j,k-1,3)

                ! extrapolate all components of velocity to right face
                urz(i,j,k,1) = utilde(i,j,k,1) - minu * slopez(i,j,k,1)
                urz(i,j,k,2) = utilde(i,j,k,2) - minu * slopez(i,j,k,2)
                urz(i,j,k,3) = utilde(i,j,k,3) - minu * slopez(i,j,k,3)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                ! extrapolate all components of velocity to left face
                ulz(i,j,k,1) = Ipu(i,j,k-1,3)
                ulz(i,j,k,2) = Ipv(i,j,k-1,3)
                ulz(i,j,k,3) = Ipw(i,j,k-1,3)

                ! extrapolate all components of velocity to right face
                urz(i,j,k,1) = Imu(i,j,k,3)
                urz(i,j,k,2) = Imv(i,j,k,3)
                urz(i,j,k,3) = Imw(i,j,k,3)
             end do
          end do
       end do
    end if

    call bl_deallocate(slopez)
    call bl_deallocate(Ipu)
    call bl_deallocate(Imu)
    call bl_deallocate(Ipv)
    call bl_deallocate(Imv)
    call bl_deallocate(Ipw)
    call bl_deallocate(Imw)

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       select case(phys_bc(3,1))
       case (Inflow)
          ulz(is-1:ie+1,js-1:je+1,ks,1:3) = utilde(is-1:ie+1,js-1:je+1,ks-1,1:3)
          urz(is-1:ie+1,js-1:je+1,ks,1:3) = utilde(is-1:ie+1,js-1:je+1,ks-1,1:3)
       case (SlipWall, Symmetry)
          ulz(is-1:ie+1,js-1:je+1,ks,1) = urz(is-1:ie+1,js-1:je+1,ks,1)
          ulz(is-1:ie+1,js-1:je+1,ks,2) = urz(is-1:ie+1,js-1:je+1,ks,2)
          ulz(is-1:ie+1,js-1:je+1,ks,3) = ZERO
          urz(is-1:ie+1,js-1:je+1,ks,3) = ZERO
       case (NoSlipWall)
          ulz(is-1:ie+1,js-1:je+1,ks,1:3) = ZERO
          urz(is-1:ie+1,js-1:je+1,ks,1:3) = ZERO
       case (Outflow)
          urz(is-1:ie+1,js-1:je+1,ks,3) = min(urz(is-1:ie+1,js-1:je+1,ks,3),ZERO)
          ulz(is-1:ie+1,js-1:je+1,ks,1:3) = urz(is-1:ie+1,js-1:je+1,ks,1:3)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       select case(phys_bc(3,2))
       case (Inflow)
          ulz(is-1:ie+1,js-1:je+1,ke+1,1:3) = utilde(is-1:ie+1,js-1:je+1,ke+1,1:3)
          urz(is-1:ie+1,js-1:je+1,ke+1,1:3) = utilde(is-1:ie+1,js-1:je+1,ke+1,1:3)
       case (SlipWall, Symmetry)
          urz(is-1:ie+1,js-1:je+1,ke+1,1) = ulz(is-1:ie+1,js-1:je+1,ke+1,1)
          urz(is-1:ie+1,js-1:je+1,ke+1,2) = ulz(is-1:ie+1,js-1:je+1,ke+1,2)
          ulz(is-1:ie+1,js-1:je+1,ke+1,3) = ZERO
          urz(is-1:ie+1,js-1:je+1,ke+1,3) = ZERO
       case (NoSlipWall)
          ulz(is-1:ie+1,js-1:je+1,ke+1,1:3) = ZERO
          urz(is-1:ie+1,js-1:je+1,ke+1,1:3) = ZERO
       case (Outflow)
          ulz(is-1:ie+1,js-1:je+1,ke+1,3) = max(ulz(is-1:ie+1,js-1:je+1,ke+1,3),ZERO)
          urz(is-1:ie+1,js-1:je+1,ke+1,1:3) = ulz(is-1:ie+1,js-1:je+1,ke+1,1:3)
       case (Interior)
       case default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
       end select
    end if

    call bl_allocate(uimhz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1,1,3)

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             ! No need to compute uimhz(:,:,:,3) since it's equal to wtrans-w0
             ! upwind using full velocity to get transverse components of uimhz
             ! Note: wtrans already contains w0
             uimhz(i,j,k,1) = merge(ulz(i,j,k,1),urz(i,j,k,1),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,1)+urz(i,j,k,1))
             uimhz(i,j,k,1) = merge(uavg,uimhz(i,j,k,1),abs(wtrans(i,j,k)).lt.rel_eps)

             uimhz(i,j,k,2) = merge(ulz(i,j,k,2),urz(i,j,k,2),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,2)+urz(i,j,k,2))
             uimhz(i,j,k,2) = merge(uavg,uimhz(i,j,k,2),abs(wtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !******************************************************************
    ! Create u_{\i-\half\e_y}^{y|z}, etc.
    !******************************************************************

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    call bl_allocate(ulyz,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3))
    call bl_allocate(uryz,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3))
    call bl_allocate(uimhyz,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3))

    ! uimhyz loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate to faces
             ulyz(i,j,k) = uly(i,j,k,1) - (dt6/hz)*(wtrans(i,j-1,k+1)+wtrans(i,j-1,k)) &
                  * (uimhz(i,j-1,k+1,1)-uimhz(i,j-1,k,1))
             uryz(i,j,k) = ury(i,j,k,1) - (dt6/hz)*(wtrans(i,j  ,k+1)+wtrans(i,j  ,k)) &
                  * (uimhz(i,j  ,k+1,1)-uimhz(i,j  ,k,1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          ulyz(is-1:ie+1,js,ks:ke) = utilde(is-1:ie+1,js-1,ks:ke,1)
          uryz(is-1:ie+1,js,ks:ke) = utilde(is-1:ie+1,js-1,ks:ke,1)
       case (SlipWall, Symmetry, Outflow)
          ulyz(is-1:ie+1,js,ks:ke) = uryz(is-1:ie+1,js,ks:ke)
       case (NoSlipWall)
          ulyz(is-1:ie+1,js,ks:ke) = ZERO
          uryz(is-1:ie+1,js,ks:ke) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          ulyz(is-1:ie+1,je+1,ks:ke) = utilde(is-1:ie+1,je+1,ks:ke,1)
          uryz(is-1:ie+1,je+1,ks:ke) = utilde(is-1:ie+1,je+1,ks:ke,1)
       case (SlipWall, Symmetry, Outflow)
          uryz(is-1:ie+1,je+1,ks:ke) = ulyz(is-1:ie+1,je+1,ks:ke)
       case (NoSlipWall)
          ulyz(is-1:ie+1,je+1,ks:ke) = ZERO
          uryz(is-1:ie+1,je+1,ks:ke) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! upwind using full velocity
             uimhyz(i,j,k) = merge(ulyz(i,j,k),uryz(i,j,k),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulyz(i,j,k)+uryz(i,j,k))
             uimhyz(i,j,k) = merge(uavg,uimhyz(i,j,k),abs(vtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(ulyz)
    call bl_deallocate(uryz)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    call bl_allocate(ulzy,lo(1)-1,hi(1)+1,lo(2),hi(2),lo(3),hi(3)+1)
    call bl_allocate(urzy,lo(1)-1,hi(1)+1,lo(2),hi(2),lo(3),hi(3)+1)
    call bl_allocate(uimhzy,lo(1)-1,hi(1)+1,lo(2),hi(2),lo(3),hi(3)+1)

    ! uimhzy loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! extrapolate to faces
             ulzy(i,j,k) = ulz(i,j,k,1) - (dt6/hy)*(vtrans(i,j+1,k-1)+vtrans(i,j,k-1)) &
                  * (uimhy(i,j+1,k-1,1)-uimhy(i,j,k-1,1))
             urzy(i,j,k) = urz(i,j,k,1) - (dt6/hy)*(vtrans(i,j+1,k  )+vtrans(i,j,k  )) &
                  * (uimhy(i,j+1,k  ,1)-uimhy(i,j,k  ,1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       select case(phys_bc(3,1))
       case (Inflow)
          ulzy(is-1:ie+1,js:je,ks) = utilde(is-1:ie+1,js:je,ks-1,1)
          urzy(is-1:ie+1,js:je,ks) = utilde(is-1:ie+1,js:je,ks-1,1)
       case (SlipWall, Symmetry, Outflow)
          ulzy(is-1:ie+1,js:je,ks) = urzy(is-1:ie+1,js:je,ks)
       case (NoSlipWall)
          ulzy(is-1:ie+1,js:je,ks) = ZERO
          urzy(is-1:ie+1,js:je,ks) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       select case(phys_bc(3,2))
       case (Inflow)
          ulzy(is-1:ie+1,js:je,ke+1) = utilde(is-1:ie+1,js:je,ke+1,1)
          urzy(is-1:ie+1,js:je,ke+1) = utilde(is-1:ie+1,js:je,ke+1,1)
       case (SlipWall, Symmetry, Outflow)
          urzy(is-1:ie+1,js:je,ke+1) = ulzy(is-1:ie+1,js:je,ke+1)
       case (NoSlipWall)
          ulzy(is-1:ie+1,js:je,ke+1) = ZERO
          urzy(is-1:ie+1,js:je,ke+1) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
       end select
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! upwind using full velocity
             uimhzy(i,j,k) = merge(ulzy(i,j,k),urzy(i,j,k),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulzy(i,j,k)+urzy(i,j,k))
             uimhzy(i,j,k) = merge(uavg,uimhzy(i,j,k),abs(wtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(ulzy)
    call bl_deallocate(urzy)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    call bl_allocate(vlxz,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3))
    call bl_allocate(vrxz,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3))

    ! vimhxz loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate to faces
             vlxz(i,j,k) = ulx(i,j,k,2) - (dt6/hz)*(wtrans(i-1,j,k+1)+wtrans(i-1,j,k)) &
                  * (uimhz(i-1,j,k+1,2)-uimhz(i-1,j,k,2))
             vrxz(i,j,k) = urx(i,j,k,2) - (dt6/hz)*(wtrans(i  ,j,k+1)+wtrans(i  ,j,k)) &
                  * (uimhz(i  ,j,k+1,2)-uimhz(i  ,j,k,2))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(uimhz)

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          vlxz(is,js-1:je+1,ks:ke) = utilde(is-1,js-1:je+1,ks:ke,2)
          vrxz(is,js-1:je+1,ks:ke) = utilde(is-1,js-1:je+1,ks:ke,2)
       case (SlipWall, Symmetry, Outflow)
          vlxz(is,js-1:je+1,ks:ke) = vrxz(is,js-1:je+1,ks:ke)
       case (NoSlipWall)
          vlxz(is,js-1:je+1,ks:ke) = ZERO
          vrxz(is,js-1:je+1,ks:ke) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          vlxz(ie+1,js-1:je+1,ks:ke) = utilde(ie+1,js-1:je+1,ks:ke,2)
          vrxz(ie+1,js-1:je+1,ks:ke) = utilde(ie+1,js-1:je+1,ks:ke,2)
       case (SlipWall, Symmetry, Outflow)
          vrxz(ie+1,js-1:je+1,ks:ke) = vlxz(ie+1,js-1:je+1,ks:ke)
       case (NoSlipWall)
          vlxz(ie+1,js-1:je+1,ks:ke) = ZERO
          vrxz(ie+1,js-1:je+1,ks:ke) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    call bl_allocate(vimhxz,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! upwind using full velocity
             vimhxz(i,j,k) = merge(vlxz(i,j,k),vrxz(i,j,k),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(vlxz(i,j,k)+vrxz(i,j,k))
             vimhxz(i,j,k) = merge(uavg,vimhxz(i,j,k),abs(utrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(vlxz)
    call bl_deallocate(vrxz)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    call bl_allocate(vlzx,lo(1),hi(1),lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(vrzx,lo(1),hi(1),lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(vimhzx,lo(1),hi(1),lo(2)-1,hi(2)+1,lo(3),hi(3)+1)

    ! vimhzx loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! extrapolate to faces
             vlzx(i,j,k) = ulz(i,j,k,2) - (dt6/hx)*(utrans(i+1,j,k-1)+utrans(i,j,k-1)) &
                  * (uimhx(i+1,j,k-1,2)-uimhx(i,j,k-1,2))
             vrzx(i,j,k) = urz(i,j,k,2) - (dt6/hx)*(utrans(i+1,j,k  )+utrans(i,j,k  )) &
                  * (uimhx(i+1,j,k  ,2)-uimhx(i,j,k  ,2))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       select case(phys_bc(3,1))
       case (Inflow)
          vlzx(is:ie,js-1:je+1,ks) = utilde(is:ie,js-1:je+1,ks-1,2)
          vrzx(is:ie,js-1:je+1,ks) = utilde(is:ie,js-1:je+1,ks-1,2)
       case (SlipWall, Symmetry, Outflow)
          vlzx(is:ie,js-1:je+1,ks) = vrzx(is:ie,js-1:je+1,ks)
       case (NoSlipWall)
          vlzx(is:ie,js-1:je+1,ks) = ZERO
          vrzx(is:ie,js-1:je+1,ks) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       select case(phys_bc(3,2))
       case (Inflow)
          vlzx(is:ie,js-1:je+1,ke+1) = utilde(is:ie,js-1:je+1,ke+1,2)
          vrzx(is:ie,js-1:je+1,ke+1) = utilde(is:ie,js-1:je+1,ke+1,2)
       case (SlipWall, Symmetry, Outflow)
          vrzx(is:ie,js-1:je+1,ke+1) = vlzx(is:ie,js-1:je+1,ke+1)
       case (NoSlipWall)
          vlzx(is:ie,js-1:je+1,ke+1) = ZERO
          vrzx(is:ie,js-1:je+1,ke+1) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
       end select
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! upwind using full velocity
             vimhzx(i,j,k) = merge(vlzx(i,j,k),vrzx(i,j,k),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(vlzx(i,j,k)+vrzx(i,j,k))
             vimhzx(i,j,k) = merge(uavg,vimhzx(i,j,k),abs(wtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(vlzx)
    call bl_deallocate(vrzx)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    call bl_allocate(wlxy,lo(1),hi(1)+1,lo(2),hi(2),lo(3)-1,hi(3)+1)
    call bl_allocate(wrxy,lo(1),hi(1)+1,lo(2),hi(2),lo(3)-1,hi(3)+1)

    ! wimhxy loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! extrapolate to faces
             wlxy(i,j,k) = ulx(i,j,k,3) - (dt6/hy)*(vtrans(i-1,j+1,k)+vtrans(i-1,j,k)) &
                  * (uimhy(i-1,j+1,k,3)-uimhy(i-1,j,k,3))
             wrxy(i,j,k) = urx(i,j,k,3) - (dt6/hy)*(vtrans(i  ,j+1,k)+vtrans(i  ,j,k)) &
                  * (uimhy(i  ,j+1,k,3)-uimhy(i  ,j,k,3))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(uimhy)

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          wlxy(is,js:je,ks-1:ke+1) = utilde(is-1,js:je,ks-1:ke+1,3)
          wrxy(is,js:je,ks-1:ke+1) = utilde(is-1,js:je,ks-1:ke+1,3)
       case (SlipWall, Symmetry, Outflow)
          wlxy(is,js:je,ks-1:ke+1) = wrxy(is,js:je,ks-1:ke+1)
       case (NoSlipWall)
          wlxy(is,js:je,ks-1:ke+1) = ZERO
          wrxy(is,js:je,ks-1:ke+1) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          wlxy(ie+1,js:je,ks-1:ke+1) = utilde(ie+1,js:je,ks-1:ke+1,3)
          wrxy(ie+1,js:je,ks-1:ke+1) = utilde(ie+1,js:je,ks-1:ke+1,3)
       case (SlipWall, Symmetry, Outflow)
          wrxy(ie+1,js:je,ks-1:ke+1) = wlxy(ie+1,js:je,ks-1:ke+1)
       case (NoSlipWall)
          wlxy(ie+1,js:je,ks-1:ke+1) = ZERO
          wrxy(ie+1,js:je,ks-1:ke+1) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    call bl_allocate(wimhxy,lo(1),hi(1)+1,lo(2),hi(2),lo(3)-1,hi(3)+1)

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! upwind using full velocity
             wimhxy(i,j,k) = merge(wlxy(i,j,k),wrxy(i,j,k),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(wlxy(i,j,k)+wrxy(i,j,k))
             wimhxy(i,j,k) = merge(uavg,wimhxy(i,j,k),abs(utrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(wlxy)
    call bl_deallocate(wrxy)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    call bl_allocate(wlyx,lo(1),hi(1),lo(2),hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(wryx,lo(1),hi(1),lo(2),hi(2)+1,lo(3)-1,hi(3)+1)

    ! wimhyx loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! extrapolate to faces
             wlyx(i,j,k) = uly(i,j,k,3) - (dt6/hx)*(utrans(i+1,j-1,k)+utrans(i,j-1,k)) &
                  * (uimhx(i+1,j-1,k,3)-uimhx(i,j-1,k,3))
             wryx(i,j,k) = ury(i,j,k,3) - (dt6/hx)*(utrans(i+1,j  ,k)+utrans(i,j  ,k)) &
                  * (uimhx(i+1,j  ,k,3)-uimhx(i,j  ,k,3))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(uimhx)

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          wlyx(is:ie,js,ks-1:ke+1) = utilde(is:ie,js-1,ks-1:ke+1,3)
          wryx(is:ie,js,ks-1:ke+1) = utilde(is:ie,js-1,ks-1:ke+1,3)
       case (SlipWall, Symmetry, Outflow)
          wlyx(is:ie,js,ks-1:ke+1) = wryx(is:ie,js,ks-1:ke+1)
       case (NoSlipWall)
          wlyx(is:ie,js,ks-1:ke+1) = ZERO
          wryx(is:ie,js,ks-1:ke+1) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          wlyx(is:ie,je+1,ks-1:ke+1) = utilde(is:ie,je+1,ks-1:ke+1,3)
          wryx(is:ie,je+1,ks-1:ke+1) = utilde(is:ie,je+1,ks-1:ke+1,3)
       case (SlipWall, Symmetry, Outflow)
          wryx(is:ie,je+1,ks-1:ke+1) = wlyx(is:ie,je+1,ks-1:ke+1)
       case (NoSlipWall)
          wlyx(is:ie,je+1,ks-1:ke+1) = ZERO
          wryx(is:ie,je+1,ks-1:ke+1) = ZERO
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    call bl_allocate(wimhyx,lo(1),hi(1),lo(2),hi(2)+1,lo(3)-1,hi(3)+1)

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! upwind using full velocity
             wimhyx(i,j,k) = merge(wlyx(i,j,k),wryx(i,j,k),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(wlyx(i,j,k)+wryx(i,j,k))
             wimhyx(i,j,k) = merge(uavg,wimhyx(i,j,k),abs(vtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(wlyx)
    call bl_deallocate(wryx)

    !******************************************************************
    ! Create umac, etc.
    !******************************************************************

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    call bl_allocate(umacl,lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))
    call bl_allocate(umacr,lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))

    !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! use the traced force if ppm_trace_forces = 1
             fl = merge(force(i-1,j,k,1), Ipfx(i-1,j,k,1), ppm_trace_forces == 0)
             fr = merge(force(i  ,j,k,1), Imfx(i  ,j,k,1), ppm_trace_forces == 0)

             ! extrapolate to edges
             umacl(i,j,k) = ulx(i,j,k,1) &
                  - (dt4/hy)*(vtrans(i-1,j+1,k  )+vtrans(i-1,j,k)) &
                  * (uimhyz(i-1,j+1,k  )-uimhyz(i-1,j,k)) &
                  - (dt4/hz)*(wtrans(i-1,j  ,k+1)+wtrans(i-1,j,k)) &
                  * (uimhzy(i-1,j  ,k+1)-uimhzy(i-1,j,k)) &
                  + dt2*fl
             umacr(i,j,k) = urx(i,j,k,1) &
                  - (dt4/hy)*(vtrans(i  ,j+1,k  )+vtrans(i  ,j,k)) &
                  * (uimhyz(i  ,j+1,k  )-uimhyz(i  ,j,k)) &
                  - (dt4/hz)*(wtrans(i  ,j  ,k+1)+wtrans(i  ,j,k)) &
                  * (uimhzy(i  ,j  ,k+1)-uimhzy(i  ,j,k)) &
                  + dt2*fr
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(ulx)
    call bl_deallocate(urx)
    call bl_deallocate(uimhyz)
    call bl_deallocate(uimhzy)

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! solve Riemann problem using full velocity
                uavg = HALF*(umacl(i,j,k)+umacr(i,j,k))
                test = ((umacl(i,j,k)+w0macx(i,j,k) .le. ZERO .and. &
                     umacr(i,j,k)+w0macx(i,j,k) .ge. ZERO) .or. &
                     (abs(umacl(i,j,k)+umacr(i,j,k)+TWO*w0macx(i,j,k)) .lt. rel_eps))
                umac(i,j,k) = merge(umacl(i,j,k),umacr(i,j,k),uavg+w0macx(i,j,k) .gt. ZERO)
                umac(i,j,k) = merge(ZERO,umac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! solve Riemann problem using full velocity
                uavg = HALF*(umacl(i,j,k)+umacr(i,j,k))
                test = ((umacl(i,j,k) .le. ZERO .and. umacr(i,j,k) .ge. ZERO) .or. &
                     (abs(umacl(i,j,k)+umacr(i,j,k)) .lt. rel_eps))
                umac(i,j,k) = merge(umacl(i,j,k),umacr(i,j,k),uavg .gt. ZERO)
                umac(i,j,k) = merge(ZERO,umac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          umac(is,js:je,ks:ke) = utilde(is-1,js:je,ks:ke,1)
       case (SlipWall, NoSlipWall, Symmetry)
          umac(is,js:je,ks:ke) = ZERO
       case (Outflow)
          umac(is,js:je,ks:ke) = min(umacr(is,js:je,ks:ke),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          umac(ie+1,js:je,ks:ke) = utilde(ie+1,js:je,ks:ke,1)
       case (SlipWall, Symmetry, NoSlipWall)
          umac(ie+1,js:je,ks:ke) = ZERO
       case (Outflow)
          umac(ie+1,js:je,ks:ke) = max(umacl(ie+1,js:je,ks:ke),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    call bl_deallocate(umacl)
    call bl_deallocate(umacr)

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    call bl_allocate(vmacl,lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))
    call bl_allocate(vmacr,lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))

    !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! use the traced force if ppm_trace_forces = 1
             fl = merge(force(i,j-1,k,2), Ipfy(i,j-1,k,2), ppm_trace_forces == 0)
             fr = merge(force(i,j  ,k,2), Imfy(i,j  ,k,2), ppm_trace_forces == 0)

             ! extrapolate to edges
             vmacl(i,j,k) = uly(i,j,k,2) &
                  - (dt4/hx)*(utrans(i+1,j-1,k  )+utrans(i,j-1,k)) &
                  * (vimhxz(i+1,j-1,k  )-vimhxz(i,j-1,k)) &
                  - (dt4/hz)*(wtrans(i  ,j-1,k+1)+wtrans(i,j-1,k)) &
                  * (vimhzx(i  ,j-1,k+1)-vimhzx(i,j-1,k)) &
                  + dt2*fl
             vmacr(i,j,k) = ury(i,j,k,2) &
                  - (dt4/hx)*(utrans(i+1,j  ,k  )+utrans(i,j  ,k)) &
                  * (vimhxz(i+1,j  ,k  )-vimhxz(i,j  ,k)) &
                  - (dt4/hz)*(wtrans(i  ,j  ,k+1)+wtrans(i,j  ,k)) &
                  * (vimhzx(i  ,j  ,k+1)-vimhzx(i,j  ,k)) &
                  + dt2*fr
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(uly)
    call bl_deallocate(ury)
    call bl_deallocate(vimhxz)
    call bl_deallocate(vimhzx)

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(vmacl(i,j,k)+vmacr(i,j,k))
                test = ((vmacl(i,j,k)+w0macy(i,j,k) .le. ZERO .and. &
                     vmacr(i,j,k)+w0macy(i,j,k) .ge. ZERO) .or. &
                     (abs(vmacl(i,j,k)+vmacr(i,j,k)+TWO*w0macy(i,j,k)) .lt. rel_eps))
                vmac(i,j,k) = merge(vmacl(i,j,k),vmacr(i,j,k),uavg+w0macy(i,j,k) .gt. ZERO)
                vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(vmacl(i,j,k)+vmacr(i,j,k))
                test = ((vmacl(i,j,k) .le. ZERO .and. vmacr(i,j,k) .ge. ZERO) .or. &
                     (abs(vmacl(i,j,k)+vmacr(i,j,k)) .lt. rel_eps))
                vmac(i,j,k) = merge(vmacl(i,j,k),vmacr(i,j,k),uavg .gt. ZERO)
                vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          vmac(is:ie,js,ks:ke) = utilde(is:ie,js-1,ks:ke,2)
       case (SlipWall, Symmetry, NoSlipWall)
          vmac(is:ie,js,ks:ke) = ZERO
       case (Outflow)
          vmac(is:ie,js,ks:ke) = min(vmacr(is:ie,js,ks:ke),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          vmac(is:ie,je+1,ks:ke) = utilde(is:ie,je+1,ks:ke,2)
       case (SlipWall, Symmetry, NoSlipWall)
          vmac(is:ie,je+1,ks:ke) = ZERO
       case (Outflow)
          vmac(is:ie,je+1,ks:ke) = max(vmacl(is:ie,je+1,ks:ke),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    call bl_deallocate(vmacl)
    call bl_deallocate(vmacr)

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    call bl_allocate(wmacl,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)
    call bl_allocate(wmacr,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)

    !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! use the traced force if ppm_trace_forces = 1
             fl = merge(force(i,j,k-1,3), Ipfz(i,j,k-1,3), ppm_trace_forces == 0)
             fr = merge(force(i,j,k  ,3), Imfz(i,j,k  ,3), ppm_trace_forces == 0)

             ! extrapolate to edges
             wmacl(i,j,k) = ulz(i,j,k,3) &
                  - (dt4/hx)*(utrans(i+1,j  ,k-1)+utrans(i,j,k-1)) &
                  * (wimhxy(i+1,j  ,k-1)-wimhxy(i,j,k-1)) &
                  - (dt4/hy)*(vtrans(i  ,j+1,k-1)+vtrans(i,j,k-1)) &
                  * (wimhyx(i  ,j+1,k-1)-wimhyx(i,j,k-1)) &
                  + dt2*fl
             wmacr(i,j,k) = urz(i,j,k,3) &
                  - (dt4/hx)*(utrans(i+1,j  ,k  )+utrans(i,j,k  )) &
                  * (wimhxy(i+1,j  ,k  )-wimhxy(i,j,k  )) &
                  - (dt4/hy)*(vtrans(i  ,j+1,k  )+vtrans(i,j,k  )) &
                  * (wimhyx(i  ,j+1,k  )-wimhyx(i,j,k  )) &
                  + dt2*fr
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(ulz)
    call bl_deallocate(urz)
    call bl_deallocate(wimhxy)
    call bl_deallocate(wimhyx)

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(wmacl(i,j,k)+wmacr(i,j,k))
                test = ((wmacl(i,j,k)+w0macz(i,j,k) .le. ZERO .and. &
                     wmacr(i,j,k)+w0macz(i,j,k) .ge. ZERO) .or. &
                     (abs(wmacl(i,j,k)+wmacr(i,j,k)+TWO*w0macz(i,j,k)) .lt. rel_eps))
                wmac(i,j,k) = merge(wmacl(i,j,k),wmacr(i,j,k),uavg+w0macz(i,j,k) .gt. ZERO)
                wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(wmacl(i,j,k)+wmacr(i,j,k))
                test = ((wmacl(i,j,k)+w0(lev,k) .le. ZERO .and. &
                     wmacr(i,j,k)+w0(lev,k) .ge. ZERO) .or. &
                     (abs(wmacl(i,j,k)+wmacr(i,j,k)+TWO*w0(lev,k)) .lt. rel_eps))
                wmac(i,j,k) = merge(wmacl(i,j,k),wmacr(i,j,k),uavg+w0(lev,k) .gt. ZERO)
                wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    ! impose hi side bc's
    if (lo(3) .eq. domlo(3)) then
       select case(phys_bc(3,1))
       case (Inflow)
          wmac(is:ie,js:je,ks) = utilde(is:ie,js:je,ks-1,3)
       case (SlipWall, Symmetry, NoSlipWall)
          wmac(is:ie,js:je,ks) = ZERO
       case (Outflow)
          wmac(is:ie,js:je,ks) = min(wmacr(is:ie,js:je,ks),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
       end select
    end if

    ! impose lo side bc's
    if (hi(3) .eq. domhi(3)) then
       select case(phys_bc(3,2))
       case (Inflow)
          wmac(is:ie,js:je,ke+1) = utilde(is:ie,js:je,ke+1,3)
       case (SlipWall, Symmetry, NoSlipWall)
          wmac(is:ie,js:je,ke+1) = ZERO
       case (Outflow)
          wmac(is:ie,js:je,ke+1) = max(wmacl(is:ie,js:je,ke+1),ZERO)
       case (Interior)
       case  default
          call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
       end select
    end if

    call bl_deallocate(wmacl)
    call bl_deallocate(wmacr)

    call bl_deallocate(Ipfx)
    call bl_deallocate(Imfx)
    call bl_deallocate(Ipfy)
    call bl_deallocate(Imfy)
    call bl_deallocate(Ipfz)
    call bl_deallocate(Imfz)

  end subroutine velpred_3d
#endif

end module velpred_module
