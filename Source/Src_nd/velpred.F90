
#include "AMReX_BC_TYPES.H"

module velpred_module
  ! velpred is called by advance_premac -- it is used to predict the
  ! normal velocities to the interfaces.  We don't care about the
  ! transverse velocities here.  The prediction is done piecewise linear (for now)

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

subroutine velpred_interface_2d(lo, hi, idir, domlo, domhi, &
     utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
     ufull, uf_lo, uf_hi, nc_uf, ng_uf, &
     utrans, uu_lo, uu_hi, &
     Imu, imu_lo, imu_hi, &
     Ipu, ipu_lo, ipu_hi, &
     Imv, imv_lo, imv_hi, &
     Ipv, ipv_lo, ipv_hi, &
     ul, ul_lo, ul_hi, &
     ur, ur_lo, ur_hi, &
     uimh, ui_lo, ui_hi, &
     dx,dt,adv_bc,phys_bc) bind(C,name="velpred_interface_2d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
  integer, value,   intent(in   ) :: idir, ng_ut, nc_ut
  integer         , intent(in   ) :: uf_lo(3), uf_hi(3)
  integer, value,   intent(in   ) :: ng_uf, nc_uf
  integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
  integer         , intent(in   ) :: ipu_lo(3), ipu_hi(3)
  integer         , intent(in   ) :: ipv_lo(3), ipv_hi(3)
  integer         , intent(in   ) :: imu_lo(3), imu_hi(3)
  integer         , intent(in   ) :: imv_lo(3), imv_hi(3)
  integer         , intent(in   ) :: ul_lo(3), ul_hi(3)
  integer         , intent(in   ) :: ur_lo(3), ur_hi(3)
  integer         , intent(in   ) :: ui_lo(3), ui_hi(3)
  double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
  double precision, intent(in   ) :: ufull(uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
  double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
  double precision, intent(in   ) :: Imu (imu_lo(1):imu_hi(1),imu_lo(2):imu_hi(2),imu_lo(3):imu_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipu (ipu_lo(1):ipu_hi(1),ipu_lo(2):ipu_hi(2),ipu_lo(3):ipu_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Imv (imv_lo(1):imv_hi(1),imv_lo(2):imv_hi(2),imv_lo(3):imv_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipv (ipv_lo(1):ipv_hi(1),ipv_lo(2):ipv_hi(2),ipv_lo(3):ipv_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: ul (ul_lo(1):ul_hi(1),ul_lo(2):ul_hi(2),ul_lo(3):ul_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: ur (ur_lo(1):ur_hi(1),ur_lo(2):ur_hi(2),ur_lo(3):ur_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: uimh (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer         , intent(in   ) :: adv_bc(2,2,2), phys_bc(2,2) ! dim, lohi, (comp)


  double precision :: hx, hy, dt2, dt4, uavg, maxu, minu
  double precision :: fl, fr

  integer :: i,j,k

  logical :: test

  !$gpu

  k = lo(3)

  dt2 = HALF*dt
  dt4 = dt/4.0d0

  hx = dx(1)
  hy = dx(2)

  ! NOTE: for ppm_type == 0, slopex == Ipu, slopey == Imv

  !******************************************************************
  ! Create u_{\i-\half\e_x}^x, etc.
  !******************************************************************

  if (idir == 1) then

     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           if (ppm_type .eq. 0) then
              maxu = max(ZERO,ufull(i-1,j,k,1))
              minu = min(ZERO,ufull(i  ,j,k,1))
              ! extrapolate both components of velocity to left face
              ul(i,j,k,1) = utilde(i-1,j,k,1) + (HALF - (dt2/hx)*maxu)*Ipu(i-1,j,k,1)
              ul(i,j,k,2) = utilde(i-1,j,k,2) + (HALF - (dt2/hx)*maxu)*Ipu(i-1,j,k,2)
              ! extrapolate both components of velocity to right face
              ur(i,j,k,1) = utilde(i  ,j,k,1) - (HALF + (dt2/hx)*minu)*Ipu(i  ,j,k,1)
              ur(i,j,k,2) = utilde(i  ,j,k,2) - (HALF + (dt2/hx)*minu)*Ipu(i  ,j,k,2)
           else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
              ! extrapolate both components of velocity to left face
              ul(i,j,k,1) = Ipu(i-1,j,k,1)
              ul(i,j,k,2) = Ipv(i-1,j,k,1)
              ! extrapolate both components of velocity to right face
              ur(i,j,k,1) = Imu(i,j,k,1)
              ur(i,j,k,2) = Imv(i,j,k,1)
           end if

           ! impose lo side bc's
           if (i .eq. domlo(1)) then
              select case(phys_bc(1,1))
              case (Inflow)
                 ul(i,j,k,1:2) = utilde(i-1,j,k,1:2)
                 ur(i,j,k,1:2) = utilde(i-1,j,k,1:2)
              case (SlipWall, Symmetry)
                 ul(i,j,k,1) = ZERO
                 ur(i,j,k,1) = ZERO
                 ul(i,j,k,2) = ur(i,j,k,2)
              case (NoSlipWall)
                 ul(i,j,k,1:2) = ZERO
                 ur(i,j,k,1:2) = ZERO
              case (Outflow)
                 ur(i,j,k,1) = min(ur(i,j,k,1),ZERO)
                 ur(i,j,k,1:2) = ul(i,j,k,1:2)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(1,1)")
#endif
              end select
           end if

           ! impose hi side bc's
           if (i .eq. domhi(1)+1) then
              select case(phys_bc(1,2))
              case (Inflow)
                 ul(i,j,k,1:2) = utilde(i,j,k,1:2)
                 ur(i,j,k,1:2) = utilde(i,j,k,1:2)
              case (SlipWall, Symmetry)
                 ul(i,j,k,1) = ZERO
                 ur(i,j,k,1) = ZERO
                 ur(i,j,k,2) = ul(i,j,k,2)
              case (NoSlipWall)
                 ul(i,j,k,1:2) = ZERO
                 ur(i,j,k,1:2) = ZERO
              case (Outflow)
                 ul(i,j,k,1) = max(ul(i,j,k,1),ZERO)
                 ur(i,j,k,1:2) = ul(i,j,k,1:2)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(1,2)")
#endif
              end select
           end if

           ! No need to compute uimh(:,:,1) since it's equal to utrans-w0
           ! upwind using full velocity to get transverse component of uimhx
           ! Note: utrans already contains w0
           uimh(i,j,k,2) = merge(ul(i,j,k,2),ur(i,j,k,2),utrans(i,j,k).gt.ZERO)
           uavg = HALF*(ul(i,j,k,2)+ur(i,j,k,2))
           uimh(i,j,k,2) = merge(uavg,uimh(i,j,k,2),abs(utrans(i,j,k)).lt.rel_eps)
        enddo
     enddo

  else ! idir == 2

     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           if (ppm_type .eq. 0) then

              maxu = max(ZERO,ufull(i,j-1,k,2))
              minu = min(ZERO,ufull(i,j,k,2))
              ! extrapolate both components of velocity to left face
              ul(i,j,k,1) = utilde(i,j-1,k,1) + (HALF-(dt2/hy)*maxu)*Imv(i,j-1,k,1)
              ul(i,j,k,2) = utilde(i,j-1,k,2) + (HALF-(dt2/hy)*maxu)*Imv(i,j-1,k,2)
              ! extrapolate both components of velocity to right face
              ur(i,j,k,1) = utilde(i,j,k,1) - (HALF+(dt2/hy)*minu)*Imv(i,j,k,1)
              ur(i,j,k,2) = utilde(i,j,k,2) - (HALF+(dt2/hy)*minu)*Imv(i,j,k,2)
           else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
              ! extrapolate both components of velocity to left face
              ul(i,j,k,1) = Ipu(i,j-1,k,2)
              ul(i,j,k,2) = Ipv(i,j-1,k,2)
              ! extrapolate both components of velocity to right face
              ur(i,j,k,1) = Imu(i,j,k,2)
              ur(i,j,k,2) = Imv(i,j,k,2)
           end if

           ! impose lo side bc's
           if (j .eq. domlo(2)) then
              select case(phys_bc(2,1))
              case (Inflow)
                 ul(i,j,k,1:2) = utilde(i,j-1,k,1:2)
                 ur(i,j,k,1:2) = utilde(i,j-1,k,1:2)
              case (SlipWall, Symmetry)
                 ul(i,j,k,1) = ur(i,j,k,1)
                 ul(i,j,k,2) = ZERO
                 ur(i,j,k,2) = ZERO
              case (NoSlipWall)
                 ul(i,j,k,1:2) = ZERO
                 ur(i,j,k,1:2) = ZERO
              case (Outflow)
                 ur(i,j,k,2) = min(ur(i,j,k,2),ZERO)
                 ul(i,j,k,1:2) = ur(i,j,k,1:2)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(2,1)")
#endif
              end select
           end if

           ! impose hi side bc's
           if (j .eq. domhi(2)+1) then
              select case(phys_bc(2,2))
              case (Inflow)
                 ul(i,j,k,1:2) = utilde(i,j,k,1:2)
                 ur(i,j,k,1:2) = utilde(i,j,k,1:2)
              case (SlipWall, Symmetry)
                 ur(i,j,k,1) = ul(i,j,k,1)
                 ul(i,j,k,2) = ZERO
                 ur(i,j,k,2) = ZERO
              case (NoSlipWall)
                 ul(i,j,k,1:2) = ZERO
                 ur(i,j,k,1:2) = ZERO
              case (Outflow)
                 ul(i,j,k,2)   = max(ul(i,j,k,2),ZERO)
                 ur(i,j,k,1:2) = ul(i,j,k,1:2)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(2,2)")
#endif
              end select
           end if
           ! No need to compute uimh(:,:,2) since it's equal to utrans-w0
           ! upwind using full velocity to get transverse component of uimhy
           ! Note: utrans already contains w0
           uimh(i,j,k,1) = merge(ul(i,j,k,1),ur(i,j,k,1),utrans(i,j,k).gt.ZERO)
           uavg = HALF*(ul(i,j,k,1)+ur(i,j,k,1))
           uimh(i,j,k,1) = merge(uavg,uimh(i,j,k,1),abs(utrans(i,j,k)).lt.rel_eps)
        enddo
     enddo

  endif

end subroutine velpred_interface_2d


subroutine velpred_2d(lo, hi, lev, idir, domlo, domhi, &
     utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
     utrans, uu_lo, uu_hi, &
     vtrans, uv_lo, uv_hi, &
     umac  , mu_lo, mu_hi, &
     vmac  , mv_lo, mv_hi, &
     Imf, imf_lo, imf_hi, &
     Ipf, ipf_lo, ipf_hi, &
     ulx, ulx_lo, ulx_hi, &
     urx, urx_lo, urx_hi, &
     uimhx, uix_lo, uix_hi, &
     uly, uly_lo, uly_hi, &
     ury, ury_lo, ury_hi, &
     uimhy, uiy_lo, uiy_hi, &
     force,   f_lo,  f_hi, nc_f, ng_f, &
     w0_cart, w_lo, w_hi, &
     dx,dt,adv_bc,phys_bc) bind(C,name="velpred_2d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
  integer, value,   intent(in   ) :: lev, idir, ng_ut, nc_ut
  integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
  integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
  integer         , intent(in   ) :: mu_lo(3), mu_hi(3)
  integer         , intent(in   ) :: mv_lo(3), mv_hi(3)
  integer         , intent(in   ) :: ipf_lo(3), ipf_hi(3)
  integer         , intent(in   ) :: imf_lo(3), imf_hi(3)
  integer         , intent(in   ) :: ulx_lo(3), ulx_hi(3)
  integer         , intent(in   ) :: urx_lo(3), urx_hi(3)
  integer         , intent(in   ) :: uix_lo(3), uix_hi(3)
  integer         , intent(in   ) :: uly_lo(3), uly_hi(3)
  integer         , intent(in   ) :: ury_lo(3), ury_hi(3)
  integer         , intent(in   ) :: uiy_lo(3), uiy_hi(3)
  integer         , intent(in   ) :: f_lo(3),  f_hi(3)
  integer         , intent(in   ) :: w_lo(3),  w_hi(3)
  integer, value,   intent(in   ) :: ng_f, nc_f
  double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
  double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
  double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
  double precision, intent(inout) :: umac  (mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))
  double precision, intent(inout) :: vmac  (mv_lo(1):mv_hi(1),mv_lo(2):mv_hi(2),mv_lo(3):mv_hi(3))
  double precision, intent(in   ) :: Imf (imf_lo(1):imf_hi(1),imf_lo(2):imf_hi(2),imf_lo(3):imf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipf (ipf_lo(1):ipf_hi(1),ipf_lo(2):ipf_hi(2),ipf_lo(3):ipf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: ulx (ulx_lo(1):ulx_hi(1),ulx_lo(2):ulx_hi(2),ulx_lo(3):ulx_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: urx (urx_lo(1):urx_hi(1),urx_lo(2):urx_hi(2),urx_lo(3):urx_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uimhx (uix_lo(1):uix_hi(1),uix_lo(2):uix_hi(2),uix_lo(3):uix_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uly (uly_lo(1):uly_hi(1),uly_lo(2):uly_hi(2),uly_lo(3):uly_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: ury (ury_lo(1):ury_hi(1),ury_lo(2):ury_hi(2),ury_lo(3):ury_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uimhy (uiy_lo(1):uiy_hi(1),uiy_lo(2):uiy_hi(2),uiy_lo(3):uiy_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: force ( f_lo(1): f_hi(1), f_lo(2): f_hi(2),f_lo(3): f_hi(3),nc_f)
  double precision, intent(in   ) :: w0_cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3): w_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer         , intent(in   ) :: adv_bc(2,2,2), phys_bc(2,2) ! dim, lohi, (comp)

  ! these correspond to umac_L, etc.
  double precision :: umacl,umacr
  double precision :: vmacl,vmacr

  double precision :: hx, hy, dt2, dt4, uavg
  double precision :: fl, fr

  integer :: i,j,k

  logical :: test

  !$gpu

  k = lo(3)

  dt2 = HALF*dt
  dt4 = dt/4.0d0

  hx = dx(1)
  hy = dx(2)


  !******************************************************************
  ! Create umac and vmac
  !******************************************************************

  if (idir == 1) then

     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           ! use the traced force if ppm_trace_forces = 1
           fl = merge(force(i-1,j,k,1), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
           fr = merge(force(i,j,k,1), Imf(i,  j,k,1), ppm_trace_forces == 0)

           ! extrapolate to edges
           umacl = ulx(i,j,k,1) &
                - (dt4/hy)*(vtrans(i-1,j+1,k)+vtrans(i-1,j,k)) &
                * (uimhy(i-1,j+1,k,1)-uimhy(i-1,j,k,1)) + dt2*fl
           umacr = urx(i,j,k,1) &
                - (dt4/hy)*(vtrans(i  ,j+1,k)+vtrans(i  ,j,k)) &
                * (uimhy(i  ,j+1,k,1)-uimhy(i  ,j,k,1)) + dt2*fr

           ! solve Riemann problem using full velocity
           uavg = HALF*(umacl+umacr)
           test = ((umacl .le. ZERO .and. umacr .ge. ZERO) .or. &
                (abs(umacl+umacr) .lt. rel_eps))
           umac(i,j,k) = merge(umacl,umacr,uavg .gt. ZERO)
           umac(i,j,k) = merge(ZERO,umac(i,j,k),test)

           ! impose lo side bc's
           if (i .eq. domlo(1)) then
              select case(phys_bc(1,1))
              case (Inflow)
                 umac(i,j,k) = utilde(i-1,j,k,1)
              case (SlipWall, NoSlipWall, Symmetry)
                 umac(i,j,k) = ZERO
              case (Outflow)
                 umac(i,j,k) = min(umacr,ZERO)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(1,1)")
#endif
              end select
           end if

           ! impose hi side bc's
           if (i .eq. domhi(1)+1) then
              select case(phys_bc(1,2))
              case (Inflow)
                 umac(i,j,k) = utilde(i,j,k,1)
              case (SlipWall, NoSlipWall, Symmetry)
                 umac(i,j,k) = ZERO
              case (Outflow)
                 umac(i,j,k) = max(umacl,ZERO)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(1,2)")
#endif
              end select
           end if
        enddo
     enddo

  else ! idir == 2

     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           ! use the traced force if ppm_trace_forces = 1
           fl = merge(force(i,j-1,k,2), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
           fr = merge(force(i,j,k,2), Imf(i,j,k,2), ppm_trace_forces == 0)

           ! extrapolate to edges
           vmacl = uly(i,j,k,2) &
                - (dt4/hx)*(utrans(i+1,j-1,k)+utrans(i,j-1,k)) &
                * (uimhx(i+1,j-1,k,2)-uimhx(i,j-1,k,2)) + dt2*fl
           vmacr = ury(i,j,k,2) &
                - (dt4/hx)*(utrans(i+1,j,k)+utrans(i,j,k)) &
                * (uimhx(i+1,j,k,2)-uimhx(i,j,k,2)) + dt2*fr

           ! solve Riemann problem using full velocity
           uavg = HALF*(vmacl+vmacr)
           test = ((vmacl+w0_cart(i,j,k,AMREX_SPACEDIM) .le. ZERO .and. vmacr+w0_cart(i,j,k,AMREX_SPACEDIM) .ge. ZERO) .or. &
                (abs(vmacl+vmacr+TWO*w0_cart(i,j,k,AMREX_SPACEDIM)) .lt. rel_eps))
           vmac(i,j,k) = merge(vmacl,vmacr,uavg+w0_cart(i,j,k,AMREX_SPACEDIM) .gt. ZERO)
           vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)

           ! impose lo side bc's
           if (j .eq. domlo(2)) then
              select case(phys_bc(2,1))
              case (Inflow)
                 vmac(i,j,k) = utilde(i,j-1,k,2)
              case (SlipWall, NoSlipWall, Symmetry)
                 vmac(i,j,k) = ZERO
              case (Outflow)
                 vmac(i,j,k) = min(vmacr,ZERO)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(2,1)")
#endif
              end select
           end if

           ! impose hi side bc's
           if (j .eq. domhi(2)+1) then
              select case(phys_bc(2,2))
              case (Inflow)
                 vmac(i,j,k) = utilde(i,j,k,2)
              case (SlipWall, NoSlipWall, Symmetry)
                 vmac(i,j,k) = ZERO
              case (Outflow)
                 vmac(i,j,k) = max(vmacl,ZERO)
              case (Interior)
              case  default
#ifndef AMREX_USE_CUDA
                 call amrex_error("velpred_2d: invalid boundary type phys_bc(2,2)")
#endif
              end select
           end if

        enddo
     enddo

  endif

end subroutine velpred_2d
#endif


#if (AMREX_SPACEDIM == 3)
subroutine velpred_interface_3d(lo, hi, idir, domlo, domhi, &
     utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
     ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
     utrans, uu_lo, uu_hi, &
     Imu, imu_lo, imu_hi, &
     Ipu, ipu_lo, ipu_hi, &
     Imv, imv_lo, imv_hi, &
     Ipv, ipv_lo, ipv_hi, &
     Imw, imw_lo, imw_hi, &
     Ipw, ipw_lo, ipw_hi, &
     ul, ul_lo, ul_hi, &
     ur, ur_lo, ur_hi, &
     uimh, ui_lo, ui_hi, &
     dx,dt,phys_bc) bind(C,name="velpred_interface_3d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
  integer, value,   intent(in   ) :: idir, ng_ut, nc_ut
  integer         , intent(in   ) :: uf_lo(3), uf_hi(3)
  integer, value,   intent(in   ) :: ng_uf, nc_uf
  integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
  integer         , intent(in   ) :: ipu_lo(3), ipu_hi(3)
  integer         , intent(in   ) :: ipv_lo(3), ipv_hi(3)
  integer         , intent(in   ) :: ipw_lo(3), ipw_hi(3)
  integer         , intent(in   ) :: imu_lo(3), imu_hi(3)
  integer         , intent(in   ) :: imv_lo(3), imv_hi(3)
  integer         , intent(in   ) :: imw_lo(3), imw_hi(3)
  integer         , intent(in   ) :: ul_lo(3), ul_hi(3)
  integer         , intent(in   ) :: ur_lo(3), ur_hi(3)
  integer         , intent(in   ) :: ui_lo(3), ui_hi(3)
  double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
  double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
  double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
  double precision, intent(in   ) :: Imu (imu_lo(1):imu_hi(1),imu_lo(2):imu_hi(2),imu_lo(3):imu_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipu (ipu_lo(1):ipu_hi(1),ipu_lo(2):ipu_hi(2),ipu_lo(3):ipu_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Imv (imv_lo(1):imv_hi(1),imv_lo(2):imv_hi(2),imv_lo(3):imv_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipv (ipv_lo(1):ipv_hi(1),ipv_lo(2):ipv_hi(2),ipv_lo(3):ipv_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Imw (imw_lo(1):imw_hi(1),imw_lo(2):imw_hi(2),imw_lo(3):imw_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipw (ipw_lo(1):ipw_hi(1),ipw_lo(2):ipw_hi(2),ipw_lo(3):ipw_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: ul (ul_lo(1):ul_hi(1),ul_lo(2):ul_hi(2),ul_lo(3):ul_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: ur (ur_lo(1):ur_hi(1),ur_lo(2):ur_hi(2),ur_lo(3):ur_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: uimh (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer         , intent(in   ) :: phys_bc(3,2) ! dim, lohi, (comp)

  ! local variables

  double precision :: hx, hy, hz, dt2, dt4, dt6, uavg, maxu, minu

  integer :: i,j,k

  !$gpu

  dt2 = HALF*dt
  dt4 = dt/4.0d0
  dt6 = dt/6.0d0

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  !******************************************************************
  ! Create u_{\i-\half\e_x}^x, etc.
  !******************************************************************

  if (idir == 1) then

     ! normal predictor states
     ! Allocated from lo:hi+1 in the normal direction
     ! lo-1:hi+1 in the transverse directions

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (ppm_type .eq. 0) then
                 maxu = (HALF - dt2*max(ZERO,ufull(i-1,j,k,1))/hx)
                 minu = (HALF + dt2*min(ZERO,ufull(i  ,j,k,1))/hx)

                 ! extrapolate all components of velocity to left face
                 ul(i,j,k,1) = utilde(i-1,j,k,1) + maxu * Ipu(i-1,j,k,1)
                 ul(i,j,k,2) = utilde(i-1,j,k,2) + maxu * Ipu(i-1,j,k,2)
                 ul(i,j,k,3) = utilde(i-1,j,k,3) + maxu * Ipu(i-1,j,k,3)

                 ! extrapolate all components of velocity to right face
                 ur(i,j,k,1) = utilde(i,j,k,1) - minu * Ipu(i,j,k,1)
                 ur(i,j,k,2) = utilde(i,j,k,2) - minu * Ipu(i,j,k,2)
                 ur(i,j,k,3) = utilde(i,j,k,3) - minu * Ipu(i,j,k,3)
              else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                 ! extrapolate all components of velocity to left face
                 ul(i,j,k,1) = Ipu(i-1,j,k,1)
                 ul(i,j,k,2) = Ipv(i-1,j,k,1)
                 ul(i,j,k,3) = Ipw(i-1,j,k,1)

                 ! extrapolate all components of velocity to right face
                 ur(i,j,k,1) = Imu(i,j,k,1)
                 ur(i,j,k,2) = Imv(i,j,k,1)
                 ur(i,j,k,3) = Imw(i,j,k,1)
              end if

              ! impose lo side bc's
              if (i .eq. domlo(1)) then
                 select case(phys_bc(1,1))
                 case (Inflow)
                    ul(i,j,k,1:3) = utilde(i-1,j,k,1:3)
                    ur(i,j,k,1:3) = utilde(i-1,j,k,1:3)
                 case (SlipWall, Symmetry)
                    ul(i,j,k,1) = ZERO
                    ur(i,j,k,1) = ZERO
                    ul(i,j,k,2) = ur(i,j,k,2)
                    ul(i,j,k,3) = ur(i,j,k,3)
                 case (NoSlipWall)
                    ul(i,j,k,1:3) = ZERO
                    ur(i,j,k,1:3) = ZERO
                 case (Outflow)
                    ur(i,j,k,1) = min(ur(i,j,k,1),ZERO)
                    ul(i,j,k,1:3) = ul(i,j,k,1:3)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (i .eq. domhi(1)+1) then
                 select case(phys_bc(1,2))
                 case (Inflow)
                    ul(i,j,k,1:3) = utilde(i,j,k,1:)
                    ur(i,j,k,1:3) = utilde(i,j,k,1:3)
                 case (SlipWall, Symmetry)
                    ul(i,j,k,1) = ZERO
                    ur(i,j,k,1) = ZERO
                    ur(i,j,k,2) = ul(i,j,k,2)
                    ur(i,j,k,3) = ul(i,j,k,3)
                 case (NoSlipWall)
                    ul(i,j,k,1:3) = ZERO
                    ur(i,j,k,1:3) = ZERO
                 case (Outflow)
                    ul(i,j,k,1) = max(ul(i,j,k,1),ZERO)
                    ur(i,j,k,1:3) = ul(i,j,k,1:3)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
#endif
                 end select
              end if

              ! No need to compute uimhx(:,:,:,1) since it's equal to utrans-w0
              ! upwind using full velocity to get transverse components of uimhx
              ! Note: utrans already contains w0
              uimh(i,j,k,2) = merge(ul(i,j,k,2),ur(i,j,k,2),utrans(i,j,k).gt.ZERO)
              uavg = HALF*(ul(i,j,k,2)+ur(i,j,k,2))
              uimh(i,j,k,2) = merge(uavg,uimh(i,j,k,2),abs(utrans(i,j,k)).lt.rel_eps)

              uimh(i,j,k,3) = merge(ul(i,j,k,3),ur(i,j,k,3),utrans(i,j,k).gt.ZERO)
              uavg = HALF*(ul(i,j,k,3)+ur(i,j,k,3))
              uimh(i,j,k,3) = merge(uavg,uimh(i,j,k,3),abs(utrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else if (idir == 2) then

     ! normal predictor states

     ! Allocated from lo:hi+1 in the normal direction
     ! lo-1:hi+1 in the transverse directions

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (ppm_type .eq. 0) then
                 maxu = (HALF - dt2*max(ZERO,ufull(i,j-1,k,2))/hy)
                 minu = (HALF + dt2*min(ZERO,ufull(i,j  ,k,2))/hy)

                 ! extrapolate all components of velocity to left face
                 ul(i,j,k,1) = utilde(i,j-1,k,1) + maxu * Imv(i,j-1,k,1)
                 ul(i,j,k,2) = utilde(i,j-1,k,2) + maxu * Imv(i,j-1,k,2)
                 ul(i,j,k,3) = utilde(i,j-1,k,3) + maxu * Imv(i,j-1,k,3)

                 ! extrapolate all components of velocity to right face
                 ur(i,j,k,1) = utilde(i,j,k,1) - minu * Imv(i,j,k,1)
                 ur(i,j,k,2) = utilde(i,j,k,2) - minu * Imv(i,j,k,2)
                 ur(i,j,k,3) = utilde(i,j,k,3) - minu * Imv(i,j,k,3)

              else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                 ! extrapolate all components of velocity to left face
                 ul(i,j,k,1) = Ipu(i,j-1,k,2)
                 ul(i,j,k,2) = Ipv(i,j-1,k,2)
                 ul(i,j,k,3) = Ipw(i,j-1,k,2)

                 ! extrapolate all components of velocity to right face
                 ur(i,j,k,1) = Imu(i,j,k,2)
                 ur(i,j,k,2) = Imv(i,j,k,2)
                 ur(i,j,k,3) = Imw(i,j,k,2)
              end if

              ! impose lo side bc's
              if (j .eq. domlo(2)) then
                 select case(phys_bc(2,1))
                 case (Inflow)
                    ul(i,j,k,1:3) = utilde(i,j-1,k,1:3)
                    ur(i,j,k,1:3) = utilde(i,j-1,k,1:3)
                 case (SlipWall, Symmetry)
                    ul(i,j,k,1) = ur(i,j,k,1)
                    ul(i,j,k,2) = ZERO
                    ur(i,j,k,2) = ZERO
                    ul(i,j,k,3) = ur(i,j,k,3)
                 case (NoSlipWall)
                    ul(i,j,k,1:3) = ZERO
                    ur(i,j,k,1:3) = ZERO
                 case (Outflow)
                    ur(i,j,k,2) = min(ur(i,j,k,2),ZERO)
                    ul(i,j,k,1:3) = ur(i,j,k,1:3)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (j .eq. domhi(2)+1) then
                 select case(phys_bc(2,2))
                 case (Inflow)
                    ul(i,j,k,1:3) = utilde(i,j,k,1:3)
                    ur(i,j,k,1:3) = utilde(i,j,k,1:3)
                 case (SlipWall, Symmetry)
                    ur(i,j,k,1) = ul(i,j,k,1)
                    ul(i,j,k,2) = ZERO
                    ur(i,j,k,2) = ZERO
                    ur(i,j,k,3) = ul(i,j,k,3)
                 case (NoSlipWall)
                    ul(i,j,k,1:3) = ZERO
                    ur(i,j,k,1:3) = ZERO
                 case (Outflow)
                    ul(i,j,k,2) = max(ul(i,j,k,2),ZERO)
                    ur(i,j,k,1:3) = ul(i,j,k,1:3)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
#endif
                 end select
              end if

              ! No need to compute uimhy(:,:,:,2) since it's equal to vtrans-w0
              ! upwind using full velocity to get transverse components of uimhy
              ! Note: vtrans already contains w0
              uimh(i,j,k,1) = merge(ul(i,j,k,1),ur(i,j,k,1),utrans(i,j,k).gt.ZERO)
              uavg = HALF*(ul(i,j,k,1)+ur(i,j,k,1))
              uimh(i,j,k,1) = merge(uavg,uimh(i,j,k,1),abs(utrans(i,j,k)).lt.rel_eps)

              uimh(i,j,k,3) = merge(ul(i,j,k,3),ur(i,j,k,3),utrans(i,j,k).gt.ZERO)
              uavg = HALF*(ul(i,j,k,3)+ur(i,j,k,3))
              uimh(i,j,k,3) = merge(uavg,uimh(i,j,k,3),abs(utrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else

     ! normal predictor states
     ! Allocated from lo:hi+1 in the normal direction
     ! lo-1:hi+1 in the transverse directions

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (ppm_type .eq. 0) then
                 maxu = (HALF - dt2*max(ZERO,ufull(i,j,k-1,3))/hz)
                 minu = (HALF + dt2*min(ZERO,ufull(i,j,k  ,3))/hz)

                 ! extrapolate all components of velocity to left face
                 ul(i,j,k,1) = utilde(i,j,k-1,1) + maxu * Imw(i,j,k-1,1)
                 ul(i,j,k,2) = utilde(i,j,k-1,2) + maxu * Imw(i,j,k-1,2)
                 ul(i,j,k,3) = utilde(i,j,k-1,3) + maxu * Imw(i,j,k-1,3)

                 ! extrapolate all components of velocity to right face
                 ur(i,j,k,1) = utilde(i,j,k,1) - minu * Imw(i,j,k,1)
                 ur(i,j,k,2) = utilde(i,j,k,2) - minu * Imw(i,j,k,2)
                 ur(i,j,k,3) = utilde(i,j,k,3) - minu * Imw(i,j,k,3)

              else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

                 ! extrapolate all components of velocity to left face
                 ul(i,j,k,1) = Ipu(i,j,k-1,3)
                 ul(i,j,k,2) = Ipv(i,j,k-1,3)
                 ul(i,j,k,3) = Ipw(i,j,k-1,3)

                 ! extrapolate all components of velocity to right face
                 ur(i,j,k,1) = Imu(i,j,k,3)
                 ur(i,j,k,2) = Imv(i,j,k,3)
                 ur(i,j,k,3) = Imw(i,j,k,3)
              end if

              ! impose lo side bc's
              if (k .eq. domlo(3)) then
                 select case(phys_bc(3,1))
                 case (Inflow)
                    ul(i,j,k,1:3) = utilde(i,j,k-1,1:3)
                    ur(i,j,k,1:3) = utilde(i,j,k-1,1:3)
                 case (SlipWall, Symmetry)
                    ul(i,j,k,1) = ur(i,j,k,1)
                    ul(i,j,k,2) = ur(i,j,k,2)
                    ul(i,j,k,3) = ZERO
                    ur(i,j,k,3) = ZERO
                 case (NoSlipWall)
                    ul(i,j,k,1:3) = ZERO
                    ur(i,j,k,1:3) = ZERO
                 case (Outflow)
                    ur(i,j,k,3) = min(ur(i,j,k,3),ZERO)
                    ul(i,j,k,1:3) = ur(i,j,k,1:3)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (k .eq. domhi(3)+1) then
                 select case(phys_bc(3,2))
                 case (Inflow)
                    ul(i,j,k,1:3) = utilde(i,j,k,1:3)
                    ur(i,j,k,1:3) = utilde(i,j,k,1:3)
                 case (SlipWall, Symmetry)
                    ur(i,j,k,1) = ul(i,j,k,1)
                    ur(i,j,k,2) = ul(i,j,k,2)
                    ul(i,j,k,3) = ZERO
                    ur(i,j,k,3) = ZERO
                 case (NoSlipWall)
                    ul(i,j,k,1:3) = ZERO
                    ur(i,j,k,1:3) = ZERO
                 case (Outflow)
                    ul(i,j,k,3) = max(ul(i,j,k,3),ZERO)
                    ur(i,j,k,1:3) = ul(i,j,k,1:3)
                 case (Interior)
                 case default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
#endif
                 end select
              end if

              ! No need to compute uimhz(:,:,:,3) since it's equal to wtrans-w0
              ! upwind using full velocity to get transverse components of uimhz
              ! Note: wtrans already contains w0
              uimh(i,j,k,1) = merge(ul(i,j,k,1),ur(i,j,k,1),utrans(i,j,k).gt.ZERO)
              uavg = HALF*(ul(i,j,k,1)+ur(i,j,k,1))
              uimh(i,j,k,1) = merge(uavg,uimh(i,j,k,1),abs(utrans(i,j,k)).lt.rel_eps)

              uimh(i,j,k,2) = merge(ul(i,j,k,2),ur(i,j,k,2),utrans(i,j,k).gt.ZERO)
              uavg = HALF*(ul(i,j,k,2)+ur(i,j,k,2))
              uimh(i,j,k,2) = merge(uavg,uimh(i,j,k,2),abs(utrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  endif

end subroutine velpred_interface_3d


subroutine velpred_transverse_3d(lo, hi, base_dir, norm_dir, &
     domlo, domhi, &
     utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
     utrans, uu_lo, uu_hi, &
     vtrans, uv_lo, uv_hi, &
     wtrans, uw_lo, uw_hi, &
     ulx, ulx_lo, ulx_hi, &
     urx, urx_lo, urx_hi, &
     uimhx, uix_lo, uix_hi, &
     uly, uly_lo, uly_hi, &
     ury, ury_lo, ury_hi, &
     uimhy, uiy_lo, uiy_hi, &
     ulz, ulz_lo, ulz_hi, &
     urz, urz_lo, urz_hi, &
     uimhz, uiz_lo, uiz_hi, &
     uimh_trans, uit_lo, uit_hi, &
     dx,dt,phys_bc) bind(C,name="velpred_transverse_3d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
  integer, value,   intent(in   ) :: base_dir, norm_dir, ng_ut, nc_ut
  integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
  integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
  integer         , intent(in   ) :: uw_lo(3), uw_hi(3)
  integer         , intent(in   ) :: ulx_lo(3), ulx_hi(3)
  integer         , intent(in   ) :: urx_lo(3), urx_hi(3)
  integer         , intent(in   ) :: uix_lo(3), uix_hi(3)
  integer         , intent(in   ) :: uly_lo(3), uly_hi(3)
  integer         , intent(in   ) :: ury_lo(3), ury_hi(3)
  integer         , intent(in   ) :: uiy_lo(3), uiy_hi(3)
  integer         , intent(in   ) :: ulz_lo(3), ulz_hi(3)
  integer         , intent(in   ) :: urz_lo(3), urz_hi(3)
  integer         , intent(in   ) :: uiz_lo(3), uiz_hi(3)
  integer         , intent(in   ) :: uit_lo(3), uit_hi(3)
  double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
  double precision, intent(in   ) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
  double precision, intent(in   ) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
  double precision, intent(in   ) :: wtrans(uw_lo(1):uw_hi(1),uw_lo(2):uw_hi(2),uw_lo(3):uw_hi(3))
  double precision, intent(in   ) :: ulx (ulx_lo(1):ulx_hi(1),ulx_lo(2):ulx_hi(2),ulx_lo(3):ulx_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: urx (urx_lo(1):urx_hi(1),urx_lo(2):urx_hi(2),urx_lo(3):urx_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uimhx (uix_lo(1):uix_hi(1),uix_lo(2):uix_hi(2),uix_lo(3):uix_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uly (uly_lo(1):uly_hi(1),uly_lo(2):uly_hi(2),uly_lo(3):uly_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: ury (ury_lo(1):ury_hi(1),ury_lo(2):ury_hi(2),ury_lo(3):ury_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uimhy (uiy_lo(1):uiy_hi(1),uiy_lo(2):uiy_hi(2),uiy_lo(3):uiy_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: ulz (ulz_lo(1):ulz_hi(1),ulz_lo(2):ulz_hi(2),ulz_lo(3):ulz_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: urz (urz_lo(1):urz_hi(1),urz_lo(2):urz_hi(2),urz_lo(3):urz_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uimhz (uiz_lo(1):uiz_hi(1),uiz_lo(2):uiz_hi(2),uiz_lo(3):uiz_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: uimh_trans (uit_lo(1):uit_hi(1),uit_lo(2):uit_hi(2),uit_lo(3):uit_hi(3))
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer         , intent(in   ) :: phys_bc(3,2) ! dim, lohi, (comp)

  ! local variables

  ! these correspond to u_L^{y|z}, etc.
  double precision :: ulyz, uryz
  double precision :: ulzy, urzy
  double precision :: vlxz, vrxz
  double precision :: vlzx, vrzx
  double precision :: wlxy, wrxy
  double precision :: wlyx, wryx

  double precision :: hx, hy, hz, dt2, dt4, dt6, uavg

  integer :: i,j,k

  !$gpu

  dt2 = HALF*dt
  dt4 = dt/4.0d0
  dt6 = dt/6.0d0

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  !******************************************************************
  ! Create u_{\i-\half\e_y}^{y|z}, etc.
  !******************************************************************

  if (base_dir == 1 .and. norm_dir == 2) then

     ! transverse states
     ! lo-1:hi+1 in base direction
     ! lo:hi+1 in normal direction
     ! lo:hi in transverse direction
     ! uimhyz loop
     ! do k=lo(3),hi(3)
     !    do j=lo(2),hi(2)+1
     !       do i=lo(1)-1,hi(1)+1
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! extrapolate to faces
              ulyz = uly(i,j,k,1) - (dt6/hz)*(wtrans(i,j-1,k+1)+wtrans(i,j-1,k)) &
                   * (uimhz(i,j-1,k+1,1)-uimhz(i,j-1,k,1))
              uryz = ury(i,j,k,1) - (dt6/hz)*(wtrans(i,j  ,k+1)+wtrans(i,j  ,k)) &
                   * (uimhz(i,j  ,k+1,1)-uimhz(i,j  ,k,1))

              ! impose lo side bc's
              if (j .eq. domlo(2)) then
                 select case(phys_bc(2,1))
                 case (Inflow)
                    ulyz = utilde(i,j-1,k,1)
                    uryz = utilde(i,j-1,k,1)
                 case (SlipWall, Symmetry, Outflow)
                    ulyz = uryz
                 case (NoSlipWall)
                    ulyz = ZERO
                    uryz = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (j .eq. domhi(2)+1) then
                 select case(phys_bc(2,2))
                 case (Inflow)
                    ulyz = utilde(i,j,k,1)
                    uryz = utilde(i,j,k,1)
                 case (SlipWall, Symmetry, Outflow)
                    uryz = ulyz
                 case (NoSlipWall)
                    ulyz = ZERO
                    uryz = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
#endif
                 end select
              end if

              ! upwind using full velocity
              uimh_trans(i,j,k) = merge(ulyz,uryz,vtrans(i,j,k).gt.ZERO)
              uavg = HALF*(ulyz+uryz)
              uimh_trans(i,j,k) = merge(uavg,uimh_trans(i,j,k),abs(vtrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else if (base_dir == 1 .and. norm_dir == 3) then


     ! transverse states
     ! lo-1:hi+1 in base direction
     ! lo:hi+1 in normal direction
     ! lo:hi in transverse direction

     ! uimhzy loop
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! extrapolate to faces
              ulzy = ulz(i,j,k,1) - (dt6/hy)*(vtrans(i,j+1,k-1)+vtrans(i,j,k-1)) &
                   * (uimhy(i,j+1,k-1,1)-uimhy(i,j,k-1,1))
              urzy = urz(i,j,k,1) - (dt6/hy)*(vtrans(i,j+1,k  )+vtrans(i,j,k  )) &
                   * (uimhy(i,j+1,k  ,1)-uimhy(i,j,k  ,1))

              ! impose lo side bc's
              if (k .eq. domlo(3)) then
                 select case(phys_bc(3,1))
                 case (Inflow)
                    ulzy = utilde(i,j,k-1,1)
                    urzy = utilde(i,j,k-1,1)
                 case (SlipWall, Symmetry, Outflow)
                    ulzy = urzy
                 case (NoSlipWall)
                    ulzy = ZERO
                    urzy = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (k .eq. domhi(3)+1) then
                 select case(phys_bc(3,2))
                 case (Inflow)
                    ulzy = utilde(i,j,k,1)
                    urzy = utilde(i,j,k,1)
                 case (SlipWall, Symmetry, Outflow)
                    urzy = ulzy
                 case (NoSlipWall)
                    ulzy = ZERO
                    urzy = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
#endif
                 end select
              end if

              ! upwind using full velocity
              uimh_trans(i,j,k) = merge(ulzy,urzy,wtrans(i,j,k).gt.ZERO)
              uavg = HALF*(ulzy+urzy)
              uimh_trans(i,j,k) = merge(uavg,uimh_trans(i,j,k),abs(wtrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else if (base_dir == 2 .and. norm_dir == 1) then

     ! transverse states
     ! lo-1:hi+1 in base direction
     ! lo:hi+1 in normal direction
     ! lo:hi in transverse direction

     ! vimhxz loop
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! extrapolate to faces
              vlxz = ulx(i,j,k,2) - (dt6/hz)*(wtrans(i-1,j,k+1)+wtrans(i-1,j,k)) &
                   * (uimhz(i-1,j,k+1,2)-uimhz(i-1,j,k,2))
              vrxz = urx(i,j,k,2) - (dt6/hz)*(wtrans(i  ,j,k+1)+wtrans(i  ,j,k)) &
                   * (uimhz(i  ,j,k+1,2)-uimhz(i  ,j,k,2))


              ! impose lo side bc's
              if (i .eq. domlo(1)) then
                 select case(phys_bc(1,1))
                 case (Inflow)
                    vlxz = utilde(i-1,j,k,2)
                    vrxz = utilde(i-1,j,k,2)
                 case (SlipWall, Symmetry, Outflow)
                    vlxz = vrxz
                 case (NoSlipWall)
                    vlxz = ZERO
                    vrxz = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (i .eq. domhi(1)+1) then
                 select case(phys_bc(1,2))
                 case (Inflow)
                    vlxz = utilde(i,j,k,2)
                    vrxz = utilde(i,j,k,2)
                 case (SlipWall, Symmetry, Outflow)
                    vrxz = vlxz
                 case (NoSlipWall)
                    vlxz = ZERO
                    vrxz = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
#endif
                 end select
              end if

              ! upwind using full velocity
              uimh_trans(i,j,k) = merge(vlxz,vrxz,utrans(i,j,k).gt.ZERO)
              uavg = HALF*(vlxz+vrxz)
              uimh_trans(i,j,k) = merge(uavg,uimh_trans(i,j,k),abs(utrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else if (base_dir == 2 .and. norm_dir == 3) then

     ! transverse states
     ! lo-1:hi+1 in base direction
     ! lo:hi+1 in normal direction
     ! lo:hi in transverse direction

     ! vimhzx loop
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! extrapolate to faces
              vlzx = ulz(i,j,k,2) - (dt6/hx)*(utrans(i+1,j,k-1)+utrans(i,j,k-1)) &
                   * (uimhx(i+1,j,k-1,2)-uimhx(i,j,k-1,2))
              vrzx = urz(i,j,k,2) - (dt6/hx)*(utrans(i+1,j,k  )+utrans(i,j,k  )) &
                   * (uimhx(i+1,j,k  ,2)-uimhx(i,j,k  ,2))

              ! impose lo side bc's
              if (k .eq. domlo(3)) then
                 select case(phys_bc(3,1))
                 case (Inflow)
                    vlzx = utilde(i,j,k-1,2)
                    vrzx = utilde(i,j,k-1,2)
                 case (SlipWall, Symmetry, Outflow)
                    vlzx = vrzx
                 case (NoSlipWall)
                    vlzx = ZERO
                    vrzx = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (k .eq. domhi(3)+1) then
                 select case(phys_bc(3,2))
                 case (Inflow)
                    vlzx = utilde(i,j,k,2)
                    vrzx = utilde(i,j,k,2)
                 case (SlipWall, Symmetry, Outflow)
                    vrzx = vlzx
                 case (NoSlipWall)
                    vlzx = ZERO
                    vrzx = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
#endif
                 end select
              end if

              ! upwind using full velocity
              uimh_trans(i,j,k) = merge(vlzx,vrzx,wtrans(i,j,k).gt.ZERO)
              uavg = HALF*(vlzx+vrzx)
              uimh_trans(i,j,k) = merge(uavg,uimh_trans(i,j,k),abs(wtrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else if (base_dir == 3 .and. norm_dir == 1) then

     ! transverse states
     ! lo-1:hi+1 in base direction
     ! lo:hi+1 in normal direction
     ! lo:hi in transverse direction

     ! wimhxy loop
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! extrapolate to faces
              wlxy = ulx(i,j,k,3) - (dt6/hy)*(vtrans(i-1,j+1,k)+vtrans(i-1,j,k)) &
                   * (uimhy(i-1,j+1,k,3)-uimhy(i-1,j,k,3))
              wrxy = urx(i,j,k,3) - (dt6/hy)*(vtrans(i  ,j+1,k)+vtrans(i  ,j,k)) &
                   * (uimhy(i  ,j+1,k,3)-uimhy(i  ,j,k,3))

              ! impose lo side bc's
              if (i .eq. domlo(1)) then
                 select case(phys_bc(1,1))
                 case (Inflow)
                    wlxy = utilde(i-1,j,k,3)
                    wrxy = utilde(i-1,j,k,3)
                 case (SlipWall, Symmetry, Outflow)
                    wlxy = wrxy
                 case (NoSlipWall)
                    wlxy = ZERO
                    wrxy = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (i .eq. domhi(1)+1) then
                 select case(phys_bc(1,2))
                 case (Inflow)
                    wlxy = utilde(i,j,k,3)
                    wrxy = utilde(i,j,k,3)
                 case (SlipWall, Symmetry, Outflow)
                    wrxy = wlxy
                 case (NoSlipWall)
                    wlxy = ZERO
                    wrxy = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
#endif
                 end select
              end if

              ! upwind using full velocity
              uimh_trans(i,j,k) = merge(wlxy,wrxy,utrans(i,j,k).gt.ZERO)
              uavg = HALF*(wlxy+wrxy)
              uimh_trans(i,j,k) = merge(uavg,uimh_trans(i,j,k),abs(utrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  else if (base_dir == 3 .and. norm_dir == 2) then

     ! transverse states
     ! lo-1:hi+1 in base direction
     ! lo:hi+1 in normal direction
     ! lo:hi in transverse direction

     ! wimhyx loop
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! extrapolate to faces
              wlyx = uly(i,j,k,3) - (dt6/hx)*(utrans(i+1,j-1,k)+utrans(i,j-1,k)) &
                   * (uimhx(i+1,j-1,k,3)-uimhx(i,j-1,k,3))
              wryx = ury(i,j,k,3) - (dt6/hx)*(utrans(i+1,j  ,k)+utrans(i,j  ,k)) &
                   * (uimhx(i+1,j  ,k,3)-uimhx(i,j  ,k,3))

              ! impose lo side bc's
              if (j .eq. domlo(2)) then
                 select case(phys_bc(2,1))
                 case (Inflow)
                    wlyx = utilde(i,j-1,k,3)
                    wryx = utilde(i,j-1,k,3)
                 case (SlipWall, Symmetry, Outflow)
                    wlyx = wryx
                 case (NoSlipWall)
                    wlyx = ZERO
                    wryx = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (j .eq. domhi(2)+1) then
                 select case(phys_bc(2,2))
                 case (Inflow)
                    wlyx = utilde(i,j,k,3)
                    wryx = utilde(i,j,k,3)
                 case (SlipWall, Symmetry, Outflow)
                    wryx = wlyx
                 case (NoSlipWall)
                    wlyx = ZERO
                    wryx = ZERO
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
#endif
                 end select
              end if

              ! upwind using full velocity
              uimh_trans(i,j,k) = merge(wlyx,wryx,vtrans(i,j,k).gt.ZERO)
              uavg = HALF*(wlyx+wryx)
              uimh_trans(i,j,k) = merge(uavg,uimh_trans(i,j,k),abs(vtrans(i,j,k)).lt.rel_eps)
           enddo
        enddo
     enddo

  endif
end subroutine velpred_transverse_3d


subroutine velpred_3d(lo, hi, lev, idir, domlo, domhi, &
     utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
     utrans, uu_lo, uu_hi, &
     vtrans, uv_lo, uv_hi, &
     wtrans, uw_lo, uw_hi, &
     umac,   mu_lo, mu_hi, &
     vmac,   mv_lo, mv_hi, &
     wmac,   mw_lo, mw_hi, &
     w0macx, wx_lo, wx_hi, &
     w0macy, wy_lo, wy_hi, &
     w0macz, wz_lo, wz_hi, &
     Imf, imf_lo, imf_hi, &
     Ipf, ipf_lo, ipf_hi, &
     ul, ul_lo, ul_hi, &
     ur, ur_lo, ur_hi, &
     uimhyz, uyz_lo, uyz_hi, &
     uimhzy, uzy_lo, uzy_hi, &
     vimhxz, vxz_lo, vxz_hi, &
     vimhzx, vzx_lo, vzx_hi, &
     wimhxy, wxy_lo, wxy_hi, &
     wimhyx, wyx_lo, wyx_hi, &
     force,   f_lo,  f_hi, nc_f, ng_f, &
     w0_cart, w_lo, w_hi, &
     dx,dt,phys_bc) bind(C,name="velpred_3d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: ut_lo(3), ut_hi(3)
  integer, value,   intent(in   ) :: lev, idir, ng_ut, nc_ut
  integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
  integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
  integer         , intent(in   ) :: uw_lo(3), uw_hi(3)
  integer         , intent(in   ) :: mu_lo(3), mu_hi(3)
  integer         , intent(in   ) :: mv_lo(3), mv_hi(3)
  integer         , intent(in   ) :: mw_lo(3), mw_hi(3)
  integer         , intent(in   ) :: ipf_lo(3), ipf_hi(3)
  integer         , intent(in   ) :: imf_lo(3), imf_hi(3)
  integer         , intent(in   ) :: ul_lo(3), ul_hi(3)
  integer         , intent(in   ) :: ur_lo(3), ur_hi(3)
  integer         , intent(in   ) :: uyz_lo(3), uyz_hi(3)
  integer         , intent(in   ) :: uzy_lo(3), uzy_hi(3)
  integer         , intent(in   ) :: vxz_lo(3), vxz_hi(3)
  integer         , intent(in   ) :: vzx_lo(3), vzx_hi(3)
  integer         , intent(in   ) :: wxy_lo(3), wxy_hi(3)
  integer         , intent(in   ) :: wyx_lo(3), wyx_hi(3)
  integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
  integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
  integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
  integer         , intent(in   ) ::  f_lo(3),  f_hi(3)
  integer         , intent(in   ) :: w_lo(3),  w_hi(3)
  integer, value,   intent(in   ) :: ng_f, nc_f
  double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
  double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
  double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
  double precision, intent(inout) :: wtrans(uw_lo(1):uw_hi(1),uw_lo(2):uw_hi(2),uw_lo(3):uw_hi(3))
  double precision, intent(inout) :: umac  (mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))
  double precision, intent(inout) :: vmac  (mv_lo(1):mv_hi(1),mv_lo(2):mv_hi(2),mv_lo(3):mv_hi(3))
  double precision, intent(inout) :: wmac  (mw_lo(1):mw_hi(1),mw_lo(2):mw_hi(2),mw_lo(3):mw_hi(3))
  double precision, intent(in   ) :: Imf (imf_lo(1):imf_hi(1),imf_lo(2):imf_hi(2),imf_lo(3):imf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: Ipf (ipf_lo(1):ipf_hi(1),ipf_lo(2):ipf_hi(2),ipf_lo(3):ipf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: ul (ul_lo(1):ul_hi(1),ul_lo(2):ul_hi(2),ul_lo(3):ul_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: ur (ur_lo(1):ur_hi(1),ur_lo(2):ur_hi(2),ur_lo(3):ur_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: uimhyz (uyz_lo(1):uyz_hi(1),uyz_lo(2):uyz_hi(2),uyz_lo(3):uyz_hi(3))
  double precision, intent(in   ) :: uimhzy (uzy_lo(1):uzy_hi(1),uzy_lo(2):uzy_hi(2),uzy_lo(3):uzy_hi(3))
  double precision, intent(in   ) :: vimhxz (vxz_lo(1):vxz_hi(1),vxz_lo(2):vxz_hi(2),vxz_lo(3):vxz_hi(3))
  double precision, intent(in   ) :: vimhzx (vzx_lo(1):vzx_hi(1),vzx_lo(2):vzx_hi(2),vzx_lo(3):vzx_hi(3))
  double precision, intent(in   ) :: wimhxy (wxy_lo(1):wxy_hi(1),wxy_lo(2):wxy_hi(2),wxy_lo(3):wxy_hi(3))
  double precision, intent(in   ) :: wimhyx (wyx_lo(1):wyx_hi(1),wyx_lo(2):wyx_hi(2),wyx_lo(3):wyx_hi(3))
  double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
  double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
  double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
  double precision, intent(in   ) :: force ( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3),nc_f)
  double precision, intent(in   ) :: w0_cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3): w_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer         , intent(in   ) :: phys_bc(3,2) ! dim, lohi, (comp)

  ! local variables

  ! these correspond to umac_L, etc.
  double precision :: umacl,umacr
  double precision :: vmacl,vmacr
  double precision :: wmacl,wmacr

  double precision :: hx, hy, hz, dt2, dt4, dt6, uavg
  double precision :: fl, fr

  integer :: i,j,k

  logical :: test

  !$gpu

  dt2 = HALF*dt
  dt4 = dt/4.0d0
  dt6 = dt/6.0d0

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  !******************************************************************
  ! Create umac, etc.
  !******************************************************************

  ! mac states
  ! Allocated from lo:hi+1 in the normal direction
  ! lo:hi in the transverse direction

  if (idir == 1) then

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! use the traced force if ppm_trace_forces = 1
              fl = merge(force(i-1,j,k,1), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
              fr = merge(force(i  ,j,k,1), Imf(i  ,j,k,1), ppm_trace_forces == 0)

              ! extrapolate to edges
              umacl = ul(i,j,k,1) &
                   - (dt4/hy)*(vtrans(i-1,j+1,k  )+vtrans(i-1,j,k)) &
                   * (uimhyz(i-1,j+1,k  )-uimhyz(i-1,j,k)) &
                   - (dt4/hz)*(wtrans(i-1,j  ,k+1)+wtrans(i-1,j,k)) &
                   * (uimhzy(i-1,j  ,k+1)-uimhzy(i-1,j,k)) &
                   + dt2*fl
              umacr = ur(i,j,k,1) &
                   - (dt4/hy)*(vtrans(i  ,j+1,k  )+vtrans(i  ,j,k)) &
                   * (uimhyz(i  ,j+1,k  )-uimhyz(i  ,j,k)) &
                   - (dt4/hz)*(wtrans(i  ,j  ,k+1)+wtrans(i  ,j,k)) &
                   * (uimhzy(i  ,j  ,k+1)-uimhzy(i  ,j,k)) &
                   + dt2*fr

              if (spherical .eq. 1) then
                 ! solve Riemann problem using full velocity
                 uavg = HALF*(umacl+umacr)
                 test = ((umacl+w0macx(i,j,k) .le. ZERO .and. &
                      umacr+w0macx(i,j,k) .ge. ZERO) .or. &
                      (abs(umacl+umacr+TWO*w0macx(i,j,k)) .lt. rel_eps))
                 umac(i,j,k) = merge(umacl,umacr,uavg+w0macx(i,j,k) .gt. ZERO)
                 umac(i,j,k) = merge(ZERO,umac(i,j,k),test)
              else
                 ! solve Riemann problem using full velocity
                 uavg = HALF*(umacl+umacr)
                 test = ((umacl .le. ZERO .and. umacr .ge. ZERO) .or. &
                      (abs(umacl+umacr) .lt. rel_eps))
                 umac(i,j,k) = merge(umacl,umacr,uavg .gt. ZERO)
                 umac(i,j,k) = merge(ZERO,umac(i,j,k),test)

              end if

              ! impose lo side bc's
              if (i .eq. domlo(1)) then
                 select case(phys_bc(1,1))
                 case (Inflow)
                    umac(i,j,k) = utilde(i-1,j,k,1)
                 case (SlipWall, NoSlipWall, Symmetry)
                    umac(i,j,k) = ZERO
                 case (Outflow)
                    umac(i,j,k) = min(umacr,ZERO)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (i .eq. domhi(1)+1) then
                 select case(phys_bc(1,2))
                 case (Inflow)
                    umac(i,j,k) = utilde(i,j,k,1)
                 case (SlipWall, Symmetry, NoSlipWall)
                    umac(i,j,k) = ZERO
                 case (Outflow)
                    umac(i,j,k) = max(umacl,ZERO)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(1,2)")
#endif
                 end select
              end if
           enddo
        enddo
     enddo

  else if (idir == 2) then

     ! mac states
     ! Allocated from lo:hi+1 in the normal direction
     ! lo:hi in the transverse direction

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! use the traced force if ppm_trace_forces = 1
              fl = merge(force(i,j-1,k,2), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
              fr = merge(force(i,j  ,k,2), Imf(i,j  ,k,2), ppm_trace_forces == 0)

              ! extrapolate to edges
              vmacl = ul(i,j,k,2) &
                   - (dt4/hx)*(utrans(i+1,j-1,k  )+utrans(i,j-1,k)) &
                   * (vimhxz(i+1,j-1,k  )-vimhxz(i,j-1,k)) &
                   - (dt4/hz)*(wtrans(i  ,j-1,k+1)+wtrans(i,j-1,k)) &
                   * (vimhzx(i  ,j-1,k+1)-vimhzx(i,j-1,k)) &
                   + dt2*fl
              vmacr = ur(i,j,k,2) &
                   - (dt4/hx)*(utrans(i+1,j  ,k  )+utrans(i,j  ,k)) &
                   * (vimhxz(i+1,j  ,k  )-vimhxz(i,j  ,k)) &
                   - (dt4/hz)*(wtrans(i  ,j  ,k+1)+wtrans(i,j  ,k)) &
                   * (vimhzx(i  ,j  ,k+1)-vimhzx(i,j  ,k)) &
                   + dt2*fr

              if (spherical .eq. 1) then
                 ! solve Riemann problem using full velocity
                 uavg = HALF*(vmacl+vmacr)
                 test = ((vmacl+w0macy(i,j,k) .le. ZERO .and. &
                      vmacr+w0macy(i,j,k) .ge. ZERO) .or. &
                      (abs(vmacl+vmacr+TWO*w0macy(i,j,k)) .lt. rel_eps))
                 vmac(i,j,k) = merge(vmacl,vmacr,uavg+w0macy(i,j,k) .gt. ZERO)
                 vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
              else
                 ! solve Riemann problem using full velocity
                 uavg = HALF*(vmacl+vmacr)
                 test = ((vmacl .le. ZERO .and. vmacr .ge. ZERO) .or. &
                      (abs(vmacl+vmacr) .lt. rel_eps))
                 vmac(i,j,k) = merge(vmacl,vmacr,uavg .gt. ZERO)
                 vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
              end if

              ! impose lo side bc's
              if (j .eq. domlo(2)) then
                 select case(phys_bc(2,1))
                 case (Inflow)
                    vmac(i,j,k) = utilde(i,j-1,k,2)
                 case (SlipWall, Symmetry, NoSlipWall)
                    vmac(i,j,k) = ZERO
                 case (Outflow)
                    vmac(i,j,k) = min(vmacr,ZERO)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,1)")
#endif
                 end select
              end if

              ! impose hi side bc's
              if (j .eq. domhi(2)+1) then
                 select case(phys_bc(2,2))
                 case (Inflow)
                    vmac(i,j,k) = utilde(i,j,k,2)
                 case (SlipWall, Symmetry, NoSlipWall)
                    vmac(i,j,k) = ZERO
                 case (Outflow)
                    vmac(i,j,k) = max(vmacl,ZERO)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(2,2)")
#endif
                 end select
              end if
           enddo
        enddo
     enddo

  else ! idir == 3

     ! mac states
     ! Allocated from lo:hi+1 in the normal direction
     ! lo:hi in the transverse direction

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              ! use the traced force if ppm_trace_forces = 1
              fl = merge(force(i,j,k-1,3), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
              fr = merge(force(i,j,k  ,3), Imf(i,j,k  ,3), ppm_trace_forces == 0)

              ! extrapolate to edges
              wmacl = ul(i,j,k,3) &
                   - (dt4/hx)*(utrans(i+1,j  ,k-1)+utrans(i,j,k-1)) &
                   * (wimhxy(i+1,j  ,k-1)-wimhxy(i,j,k-1)) &
                   - (dt4/hy)*(vtrans(i  ,j+1,k-1)+vtrans(i,j,k-1)) &
                   * (wimhyx(i  ,j+1,k-1)-wimhyx(i,j,k-1)) &
                   + dt2*fl
              wmacr = ur(i,j,k,3) &
                   - (dt4/hx)*(utrans(i+1,j  ,k  )+utrans(i,j,k  )) &
                   * (wimhxy(i+1,j  ,k  )-wimhxy(i,j,k  )) &
                   - (dt4/hy)*(vtrans(i  ,j+1,k  )+vtrans(i,j,k  )) &
                   * (wimhyx(i  ,j+1,k  )-wimhyx(i,j,k  )) &
                   + dt2*fr

              if (spherical .eq. 1) then
                 ! solve Riemann problem using full velocity
                 uavg = HALF*(wmacl+wmacr)
                 test = ((wmacl+w0macz(i,j,k) .le. ZERO .and. &
                      wmacr+w0macz(i,j,k) .ge. ZERO) .or. &
                      (abs(wmacl+wmacr+TWO*w0macz(i,j,k)) .lt. rel_eps))
                 wmac(i,j,k) = merge(wmacl,wmacr,uavg+w0macz(i,j,k) .gt. ZERO)
                 wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)
              else
                 ! solve Riemann problem using full velocity
                 uavg = HALF*(wmacl+wmacr)
                 test = ((wmacl+w0_cart(i,j,k,AMREX_SPACEDIM) .le. ZERO .and. &
                      wmacr+w0_cart(i,j,k,AMREX_SPACEDIM) .ge. ZERO) .or. &
                      (abs(wmacl+wmacr+TWO*w0_cart(i,j,k,AMREX_SPACEDIM)) .lt. rel_eps))
                 wmac(i,j,k) = merge(wmacl,wmacr,uavg+w0_cart(i,j,k,AMREX_SPACEDIM) .gt. ZERO)
                 wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)

              end if

              ! impose hi side bc's
              if (k .eq. domlo(3)) then
                 select case(phys_bc(3,1))
                 case (Inflow)
                    wmac(i,j,k) = utilde(i,j,k-1,3)
                 case (SlipWall, Symmetry, NoSlipWall)
                    wmac(i,j,k) = ZERO
                 case (Outflow)
                    wmac(i,j,k) = min(wmacr,ZERO)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,1)")
#endif
                 end select
              end if

              ! impose lo side bc's
              if (k .eq. domhi(3)+1) then
                 select case(phys_bc(3,2))
                 case (Inflow)
                    wmac(i,j,k) = utilde(i,j,k,3)
                 case (SlipWall, Symmetry, NoSlipWall)
                    wmac(i,j,k) = ZERO
                 case (Outflow)
                    wmac(i,j,k) = max(wmacl,ZERO)
                 case (Interior)
                 case  default
#ifndef AMREX_USE_CUDA
                    call amrex_error("velpred_3d: invalid boundary type phys_bc(3,2)")
#endif
                 end select
              end if

           enddo
        enddo
     enddo

  endif

end subroutine velpred_3d
#endif

end module velpred_module
