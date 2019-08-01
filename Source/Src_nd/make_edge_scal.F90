! make_edge_scal constructs the edge state of a scalar, using a
! second-order Taylor expansion in space (through dx/2) and time
! (though dt/2) (if ppm_type = 0) or using PPM (for ppm_type = 1,2).
!
! We use only MAC-projected edge velocities in this prediction.
!
! We are computing all edge states for each variable.  This is what is
! done for the final updates of the state variables and velocity.  For
! velocity, we should set is_vel = .true.

#include "AMReX_BC_TYPES.H"

module make_edge_scal_module

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use slope_module
  use ppm_module
  use meth_params_module, only: rel_eps, ppm_type, ppm_trace_forces

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 1)
  subroutine make_edge_scal_1d(domlo, domhi, lo, hi, &
       s,      s_lo, s_hi, nc_s, ng_s, &
       sedgex, x_lo, x_hi, nc_x, &
       umac,   u_lo, u_hi, ng_um, &
       force,  f_lo, f_hi, nc_f, &
       dx, dt, is_vel, adv_bc, nbccomp, &
       comp, bccomp, is_conservative) bind(C,name="make_edge_scal_1d")

    integer         , intent(in   ) :: domlo(1), domhi(1), lo(1), hi(1)
    integer         , intent(in   ) :: s_lo(1), s_hi(1)
    integer, value,   intent(in   ) :: nc_s, ng_s
    integer         , intent(in   ) :: x_lo(1), x_hi(1)
    integer, value,   intent(in   ) :: nc_x
    integer         , intent(in   ) :: u_lo(1), u_hi(1)
    integer, value,   intent(in   ) :: ng_um
    integer         , intent(in   ) :: f_lo(1), f_hi(1)
    integer, value,   intent(in   ) :: nc_f
    double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),nc_s)
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),nc_x)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1))
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),nc_f)
    double precision, intent(in   ) :: dx(1)
    double precision, value, intent(in   ) :: dt
    integer, value, intent(in   ) :: is_vel, nbccomp, comp, bccomp, is_conservative
    integer         , intent(in   ) :: adv_bc(1,2,nbccomp)

    ! Local variables
    double precision, pointer :: slopex(:,:)

    double precision :: hx,dt2,dt4,savg,fl,fr

    integer :: i,is,ie

    double precision, pointer :: Ip(:), Ipf(:)
    double precision, pointer :: Im(:), Imf(:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    double precision, pointer:: sedgelx(:),sedgerx(:)

    allocate(Ip(lo(1)-1:hi(1)+1))
    allocate(Im(lo(1)-1:hi(1)+1))

    allocate(Ipf(lo(1)-1:hi(1)+1))
    allocate(Imf(lo(1)-1:hi(1)+1))

    allocate(slopex(lo(1)-1:hi(1)+1,1))

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(sedgelx(lo(1):hi(1)+1))
    allocate(sedgerx(lo(1):hi(1)+1))

    is = lo(1)
    ie = hi(1)

    if (ppm_type .eq. 0) then
       call slopex_1d(s(:,comp:),slopex,domlo,domhi,lo,hi,ng_s,1,adv_bc(:,:,bccomp:))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_1d(s(:,comp),ng_s,umac,ng_um,Ip,Im, &
                    domlo,domhi,lo,hi,adv_bc(:,:,bccomp),dx,dt,.true.)
       if (ppm_trace_forces .eq. 1) then
          call ppm_1d(force(:,comp),ng_s,umac,ng_um,Ipf,Imf, &
                       domlo,domhi,lo,hi,adv_bc(:,:,bccomp),dx,dt,.true.)
       endif
    end if

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! loop over appropriate x-faces
    if (ppm_type .eq. 0) then
       do i=lo(1),hi(1)+1
          ! make sedgelx, sedgerx with 1D extrapolation
          sedgelx(i) = s(i-1,comp) + (HALF - dt2*umac(i)/hx)*slopex(i-1,1)
          sedgerx(i) = s(i  ,comp) - (HALF + dt2*umac(i)/hx)*slopex(i  ,1)
       enddo
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do i=lo(1),hi(1)+1
          ! make sedgelx, sedgerx with 1D extrapolation
          sedgelx(i) = Ip(i-1)
          sedgerx(i) = Im(i  )
       end do
    end if

    ! loop over appropriate x-faces
    do i=lo(1),hi(1)+1
       ! make sedgelx, sedgerx
       fl = merge(force(i-1,comp), Ipf(i-1), ppm_trace_forces == 0)
       fr = merge(force(i  ,comp), Imf(i  ), ppm_trace_forces == 0)

       if(is_conservative .eq. 1) then
          sedgelx(i) = sedgelx(i) &
               - (dt2/hx)*s(i-1,comp)*(umac(i  )-umac(i-1)) &
               + dt2*fl
          sedgerx(i) = sedgerx(i) &
               - (dt2/hx)*s(i  ,comp)*(umac(i+1)-umac(i  )) &
               + dt2*fr
       else
          sedgelx(i) = sedgelx(i) + dt2*fl
          sedgerx(i) = sedgerx(i) + dt2*fr
       end if

       ! make sedgex by solving Riemann problem
       ! boundary conditions enforced outside of i loop
       sedgex(i,comp) = merge(sedgelx(i),sedgerx(i),umac(i) .gt. 0.d0)
       savg = HALF*(sedgelx(i)+sedgerx(i))
       sedgex(i,comp) = merge(sedgex(i,comp),savg,abs(umac(i)) .gt. rel_eps)
    enddo

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          sedgex(lo(1),comp) = s(lo(1)-1,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1) then
             sedgex(lo(1),comp) = min(sedgerx(lo(1)),0.d0)
          else
             sedgex(lo(1),comp) = sedgerx(lo(1))
          end if
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          sedgex(lo(1),comp) = sedgerx(lo(1))
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          sedgex(lo(1),comp) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
#ifndef AMREX_USE_GPU
          call amrex_error("make_edge_scal_1d: invalid boundary type adv_bc(1,1)")
#endif
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          sedgex(hi(1)+1,comp) = s(hi(1)+1,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1) then
             sedgex(hi(1)+1,comp) = max(sedgelx(hi(1)+1),0.d0)
          else
             sedgex(hi(1)+1,comp) = sedgelx(hi(1)+1)
          end if
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          sedgex(hi(1)+1,comp) = sedgelx(hi(1)+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          sedgex(hi(1)+1,comp) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
#ifndef AMREX_USE_GPU
          call amrex_error("make_edge_scal_1d: invalid boundary type adv_bc(1,2)")
#endif
       end if
    end if

    deallocate(Ip)
    deallocate(Im)
    deallocate(Ipf)
    deallocate(Imf)
    deallocate(slopex)
    deallocate(sedgelx)
    deallocate(sedgerx)

  end subroutine make_edge_scal_1d
#endif


#if (AMREX_SPACEDIM == 2)
subroutine make_edge_scal_predictor_2d(lo, hi, idir, domlo, domhi, &
     s,      s_lo, s_hi, nc_s, ng_s, &
     umac,   u_lo, u_hi, &
     vmac,   v_lo, v_hi,&
     Ip, ip_lo, ip_hi, &
     Im, im_lo, im_hi, &
     sl, sl_lo, sl_hi, &
     sr, sr_lo, sr_hi, &
     simh, si_lo, si_hi, &
     dx, dt, is_vel, adv_bc, nbccomp, &
     comp, bccomp) bind(C,name="make_edge_scal_predictor_2d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)
  integer, value,   intent(in   ) :: idir, nc_s, ng_s
  integer         , intent(in   ) :: u_lo(3), u_hi(3)
  integer         , intent(in   ) :: v_lo(3), v_hi(3)
  integer         , intent(in   ) :: ip_lo(3), ip_hi(3)
  integer         , intent(in   ) :: im_lo(3), im_hi(3)
  integer         , intent(in   ) :: sl_lo(3), sl_hi(3)
  integer         , intent(in   ) :: sr_lo(3), sr_hi(3)
  integer         , intent(in   ) :: si_lo(3), si_hi(3)
  double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),1:nc_s)
  double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: sl     (sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3))
  double precision, intent(inout) :: sr     (sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3))
  double precision, intent(inout) :: simh     (si_lo(1):si_hi(1),si_lo(2):si_hi(2),si_lo(3):si_hi(3))
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer, value, intent(in   ) :: is_vel, nbccomp, comp, bccomp
  integer         , intent(in   ) :: adv_bc(2,2,nbccomp)

  ! Local variables

  double precision :: hx,hy,dt2,dt4,savg

  integer :: i,j,k

  !$gpu

  k = s_lo(3)

  dt2 = HALF*dt
  dt4 = dt/4.0d0

  hx = dx(1)
  hy = dx(2)

  !******************************************************************
  ! Create s_{\i-\half\e_x}^x, etc.
  !******************************************************************

  if (idir == 1) then

     ! loop over appropriate x-faces
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           if (ppm_type .eq. 0) then
              ! make slx, srx with 1D extrapolation
              sl(i,j,k) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*Ip(i-1,j,k,1)
              sr(i,j,k) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*Ip(i  ,j,k,1)
           else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
              ! make slx, srx with 1D extrapolation
              sl(i,j,k) = Ip(i-1,j,k,1)
              sr(i,j,k) = Im(i  ,j,k,1)
           end if

           ! impose lo side bc's
           if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
              if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
                 sl(i,j,k) = s(i-1,j,k,comp)
                 sr(i,j,k) = s(i-1,j,k,comp)
              else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 1) then
                    sr(i,j,k) = min(sr(i,j,k),0.d0)
                 end if
                 sl(i,j,k) = sr(i,j,k)
              else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
                 sl(i,j,k) = sr(i,j,k)
              else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
                 sl(i,j,k) = 0.d0
                 sr(i,j,k) = 0.d0
              else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,1)")
#endif
              end if
           end if

           ! impose hi side bc's
           if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
              if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
                 sl(i,j,k) = s(i,j,k,comp)
                 sr(i,j,k) = s(i,j,k,comp)
              else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 1) then
                    sl(i,j,k) = max(sl(i,j,k),0.d0)
                 end if
                 sr(i,j,k) = sl(i,j,k)
              else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
                 sr(i,j,k) = sl(i,j,k)
              else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
                 sl(i,j,k) = 0.d0
                 sr(i,j,k) = 0.d0
              else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,2)")
#endif
              end if
           end if

           ! make simhx by solving Riemann problem
           simh(i,j,k) = merge(sl(i,j,k),sr(i,j,k),umac(i,j,k) .gt. 0.d0)
           savg = HALF*(sl(i,j,k)+sr(i,j,k))
           simh(i,j,k) = merge(simh(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
        enddo
     enddo

  else

     ! loop over appropriate y-faces
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           if (ppm_type .eq. 0) then
              ! make sly, sry with 1D extrapolation
              sl(i,j,k) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*Im(i,j-1,k,1)
              sr(i,j,k) = s(i,j,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*Im(i,j,k,1)
           else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
              ! make sly, sry with 1D extrapolation
              sl(i,j,k) = Ip(i,j-1,k,2)
              sr(i,j,k) = Im(i,j,k,2)
           end if

           ! impose lo side bc's
           if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
              if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
                 sl(i,j,k) = s(i,j-1,k,comp)
                 sr(i,j,k) = s(i,j-1,k,comp)
              else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 2) then
                    sr(i,j,k) = min(sr(i,j,k),0.d0)
                 end if
                 sl(i,j,k) = sr(i,j,k)
              else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
                 sl(i,j,k) = sr(i,j,k)
              else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
                 sl(i,j,k) = 0.d0
                 sr(i,j,k) = 0.d0
              else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,1)")
#endif
              end if
           end if

           ! impose hi side bc's
           if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
              if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
                 sl(i,j,k) = s(i,j,k,comp)
                 sr(i,j,k) = s(i,j,k,comp)
              else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 2) then
                    sl(i,j,k) = max(sl(i,j,k),0.d0)
                 end if
                 sr(i,j,k) = sl(i,j,k)
              else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
                 sr(i,j,k) = sl(i,j,k)
              else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
                 sl(i,j,k) = 0.d0
                 sr(i,j,k) = 0.d0
              else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,2)")
#endif
              end if
           end if

           ! make simhy by solving Riemann problem
           simh(i,j,k) = merge(sl(i,j,k),sr(i,j,k),vmac(i,j,k) .gt. 0.d0)
           savg = HALF*(sl(i,j,k)+sr(i,j,k))
           simh(i,j,k) = merge(simh(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
        enddo
     enddo

  endif

end subroutine make_edge_scal_predictor_2d

subroutine make_edge_scal_2d(lo, hi, idir, domlo, domhi, &
     s,      s_lo, s_hi, nc_s, ng_s, &
     sedge, x_lo, x_hi, nc_x, &
     umac,   u_lo, u_hi, &
     vmac,   v_lo, v_hi,&
     Ipf, ipf_lo, ipf_hi, &
     Imf, imf_lo, imf_hi, &
     sl, sl_lo, sl_hi, &
     sr, sr_lo, sr_hi, &
     simh, si_lo, si_hi, &
     force,  f_lo, f_hi, nc_f, &
     dx, dt, is_vel, adv_bc, nbccomp, &
     comp, bccomp, is_conservative) bind(C,name="make_edge_scal_2d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)
  integer, value,   intent(in   ) :: idir, nc_s, ng_s
  integer         , intent(in   ) :: x_lo(3), x_hi(3)
  integer, value,   intent(in   ) :: nc_x
  integer         , intent(in   ) :: u_lo(3), u_hi(3)
  integer         , intent(in   ) :: v_lo(3), v_hi(3)
  integer         , intent(in   ) :: ipf_lo(3), ipf_hi(3)
  integer         , intent(in   ) :: imf_lo(3), imf_hi(3)
  integer         , intent(in   ) :: sl_lo(3), sl_hi(3)
  integer         , intent(in   ) :: sr_lo(3), sr_hi(3)
  integer         , intent(in   ) :: si_lo(3), si_hi(3)
  integer         , intent(in   ) :: f_lo(3), f_hi(3)
  integer, value,   intent(in   ) :: nc_f
  double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),1:nc_s)
  double precision, intent(inout) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),1:nc_x)
  double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(in) :: Ipf(ipf_lo(1):ipf_hi(1),ipf_lo(2):ipf_hi(2),ipf_lo(3):ipf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in) :: Imf(imf_lo(1):imf_hi(1),imf_lo(2):imf_hi(2),imf_lo(3):imf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: sl     (sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3))
  double precision, intent(in   ) :: sr     (sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3))
  double precision, intent(in   ) :: simh     (si_lo(1):si_hi(1),si_lo(2):si_hi(2),si_lo(3):si_hi(3))
  double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer, value, intent(in   ) :: is_vel, nbccomp, comp, bccomp, is_conservative
  integer         , intent(in   ) :: adv_bc(2,2,nbccomp)

  ! Local variables

  double precision :: hx,hy,dt2,dt4,savg,fl,fr

  integer :: i,j,k

  ! these correspond to \mathrm{sedge}_L^x, etc.
  double precision :: sedgel,sedger

  !$gpu

  k = lo(3)

  dt2 = HALF*dt
  dt4 = dt/4.0d0

  hx = dx(1)
  hy = dx(2)

  !******************************************************************
  ! Create sedgel, etc.
  !******************************************************************

  if (idir == 1) then

     ! loop over appropriate x-faces
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           ! make sedgel, sedger
           fl = merge(force(i-1,j,k,comp), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
           fr = merge(force(i  ,j,k,comp), Imf(i  ,j,k,1), ppm_trace_forces == 0)

           if(is_conservative .eq. 1) then
              sedgel = sl(i,j,k) &
                   - (dt2/hy)*(simh(i-1,j+1,k)*vmac(i-1,j+1,k) - simh(i-1,j,k)*vmac(i-1,j,k)) &
                   - (dt2/hx)*s(i-1,j,k,comp)*(umac(i  ,j,k)-umac(i-1,j,k)) &
                   + dt2*fl
              sedger = sr(i,j,k) &
                   - (dt2/hy)*(simh(i  ,j+1,k)*vmac(i  ,j+1,k) - simh(i  ,j,k)*vmac(i  ,j,k)) &
                   - (dt2/hx)*s(i  ,j,k,comp)*(umac(i+1,j,k)-umac(i  ,j,k)) &
                   + dt2*fr
           else
              sedgel = sl(i,j,k) &
                   - (dt4/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simh(i-1,j+1,k)-simh(i-1,j,k)) &
                   + dt2*fl
              sedger = sr(i,j,k) &
                   - (dt4/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k))*(simh(i  ,j+1,k)-simh(i  ,j,k)) &
                   + dt2*fr
           end if

           ! make sedgex by solving Riemann problem
           ! boundary conditions enforced outside of i,j loop
           sedge(i,j,k,comp) = merge(sedgel,sedger,umac(i,j,k) .gt. 0.d0)
           savg = HALF*(sedgel+sedger)
           sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(umac(i,j,k)) .gt. rel_eps)

           ! impose lo side bc's
           if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
              if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
                 sedge(i,j,k,comp) = s(i-1,j,k,comp)
              else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 1) then
                    sedge(i,j,k,comp) = min(sedger,0.d0)
                 else
                    sedge(i,j,k,comp) = sedger
                 end if
              else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
                 sedge(i,j,k,comp) = sedger
              else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
                 sedge(i,j,k,comp) = 0.d0
              else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,1)")
#endif
              end if
           end if

           ! impose hi side bc's
           if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
              if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
                 sedge(i,j,k,comp) = s(i,j,k,comp)
              else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 1) then
                    sedge(i,j,k,comp) = max(sedgel,0.d0)
                 else
                    sedge(i,j,k,comp) = sedgel
                 end if
              else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
                 sedge(i,j,k,comp) = sedgel
              else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
                 sedge(i,j,k,comp) = 0.d0
              else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,2)")
#endif
              end if
           end if
        enddo
     enddo

  else

     ! loop over appropriate y-faces
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           ! make sedgel, sedger
           fl = merge(force(i,j-1,k,comp), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
           fr = merge(force(i,j,k,comp), Imf(i,j,k,2), ppm_trace_forces == 0)

           if(is_conservative .eq. 1) then
              sedgel = sl(i,j,k) &
                   - (dt2/hx)*(simh(i+1,j-1,k)*umac(i+1,j-1,k) - simh(i,j-1,k)*umac(i,j-1,k)) &
                   - (dt2/hy)*s(i,j-1,k,comp)*(vmac(i,j,k)-vmac(i,j-1,k)) &
                   + dt2*fl
              sedger = sr(i,j,k) &
                   - (dt2/hx)*(simh(i+1,j,k)*umac(i+1,j,k) - simh(i,j,k)*umac(i,j,k)) &
                   - (dt2/hy)*s(i,j,k,comp)*(vmac(i,j+1,k)-vmac(i,j,k)) &
                   + dt2*fr
           else
              sedgel = sl(i,j,k) &
                   - (dt4/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simh(i+1,j-1,k)-simh(i,j-1,k)) &
                   + dt2*fl
              sedger = sr(i,j,k) &
                   - (dt4/hx)*(umac(i+1,j,k)+umac(i,j,k))*(simh(i+1,j,k)-simh(i,j,k)) &
                   + dt2*fr
           end if

           ! make sedgey by solving Riemann problem
           ! boundary conditions enforced outside of i,j loop
           sedge(i,j,k,comp) = merge(sedgel,sedger,vmac(i,j,k) .gt. 0.d0)
           savg = HALF*(sedgel+sedger)
           sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(vmac(i,j,k)) .gt. rel_eps)

           ! impose lo side bc's
           if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
              if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
                 sedge(i,j,k,comp) = s(i,j-1,k,comp)
              else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 2) then
                    sedge(i,j,k,comp) = min(sedger,0.d0)
                 else
                    sedge(i,j,k,comp) = sedger
                 end if
              else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
                 sedge(i,j,k,comp) = sedger
              else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
                 sedge(i,j,k,comp) = 0.d0
              else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,1)")
#endif
              end if
           end if

           ! impose hi side bc's
           if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
              if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
                 sedge(i,j,k,comp) = s(i,j,k,comp)
              else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
                   adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                 if (is_vel .eq. 1 .and. comp .eq. 2) then
                    sedge(i,j,k,comp) = max(sedgel,0.d0)
                 else
                    sedge(i,j,k,comp) = sedgel
                 end if
              else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
                 sedge(i,j,k,comp) = sedgel
              else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
                 sedge(i,j,k,comp) = 0.d0
              else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
              else
#ifndef AMREX_USE_GPU
                 call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,2)")
#endif
              end if
           end if
        enddo
     enddo

  endif

end subroutine make_edge_scal_2d
#endif

#if (AMREX_SPACEDIM == 3)
subroutine make_divu(lo, hi, &
     divu, d_lo, d_hi, &
     umac,   u_lo, u_hi, &
     vmac,   v_lo, v_hi, &
     wmac,   w_lo, w_hi, &
     dx, is_conservative) bind(C,name="make_divu")

  integer         , intent(in   ) :: lo(3), hi(3)
  integer         , intent(in   ) :: d_lo(3), d_hi(3)
  integer         , intent(in   ) :: u_lo(3), u_hi(3)
  integer         , intent(in   ) :: v_lo(3), v_hi(3)
  integer         , intent(in   ) :: w_lo(3), w_hi(3)
  double precision, intent(inout) :: divu(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
  double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
  double precision, intent(in   ) :: dx(3)
  integer, value, intent(in   ) :: is_conservative

  integer :: i,j,k

  !$gpu

  if (is_conservative .eq. 1) then
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              divu(i,j,k) = (  umac(i+1,j,k)-umac(i,j,k) &
                   + vmac(i,j+1,k)-vmac(i,j,k) &
                   + wmac(i,j,k+1)-wmac(i,j,k) ) / dx(1)
           end do
        end do
     end do
  end if

end subroutine make_divu


subroutine make_edge_scal_predictor_3d(lo, hi, idir, domlo, domhi, &
     s,      s_lo, s_hi, nc_s, &
     umac,   u_lo, u_hi, &
     vmac,   v_lo, v_hi, &
     wmac,   w_lo, w_hi, &
     Ip, ip_lo, ip_hi, &
     Im, im_lo, im_hi, &
     slopez, slo_lo, slo_hi, &
     sl, sl_lo, sl_hi, &
     sr, sr_lo, sr_hi, &
     simh, si_lo, si_hi, &
     dx, dt, is_vel, adv_bc, nbccomp, &
     comp, bccomp) bind(C,name="make_edge_scal_predictor_3d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)
  integer, value,   intent(in   ) :: idir, nc_s
  integer         , intent(in   ) :: u_lo(3), u_hi(3)
  integer         , intent(in   ) :: v_lo(3), v_hi(3)
  integer         , intent(in   ) :: w_lo(3), w_hi(3)
  integer         , intent(in   ) :: ip_lo(3), ip_hi(3)
  integer         , intent(in   ) :: im_lo(3), im_hi(3)
  integer         , intent(in   ) :: slo_lo(3), slo_hi(3)
  integer         , intent(in   ) :: sl_lo(3), sl_hi(3)
  integer         , intent(in   ) :: sr_lo(3), sr_hi(3)
  integer         , intent(in   ) :: si_lo(3), si_hi(3)
  double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
  double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
  double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),AMREX_SPACEDIM)
  double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: slopez     (slo_lo(1):slo_hi(1),slo_lo(2):slo_hi(2),slo_lo(3):slo_hi(3))
  double precision, intent(inout) :: sl     (sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3))
  double precision, intent(inout) :: sr     (sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3))
  double precision, intent(inout) :: simh     (si_lo(1):si_hi(1),si_lo(2):si_hi(2),si_lo(3):si_hi(3))
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer, value, intent(in   ) :: is_vel, nbccomp, comp, bccomp
  integer         , intent(in   ) :: adv_bc(3,2,nbccomp)

  ! Local variables

  double precision :: hx,hy,hz,dt2,dt3,dt4,dt6,fl,fr
  double precision :: savg

  integer :: i,j,k

  !$gpu

  dt2 = HALF*dt
  dt3 = dt/3.0d0
  dt4 = dt/4.0d0
  dt6 = dt/6.0d0

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  !******************************************************************
  ! Create s_{\i-\half\e_x}^x, etc.
  !******************************************************************

  if (idir == 1) then
     ! Normal predictor states.
     ! call bl_allocated from lo:hi+1 in the normal direction
     ! lo-1:hi+1 in the transverse directions

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate x-faces
              if (ppm_type .eq. 0) then
                 ! mahi(3) slx, srx with 1D extrapolation
                 sl(i,j,k) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*Ip(i-1,j,k,1)
                 sr(i,j,k) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*Ip(i  ,j,k,1)
              else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                 ! make slx, srx with 1D extrapolation
                 sl(i,j,k) = Ip(i-1,j,k,1)
                 sr(i,j,k) = Im(i  ,j,k,1)
              end if

              ! impose lo side bc's
              if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                 if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
                    sl(i,j,k) = s(i-1,j,k,comp)
                    sr(i,j,k) = s(i-1,j,k,comp)
                 else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       sr(i,j,k) = min(sr(i,j,k),0.d0)
                    end if
                    sl(i,j,k) = sr(i,j,k)
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
                    sl(i,j,k) = sr(i,j,k)
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
                    sl(i,j,k) = 0.d0
                    sr(i,j,k) = 0.d0
                 else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                 if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
                    sl(i,j,k) = s(i,j,k,comp)
                    sr(i,j,k) = s(i,j,k,comp)
                 else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       sl(i,j,k) = max(sl(i,j,k),0.d0)
                    end if
                    sr(i,j,k) = sl(i,j,k)
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
                    sr(i,j,k) = sl(i+1,j,k)
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
                    sl(i,j,k) = 0.d0
                    sr(i,j,k) = 0.d0
                 else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
#endif
                 end if
              end if

              ! make simhx by solving Riemann problem
              simh(i,j,k) = merge(sl(i,j,k),sr(i,j,k),umac(i,j,k) .gt. 0.d0)
              savg = HALF*(sl(i,j,k)+sr(i,j,k))
              simh(i,j,k) = merge(simh(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else if (idir == 2) then

     ! Normal predictor states.
     ! call bl_allocated from lo:hi+1 in the normal direction
     ! lo-1:hi+1 in the transverse directions

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate y-faces
              if (ppm_type .eq. 0) then

                 ! make sly, sry with 1D extrapolation
                 sl(i,j,k) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*Im(i,j-1,k,1)
                 sr(i,j,k) = s(i,j  ,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*Im(i,j  ,k,1)
              else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                 ! make sly, sry with 1D extrapolation
                 sl(i,j,k) = Ip(i,j-1,k,2)
                 sr(i,j,k) = Im(i,j  ,k,2)
              end if

              ! impose lo side bc's
              if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                 if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
                    sl(i,j,k) = s(i,j-1,k,comp)
                    sr(i,j,k) = s(i,j-1,k,comp)
                 else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       sr(i,j,k) = min(sr(i,j,k),0.d0)
                    end if
                    sl(i,j,k) = sr(i,j,k)
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
                    sl(i,j,k) = sr(i,j,k)
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
                    sl(i,j,k) = 0.d0
                    sr(i,j,k) = 0.d0
                 else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                 if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
                    sl(i,j,k) = s(i,j,k,comp)
                    sr(i,j,k) = s(i,j,k,comp)
                 else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       sl(i,j,k) = max(sl(i,j,k),0.d0)
                    end if
                    sr(i,j,k) = sl(i,j,k)
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
                    sr(i,j,k) = sl(i,j,k)
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
                    sl(i,j,k) = 0.d0
                    sr(i,j,k) = 0.d0
                 else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
#endif
                 end if
              end if

              ! make simhy by solving Riemann problem
              simh(i,j,k) = merge(sl(i,j,k),sr(i,j,k),vmac(i,j,k) .gt. 0.d0)
              savg = HALF*(sl(i,j,k)+sr(i,j,k))
              simh(i,j,k) = merge(simh(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else ! idir == 3
     ! Normal predictor states.
     ! call bl_allocated from lo:hi+1 in the normal direction
     ! lo-1:hi+1 in the transverse directions

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate z-faces
              if (ppm_type .eq. 0) then
                 ! make slz, srz with 1D extrapolation
                 sl(i,j,k) = s(i,j,k-1,comp) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1)
                 sr(i,j,k) = s(i,j,k  ,comp) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k)
              else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
                 ! make slz, srz with 1D extrapolation
                 sl(i,j,k) = Ip(i,j,k-1,3)
                 sr(i,j,k) = Im(i,j,k  ,3)
              end if

              ! impose lo side bc's
              if (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                 if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
                    sl(i,j,k) = s(i,j,k,comp)
                    sr(i,j,k) = s(i,j,k,comp)
                 else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       sr(i,j,k) = min(sr(i,j,k),0.d0)
                    end if
                    sl(i,j,k) = sr(i,j,k)
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
                    sl(i,j,k) = sr(i,j,k)
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
                    sl(i,j,k) = 0.d0
                    sr(i,j,k) = 0.d0
                 else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (k .eq. hi(3) .and. hi(3)-1 .eq. domhi(3)) then
                 if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
                    sl(i,j,k) = s(i,j,k,comp)
                    sr(i,j,k) = s(i,j,k,comp)
                 else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       sl(i,j,k) = max(sl(i,j,k),0.d0)
                    end if
                    sr(i,j,k) = sl(i,j,k)
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
                    sr(i,j,k) = sl(i,j,k)
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
                    sl(i,j,k) = 0.d0
                    sr(i,j,k) = 0.d0
                 else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
#endif
                 end if
              end if

              ! make simhz by solving Riemann problem
              simh(i,j,k) = merge(sl(i,j,k),sr(i,j,k),wmac(i,j,k) .gt. 0.d0)
              savg = HALF*(sl(i,j,k)+sr(i,j,k))
              simh(i,j,k) = merge(simh(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  endif

end subroutine make_edge_scal_predictor_3d


subroutine make_edge_scal_transverse_3d(lo, hi, norm_dir, trans_dir, domlo, domhi, &
     s,      s_lo, s_hi, nc_s, &
     umac,   u_lo, u_hi, &
     vmac,   v_lo, v_hi, &
     wmac,   w_lo, w_hi, &
     divu, d_lo, d_hi, &
     slx, slx_lo, slx_hi, &
     srx, srx_lo, srx_hi, &
     simhx, six_lo, six_hi, &
     sly, sly_lo, sly_hi, &
     sry, sry_lo, sry_hi, &
     simhy, siy_lo, siy_hi, &
     slz, slz_lo, slz_hi, &
     srz, srz_lo, srz_hi, &
     simhz, siz_lo, siz_hi, &
     simh_trans, sit_lo, sit_hi, &
     dx, dt, is_vel, adv_bc, nbccomp, &
     comp, bccomp, is_conservative) bind(C,name="make_edge_scal_transverse_3d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)
  integer, value,   intent(in   ) :: norm_dir, trans_dir, nc_s
  integer         , intent(in   ) :: u_lo(3), u_hi(3)
  integer         , intent(in   ) :: v_lo(3), v_hi(3)
  integer         , intent(in   ) :: w_lo(3), w_hi(3)
  integer         , intent(in   ) :: d_lo(3), d_hi(3)
  integer         , intent(in   ) :: slx_lo(3), slx_hi(3)
  integer         , intent(in   ) :: srx_lo(3), srx_hi(3)
  integer         , intent(in   ) :: six_lo(3), six_hi(3)
  integer         , intent(in   ) :: sly_lo(3), sly_hi(3)
  integer         , intent(in   ) :: sry_lo(3), sry_hi(3)
  integer         , intent(in   ) :: siy_lo(3), siy_hi(3)
  integer         , intent(in   ) :: slz_lo(3), slz_hi(3)
  integer         , intent(in   ) :: srz_lo(3), srz_hi(3)
  integer         , intent(in   ) :: siz_lo(3), siz_hi(3)
  integer         , intent(in   ) :: sit_lo(3), sit_hi(3)
  double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
  double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
  double precision, intent(in   ) :: divu  (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
  double precision, intent(in   ) :: slx   (slx_lo(1):slx_hi(1),slx_lo(2):slx_hi(2),slx_lo(3):slx_hi(3))
  double precision, intent(in   ) :: srx   (srx_lo(1):srx_hi(1),srx_lo(2):srx_hi(2),srx_lo(3):srx_hi(3))
  double precision, intent(in   ) :: simhx (six_lo(1):six_hi(1),six_lo(2):six_hi(2),six_lo(3):six_hi(3))
  double precision, intent(in   ) :: sly   (sly_lo(1):sly_hi(1),sly_lo(2):sly_hi(2),sly_lo(3):sly_hi(3))
  double precision, intent(in   ) :: sry   (sry_lo(1):sry_hi(1),sry_lo(2):sry_hi(2),sry_lo(3):sry_hi(3))
  double precision, intent(in   ) :: simhy (siy_lo(1):siy_hi(1),siy_lo(2):siy_hi(2),siy_lo(3):siy_hi(3))
  double precision, intent(in   ) :: slz   (slz_lo(1):slz_hi(1),slz_lo(2):slz_hi(2),slz_lo(3):slz_hi(3))
  double precision, intent(in   ) :: srz   (srz_lo(1):srz_hi(1),srz_lo(2):srz_hi(2),srz_lo(3):srz_hi(3))
  double precision, intent(in   ) :: simhz (siz_lo(1):siz_hi(1),siz_lo(2):siz_hi(2),siz_lo(3):siz_hi(3))
  double precision, intent(inout) :: simh_trans(sit_lo(1):sit_hi(1),sit_lo(2):sit_hi(2),sit_lo(3):sit_hi(3))
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer,   value, intent(in   ) :: is_vel, nbccomp, comp, bccomp, is_conservative
  integer         , intent(in   ) :: adv_bc(3,2,nbccomp)

  ! Local variables

  double precision :: hx,hy,hz,dt2,dt3,dt4,dt6,fl,fr
  double precision :: savg

  integer :: i,j,k

  ! these correspond to s_L^{x|y}, etc.
  double precision :: slxy,srxy,slxz,srxz
  double precision :: slyx,sryx,slyz,sryz
  double precision :: slzx,srzx,slzy,srzy

  !$gpu

  dt2 = HALF*dt
  dt3 = dt/3.0d0
  dt4 = dt/4.0d0
  dt6 = dt/6.0d0

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  ! These are transverse terms.
  ! lo:hi+1 in normal direction
  ! lo:hi in transverse direction
  ! lo-1:hi+1 in unused direction

  !******************************************************************
  ! Create s_{\i-\half\e_x}^{x|y}, etc.
  !******************************************************************

  if (norm_dir == 1 .and. trans_dir == 2) then
     ! simhxy
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate xy faces
              if (is_conservative .eq. 1) then
                 ! make slxy, srxy by updating 1D extrapolation
                 slxy = slx(i,j,k) &
                      - (dt3/hy)*(simhy(i-1,j+1,k)*vmac(i-1,j+1,k) &
                      - simhy(i-1,j,k)*vmac(i-1,j,k)) &
                      - dt3*s(i-1,j,k,comp)*divu(i-1,j,k) &
                      + (dt3/hy)*s(i-1,j,k,comp)*(vmac(i-1,j+1,k)-vmac(i-1,j,k))
                 srxy = srx(i,j,k) &
                      - (dt3/hy)*(simhy(i  ,j+1,k)*vmac(i  ,j+1,k) &
                      - simhy(i  ,j,k)*vmac(i  ,j,k)) &
                      - dt3*s(i,j,k,comp)*divu(i,j,k) &
                      + (dt3/hy)*s(i,j,k,comp)*(vmac(i,j+1,k)-vmac(i,j,k))

              else

                 ! make slxy, srxy by updating 1D extrapolation
                 slxy = slx(i,j,k) &
                      - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k)) &
                      *(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                 srxy = srx(i,j,k) &
                      - (dt6/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k)) &
                      *(simhy(i  ,j+1,k)-simhy(i  ,j,k))

              end if

              ! impose lo side bc's
              if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                 if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
                    slxy = s(i-1,j,k,comp)
                    srxy = s(i-1,j,k,comp)
                 else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       srxy = min(srxy,0.d0)
                    end if
                    slxy = srxy
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
                    slxy = srxy
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
                    slxy = 0.d0
                    srxy = 0.d0
                 else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                 if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
                    slxy = s(i,j,k,comp)
                    srxy = s(i,j,k,comp)
                 else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       slxy = max(slxy,0.d0)
                    end if
                    srxy = slxy
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
                    srxy = slxy
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
                    slxy = 0.d0
                    srxy = 0.d0
                 else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
#endif
                 end if
              end if

              ! make simhxy by solving Riemann problem
              simh_trans(i,j,k) = merge(slxy,srxy,umac(i,j,k) .gt. 0.d0)
              savg = HALF*(slxy+srxy)
              simh_trans(i,j,k) = merge(simh_trans(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else if (norm_dir == 1 .and. trans_dir == 3) then
     ! loop over appropriate xz faces
     ! simhxz
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (is_conservative .eq. 1) then
                 ! make slxz, srxz by updating 1D extrapolation
                 slxz = slx(i,j,k) &
                      - (dt3/hz)*(simhz(i-1,j,k+1)*wmac(i-1,j,k+1) &
                      - simhz(i-1,j,k)*wmac(i-1,j,k)) &
                      - dt3*s(i-1,j,k,comp)*divu(i-1,j,k) &
                      + (dt3/hz)*s(i-1,j,k,comp)*(wmac(i-1,j,k+1)-wmac(i-1,j,k))
                 srxz = srx(i,j,k) &
                      - (dt3/hz)*(simhz(i  ,j,k+1)*wmac(i  ,j,k+1) &
                      - simhz(i  ,j,k)*wmac(i  ,j,k)) &
                      - dt3*s(i,j,k,comp)*divu(i,j,k) &
                      + (dt3/hz)*s(i,j,k,comp)*(wmac(i,j,k+1)-wmac(i,j,k))
              else
                 ! make slxz, srxz by updating 1D extrapolation
                 slxz = slx(i,j,k) &
                      - (dt6/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k)) &
                      *(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                 srxz = srx(i,j,k) &
                      - (dt6/hz)*(wmac(i  ,j,k+1)+wmac(i  ,j,k)) &
                      *(simhz(i  ,j,k+1)-simhz(i  ,j,k))
              end if

              ! impose lo side bc's
              if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                 if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
                    slxz = s(i-1,j,k,comp)
                    srxz = s(i-1,j,k,comp)
                 else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       srxz = min(srxz,0.d0)
                    end if
                    slxz = srxz
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
                    slxz = srxz
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
                    slxz = 0.d0
                    srxz = 0.d0
                 else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                 if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
                    slxz = s(i,j,k,comp)
                    srxz = s(i,j,k,comp)
                 else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       slxz = max(slxz,0.d0)
                    end if
                    srxz = slxz
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
                    srxz = slxz
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
                    slxz = 0.d0
                    srxz = 0.d0
                 else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
#endif
                 end if
              end if

              ! make simhxz by solving Riemann problem
              simh_trans(i,j,k) = merge(slxz,srxz,umac(i,j,k) .gt. 0.d0)
              savg = HALF*(slxz+srxz)
              simh_trans(i,j,k) = merge(simh_trans(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else if (norm_dir == 2 .and. trans_dir == 1) then
     ! simhyx
     ! loop over appropriate yx faces

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (is_conservative .eq. 1) then
                 ! make slyx, sryx by updating 1D extrapolation
                 slyx = sly(i,j,k) &
                      - (dt3/hx)*(simhx(i+1,j-1,k)*umac(i+1,j-1,k) &
                      - simhx(i,j-1,k)*umac(i,j-1,k)) &
                      - dt3*s(i,j-1,k,comp)*divu(i,j-1,k) &
                      + (dt3/hx)*s(i,j-1,k,comp)*(umac(i+1,j-1,k)-umac(i,j-1,k))
                 sryx = sry(i,j,k) &
                      - (dt3/hx)*(simhx(i+1,j  ,k)*umac(i+1,j  ,k) &
                      - simhx(i,j  ,k)*umac(i,j  ,k)) &
                      - dt3*s(i,j,k,comp)*divu(i,j,k) &
                      + (dt3/hx)*s(i,j,k,comp)*(umac(i+1,j,k)-umac(i,j,k))
              else
                 ! make slyx, sryx by updating 1D extrapolation
                 slyx = sly(i,j,k) &
                      - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k)) &
                      *(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                 sryx = sry(i,j,k) &
                      - (dt6/hx)*(umac(i+1,j  ,k)+umac(i,j  ,k)) &
                      *(simhx(i+1,j  ,k)-simhx(i,j  ,k))
              end if

              ! impose lo side bc's
              if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                 if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
                    slyx = s(i,j-1,k,comp)
                    sryx = s(i,j-1,k,comp)
                 else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       sryx = min(sryx,0.d0)
                    end if
                    slyx = sryx
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
                    slyx = sryx
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
                    slyx = 0.d0
                    sryx = 0.d0
                 else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                 if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
                    slyx = s(i,j,k,comp)
                    sryx = s(i,j,k,comp)
                 else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       slyx = max(slyx,0.d0)
                    end if
                    sryx = slyx
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
                    sryx = slyx
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
                    slyx = 0.d0
                    sryx = 0.d0
                 else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
#endif
                 end if
              end if

              ! make simhyx by solving Riemann problem
              simh_trans(i,j,k) = merge(slyx,sryx,vmac(i,j,k) .gt. 0.d0)
              savg = HALF*(slyx+sryx)
              simh_trans(i,j,k) = merge(simh_trans(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else if (norm_dir == 2 .and. trans_dir == 3) then
     ! simhyz
     ! loop over appropriate yz faces

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (is_conservative .eq. 1) then
                 ! make slyz, sryz by updating 1D extrapolation
                 slyz = sly(i,j,k) &
                      - (dt3/hz)*(simhz(i,j-1,k+1)*wmac(i,j-1,k+1) &
                      - simhz(i,j-1,k)*wmac(i,j-1,k)) &
                      - dt3*s(i,j-1,k,comp)*divu(i,j-1,k) &
                      + (dt3/hz)*s(i,j-1,k,comp)*(wmac(i,j-1,k+1)-wmac(i,j-1,k))
                 sryz = sry(i,j,k) &
                      - (dt3/hz)*(simhz(i,j  ,k+1)*wmac(i,j  ,k+1) &
                      - simhz(i,j  ,k)*wmac(i,j  ,k)) &
                      - dt3*s(i,j,k,comp)*divu(i,j,k) &
                      + (dt3/hz)*s(i,j,k,comp)*(wmac(i,j,k+1)-wmac(i,j,k))
              else
                 ! make slyz, sryz by updating 1D extrapolation
                 slyz = sly(i,j,k) &
                      - (dt6/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k)) &
                      *(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                 sryz = sry(i,j,k) &
                      - (dt6/hz)*(wmac(i,j  ,k+1)+wmac(i,j  ,k)) &
                      *(simhz(i,j  ,k+1)-simhz(i,j  ,k))
              end if

              ! impose lo side bc's
              if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                 if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
                    slyz = s(i,j-1,k,comp)
                    sryz = s(i,j-1,k,comp)
                 else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       sryz = min(sryz,0.d0)
                    end if
                    slyz = sryz
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
                    slyz = sryz
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
                    slyz = 0.d0
                    sryz = 0.d0
                 else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                 if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
                    slyz = s(i,j,k,comp)
                    sryz = s(i,j,k,comp)
                 else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       slyz = max(slyz,0.d0)
                    end if
                    sryz = slyz
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
                    sryz = slyz
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
                    slyz = 0.d0
                    sryz = 0.d0
                 else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
#endif
                 end if
              end if

              ! make simhyz by solving Riemann problem
              simh_trans(i,j,k) = merge(slyz,sryz,vmac(i,j,k) .gt. 0.d0)
              savg = HALF*(slyz+sryz)
              simh_trans(i,j,k) = merge(simh_trans(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else if (norm_dir == 3 .and. trans_dir == 1) then
     ! simhzx
     ! loop over appropriate zx faces

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (is_conservative .eq. 1) then
                 ! make slzx, srzx by updating 1D extrapolation
                 slzx = slz(i,j,k) &
                      - (dt3/hx)*(simhx(i+1,j,k-1)*umac(i+1,j,k-1) &
                      - simhx(i,j,k-1)*umac(i,j,k-1)) &
                      - dt3*s(i,j,k-1,comp)*divu(i,j,k-1) &
                      + (dt3/hx)*s(i,j,k-1,comp)*(umac(i+1,j,k-1)-umac(i,j,k-1))
                 srzx = srz(i,j,k) &
                      - (dt3/hx)*(simhx(i+1,j,k  )*umac(i+1,j,k  ) &
                      - simhx(i,j,k  )*umac(i,j,k  )) &
                      - dt3*s(i,j,k,comp)*divu(i,j,k) &
                      + (dt3/hx)*s(i,j,k,comp)*(umac(i+1,j,k)-umac(i,j,k))
              else
                 ! make slzx, srzx by updating 1D extrapolation
                 slzx = slz(i,j,k) &
                      - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1)) &
                      *(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                 srzx = srz(i,j,k) &
                      - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  )) &
                      *(simhx(i+1,j,k  )-simhx(i,j,k  ))
              end if

              ! impose lo side bc's
              if (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                 if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
                    slzx = s(i,j,k-1,comp)
                    srzx = s(i,j,k-1,comp)
                 else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       srzx = min(srzx,0.d0)
                    end if
                    slzx = srzx
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
                    slzx = srzx
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
                    slzx = 0.d0
                    srzx = 0.d0
                 else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (k .eq. hi(3) .and. hi(3)-1 .eq. domhi(3)) then
                 if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
                    slzx = s(i,j,k,comp)
                    srzx = s(i,j,k,comp)
                 else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       slzx = max(slzx,0.d0)
                    end if
                    srzx = slzx
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
                    srzx = slzx
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
                    slzx = 0.d0
                    srzx = 0.d0
                 else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
#endif
                 end if
              end if

              ! make simhzx by solving Riemann problem
              simh_trans(i,j,k) = merge(slzx,srzx,wmac(i,j,k) .gt. 0.d0)
              savg = HALF*(slzx+srzx)
              simh_trans(i,j,k) = merge(simh_trans(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  else if (norm_dir == 3 .and. trans_dir == 2) then
     ! simhzy
     ! loop over appropriate zy faces

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              if (is_conservative .eq. 1) then
                 ! make slzy, srzy by updating 1D extrapolation
                 slzy = slz(i,j,k) &
                      - (dt3/hy)*(simhy(i,j+1,k-1)*vmac(i,j+1,k-1) &
                      - simhy(i,j,k-1)*vmac(i,j,k-1)) &
                      - dt3*s(i,j,k-1,comp)*divu(i,j,k-1) &
                      + (dt3/hy)*s(i,j,k-1,comp)*(vmac(i,j+1,k-1)-vmac(i,j,k-1))
                 srzy = srz(i,j,k) &
                      - (dt3/hy)*(simhy(i,j+1,k  )*vmac(i,j+1,k  ) &
                      - simhy(i,j,k  )*vmac(i,j,k  )) &
                      - dt3*s(i,j,k,comp)*divu(i,j,k) &
                      + (dt3/hy)*s(i,j,k,comp)*(vmac(i,j+1,k)-vmac(i,j,k))
              else
                 ! make slzy, srzy by updating 1D extrapolation
                 slzy = slz(i,j,k) &
                      - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1)) &
                      *(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                 srzy = srz(i,j,k) &
                      - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  )) &
                      *(simhy(i,j+1,k  )-simhy(i,j,k  ))
              end if

              ! impose lo side bc's
              if (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                 if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
                    slzy = s(i,j,k-1,comp)
                    srzy = s(i,j,k-1,comp)
                 else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       srzy = min(srzy,0.d0)
                    end if
                    slzy = srzy
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
                    slzy = srzy
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
                    slzy = 0.d0
                    srzy = 0.d0
                 else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (k .eq. hi(3) .and. hi(3)-1 .eq. domhi(3)) then
                 if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
                    slzy = s(i,j,k,comp)
                    srzy = s(i,j,k,comp)
                 else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       slzy = max(slzy,0.d0)
                    end if
                    srzy = slzy
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
                    srzy = slzy
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
                    slzy = 0.d0
                    srzy = 0.d0
                 else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
#endif
                 end if
              end if

              ! make simhzy by solving Riemann problem
              simh_trans(i,j,k) = merge(slzy,srzy,wmac(i,j,k) .gt. 0.d0)
              savg = HALF*(slzy+srzy)
              simh_trans(i,j,k) = merge(simh_trans(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
           enddo
        enddo
     enddo

  endif

end subroutine make_edge_scal_transverse_3d



subroutine make_edge_scal_3d(lo, hi, idir, domlo, domhi, &
     s,      s_lo, s_hi, nc_s, &
     sedge, x_lo, x_hi, nc_x, &
     umac,   u_lo, u_hi, &
     vmac,   v_lo, v_hi, &
     wmac,   w_lo, w_hi, &
     Ipf, ipf_lo, ipf_hi, &
     Imf, imf_lo, imf_hi, &
     sl, sl_lo, sl_hi, &
     sr, sr_lo, sr_hi, &
     simhxy, xy_lo, xy_hi, &
     simhxz, xz_lo, xz_hi, &
     simhyx, yx_lo, yx_hi, &
     simhyz, yz_lo, yz_hi, &
     simhzx, zx_lo, zx_hi, &
     simhzy, zy_lo, zy_hi, &
     force,  f_lo, f_hi, nc_f, &
     dx, dt, is_vel, adv_bc, nbccomp, &
     comp, bccomp, is_conservative) bind(C,name="make_edge_scal_3d")

  integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)
  integer, value,   intent(in   ) :: idir, nc_s
  integer         , intent(in   ) :: x_lo(3), x_hi(3)
  integer, value,   intent(in   ) :: nc_x
  integer         , intent(in   ) :: u_lo(3), u_hi(3)
  integer         , intent(in   ) :: v_lo(3), v_hi(3)
  integer         , intent(in   ) :: w_lo(3), w_hi(3)
  integer         , intent(in   ) :: ipf_lo(3), ipf_hi(3)
  integer         , intent(in   ) :: imf_lo(3), imf_hi(3)
  integer         , intent(in   ) :: sl_lo(3), sl_hi(3)
  integer         , intent(in   ) :: sr_lo(3), sr_hi(3)
  integer         , intent(in   ) :: xy_lo(3), xy_hi(3)
  integer         , intent(in   ) :: xz_lo(3), xz_hi(3)
  integer         , intent(in   ) :: yx_lo(3), yx_hi(3)
  integer         , intent(in   ) :: yz_lo(3), yz_hi(3)
  integer         , intent(in   ) :: zx_lo(3), zx_hi(3)
  integer         , intent(in   ) :: zy_lo(3), zy_hi(3)
  integer         , intent(in   ) :: f_lo(3), f_hi(3)
  integer, value,   intent(in   ) :: nc_f
  double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
  double precision, intent(inout) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
  double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
  double precision, intent(in) :: Ipf(ipf_lo(1):ipf_hi(1),ipf_lo(2):ipf_hi(2),ipf_lo(3):ipf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in) :: Imf(imf_lo(1):imf_hi(1),imf_lo(2):imf_hi(2),imf_lo(3):imf_hi(3),AMREX_SPACEDIM)
  double precision, intent(in   ) :: sl    (sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3))
  double precision, intent(in   ) :: sr    (sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3))
  double precision, intent(in   ) :: simhxy(xy_lo(1):xy_hi(1),xy_lo(2):xy_hi(2),xy_lo(3):xy_hi(3))
  double precision, intent(in   ) :: simhxz(xz_lo(1):xz_hi(1),xz_lo(2):xz_hi(2),xz_lo(3):xz_hi(3))
  double precision, intent(in   ) :: simhyx(yx_lo(1):yx_hi(1),yx_lo(2):yx_hi(2),yx_lo(3):yx_hi(3))
  double precision, intent(in   ) :: simhyz(yz_lo(1):yz_hi(1),yz_lo(2):yz_hi(2),yz_lo(3):yz_hi(3))
  double precision, intent(in   ) :: simhzx(zx_lo(1):zx_hi(1),zx_lo(2):zx_hi(2),zx_lo(3):zx_hi(3))
  double precision, intent(in   ) :: simhzy(zy_lo(1):zy_hi(1),zy_lo(2):zy_hi(2),zy_lo(3):zy_hi(3))
  double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
  double precision, intent(in   ) :: dx(3)
  double precision, value, intent(in   ) :: dt
  integer,   value, intent(in   ) :: is_vel, nbccomp, comp, bccomp, is_conservative
  integer         , intent(in   ) :: adv_bc(3,2,nbccomp)

  ! Local variables

  double precision :: hx,hy,hz,dt2,dt3,dt4,dt6,fl,fr
  double precision :: savg

  integer :: i,j,k

  ! these correspond to \mathrm{sedge}_L^x, etc.
  double precision :: sedgelx,sedgerx
  double precision :: sedgely,sedgery
  double precision :: sedgelz,sedgerz

  !$gpu

  dt2 = HALF*dt
  dt3 = dt/3.0d0
  dt4 = dt/4.0d0
  dt6 = dt/6.0d0

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  !******************************************************************
  ! Create sedgelx, etc.
  !******************************************************************

  if (idir == 1) then
     ! Final edge states.
     ! lo:hi+1 in the normal direction
     ! lo:hi in the transverse directions
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate x-faces
              if (is_conservative .eq. 1) then
                 ! make sedgelx, sedgerx
                 fl = merge(force(i-1,j,k,comp), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
                 fr = merge(force(i  ,j,k,comp), Imf(i  ,j,k,1), ppm_trace_forces == 0)

                 sedgelx = sl(i,j,k) &
                      - (dt2/hy)*(simhyz(i-1,j+1,k  )*vmac(i-1,j+1,k  ) &
                      - simhyz(i-1,j,k)*vmac(i-1,j,k)) &
                      - (dt2/hz)*(simhzy(i-1,j  ,k+1)*wmac(i-1,j  ,k+1) &
                      - simhzy(i-1,j,k)*wmac(i-1,j,k)) &
                      - (dt2/hx)*s(i-1,j,k,comp)*(umac(i  ,j,k)-umac(i-1,j,k)) &
                      + dt2*fl

                 sedgerx = sr(i,j,k) &
                      - (dt2/hy)*(simhyz(i  ,j+1,k  )*vmac(i  ,j+1,  k) &
                      - simhyz(i  ,j,k)*vmac(i  ,j,k)) &
                      - (dt2/hz)*(simhzy(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) &
                      - simhzy(i  ,j,k)*wmac(i  ,j,k)) &
                      - (dt2/hx)*s(i  ,j,k,comp)*(umac(i+1,j,k)-umac(i  ,j,k)) &
                      + dt2*fr
              else
                 ! make sedgelx, sedgerx
                 fl = merge(force(i-1,j,k,comp), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
                 fr = merge(force(i  ,j,k,comp), Ipf(i  ,j,k,1), ppm_trace_forces == 0)

                 sedgelx = sl(i,j,k) &
                      - (dt4/hy)*(vmac(i-1,j+1,k  )+vmac(i-1,j,k))* &
                      (simhyz(i-1,j+1,k  )-simhyz(i-1,j,k)) &
                      - (dt4/hz)*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k))* &
                      (simhzy(i-1,j  ,k+1)-simhzy(i-1,j,k)) &
                      + dt2*fl

                 sedgerx = sr(i,j,k) &
                      - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i  ,j,k))* &
                      (simhyz(i  ,j+1,k  )-simhyz(i  ,j,k)) &
                      - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k))* &
                      (simhzy(i  ,j  ,k+1)-simhzy(i  ,j,k)) &
                      + dt2*fr
              end if

              ! make sedgex by solving Riemann problem
              ! boundary conditions enforced outside of i,j,k loop
              sedge(i,j,k,comp) = merge(sedgelx,sedgerx,umac(i,j,k) .gt. 0.d0)
              savg = HALF*(sedgelx+sedgerx)
              sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(umac(i,j,k)).gt.rel_eps)

              ! impose lo side bc's
              if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                 if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
                    sedge(i,j,k,comp) = s(i-1,j,k,comp)
                 else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       sedge(i,j,k,comp) = min(sedgerx,0.d0)
                    else
                       sedge(i,j,k,comp) = sedgerx
                    end if
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
                    sedge(i,j,k,comp) = sedgerx
                 else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
                    sedge(i,j,k,comp) = 0.d0
                 else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                 if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
                    sedge(i,j,k,comp) = s(i,j,k,comp)
                 else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 1) then
                       sedge(i,j,k,comp) = max(sedgelx,0.d0)
                    else
                       sedge(i,j,k,comp) = sedgelx
                    end if
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
                    sedge(i,j,k,comp) = sedgelx
                 else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
                    sedge(i,j,k,comp) = 0.d0
                 else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
#endif
                 end if
              end if
           enddo
        enddo
     enddo

  else if (idir == 2) then

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate y-faces
              if (is_conservative .eq. 1) then
                 ! make sedgely, sedgery
                 fl = merge(force(i,j-1,k,comp), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
                 fr = merge(force(i,j  ,k,comp), Imf(i,j  ,k,2), ppm_trace_forces == 0)

                 sedgely = sl(i,j,k) &
                      - (dt2/hx)*(simhxz(i+1,j-1,k  )*umac(i+1,j-1,k  ) &
                      - simhxz(i,j-1,k)*umac(i,j-1,k)) &
                      - (dt2/hz)*(simhzx(i  ,j-1,k+1)*wmac(i  ,j-1,k+1) &
                      - simhzx(i,j-1,k)*wmac(i,j-1,k)) &
                      - (dt2/hy)*s(i,j-1,k,comp)*(vmac(i,j  ,k)-vmac(i,j-1,k)) &
                      + dt2*fl

                 sedgery = sr(i,j,k) &
                      - (dt2/hx)*(simhxz(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
                      - simhxz(i,j  ,k)*umac(i,j  ,k)) &
                      - (dt2/hz)*(simhzx(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) &
                      - simhzx(i,j  ,k)*wmac(i,j  ,k)) &
                      - (dt2/hy)*s(i,j  ,k,comp)*(vmac(i,j+1,k)-vmac(i,j  ,k)) &
                      + dt2*fr
              else
                 ! make sedgely, sedgery
                 fl = merge(force(i,j-1,k,comp), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
                 fr = merge(force(i,j  ,k,comp), Imf(i,j  ,k,2), ppm_trace_forces == 0)

                 sedgely = sl(i,j,k) &
                      - (dt4/hx)*(umac(i+1,j-1,k  )+umac(i,j-1,k))* &
                      (simhxz(i+1,j-1,k  )-simhxz(i,j-1,k)) &
                      - (dt4/hz)*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k))* &
                      (simhzx(i  ,j-1,k+1)-simhzx(i,j-1,k)) &
                      + dt2*fl

                 sedgery = sr(i,j,k) &
                      - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j  ,k))* &
                      (simhxz(i+1,j  ,k  )-simhxz(i,j  ,k)) &
                      - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k))* &
                      (simhzx(i  ,j  ,k+1)-simhzx(i,j  ,k)) &
                      + dt2*fr
              end if

              ! make sedgey by solving Riemann problem
              ! boundary conditions enforced outside of i,j,k loop
              sedge(i,j,k,comp) = merge(sedgely,sedgery,vmac(i,j,k) .gt. 0.d0)
              savg = HALF*(sedgely+sedgery)
              sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(vmac(i,j,k)).gt.rel_eps)

              ! impose lo side bc's
              if (j .eq. lo(2) .and. lo(2) .eq. domlo(2)) then
                 if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
                    sedge(i,j,k,comp) = s(i,j-1,k,comp)
                 else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       sedge(i,j,k,comp) = min(sedgery,0.d0)
                    else
                       sedge(i,j,k,comp) = sedgery
                    end if
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
                    sedge(i,j,k,comp) = sedgery
                 else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
                    sedge(i,j,k,comp) = 0.d0
                 else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                 if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
                    sedge(i,j,k,comp) = s(i,j,k,comp)
                 else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 2) then
                       sedge(i,j,k,comp) = max(sedgely,0.d0)
                    else
                       sedge(i,j,k,comp) = sedgely
                    end if
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
                    sedge(i,j,k,comp) = sedgely
                 else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
                    sedge(i,j,k,comp) = 0.d0
                 else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
#endif
                 end if
              end if
           enddo
        enddo
     enddo

  else ! idir == 3

     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

              ! loop over appropriate z-faces
              if (is_conservative .eq. 1) then
                 ! make sedgelz, sedgerz
                 fl = merge(force(i,j,k-1,comp), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
                 fr = merge(force(i,j,k  ,comp), Imf(i,j,k  ,3), ppm_trace_forces == 0)

                 sedgelz = sl(i,j,k) &
                      - (dt2/hx)*(simhxy(i+1,j  ,k-1)*umac(i+1,j  ,k-1) &
                      - simhxy(i,j,k-1)*umac(i,j,k-1)) &
                      - (dt2/hy)*(simhyx(i  ,j+1,k-1)*vmac(i  ,j+1,k-1) &
                      - simhyx(i,j,k-1)*vmac(i,j,k-1)) &
                      - (dt2/hz)*s(i,j,k-1,comp)*(wmac(i,j,k  )-wmac(i,j,k-1)) &
                      + dt2*fl

                 sedgerz = sr(i,j,k) &
                      - (dt2/hx)*(simhxy(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
                      - simhxy(i,j,k  )*umac(i,j,k  )) &
                      - (dt2/hy)*(simhyx(i  ,j+1,k  )*vmac(i  ,j+1,k  ) &
                      - simhyx(i,j,k  )*vmac(i,j,k  )) &
                      - (dt2/hz)*s(i,j,k  ,comp)*(wmac(i,j,k+1)-wmac(i,j,k  )) &
                      + dt2*fr
              else
                 ! make sedgelz, sedgerz
                 fl = merge(force(i,j,k-1,comp), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
                 fr = merge(force(i,j,k  ,comp), Imf(i,j,k  ,3), ppm_trace_forces == 0)

                 sedgelz = sl(i,j,k) &
                      - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j,k-1)) &
                      *(simhxy(i+1,j  ,k-1)-simhxy(i,j,k-1)) &
                      - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1)) &
                      *(simhyx(i  ,j+1,k-1)-simhyx(i,j,k-1)) &
                      + dt2*fl

                 sedgerz = sr(i,j,k) &
                      - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j,k  )) &
                      *(simhxy(i+1,j  ,k  )-simhxy(i,j,k  )) &
                      - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i,j,k  )) &
                      *(simhyx(i  ,j+1,k  )-simhyx(i,j,k  )) &
                      + dt2*fr
              end if

              ! make sedgez by solving Riemann problem
              ! boundary conditions enforced outside of i,j,k loop
              sedge(i,j,k,comp) = merge(sedgelz,sedgerz,wmac(i,j,k) .gt. 0.d0)
              savg = HALF*(sedgelz+sedgerz)
              sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(wmac(i,j,k)).gt.rel_eps)\

              ! impose lo side bc's
              if (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                 if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
                    sedge(i,j,k,comp) = s(i,j,k-1,comp)
                 else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       sedge(i,j,k,comp) = min(sedgerz,0.d0)
                    else
                       sedge(i,j,k,comp) = sedgerz
                    end if
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
                    sedge(i,j,k,comp) = sedgerz
                 else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
                    sedge(i,j,k,comp) = 0.d0
                 else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
#endif
                 end if
              end if

              ! impose hi side bc's
              if (k .eq. hi(3) .and. hi(3)-1 .eq. domhi(3)) then
                 if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
                    sedge(i,j,k,comp) = s(i,j,k,comp)
                 else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
                      adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
                    if (is_vel .eq. 1 .and. comp .eq. 3) then
                       sedge(i,j,k,comp) = max(sedgelz,0.d0)
                    else
                       sedge(i,j,k,comp) = sedgelz
                    end if
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
                    sedge(i,j,k,comp) = sedgelz
                 else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
                    sedge(i,j,k,comp) = 0.d0
                 else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
                 else
#ifndef AMREX_USE_GPU
                    call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
#endif
                 end if
              end if
           enddo
        enddo
     enddo

  end if

end subroutine make_edge_scal_3d

#endif

end module make_edge_scal_module
