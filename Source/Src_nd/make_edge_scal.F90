
#include "AMReX_BC_TYPES.H"

module make_edge_scal_module
  ! make_edge_scal constructs the edge state of a scalar, using a
  ! second-order Taylor expansion in space (through dx/2) and time
  ! (though dt/2) (if ppm_type = 0) or using PPM (for ppm_type = 1,2).
  !
  ! We use only MAC-projected edge velocities in this prediction.
  !
  ! We are computing all edge states for each variable.  This is what is
  ! done for the final updates of the state variables and velocity.  For
  ! velocity, we should set is_vel = .true.

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use slope_module
  use ppm_module
  use meth_params_module, only: rel_eps, ppm_type, ppm_trace_forces

  implicit none

  private

contains

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
           if (i .eq. domlo(1)) then
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

            ! impose hi side bc's
            else if (i .eq. domhi(1)+1) then
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
           if (j .eq. domlo(2)) then
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

            ! impose hi side bc's
            else if (j .eq. domhi(2)+1) then
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
           if (i .eq. domlo(1)) then
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

           ! impose hi side bc's
            else if (i .eq. domhi(1)+1) then
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
           if (j .eq. domlo(2)) then
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

           ! impose hi side bc's
            else if (j .eq. domhi(2)+1) then
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
     comp, bccomp, is_conservative, &
     sll, sll_lo, sll_hi, &
     srr, srr_lo, srr_hi) bind(C,name="make_edge_scal_3d")

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
  integer         , intent(in   ) :: sll_lo(3), sll_hi(3)
  integer         , intent(in   ) :: srr_lo(3), srr_hi(3)
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
  double precision, intent(inout) :: sll    (sll_lo(1):sll_hi(1),sll_lo(2):sll_hi(2),sll_lo(3):sll_hi(3))
  double precision, intent(inout) :: srr    (srr_lo(1):srr_hi(1),srr_lo(2):srr_hi(2),srr_lo(3):srr_hi(3))
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

            sedgelx = sll(i,j,k)
            sedgerx = srr(i,j,k)

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
                 fr = merge(force(i  ,j,k,comp), Imf(i  ,j,k,1), ppm_trace_forces == 0)

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

              if (i .eq. domlo(1)+1 .and. j .eq. domlo(2)+1 .and. k .eq. domlo(3)+1) then 
                write(*,*) "fortran fl, fr  =", fl, fr
              endif

              sll(i,j,k) = sedgelx
              srr(i,j,k) = sedgerx

              ! make sedgex by solving Riemann problem
              ! boundary conditions enforced outside of i,j,k loop
              sedge(i,j,k,comp) = merge(sedgelx,sedgerx,umac(i,j,k) .gt. 0.d0)
              savg = HALF*(sedgelx+sedgerx)
              sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(umac(i,j,k)).gt.rel_eps)

              ! impose lo side bc's
              if (i .eq. domlo(1)) then
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

              ! impose hi side bc's
              else if (i .eq. domhi(1)+1) then
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
              if (j .eq. domlo(2)) then
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

              ! impose hi side bc's
              else if (j .eq. domhi(2)+1) then
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

            sedgelz = sll(i,j,k)
            sedgerz = srr(i,j,k)

              ! loop over appropriate z-faces
            !   if (is_conservative .eq. 1) then
            !      ! make sedgelz, sedgerz
            !      fl = merge(force(i,j,k-1,comp), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
            !      fr = merge(force(i,j,k  ,comp), Imf(i,j,k  ,3), ppm_trace_forces == 0)

            !      sedgelz = sl(i,j,k) &
            !           - (dt2/hx)*(simhxy(i+1,j  ,k-1)*umac(i+1,j  ,k-1) &
            !           - simhxy(i,j,k-1)*umac(i,j,k-1)) &
            !           - (dt2/hy)*(simhyx(i  ,j+1,k-1)*vmac(i  ,j+1,k-1) &
            !           - simhyx(i,j,k-1)*vmac(i,j,k-1)) &
            !           - (dt2/hz)*s(i,j,k-1,comp)*(wmac(i,j,k  )-wmac(i,j,k-1)) &
            !           + dt2*fl

            !      sedgerz = sr(i,j,k) &
            !           - (dt2/hx)*(simhxy(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
            !           - simhxy(i,j,k  )*umac(i,j,k  )) &
            !           - (dt2/hy)*(simhyx(i  ,j+1,k  )*vmac(i  ,j+1,k  ) &
            !           - simhyx(i,j,k  )*vmac(i,j,k  )) &
            !           - (dt2/hz)*s(i,j,k  ,comp)*(wmac(i,j,k+1)-wmac(i,j,k  )) &
            !           + dt2*fr
            !   else
            !      ! make sedgelz, sedgerz
            !      fl = merge(force(i,j,k-1,comp), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
            !      fr = merge(force(i,j,k  ,comp), Imf(i,j,k  ,3), ppm_trace_forces == 0)

            !      sedgelz = sl(i,j,k) &
            !           - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j,k-1)) &
            !           *(simhxy(i+1,j  ,k-1)-simhxy(i,j,k-1)) &
            !           - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1)) &
            !           *(simhyx(i  ,j+1,k-1)-simhyx(i,j,k-1)) &
            !           + dt2*fl

            !      sedgerz = sr(i,j,k) &
            !           - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j,k  )) &
            !           *(simhxy(i+1,j  ,k  )-simhxy(i,j,k  )) &
            !           - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i,j,k  )) &
            !           *(simhyx(i  ,j+1,k  )-simhyx(i,j,k  )) &
            !           + dt2*fr
            !   end if

              ! make sedgez by solving Riemann problem
              ! boundary conditions enforced outside of i,j,k loop
            !   sedge(i,j,k,comp) = merge(sedgelz,sedgerz,wmac(i,j,k) .gt. 0.d0)
            !   savg = HALF*(sedgelz+sedgerz)
            !   sedge(i,j,k,comp) = merge(sedge(i,j,k,comp),savg,abs(wmac(i,j,k)).gt.rel_eps)

!               ! impose lo side bc's
              if (k .eq. domlo(3)) then
!                  if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
!                     sedge(i,j,k,comp) = s(i,j,k-1,comp)
!                  else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
!                       adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
!                     if (is_vel .eq. 1 .and. comp .eq. 3) then
!                        sedge(i,j,k,comp) = min(sedgerz,0.d0)
!                     else
!                        sedge(i,j,k,comp) = sedgerz
!                     end if
!                  else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
!                     sedge(i,j,k,comp) = sedgerz
!                  else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
!                     sedge(i,j,k,comp) = 0.d0
!                  else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
!                  else
! #ifndef AMREX_USE_GPU
!                     call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
! #endif
!                  end if

              ! impose hi side bc's
              else if (k .eq. domhi(3)+1) then
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
