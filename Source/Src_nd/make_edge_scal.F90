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

#if (AMREX_SPACEDIM == 2)
  subroutine make_edge_scal_2d(domlo, domhi, lo, hi, &
       s,      s_lo, s_hi, nc_s, ng_s, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, ng_um, &
       force,  f_lo, f_hi, nc_f, &
       dx, dt, is_vel, adv_bc, nbccomp, &
       comp, bccomp, is_conservative) bind(C,name="make_edge_scal_2d")

    integer         , intent(in   ) :: domlo(2), domhi(2), lo(2), hi(2)
    integer         , intent(in   ) :: s_lo(2), s_hi(2), nc_s
    integer, value,   intent(in   ) :: ng_s
    integer         , intent(in   ) :: x_lo(2), x_hi(2), nc_x
    integer         , intent(in   ) :: y_lo(2), y_hi(2), nc_y
    integer         , intent(in   ) :: u_lo(2), u_hi(2)
    integer         , intent(in   ) :: v_lo(2), v_hi(2)
    integer, value,   intent(in   ) :: ng_um
    integer         , intent(in   ) :: f_lo(2), f_hi(2), nc_f
    double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),nc_s)
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),nc_x)
    double precision, intent(inout) :: sedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),nc_y)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2))
    double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2))
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),nc_f)
    double precision, intent(in   ) :: dx(2), dt
    integer         , intent(in   ) :: is_vel, nbccomp, comp, bccomp, is_conservative
    integer         , intent(in   ) :: adv_bc(2,2,nbccomp)

    ! Local variables
    double precision, pointer :: slopex(:,:,:)
    double precision, pointer :: slopey(:,:,:)

    double precision :: hx,hy,dt2,dt4,savg,fl,fr

    integer :: i,j,is,js,ie,je

    double precision, pointer :: Ip(:,:,:), Ipf(:,:,:)
    double precision, pointer :: Im(:,:,:), Imf(:,:,:)

    ! these correspond to s_L^x, etc.
    double precision, pointer:: slx(:,:),srx(:,:)
    double precision, pointer:: sly(:,:),sry(:,:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    double precision, pointer:: simhx(:,:),simhy(:,:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    double precision, pointer:: sedgelx(:,:),sedgerx(:,:)
    double precision, pointer:: sedgely(:,:),sedgery(:,:)

    call bl_allocate(Ip,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Im,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    call bl_allocate(Ipf,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Imf,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    call bl_allocate(slopex,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(slopey,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse direction
    call bl_allocate(slx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1)
    call bl_allocate(srx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1)
    call bl_allocate(simhx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1)

    call bl_allocate(sly,lo(1)-1,hi(1)+1,lo(2),hi(2)+1)
    call bl_allocate(sry,lo(1)-1,hi(1)+1,lo(2),hi(2)+1)
    call bl_allocate(simhy,lo(1)-1,hi(1)+1,lo(2),hi(2)+1)

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    call bl_allocate(sedgelx,lo(1),hi(1)+1,lo(2),hi(2))
    call bl_allocate(sedgerx,lo(1),hi(1)+1,lo(2),hi(2))
    call bl_allocate(sedgely,lo(1),hi(1),lo(2),hi(2)+1)
    call bl_allocate(sedgery,lo(1),hi(1),lo(2),hi(2)+1)

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    if (ppm_type .eq. 0) then
       call slopex_2d(s(:,:,comp:comp),slopex,domlo,domhi,lo,hi,ng_s,1,adv_bc(:,:,bccomp:bccomp))
       call slopey_2d(s(:,:,comp:comp),slopey,domlo,domhi,lo,hi,ng_s,1,adv_bc(:,:,bccomp:bccomp))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(s(:,:,comp),ng_s,umac,vmac,ng_um,Ip,Im, &
                    domlo,domhi,lo,hi,adv_bc(:,:,bccomp),dx,dt,.true.)
       if (ppm_trace_forces .eq. 1) then
          call ppm_2d(force(:,:,comp),ng_s,umac,vmac,ng_um,Ipf,Imf, &
                       domlo,domhi,lo,hi,adv_bc(:,:,bccomp),dx,dt,.true.)
       end if
    end if
    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    !******************************************************************
    ! Create s_{\i-\half\e_x}^x, etc.
    !******************************************************************

    ! loop over appropriate x-faces
    if (ppm_type .eq. 0) then
       do j=js-1,je+1
          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,j) = s(i-1,j,comp) + (HALF - dt2*umac(i,j)/hx)*slopex(i-1,j,1)
             srx(i,j) = s(i  ,j,comp) - (HALF + dt2*umac(i,j)/hx)*slopex(i  ,j,1)
          enddo
       enddo
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js-1,je+1
          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,j) = Ip(i-1,j,1)
             srx(i,j) = Im(i  ,j,1)
          end do
       end do
    end if

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          slx(is,js-1:je+1) = s(is-1,js-1:je+1,comp)
          srx(is,js-1:je+1) = s(is-1,js-1:je+1,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             srx(is,js-1:je+1) = min(srx(is,js-1:je+1),0.d0)
          end if
          slx(is,js-1:je+1) = srx(is,js-1:je+1)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          slx(is,js-1:je+1) = srx(is,js-1:je+1)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          slx(is,js-1:je+1) = 0.d0
          srx(is,js-1:je+1) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          slx(ie+1,js-1:je+1) = s(ie+1,js-1:je+1,comp)
          srx(ie+1,js-1:je+1) = s(ie+1,js-1:je+1,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             slx(ie+1,js-1:je+1) = max(slx(ie+1,js-1:je+1),0.d0)
          end if
          srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          slx(ie+1,js-1:je+1) = 0.d0
          srx(ie+1,js-1:je+1) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,2)")
       end if
    end if

    do j=js-1,je+1
       do i=is,ie+1
          ! make simhx by solving Riemann problem
          simhx(i,j) = merge(slx(i,j),srx(i,j),umac(i,j) .gt. 0.d0)
          savg = HALF*(slx(i,j)+srx(i,j))
          simhx(i,j) = merge(simhx(i,j),savg,abs(umac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! loop over appropriate y-faces
    if (ppm_type .eq. 0) then
       do j=js,je+1
          do i=is-1,ie+1
             ! make sly, sry with 1D extrapolation
             sly(i,j) = s(i,j-1,comp) + (HALF - dt2*vmac(i,j)/hy)*slopey(i,j-1,1)
             sry(i,j) = s(i,j  ,comp) - (HALF + dt2*vmac(i,j)/hy)*slopey(i,j  ,1)
          enddo
       enddo
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je+1
          do i=is-1,ie+1
             ! make sly, sry with 1D extrapolation
             sly(i,j) = Ip(i,j-1,2)
             sry(i,j) = Im(i,j  ,2)
          enddo
       enddo
    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
          sly(is-1:ie+1,js) = s(is-1:ie+1,js-1,comp)
          sry(is-1:ie+1,js) = s(is-1:ie+1,js-1,comp)
       else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sry(is-1:ie+1,js) = min(sry(is-1:ie+1,js),0.d0)
          end if
          sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
          sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
          sly(is-1:ie+1,js) = 0.d0
          sry(is-1:ie+1,js) = 0.d0
       else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
          sly(is-1:ie+1,je+1) = s(is-1:ie+1,je+1,comp)
          sry(is-1:ie+1,je+1) = s(is-1:ie+1,je+1,comp)
       else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sly(is-1:ie+1,je+1) = max(sly(is-1:ie+1,je+1),0.d0)
          end if
          sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
          sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
          sly(is-1:ie+1,je+1) = 0.d0
          sry(is-1:ie+1,je+1) = 0.d0
       else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,2)")
       end if
    end if

    do j=js,je+1
       do i=is-1,ie+1
          ! make simhy by solving Riemann problem
          simhy(i,j) = merge(sly(i,j),sry(i,j),vmac(i,j) .gt. 0.d0)
          savg = HALF*(sly(i,j)+sry(i,j))
          simhy(i,j) = merge(simhy(i,j),savg,abs(vmac(i,j)) .gt. rel_eps)
       enddo
    enddo

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! loop over appropriate x-faces
    do j=js,je
       do i=is,ie+1
          ! make sedgelx, sedgerx
          fl = merge(force(i-1,j,comp), Ipf(i-1,j,1), ppm_trace_forces == 0)
          fr = merge(force(i  ,j,comp), Imf(i  ,j,1), ppm_trace_forces == 0)

          if(is_conservative .eq. 1) then
             sedgelx(i,j) = slx(i,j) &
                  - (dt2/hy)*(simhy(i-1,j+1)*vmac(i-1,j+1) - simhy(i-1,j)*vmac(i-1,j)) &
                  - (dt2/hx)*s(i-1,j,comp)*(umac(i  ,j)-umac(i-1,j)) &
                  + dt2*fl
             sedgerx(i,j) = srx(i,j) &
                  - (dt2/hy)*(simhy(i  ,j+1)*vmac(i  ,j+1) - simhy(i  ,j)*vmac(i  ,j)) &
                  - (dt2/hx)*s(i  ,j,comp)*(umac(i+1,j)-umac(i  ,j)) &
                  + dt2*fr
          else
             sedgelx(i,j) = slx(i,j) &
                  - (dt4/hy)*(vmac(i-1,j+1)+vmac(i-1,j))*(simhy(i-1,j+1)-simhy(i-1,j)) &
                  + dt2*fl
             sedgerx(i,j) = srx(i,j) &
                  - (dt4/hy)*(vmac(i  ,j+1)+vmac(i  ,j))*(simhy(i  ,j+1)-simhy(i  ,j)) &
                  + dt2*fr
          end if

          ! make sedgex by solving Riemann problem
          ! boundary conditions enforced outside of i,j loop
          sedgex(i,j,comp) = merge(sedgelx(i,j),sedgerx(i,j),umac(i,j) .gt. 0.d0)
          savg = HALF*(sedgelx(i,j)+sedgerx(i,j))
          sedgex(i,j,comp) = merge(sedgex(i,j,comp),savg,abs(umac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          sedgex(is,js:je,comp) = s(is-1,js:je,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             sedgex(is,js:je,comp) = min(sedgerx(is,js:je),0.d0)
          else
             sedgex(is,js:je,comp) = sedgerx(is,js:je)
          end if
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          sedgex(is,js:je,comp) = sedgerx(is,js:je)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          sedgex(is,js:je,comp) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          sedgex(ie+1,js:je,comp) = s(ie+1,js:je,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             sedgex(ie+1,js:je,comp) = max(sedgelx(ie+1,js:je),0.d0)
          else
             sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
          end if
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          sedgex(ie+1,js:je,comp) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(1,2)")
       end if
    end if

    ! loop over appropriate y-faces
    do j=js,je+1
       do i=is,ie
          ! make sedgely, sedgery
          fl = merge(force(i,j-1,comp), Ipf(i,j-1,2), ppm_trace_forces == 0)
          fr = merge(force(i,j  ,comp), Imf(i,j  ,2), ppm_trace_forces == 0)

          if(is_conservative .eq. 1) then
             sedgely(i,j) = sly(i,j) &
                  - (dt2/hx)*(simhx(i+1,j-1)*umac(i+1,j-1) - simhx(i,j-1)*umac(i,j-1)) &
                  - (dt2/hy)*s(i,j-1,comp)*(vmac(i,j  )-vmac(i,j-1)) &
                  + dt2*fl
             sedgery(i,j) = sry(i,j) &
                  - (dt2/hx)*(simhx(i+1,j  )*umac(i+1,j  ) - simhx(i,j  )*umac(i,j  )) &
                  - (dt2/hy)*s(i,j  ,comp)*(vmac(i,j+1)-vmac(i,j  )) &
                  + dt2*fr
          else
             sedgely(i,j) = sly(i,j) &
                  - (dt4/hx)*(umac(i+1,j-1)+umac(i,j-1))*(simhx(i+1,j-1)-simhx(i,j-1)) &
                  + dt2*fl
             sedgery(i,j) = sry(i,j) &
                  - (dt4/hx)*(umac(i+1,j  )+umac(i,j  ))*(simhx(i+1,j  )-simhx(i,j  )) &
                  + dt2*fr
          end if

          ! make sedgey by solving Riemann problem
          ! boundary conditions enforced outside of i,j loop
          sedgey(i,j,comp) = merge(sedgely(i,j),sedgery(i,j),vmac(i,j) .gt. 0.d0)
          savg = HALF*(sedgely(i,j)+sedgery(i,j))
          sedgey(i,j,comp) = merge(sedgey(i,j,comp),savg,abs(vmac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
          sedgey(is:ie,js,comp) = s(is:ie,js-1,comp)
       else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sedgey(is:ie,js,comp) = min(sedgery(is:ie,js),0.d0)
          else
             sedgey(is:ie,js,comp) = sedgery(is:ie,js)
          end if
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
          sedgey(is:ie,js,comp) = sedgery(is:ie,js)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
          sedgey(is:ie,js,comp) = 0.d0
       else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
          sedgey(is:ie,je+1,comp) = s(is:ie,je+1,comp)
       else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sedgey(is:ie,je+1,comp) = max(sedgely(is:ie,je+1),0.d0)
          else
             sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
          end if
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
          sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
          sedgey(is:ie,je+1,comp) = 0.d0
       else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_2d: invalid boundary type adv_bc(2,2)")
       end if
    end if

    call bl_deallocate(Ip)
    call bl_deallocate(Im)
    call bl_deallocate(Ipf)
    call bl_deallocate(Imf)
    call bl_deallocate(slopex)
    call bl_deallocate(slopey)
    call bl_deallocate(slx)
    call bl_deallocate(srx)
    call bl_deallocate(simhx)
    call bl_deallocate(sly)
    call bl_deallocate(sry)
    call bl_deallocate(simhy)
    call bl_deallocate(sedgelx)
    call bl_deallocate(sedgerx)
    call bl_deallocate(sedgely)
    call bl_deallocate(sedgery)

  end subroutine make_edge_scal_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine make_edge_scal_3d(domlo, domhi, lo, hi, &
       s,      s_lo, s_hi, nc_s, ng_s, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       sedgez, z_lo, z_hi, nc_z, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, &
       wmac,   w_lo, w_hi, ng_um, &
       force,  f_lo, f_hi, nc_f, &
       dx, dt, is_vel, adv_bc, nbccomp, &
       comp, bccomp, is_conservative) bind(C,name="make_edge_scal_3d")

    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer, value,   intent(in   ) :: ng_s
    integer         , intent(in   ) :: x_lo(3), x_hi(3), nc_x
    integer         , intent(in   ) :: y_lo(3), y_hi(3), nc_y
    integer         , intent(in   ) :: z_lo(3), z_hi(3), nc_z
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    integer, value,   intent(in   ) :: ng_um
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    double precision, intent(in   ) :: s     (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    double precision, intent(inout) :: sedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nc_y)
    double precision, intent(inout) :: sedgez(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),nc_z)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: is_vel, nbccomp, comp, bccomp, is_conservative
    integer         , intent(in   ) :: adv_bc(3,2,nbccomp)

    ! Local variables
    double precision, pointer :: slopex(:,:,:,:)
    double precision, pointer :: slopey(:,:,:,:)
    double precision, pointer :: slopez(:,:,:,:)

    double precision :: hx,hy,hz,dt2,dt3,dt4,dt6,fl,fr
    double precision :: savg

    integer :: i,j,k,is,js,ks,ie,je,ke

    double precision, pointer :: Ip(:,:,:,:), Ipf(:,:,:,:)
    double precision, pointer :: Im(:,:,:,:), Imf(:,:,:,:)

    ! these correspond to s_L^x, etc.
    double precision, pointer:: slx(:,:,:),srx(:,:,:)
    double precision, pointer:: sly(:,:,:),sry(:,:,:)
    double precision, pointer:: slz(:,:,:),srz(:,:,:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    double precision, pointer:: simhx(:,:,:),simhy(:,:,:),simhz(:,:,:)

    ! these correspond to s_L^{x|y}, etc.
    double precision, pointer:: slxy(:,:,:),srxy(:,:,:),slxz(:,:,:),srxz(:,:,:)
    double precision, pointer:: slyx(:,:,:),sryx(:,:,:),slyz(:,:,:),sryz(:,:,:)
    double precision, pointer:: slzx(:,:,:),srzx(:,:,:),slzy(:,:,:),srzy(:,:,:)

    ! these correspond to s_{\i-\half\e_x}^{x|y}, etc.
    double precision, pointer:: simhxy(:,:,:),simhxz(:,:,:)
    double precision, pointer:: simhyx(:,:,:),simhyz(:,:,:)
    double precision, pointer:: simhzx(:,:,:),simhzy(:,:,:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    double precision, pointer:: sedgelx(:,:,:),sedgerx(:,:,:)
    double precision, pointer:: sedgely(:,:,:),sedgery(:,:,:)
    double precision, pointer:: sedgelz(:,:,:),sedgerz(:,:,:)

    ! used in corner coupling for conservative quantities
    double precision, pointer:: divu(:,:,:)

    call bl_allocate(slopex,lo-1,hi+1,1)
    call bl_allocate(slopey,lo-1,hi+1,1)
    call bl_allocate(slopez,lo-1,hi+1,1)

    call bl_allocate(Ip,lo-1,hi+1,3)
    call bl_allocate(Im,lo-1,hi+1,3)

    call bl_allocate(Ipf,lo-1,hi+1,3)
    call bl_allocate(Imf,lo-1,hi+1,3)

    if (is_conservative .eq. 1) then
       call bl_allocate(divu,lo-1,hi+1)
    end if

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)

    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                divu(i,j,k) = (  umac(i+1,j,k)-umac(i,j,k) &
                     + vmac(i,j+1,k)-vmac(i,j,k) &
                     + wmac(i,j,k+1)-wmac(i,j,k) ) / dx(1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (ppm_type .eq. 0) then

       do k = lo(3)-1,hi(3)+1
          call slopex_2d(s(:,:,k,comp:),slopex(:,:,k,:),domlo,domhi,lo,hi,ng_s,1,adv_bc(:,:,bccomp:))
          call slopey_2d(s(:,:,k,comp:),slopey(:,:,k,:),domlo,domhi,lo,hi,ng_s,1,adv_bc(:,:,bccomp:))
       end do
       call slopez_3d(s(:,:,:,comp:),slopez,domlo,domhi,lo,hi,ng_s,1,adv_bc(:,:,bccomp:))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(s(:,:,:,comp),ng_s,umac,vmac,wmac,ng_um,Ip,Im, &
                    domlo,domhi,lo,hi,adv_bc(:,:,bccomp),dx,dt,.true.)
       if (ppm_trace_forces .eq. 1) then
          call ppm_3d(force(:,:,:,comp),ng_s,umac,vmac,wmac,ng_um,Ipf,Imf, &
                       domlo,domhi,lo,hi,adv_bc(:,:,bccomp),dx,dt,.true.)
       endif
    end if

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

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    call bl_allocate(slx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(srx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(simhx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    ! loop over appropriate x-faces
    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! make slx, srx with 1D extrapolation
                slx(i,j,k) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,1)
                srx(i,j,k) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i  ,j,k,1)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! make slx, srx with 1D extrapolation
                slx(i,j,k) = Ip(i-1,j,k,1)
                srx(i,j,k) = Im(i  ,j,k,1)
             end do
          end do
       end do
    end if

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          slx(is,js-1:je+1,ks-1:ke+1) = s(is-1,js-1:je+1,ks-1:ke+1,comp)
          srx(is,js-1:je+1,ks-1:ke+1) = s(is-1,js-1:je+1,ks-1:ke+1,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             srx(is,js-1:je+1,ks-1:ke+1) = min(srx(is,js-1:je+1,ks-1:ke+1),0.d0)
          end if
          slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          slx(is,js-1:je+1,ks-1:ke+1) = 0.d0
          srx(is,js-1:je+1,ks-1:ke+1) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          slx(ie+1,js-1:je+1,ks-1:ke+1) = s(ie+1,js-1:je+1,ks-1:ke+1,comp)
          srx(ie+1,js-1:je+1,ks-1:ke+1) = s(ie+1,js-1:je+1,ks-1:ke+1,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             slx(ie+1,js-1:je+1,ks-1:ke+1) = max(slx(ie+1,js-1:je+1,ks-1:ke+1),0.d0)
          end if
          srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          slx(ie+1,js-1:je+1,ks-1:ke+1) = 0.d0
          srx(ie+1,js-1:je+1,ks-1:ke+1) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             ! make simhx by solving Riemann problem
             simhx(i,j,k) = merge(slx(i,j,k),srx(i,j,k),umac(i,j,k) .gt. 0.d0)
             savg = HALF*(slx(i,j,k)+srx(i,j,k))
             simhx(i,j,k) = merge(simhx(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(slopex)

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    call bl_allocate(sly,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(sry,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(simhy,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1)

    ! loop over appropriate y-faces
    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                ! make sly, sry with 1D extrapolation
                sly(i,j,k) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,1)
                sry(i,j,k) = s(i,j  ,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j  ,k,1)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                ! make sly, sry with 1D extrapolation
                sly(i,j,k) = Ip(i,j-1,k,2)
                sry(i,j,k) = Im(i,j  ,k,2)
             enddo
          enddo
       enddo
    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
          sly(is-1:ie+1,js,ks-1:ke+1) = s(is-1:ie+1,js-1,ks-1:ke+1,comp)
          sry(is-1:ie+1,js,ks-1:ke+1) = s(is-1:ie+1,js-1,ks-1:ke+1,comp)
       else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sry(is-1:ie+1,js,ks-1:ke+1) = min(sry(is-1:ie+1,js,ks-1:ke+1),0.d0)
          end if
          sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
          sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
          sly(is-1:ie+1,js,ks-1:ke+1) = 0.d0
          sry(is-1:ie+1,js,ks-1:ke+1) = 0.d0
       else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
          sly(is-1:ie+1,je+1,ks-1:ke+1) = s(is-1:ie+1,je+1,ks-1:ke+1,comp)
          sry(is-1:ie+1,je+1,ks-1:ke+1) = s(is-1:ie+1,je+1,ks-1:ke+1,comp)
       else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sly(is-1:ie+1,je+1,ks-1:ke+1) = max(sly(is-1:ie+1,je+1,ks-1:ke+1),0.d0)
          end if
          sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
          sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
          sly(is-1:ie+1,je+1,ks-1:ke+1) = 0.d0
          sry(is-1:ie+1,je+1,ks-1:ke+1) = 0.d0
       else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             ! make simhy by solving Riemann problem
             simhy(i,j,k) = merge(sly(i,j,k),sry(i,j,k),vmac(i,j,k) .gt. 0.d0)
             savg = HALF*(sly(i,j,k)+sry(i,j,k))
             simhy(i,j,k) = merge(simhy(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(slopey)

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    call bl_allocate(slz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(srz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(simhz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1)

    ! loop over appropriate z-faces
    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                ! make slz, srz with 1D extrapolation
                slz(i,j,k) = s(i,j,k-1,comp) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,1)
                srz(i,j,k) = s(i,j,k  ,comp) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k  ,1)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                ! make slz, srz with 1D extrapolation
                slz(i,j,k) = Ip(i,j,k-1,3)
                srz(i,j,k) = Im(i,j,k  ,3)
             enddo
          enddo
       enddo
    end if

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
          slz(is-1:ie+1,js-1:je+1,ks) = s(is-1:ie+1,js-1:je+1,ks,comp)
          srz(is-1:ie+1,js-1:je+1,ks) = s(is-1:ie+1,js-1:je+1,ks,comp)
       else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             srz(is-1:ie+1,js-1:je+1,ks) = min(srz(is-1:ie+1,js-1:je+1,ks),0.d0)
          end if
          slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
          slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
          slz(is-1:ie+1,js-1:je+1,ks) = 0.d0
          srz(is-1:ie+1,js-1:je+1,ks) = 0.d0
       else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
          slz(is-1:ie+1,js-1:je+1,ke+1) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
          srz(is-1:ie+1,js-1:je+1,ke+1) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
       else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             slz(is-1:ie+1,js-1:je+1,ke+1) = max(slz(is-1:ie+1,js-1:je+1,ke+1),0.d0)
          end if
          srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
          srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
          slz(is-1:ie+1,js-1:je+1,ke+1) = 0.d0
          srz(is-1:ie+1,js-1:je+1,ke+1) = 0.d0
       else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             ! make simhz by solving Riemann problem
             simhz(i,j,k) = merge(slz(i,j,k),srz(i,j,k),wmac(i,j,k) .gt. 0.d0)
             savg = HALF*(slz(i,j,k)+srz(i,j,k))
             simhz(i,j,k) = merge(simhz(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(slopez)

    !******************************************************************
    ! Create s_{\i-\half\e_x}^{x|y}, etc.
    !******************************************************************

    call bl_allocate(slxy,lo(1),hi(1)+1,lo(2),hi(2),lo(3)-1,hi(3)+1)
    call bl_allocate(srxy,lo(1),hi(1)+1,lo(2),hi(2),lo(3)-1,hi(3)+1)
    call bl_allocate(simhxy,lo(1),hi(1)+1,lo(2),hi(2),lo(3)-1,hi(3)+1)

    ! loop over appropriate xy faces
    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js,je
             do i=is,ie+1
                ! make slxy, srxy by updating 1D extrapolation
                slxy(i,j,k) = slx(i,j,k) &
                     - (dt3/hy)*(simhy(i-1,j+1,k)*vmac(i-1,j+1,k) &
                     - simhy(i-1,j,k)*vmac(i-1,j,k)) &
                     - dt3*s(i-1,j,k,comp)*divu(i-1,j,k) &
                     + (dt3/hy)*s(i-1,j,k,comp)*(vmac(i-1,j+1,k)-vmac(i-1,j,k))
                srxy(i,j,k) = srx(i,j,k) &
                     - (dt3/hy)*(simhy(i  ,j+1,k)*vmac(i  ,j+1,k) &
                     - simhy(i  ,j,k)*vmac(i  ,j,k)) &
                     - dt3*s(i,j,k,comp)*divu(i,j,k) &
                     + (dt3/hy)*s(i,j,k,comp)*(vmac(i,j+1,k)-vmac(i,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js,je
             do i=is,ie+1
                ! make slxy, srxy by updating 1D extrapolation
                slxy(i,j,k) = slx(i,j,k) &
                     - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k)) &
                     *(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                srxy(i,j,k) = srx(i,j,k) &
                     - (dt6/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k)) &
                     *(simhy(i  ,j+1,k)-simhy(i  ,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          slxy(is,js:je,ks-1:ke+1) = s(is-1,js:je,ks-1:ke+1,comp)
          srxy(is,js:je,ks-1:ke+1) = s(is-1,js:je,ks-1:ke+1,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             srxy(is,js:je,ks-1:ke+1) = min(srxy(is,js:je,ks-1:ke+1),0.d0)
          end if
          slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          slxy(is,js:je,ks-1:ke+1) = 0.d0
          srxy(is,js:je,ks-1:ke+1) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          slxy(ie+1,js:je,ks-1:ke+1) = s(ie+1,js:je,ks-1:ke+1,comp)
          srxy(ie+1,js:je,ks-1:ke+1) = s(ie+1,js:je,ks-1:ke+1,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             slxy(ie+1,js:je,ks-1:ke+1) = max(slxy(ie+1,js:je,ks-1:ke+1),0.d0)
          end if
          srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          slxy(ie+1,js:je,ks-1:ke+1) = 0.d0
          srxy(ie+1,js:je,ks-1:ke+1) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! make simhxy by solving Riemann problem
             simhxy(i,j,k) = merge(slxy(i,j,k),srxy(i,j,k),umac(i,j,k) .gt. 0.d0)
             savg = HALF*(slxy(i,j,k)+srxy(i,j,k))
             simhxy(i,j,k) = merge(simhxy(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(slxy)
    call bl_deallocate(srxy)

    ! loop over appropriate xz faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    call bl_allocate(slxz,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3))
    call bl_allocate(srxz,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3))
    call bl_allocate(simhxz,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3))

    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js-1,je+1
             do i=is,ie+1
                ! make slxz, srxz by updating 1D extrapolation
                slxz(i,j,k) = slx(i,j,k) &
                     - (dt3/hz)*(simhz(i-1,j,k+1)*wmac(i-1,j,k+1) &
                     - simhz(i-1,j,k)*wmac(i-1,j,k)) &
                     - dt3*s(i-1,j,k,comp)*divu(i-1,j,k) &
                     + (dt3/hz)*s(i-1,j,k,comp)*(wmac(i-1,j,k+1)-wmac(i-1,j,k))
                srxz(i,j,k) = srx(i,j,k) &
                     - (dt3/hz)*(simhz(i  ,j,k+1)*wmac(i  ,j,k+1) &
                     - simhz(i  ,j,k)*wmac(i  ,j,k)) &
                     - dt3*s(i,j,k,comp)*divu(i,j,k) &
                     + (dt3/hz)*s(i,j,k,comp)*(wmac(i,j,k+1)-wmac(i,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js-1,je+1
             do i=is,ie+1
                ! make slxz, srxz by updating 1D extrapolation
                slxz(i,j,k) = slx(i,j,k) &
                     - (dt6/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k)) &
                     *(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                srxz(i,j,k) = srx(i,j,k) &
                     - (dt6/hz)*(wmac(i  ,j,k+1)+wmac(i  ,j,k)) &
                     *(simhz(i  ,j,k+1)-simhz(i  ,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          slxz(is,js-1:je+1,ks:ke) = s(is-1,js-1:je+1,ks:ke,comp)
          srxz(is,js-1:je+1,ks:ke) = s(is-1,js-1:je+1,ks:ke,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             srxz(is,js-1:je+1,ks:ke) = min(srxz(is,js-1:je+1,ks:ke),0.d0)
          end if
          slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          slxz(is,js-1:je+1,ks:ke) = 0.d0
          srxz(is,js-1:je+1,ks:ke) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          slxz(ie+1,js-1:je+1,ks:ke) = s(ie+1,js-1:je+1,ks:ke,comp)
          srxz(ie+1,js-1:je+1,ks:ke) = s(ie+1,js-1:je+1,ks:ke,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             slxz(ie+1,js-1:je+1,ks:ke) = max(slxz(ie+1,js-1:je+1,ks:ke),0.d0)
          end if
          srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          slxz(ie+1,js-1:je+1,ks:ke) = 0.d0
          srxz(ie+1,js-1:je+1,ks:ke) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! make simhxz by solving Riemann problem
             simhxz(i,j,k) = merge(slxz(i,j,k),srxz(i,j,k),umac(i,j,k) .gt. 0.d0)
             savg = HALF*(slxz(i,j,k)+srxz(i,j,k))
             simhxz(i,j,k) = merge(simhxz(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(slxz)
    call bl_deallocate(srxz)

    ! loop over appropriate yx faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    call bl_allocate(slyx,lo(1),hi(1),lo(2),hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(sryx,lo(1),hi(1),lo(2),hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(simhyx,lo(1),hi(1),lo(2),hi(2)+1,lo(3)-1,hi(3)+1)

    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is,ie
                ! make slyx, sryx by updating 1D extrapolation
                slyx(i,j,k) = sly(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j-1,k)*umac(i+1,j-1,k) &
                     - simhx(i,j-1,k)*umac(i,j-1,k)) &
                     - dt3*s(i,j-1,k,comp)*divu(i,j-1,k) &
                     + (dt3/hx)*s(i,j-1,k,comp)*(umac(i+1,j-1,k)-umac(i,j-1,k))
                sryx(i,j,k) = sry(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j  ,k)*umac(i+1,j  ,k) &
                     - simhx(i,j  ,k)*umac(i,j  ,k)) &
                     - dt3*s(i,j,k,comp)*divu(i,j,k) &
                     + (dt3/hx)*s(i,j,k,comp)*(umac(i+1,j,k)-umac(i,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is,ie
                ! make slyx, sryx by updating 1D extrapolation
                slyx(i,j,k) = sly(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k)) &
                     *(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                sryx(i,j,k) = sry(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j  ,k)+umac(i,j  ,k)) &
                     *(simhx(i+1,j  ,k)-simhx(i,j  ,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
          slyx(is:ie,js,ks-1:ke+1) = s(is:ie,js-1,ks-1:ke+1,comp)
          sryx(is:ie,js,ks-1:ke+1) = s(is:ie,js-1,ks-1:ke+1,comp)
       else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sryx(is:ie,js,ks-1:ke+1) = min(sryx(is:ie,js,ks-1:ke+1),0.d0)
          end if
          slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
          slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
          slyx(is:ie,js,ks-1:ke+1) = 0.d0
          sryx(is:ie,js,ks-1:ke+1) = 0.d0
       else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
          slyx(is:ie,je+1,ks-1:ke+1) = s(is:ie,je+1,ks-1:ke+1,comp)
          sryx(is:ie,je+1,ks-1:ke+1) = s(is:ie,je+1,ks-1:ke+1,comp)
       else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             slyx(is:ie,je+1,ks-1:ke+1) = max(slyx(is:ie,je+1,ks-1:ke+1),0.d0)
          end if
          sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
          sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
          slyx(is:ie,je+1,ks-1:ke+1) = 0.d0
          sryx(is:ie,je+1,ks-1:ke+1) = 0.d0
       else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! make simhyx by solving Riemann problem
             simhyx(i,j,k) = merge(slyx(i,j,k),sryx(i,j,k),vmac(i,j,k) .gt. 0.d0)
             savg = HALF*(slyx(i,j,k)+sryx(i,j,k))
             simhyx(i,j,k) = merge(simhyx(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(slyx)
    call bl_deallocate(sryx)

    ! loop over appropriate yz faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    call bl_allocate(slyz,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3))
    call bl_allocate(sryz,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3))
    call bl_allocate(simhyz,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3))

    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is-1,ie+1
                ! make slyz, sryz by updating 1D extrapolation
                slyz(i,j,k) = sly(i,j,k) &
                     - (dt3/hz)*(simhz(i,j-1,k+1)*wmac(i,j-1,k+1) &
                     - simhz(i,j-1,k)*wmac(i,j-1,k)) &
                     - dt3*s(i,j-1,k,comp)*divu(i,j-1,k) &
                     + (dt3/hz)*s(i,j-1,k,comp)*(wmac(i,j-1,k+1)-wmac(i,j-1,k))
                sryz(i,j,k) = sry(i,j,k) &
                     - (dt3/hz)*(simhz(i,j  ,k+1)*wmac(i,j  ,k+1) &
                     - simhz(i,j  ,k)*wmac(i,j  ,k)) &
                     - dt3*s(i,j,k,comp)*divu(i,j,k) &
                     + (dt3/hz)*s(i,j,k,comp)*(wmac(i,j,k+1)-wmac(i,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is-1,ie+1
                ! make slyz, sryz by updating 1D extrapolation
                slyz(i,j,k) = sly(i,j,k) &
                     - (dt6/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k)) &
                     *(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                sryz(i,j,k) = sry(i,j,k) &
                     - (dt6/hz)*(wmac(i,j  ,k+1)+wmac(i,j  ,k)) &
                     *(simhz(i,j  ,k+1)-simhz(i,j  ,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
          slyz(is-1:ie+1,js,ks:ke) = s(is-1:ie+1,js-1,ks:ke,comp)
          sryz(is-1:ie+1,js,ks:ke) = s(is-1:ie+1,js-1,ks:ke,comp)
       else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sryz(is-1:ie+1,js,ks:ke) = min(sryz(is-1:ie+1,js,ks:ke),0.d0)
          end if
          slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
          slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
          slyz(is-1:ie+1,js,ks:ke) = 0.d0
          sryz(is-1:ie+1,js,ks:ke) = 0.d0
       else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
          slyz(is-1:ie+1,je+1,ks:ke) = s(is-1:ie+1,je+1,ks:ke,comp)
          sryz(is-1:ie+1,je+1,ks:ke) = s(is-1:ie+1,je+1,ks:ke,comp)
       else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             slyz(is-1:ie+1,je+1,ks:ke) = max(slyz(is-1:ie+1,je+1,ks:ke),0.d0)
          end if
          sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
          sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
          slyz(is-1:ie+1,je+1,ks:ke) = 0.d0
          sryz(is-1:ie+1,je+1,ks:ke) = 0.d0
       else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! make simhyz by solving Riemann problem
             simhyz(i,j,k) = merge(slyz(i,j,k),sryz(i,j,k),vmac(i,j,k) .gt. 0.d0)
             savg = HALF*(slyz(i,j,k)+sryz(i,j,k))
             simhyz(i,j,k) = merge(simhyz(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(simhz)
    call bl_deallocate(slyz)
    call bl_deallocate(sryz)

    ! loop over appropriate zx faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    call bl_allocate(slzx,lo(1),hi(1),lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(srzx,lo(1),hi(1),lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(simhzx,lo(1),hi(1),lo(2)-1,hi(2)+1,lo(3),hi(3)+1)

    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is,ie
                ! make slzx, srzx by updating 1D extrapolation
                slzx(i,j,k) = slz(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j,k-1)*umac(i+1,j,k-1) &
                     - simhx(i,j,k-1)*umac(i,j,k-1)) &
                     - dt3*s(i,j,k-1,comp)*divu(i,j,k-1) &
                     + (dt3/hx)*s(i,j,k-1,comp)*(umac(i+1,j,k-1)-umac(i,j,k-1))
                srzx(i,j,k) = srz(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j,k  )*umac(i+1,j,k  ) &
                     - simhx(i,j,k  )*umac(i,j,k  )) &
                     - dt3*s(i,j,k,comp)*divu(i,j,k) &
                     + (dt3/hx)*s(i,j,k,comp)*(umac(i+1,j,k)-umac(i,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is,ie
                ! make slzx, srzx by updating 1D extrapolation
                slzx(i,j,k) = slz(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1)) &
                     *(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                srzx(i,j,k) = srz(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  )) &
                     *(simhx(i+1,j,k  )-simhx(i,j,k  ))
             enddo
          enddo
       end do
       !$OMP END PARALLEL DO
    end if

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
          slzx(is:ie,js-1:je+1,ks) = s(is:ie,js-1:je+1,ks-1,comp)
          srzx(is:ie,js-1:je+1,ks) = s(is:ie,js-1:je+1,ks-1,comp)
       else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             srzx(is:ie,js-1:je+1,ks) = min(srzx(is:ie,js-1:je+1,ks),0.d0)
          end if
          slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
          slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
          slzx(is:ie,js-1:je+1,ks) = 0.d0
          srzx(is:ie,js-1:je+1,ks) = 0.d0
       else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
          slzx(is:ie,js-1:je+1,ke+1) = s(is:ie,js-1:je+1,ke+1,comp)
          srzx(is:ie,js-1:je+1,ke+1) = s(is:ie,js-1:je+1,ke+1,comp)
       else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             slzx(is:ie,js-1:je+1,ke+1) = max(slzx(is:ie,js-1:je+1,ke+1),0.d0)
          end if
          srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
          srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
          slzx(is:ie,js-1:je+1,ke+1) = 0.d0
          srzx(is:ie,js-1:je+1,ke+1) = 0.d0
       else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! make simhzx by solving Riemann problem
             simhzx(i,j,k) = merge(slzx(i,j,k),srzx(i,j,k),wmac(i,j,k) .gt. 0.d0)
             savg = HALF*(slzx(i,j,k)+srzx(i,j,k))
             simhzx(i,j,k) = merge(simhzx(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(simhx)
    call bl_deallocate(slzx)
    call bl_deallocate(srzx)

    ! loop over appropriate zy faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    call bl_allocate(slzy,lo(1)-1,hi(1)+1,lo(2),hi(2),lo(3),hi(3)+1)
    call bl_allocate(srzy,lo(1)-1,hi(1)+1,lo(2),hi(2),lo(3),hi(3)+1)
    call bl_allocate(simhzy,lo(1)-1,hi(1)+1,lo(2),hi(2),lo(3),hi(3)+1)

    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is-1,ie+1
                ! make slzy, srzy by updating 1D extrapolation
                slzy(i,j,k) = slz(i,j,k) &
                     - (dt3/hy)*(simhy(i,j+1,k-1)*vmac(i,j+1,k-1) &
                     - simhy(i,j,k-1)*vmac(i,j,k-1)) &
                     - dt3*s(i,j,k-1,comp)*divu(i,j,k-1) &
                     + (dt3/hy)*s(i,j,k-1,comp)*(vmac(i,j+1,k-1)-vmac(i,j,k-1))
                srzy(i,j,k) = srz(i,j,k) &
                     - (dt3/hy)*(simhy(i,j+1,k  )*vmac(i,j+1,k  ) &
                     - simhy(i,j,k  )*vmac(i,j,k  )) &
                     - dt3*s(i,j,k,comp)*divu(i,j,k) &
                     + (dt3/hy)*s(i,j,k,comp)*(vmac(i,j+1,k)-vmac(i,j,k))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is-1,ie+1
                ! make slzy, srzy by updating 1D extrapolation
                slzy(i,j,k) = slz(i,j,k) &
                     - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1)) &
                     *(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                srzy(i,j,k) = srz(i,j,k) &
                     - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  )) &
                     *(simhy(i,j+1,k  )-simhy(i,j,k  ))
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
          slzy(is-1:ie+1,js:je,ks) = s(is-1:ie+1,js:je,ks-1,comp)
          srzy(is-1:ie+1,js:je,ks) = s(is-1:ie+1,js:je,ks-1,comp)
       else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             srzy(is-1:ie+1,js:je,ks) = min(srzy(is-1:ie+1,js:je,ks),0.d0)
          end if
          slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
          slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
          slzy(is-1:ie+1,js:je,ks) = 0.d0
          srzy(is-1:ie+1,js:je,ks) = 0.d0
       else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
          slzy(is-1:ie+1,js:je,ke+1) = s(is-1:ie+1,js:je,ke+1,comp)
          srzy(is-1:ie+1,js:je,ke+1) = s(is-1:ie+1,js:je,ke+1,comp)
       else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             slzy(is-1:ie+1,js:je,ke+1) = max(slzy(is-1:ie+1,js:je,ke+1),0.d0)
          end if
          srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
          srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
          slzy(is-1:ie+1,js:je,ke+1) = 0.d0
          srzy(is-1:ie+1,js:je,ke+1) = 0.d0
       else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
       end if
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! make simhzy by solving Riemann problem
             simhzy(i,j,k) = merge(slzy(i,j,k),srzy(i,j,k),wmac(i,j,k) .gt. 0.d0)
             savg = HALF*(slzy(i,j,k)+srzy(i,j,k))
             simhzy(i,j,k) = merge(simhzy(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(simhy)
    call bl_deallocate(slzy)
    call bl_deallocate(srzy)
    if (is_conservative .eq. 1) then
       call bl_deallocate(divu)
    end if

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse directions
    call bl_allocate(sedgelx,lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))
    call bl_allocate(sedgerx,lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))

    ! loop over appropriate x-faces
    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! make sedgelx, sedgerx
                fl = merge(force(i-1,j,k,comp), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
                fr = merge(force(i  ,j,k,comp), Imf(i  ,j,k,1), ppm_trace_forces == 0)

                sedgelx(i,j,k) = slx(i,j,k) &
                     - (dt2/hy)*(simhyz(i-1,j+1,k  )*vmac(i-1,j+1,k  ) &
                     - simhyz(i-1,j,k)*vmac(i-1,j,k)) &
                     - (dt2/hz)*(simhzy(i-1,j  ,k+1)*wmac(i-1,j  ,k+1) &
                     - simhzy(i-1,j,k)*wmac(i-1,j,k)) &
                     - (dt2/hx)*s(i-1,j,k,comp)*(umac(i  ,j,k)-umac(i-1,j,k)) &
                     + dt2*fl

                sedgerx(i,j,k) = srx(i,j,k) &
                     - (dt2/hy)*(simhyz(i  ,j+1,k  )*vmac(i  ,j+1,  k) &
                     - simhyz(i  ,j,k)*vmac(i  ,j,k)) &
                     - (dt2/hz)*(simhzy(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) &
                     - simhzy(i  ,j,k)*wmac(i  ,j,k)) &
                     - (dt2/hx)*s(i  ,j,k,comp)*(umac(i+1,j,k)-umac(i  ,j,k)) &
                     + dt2*fr
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! make sedgelx, sedgerx
                fl = merge(force(i-1,j,k,comp), Ipf(i-1,j,k,1), ppm_trace_forces == 0)
                fr = merge(force(i  ,j,k,comp), Ipf(i  ,j,k,1), ppm_trace_forces == 0)

                sedgelx(i,j,k) = slx(i,j,k) &
                     - (dt4/hy)*(vmac(i-1,j+1,k  )+vmac(i-1,j,k))* &
                     (simhyz(i-1,j+1,k  )-simhyz(i-1,j,k)) &
                     - (dt4/hz)*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k))* &
                     (simhzy(i-1,j  ,k+1)-simhzy(i-1,j,k)) &
                     + dt2*fl

                sedgerx(i,j,k) = srx(i,j,k) &
                     - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i  ,j,k))* &
                     (simhyz(i  ,j+1,k  )-simhyz(i  ,j,k)) &
                     - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k))* &
                     (simhzy(i  ,j  ,k+1)-simhzy(i  ,j,k)) &
                     + dt2*fr
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! make sedgex by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgex(i,j,k,comp) = merge(sedgelx(i,j,k),sedgerx(i,j,k),umac(i,j,k) .gt. 0.d0)
             savg = HALF*(sedgelx(i,j,k)+sedgerx(i,j,k))
             sedgex(i,j,k,comp) = merge(sedgex(i,j,k,comp),savg,abs(umac(i,j,k)).gt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       if (adv_bc(1,1,bccomp) .eq. EXT_DIR) then
          sedgex(is,js:je,ks:ke,comp) = s(is-1,js:je,ks:ke,comp)
       else if (adv_bc(1,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             sedgex(is,js:je,ks:ke,comp) = min(sedgerx(is,js:je,ks:ke),0.d0)
          else
             sedgex(is,js:je,ks:ke,comp) = sedgerx(is,js:je,ks:ke)
          end if
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_EVEN) then
          sedgex(is,js:je,ks:ke,comp) = sedgerx(is,js:je,ks:ke)
       else if (adv_bc(1,1,bccomp) .eq. REFLECT_ODD) then
          sedgex(is,js:je,ks:ke,comp) = 0.d0
       else if (adv_bc(1,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       if (adv_bc(1,2,bccomp) .eq. EXT_DIR) then
          sedgex(ie+1,js:je,ks:ke,comp) = s(ie+1,js:je,ks:ke,comp)
       else if (adv_bc(1,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 1) then
             sedgex(ie+1,js:je,ks:ke,comp) = max(sedgelx(ie+1,js:je,ks:ke),0.d0)
          else
             sedgex(ie+1,js:je,ks:ke,comp) = sedgelx(ie+1,js:je,ks:ke)
          end if
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_EVEN) then
          sedgex(ie+1,js:je,ks:ke,comp) = sedgelx(ie+1,js:je,ks:ke)
       else if (adv_bc(1,2,bccomp) .eq. REFLECT_ODD) then
          sedgex(ie+1,js:je,ks:ke,comp) = 0.d0
       else if (adv_bc(1,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(1,2)")
       end if
    end if

    call bl_deallocate(slx)
    call bl_deallocate(srx)
    call bl_deallocate(sedgelx)
    call bl_deallocate(sedgerx)

    call bl_allocate(sedgely,lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))
    call bl_allocate(sedgery,lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))

    ! loop over appropriate y-faces
    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! make sedgely, sedgery
                fl = merge(force(i,j-1,k,comp), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
                fr = merge(force(i,j  ,k,comp), Imf(i,j  ,k,2), ppm_trace_forces == 0)

                sedgely(i,j,k) = sly(i,j,k) &
                     - (dt2/hx)*(simhxz(i+1,j-1,k  )*umac(i+1,j-1,k  ) &
                     - simhxz(i,j-1,k)*umac(i,j-1,k)) &
                     - (dt2/hz)*(simhzx(i  ,j-1,k+1)*wmac(i  ,j-1,k+1) &
                     - simhzx(i,j-1,k)*wmac(i,j-1,k)) &
                     - (dt2/hy)*s(i,j-1,k,comp)*(vmac(i,j  ,k)-vmac(i,j-1,k)) &
                     + dt2*fl

                sedgery(i,j,k) = sry(i,j,k) &
                     - (dt2/hx)*(simhxz(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
                     - simhxz(i,j  ,k)*umac(i,j  ,k)) &
                     - (dt2/hz)*(simhzx(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) &
                     - simhzx(i,j  ,k)*wmac(i,j  ,k)) &
                     - (dt2/hy)*s(i,j  ,k,comp)*(vmac(i,j+1,k)-vmac(i,j  ,k)) &
                     + dt2*fr
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! make sedgely, sedgery
                fl = merge(force(i,j-1,k,comp), Ipf(i,j-1,k,2), ppm_trace_forces == 0)
                fr = merge(force(i,j  ,k,comp), Imf(i,j  ,k,2), ppm_trace_forces == 0)

                sedgely(i,j,k) = sly(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j-1,k  )+umac(i,j-1,k))* &
                     (simhxz(i+1,j-1,k  )-simhxz(i,j-1,k)) &
                     - (dt4/hz)*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k))* &
                     (simhzx(i  ,j-1,k+1)-simhzx(i,j-1,k)) &
                     + dt2*fl

                sedgery(i,j,k) = sry(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j  ,k))* &
                     (simhxz(i+1,j  ,k  )-simhxz(i,j  ,k)) &
                     - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k))* &
                     (simhzx(i  ,j  ,k+1)-simhzx(i,j  ,k)) &
                     + dt2*fr
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! make sedgey by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgey(i,j,k,comp) = merge(sedgely(i,j,k),sedgery(i,j,k),vmac(i,j,k) .gt. 0.d0)
             savg = HALF*(sedgely(i,j,k)+sedgery(i,j,k))
             sedgey(i,j,k,comp) = merge(sedgey(i,j,k,comp),savg,abs(vmac(i,j,k)).gt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       if (adv_bc(2,1,bccomp) .eq. EXT_DIR) then
          sedgey(is:ie,js,ks:ke,comp) = s(is:ie,js-1,ks:ke,comp)
       else if (adv_bc(2,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sedgey(is:ie,js,ks:ke,comp) = min(sedgery(is:ie,js,ks:ke),0.d0)
          else
             sedgey(is:ie,js,ks:ke,comp) = sedgery(is:ie,js,ks:ke)
          end if
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_EVEN) then
          sedgey(is:ie,js,ks:ke,comp) = sedgery(is:ie,js,ks:ke)
       else if (adv_bc(2,1,bccomp) .eq. REFLECT_ODD) then
          sedgey(is:ie,js,ks:ke,comp) = 0.d0
       else if (adv_bc(2,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       if (adv_bc(2,2,bccomp) .eq. EXT_DIR) then
          sedgey(is:ie,je+1,ks:ke,comp) = s(is:ie,je+1,ks:ke,comp)
       else if (adv_bc(2,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 2) then
             sedgey(is:ie,je+1,ks:ke,comp) = max(sedgely(is:ie,je+1,ks:ke),0.d0)
          else
             sedgey(is:ie,je+1,ks:ke,comp) = sedgely(is:ie,je+1,ks:ke)
          end if
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_EVEN) then
          sedgey(is:ie,je+1,ks:ke,comp) = sedgely(is:ie,je+1,ks:ke)
       else if (adv_bc(2,2,bccomp) .eq. REFLECT_ODD) then
          sedgey(is:ie,je+1,ks:ke,comp) = 0.d0
       else if (adv_bc(2,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(2,2)")
       end if
    end if

    call bl_deallocate(sly)
    call bl_deallocate(sry)
    call bl_deallocate(sedgely)
    call bl_deallocate(sedgery)

    call bl_allocate(sedgelz,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)
    call bl_allocate(sedgerz,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)

    ! loop over appropriate z-faces
    if (is_conservative .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! make sedgelz, sedgerz
                fl = merge(force(i,j,k-1,comp), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
                fr = merge(force(i,j,k  ,comp), Imf(i,j,k  ,3), ppm_trace_forces == 0)

                sedgelz(i,j,k) = slz(i,j,k) &
                     - (dt2/hx)*(simhxy(i+1,j  ,k-1)*umac(i+1,j  ,k-1) &
                     - simhxy(i,j,k-1)*umac(i,j,k-1)) &
                     - (dt2/hy)*(simhyx(i  ,j+1,k-1)*vmac(i  ,j+1,k-1) &
                     - simhyx(i,j,k-1)*vmac(i,j,k-1)) &
                     - (dt2/hz)*s(i,j,k-1,comp)*(wmac(i,j,k  )-wmac(i,j,k-1)) &
                     + dt2*fl

                sedgerz(i,j,k) = srz(i,j,k) &
                     - (dt2/hx)*(simhxy(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
                     - simhxy(i,j,k  )*umac(i,j,k  )) &
                     - (dt2/hy)*(simhyx(i  ,j+1,k  )*vmac(i  ,j+1,k  ) &
                     - simhyx(i,j,k  )*vmac(i,j,k  )) &
                     - (dt2/hz)*s(i,j,k  ,comp)*(wmac(i,j,k+1)-wmac(i,j,k  )) &
                     + dt2*fr
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! make sedgelz, sedgerz
                fl = merge(force(i,j,k-1,comp), Ipf(i,j,k-1,3), ppm_trace_forces == 0)
                fr = merge(force(i,j,k  ,comp), Imf(i,j,k  ,3), ppm_trace_forces == 0)

                sedgelz(i,j,k) = slz(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j,k-1)) &
                     *(simhxy(i+1,j  ,k-1)-simhxy(i,j,k-1)) &
                     - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1)) &
                     *(simhyx(i  ,j+1,k-1)-simhyx(i,j,k-1)) &
                     + dt2*fl

                sedgerz(i,j,k) = srz(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j,k  )) &
                     *(simhxy(i+1,j  ,k  )-simhxy(i,j,k  )) &
                     - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i,j,k  )) &
                     *(simhyx(i  ,j+1,k  )-simhyx(i,j,k  )) &
                     + dt2*fr
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO PRIVATE(i,j,k,savg)
    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! make sedgez by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgez(i,j,k,comp) = merge(sedgelz(i,j,k),sedgerz(i,j,k),wmac(i,j,k) .gt. 0.d0)
             savg = HALF*(sedgelz(i,j,k)+sedgerz(i,j,k))
             sedgez(i,j,k,comp) = merge(sedgez(i,j,k,comp),savg,abs(wmac(i,j,k)).gt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       if (adv_bc(3,1,bccomp) .eq. EXT_DIR) then
          sedgez(is:ie,js:je,ks,comp) = s(is:ie,js:je,ks-1,comp)
       else if (adv_bc(3,1,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,1,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             sedgez(is:ie,js:je,ks,comp) = min(sedgerz(is:ie,js:je,ks),0.d0)
          else
             sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je,ks)
          end if
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_EVEN) then
          sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je,ks)
       else if (adv_bc(3,1,bccomp) .eq. REFLECT_ODD) then
          sedgez(is:ie,js:je,ks,comp) = 0.d0
       else if (adv_bc(3,1,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,1)")
       end if
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       if (adv_bc(3,2,bccomp) .eq. EXT_DIR) then
          sedgez(is:ie,js:je,ke+1,comp) = s(is:ie,js:je,ke+1,comp)
       else if (adv_bc(3,2,bccomp) .eq. FOEXTRAP .or. &
            adv_bc(3,2,bccomp) .eq. HOEXTRAP) then
          if (is_vel .eq. 1 .and. comp .eq. 3) then
             sedgez(is:ie,js:je,ke+1,comp) = max(sedgelz(is:ie,js:je,ke+1),0.d0)
          else
             sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je,ke+1)
          end if
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_EVEN) then
          sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je,ke+1)
       else if (adv_bc(3,2,bccomp) .eq. REFLECT_ODD) then
          sedgez(is:ie,js:je,ke+1,comp) = 0.d0
       else if (adv_bc(3,2,bccomp) .eq. INT_DIR) then
       else
          call amrex_error("make_edge_scal_3d: invalid boundary type adv_bc(3,2)")
       end if
    end if

    call bl_deallocate(Ip)
    call bl_deallocate(Im)
    call bl_deallocate(Ipf)
    call bl_deallocate(Imf)

    call bl_deallocate(slz)
    call bl_deallocate(srz)

    call bl_deallocate(simhxy)
    call bl_deallocate(simhxz)
    call bl_deallocate(simhyx)
    call bl_deallocate(simhyz)
    call bl_deallocate(simhzx)
    call bl_deallocate(simhzy)

    call bl_deallocate(sedgelz)
    call bl_deallocate(sedgerz)

  end subroutine make_edge_scal_3d

#endif

end module make_edge_scal_module
