
#include "AMReX_BC_TYPES.H"

module mkutrans_module

  use amrex_error_module
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
  subroutine mkutrans_2d(lev, domlo, domhi, lo, hi, &
       utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
       ufull,  uf_lo, uf_hi, nc_uf, ng_uf, &
       utrans, uu_lo, uu_hi, &
       vtrans, uv_lo, uv_hi, &
       w0,dx,dt,adv_bc,phys_bc) bind(C,name="mkutrans_2d")

    integer         , intent(in   ) :: lev, domlo(2), domhi(2), lo(2), hi(2)
    integer         , intent(in   ) :: ut_lo(2), ut_hi(2), nc_ut
    integer, value  , intent(in   ) :: ng_ut
    integer         , intent(in   ) :: uf_lo(2), uf_hi(2), nc_uf
    integer, value  , intent(in   ) :: ng_uf
    integer         , intent(in   ) :: uu_lo(2), uu_hi(2)
    integer         , intent(in   ) :: uv_lo(2), uv_hi(2)
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2))
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(2), dt
    integer         , intent(in   ) :: adv_bc(2,2,2), phys_bc(2,2) ! dim, lohi, (comp)

    double precision, pointer :: slopex(:,:,:)
    double precision, pointer :: slopey(:,:,:)

    double precision, pointer :: Ip(:,:,:)
    double precision, pointer :: Im(:,:,:)

    double precision, pointer :: ulx(:,:),urx(:,:)
    double precision, pointer :: vly(:,:),vry(:,:)

    double precision hx,hy,dt2,uavg

    integer :: i,j,is,js,ie,je

    logical :: test

    call bl_allocate(slopex,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,1)
    call bl_allocate(slopey,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,1)

    call bl_allocate(ulx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1)
    call bl_allocate(urx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1)

    call bl_allocate(vly,lo(1)-1,hi(1)+1,lo(2),hi(2)+1)
    call bl_allocate(vry,lo(1)-1,hi(1)+1,lo(2),hi(2)+1)

    call bl_allocate(Ip,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)
    call bl_allocate(Im,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,1,2)

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)

    dt2 = HALF*dt

    hx = dx(1)
    hy = dx(2)

    if (ppm_type .eq. 0) then
       call slopex_2d(utilde(:,:,1:1),slopex,domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,1:1))
       call slopey_2d(utilde(:,:,2:2),slopey,domlo,domhi,lo,hi,ng_ut,1,adv_bc(:,:,2:2))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(utilde(:,:,1),ng_ut, &
            ufull(:,:,1),ufull(:,:,2),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
    end if

    !******************************************************************
    ! create utrans
    !******************************************************************

    if (ppm_type .eq. 0) then
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             ulx(i,j) = utilde(i-1,j,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,1)))*slopex(i-1,j,1)
             urx(i,j) = utilde(i  ,j,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,1)))*slopex(i  ,j,1)
          end do
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             ulx(i,j) = Ip(i-1,j,1)
             urx(i,j) = Im(i  ,j,1)
          end do
       end do
    end if

    ! impose lo i side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          ulx(is,js:je) = utilde(is-1,js:je,1)
          urx(is,js:je) = utilde(is-1,js:je,1)
       case (SlipWall, NoSlipWall, Symmetry)
          ulx(is,js:je) = ZERO
          urx(is,js:je) = ZERO
       case (Outflow)
          ulx(is,js:je) = min(urx(is,js:je),ZERO)
          urx(is,js:je) = ulx(is,js:je)
       case (Interior)
       case  default
          call amrex_error("mkutrans_2d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi i side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          ulx(ie+1,js:je) = utilde(ie+1,js:je,1)
          urx(ie+1,js:je) = utilde(ie+1,js:je,1)
       case (SlipWall, NoSlipWall, Symmetry)
          ulx(ie+1,js:je) = ZERO
          urx(ie+1,js:je) = ZERO
       case (Outflow)
          ulx(ie+1,js:je) = max(ulx(ie+1,js:je),ZERO)
          urx(ie+1,js:je) = ulx(ie+1,js:je)
       case (Interior)
       case  default
          call amrex_error("mkutrans_2d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    do j=js,je
       do i=is,ie+1
          ! solve Riemann problem using full velocity
          uavg = HALF*(ulx(i,j)+urx(i,j))
          test = ((ulx(i,j) .le. ZERO .and. urx(i,j) .ge. ZERO) .or. &
               (abs(ulx(i,j)+urx(i,j)) .lt. rel_eps))
          utrans(i,j) = merge(ulx(i,j),urx(i,j),uavg .gt. ZERO)
          utrans(i,j) = merge(ZERO,utrans(i,j),test)
       end do
    end do

    !******************************************************************
    ! create vtrans
    !******************************************************************
    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(utilde(:,:,2),ng_ut, &
            ufull(:,:,1),ufull(:,:,2),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
    end if

    if (ppm_type .eq. 0) then
       do j=js,je+1
          do i=is,ie
             ! extrapolate to edges
             vly(i,j) = utilde(i,j-1,2) + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,2)))*slopey(i,j-1,1)
             vry(i,j) = utilde(i,j  ,2) - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,2)))*slopey(i,j  ,1)
          end do
       end do

    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je+1
          do i=is,ie
             ! extrapolate to edges
             vly(i,j) = Ip(i,j-1,2)
             vry(i,j) = Im(i,j  ,2)
          end do
       end do
    end if

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          vly(is:ie,js) = utilde(is:ie,js-1,2)
          vry(is:ie,js) = utilde(is:ie,js-1,2)
       case (SlipWall, NoSlipWall, Symmetry)
          vly(is:ie,js) = ZERO
          vry(is:ie,js) = ZERO
       case (Outflow)
          vly(is:ie,js) = min(vry(is:ie,js),ZERO)
          vry(is:ie,js) = vly(is:ie,js)
       case (Interior)
       case  default
          call amrex_error("mkutrans_2d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          vly(is:ie,je+1) = utilde(is:ie,je+1,2)
          vry(is:ie,je+1) = utilde(is:ie,je+1,2)
       case (SlipWall, NoSlipWall, Symmetry)
          vly(is:ie,je+1) = ZERO
          vry(is:ie,je+1) = ZERO
       case (Outflow)
          vly(is:ie,je+1) = max(vly(is:ie,je+1),ZERO)
          vry(is:ie,je+1) = vly(is:ie,je+1)
       case (Interior)
       case  default
          call amrex_error("mkutrans_2d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    do j=js,je+1
       do i=is,ie
          ! solve Riemann problem using full velocity
          uavg = HALF*(vly(i,j)+vry(i,j))
          test = ((vly(i,j)+w0(lev,j) .le. ZERO .and. vry(i,j)+w0(lev,j) .ge. ZERO) .or. &
               (abs(vly(i,j)+vry(i,j)+TWO*w0(lev,j)) .lt. rel_eps))
          vtrans(i,j) = merge(vly(i,j),vry(i,j),uavg+w0(lev,j) .gt. ZERO)
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
       enddo
    enddo

    call bl_deallocate(slopex)
    call bl_deallocate(slopey)

    call bl_deallocate(Ip)
    call bl_deallocate(Im)

    call bl_deallocate(ulx)
    call bl_deallocate(urx)
    call bl_deallocate(vly)
    call bl_deallocate(vry)

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

    call bl_allocate(slopex,lo-1,hi+1,1)
    call bl_allocate(slopey,lo-1,hi+1,1)
    call bl_allocate(slopez,lo-1,hi+1,1)

    call bl_allocate(Ip,lo-1,hi+1,3)
    call bl_allocate(Im,lo-1,hi+1,3)

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

    call bl_allocate(ulx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(urx,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! extrapolate to edges
                ulx(i,j,k) = utilde(i-1,j,k,1) &
                     + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,k,1)))*slopex(i-1,j,k,1)
                urx(i,j,k) = utilde(i  ,j,k,1) &
                     - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,k,1)))*slopex(i  ,j,k,1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! extrapolate to edges
                ulx(i,j,k) = Ip(i-1,j,k,1)
                urx(i,j,k) = Im(i  ,j,k,1)
             end do
          end do
       end do
    end if

    call bl_deallocate(slopex)

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
       select case(phys_bc(1,1))
       case (Inflow)
          ulx(is,js:je,ks:ke) = utilde(is-1,js:je,ks:ke,1)
          urx(is,js:je,ks:ke) = utilde(is-1,js:je,ks:ke,1)
       case (SlipWall, NoSlipWall, Symmetry)
          ulx(is,js:je,ks:ke) = ZERO
          urx(is,js:je,ks:ke) = ZERO
       case (Outflow)
          ulx(is,js:je,ks:ke) = min(urx(is,js:je,ks:ke),ZERO)
          urx(is,js:je,ks:ke) = ulx(is,js:je,ks:ke)
       case (Interior)
       case  default
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(1,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
       select case(phys_bc(1,2))
       case (Inflow)
          ulx(ie+1,js:je,ks:ke) = utilde(ie+1,js:je,ks:ke,1)
          urx(ie+1,js:je,ks:ke) = utilde(ie+1,js:je,ks:ke,1)
       case (SlipWall, NoSlipWall, Symmetry)
          ulx(ie+1,js:je,ks:ke) = ZERO
          urx(ie+1,js:je,ks:ke) = ZERO
       case (Outflow)
          ulx(ie+1,js:je,ks:ke) = max(ulx(ie+1,js:je,ks:ke),ZERO)
          urx(ie+1,js:je,ks:ke) = ulx(ie+1,js:je,ks:ke)
       case (Interior)
       case  default
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(1,2)")
       end select
    end if

    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! solve Riemann problem using full velocity
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k)+w0macx(i,j,k) .le. ZERO .and. &
                     urx(i,j,k)+w0macx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)+TWO*w0macx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg+w0macx(i,j,k) .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
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
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k) .le. ZERO .and. urx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    call bl_deallocate(ulx)
    call bl_deallocate(urx)

    !******************************************************************
    ! create vtrans
    !******************************************************************

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(utilde(:,:,:,2),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
    end if

    call bl_allocate(vly,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(vry,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,lo(3)-1,hi(3)+1)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! extrapolate to edges
                vly(i,j,k) = utilde(i,j-1,k,2) &
                     + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,k,2)))*slopey(i,j-1,k,1)
                vry(i,j,k) = utilde(i,j  ,k,2) &
                     - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,k,2)))*slopey(i,j  ,k,1)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! extrapolate to edges
                vly(i,j,k) = Ip(i,j-1,k,2)
                vry(i,j,k) = Im(i,j  ,k,2)
             enddo
          enddo
       enddo
    end if

    call bl_deallocate(slopey)

    ! impose lo side bc's
    if (lo(2) .eq. domlo(2)) then
       select case(phys_bc(2,1))
       case (Inflow)
          vly(is:ie,js,ks:ke) = utilde(is:ie,js-1,ks:ke,2)
          vry(is:ie,js,ks:ke) = utilde(is:ie,js-1,ks:ke,2)
       case (SlipWall, NoSlipWall, Symmetry)
          vly(is:ie,js,ks:ke) = ZERO
          vry(is:ie,js,ks:ke) = ZERO
       case (Outflow)
          vly(is:ie,js,ks:ke) = min(vry(is:ie,js,ks:ke),ZERO)
          vry(is:ie,js,ks:ke) = vly(is:ie,js,ks:ke)
       case (Interior)
       case  default
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(2,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(2) .eq. domhi(2)) then
       select case(phys_bc(2,2))
       case (Inflow)
          vly(is:ie,je+1,ks:ke) = utilde(is:ie,je+1,ks:ke,2)
          vry(is:ie,je+1,ks:ke) = utilde(is:ie,je+1,ks:ke,2)
       case (SlipWall, NoSlipWall, Symmetry)
          vly(is:ie,je+1,ks:ke) = ZERO
          vry(is:ie,je+1,ks:ke) = ZERO
       case (Outflow)
          vly(is:ie,je+1,ks:ke) = max(vly(is:ie,je+1,ks:ke),ZERO)
          vry(is:ie,je+1,ks:ke) = vly(is:ie,je+1,ks:ke)
       case (Interior)
       case  default
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(2,2)")
       end select
    end if

    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k)+w0macy(i,j,k) .le. ZERO .and. &
                     vry(i,j,k)+w0macy(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)+TWO*w0macy(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg+w0macy(i,j,k) .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
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
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k) .le. ZERO .and. vry(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    call bl_deallocate(vly)
    call bl_deallocate(vry)

    !******************************************************************
    ! create wtrans
    !******************************************************************

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(utilde(:,:,:,3),ng_ut, &
            ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
            Ip,Im,domlo,domhi,lo,hi,adv_bc(:,:,3),dx,dt,.false.)
    end if

    call bl_allocate(wlz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1)
    call bl_allocate(wrz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3),hi(3)+1)

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
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
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! extrapolate to edges
                wlz(i,j,k) = Ip(i,j,k-1,3)
                wrz(i,j,k) = Im(i,j,k  ,3)
             end do
          end do
       end do
    end if

    call bl_deallocate(slopez)
    call bl_deallocate(Ip)
    call bl_deallocate(Im)

    ! impose lo side bc's
    if (lo(3) .eq. domlo(3)) then
       select case(phys_bc(3,1))
       case (Inflow)
          wlz(is:ie,js:je,ks) = utilde(is:ie,js:je,ks-1,3)
          wrz(is:ie,js:je,ks) = utilde(is:ie,js:je,ks-1,3)
       case (SlipWall, NoSlipWall, Symmetry)
          wlz(is:ie,js:je,ks) = ZERO
          wrz(is:ie,js:je,ks) = ZERO
       case (Outflow)
          wlz(is:ie,js:je,ks) = min(wrz(is:ie,js:je,ks),ZERO)
          wrz(is:ie,js:je,ks) = wlz(is:ie,js:je,ks)
       case (Interior)
       case  default
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(3,1)")
       end select
    end if

    ! impose hi side bc's
    if (hi(3) .eq. domhi(3)) then
       select case(phys_bc(3,2))
       case (Inflow)
          wlz(is:ie,js:je,ke+1) = utilde(is:ie,js:je,ke+1,3)
          wrz(is:ie,js:je,ke+1) = utilde(is:ie,js:je,ke+1,3)
       case (SlipWall, NoSlipWall, Symmetry)
          wlz(is:ie,js:je,ke+1) = ZERO
          wrz(is:ie,js:je,ke+1) = ZERO
       case (Outflow)
          wlz(is:ie,js:je,ke+1) = max(wlz(is:ie,js:je,ke+1),ZERO)
          wrz(is:ie,js:je,ke+1) = wlz(is:ie,js:je,ke+1)
       case (Interior)
       case  default
          call amrex_error("mkutrans_3d: invalid boundary type phys_bc(3,2)")
       end select
    end if

    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
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
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
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

    call bl_deallocate(wlz)
    call bl_deallocate(wrz)

  end subroutine mkutrans_3d
#endif

end module mkutrans_module
