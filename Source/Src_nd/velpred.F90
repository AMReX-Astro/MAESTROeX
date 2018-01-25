! velpred is called by advance_premac -- it is used to predict the
! normal velocities to the interfaces.  We don't care about the
! transverse velocities here.  The prediction is done piecewise linear (for now)

#include "AMReX_BC_TYPES.H"

module velpred_module

  use amrex_constants_module
  use slope_module
  use meth_params_module, only: rel_eps
  use base_state_geometry_module, only: nr_fine, max_radial_level

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 1)
  subroutine velpred_1d(lev, domlo, domhi, lo, hi, &
                        utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
                        ufull,  uf_lo, uf_hi, nc_uf, &
                        umac,   uu_lo, uu_hi, &
                        force,   f_lo,  f_hi, nc_f, &
                        w0,dx,dt,adv_bc,phys_bc) bind(C,name="velpred_1d")

    integer         , intent(in   ) :: lev, domlo(1), domhi(1), lo(1), hi(1)
    integer         , intent(in   ) :: ut_lo(1), ut_hi(1), nc_ut, ng_ut
    integer         , intent(in   ) :: uf_lo(1), uf_hi(1), nc_uf
    integer         , intent(in   ) :: uu_lo(1), uu_hi(1)
    integer         , intent(in   ) ::  f_lo(1),  f_hi(1), nc_f
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1))
    double precision, intent(in   ) :: force ( f_lo(1): f_hi(1),nc_f)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(1), dt
    integer         , intent(in   ) :: adv_bc(1,2,1), phys_bc(1,2) ! dim, lohi, (comp)

    ! Local variables
    double precision :: slopex(lo(1)-1:hi(1)+1,1)

    ! these correspond to umac_L, etc.
    double precision, allocatable :: umacl(:),umacr(:)

    double precision :: hx, dt2, dt4, uavg

    integer :: i,is,ie

    logical :: test

    allocate(umacl(lo(1):hi(1)+1))
    allocate(umacr(lo(1):hi(1)+1))

    is = lo(1)
    ie = hi(1)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)

    call slopex_1d(utilde,slopex,domlo,domhi,lo,hi,ng_u,1,adv_bc)

    !******************************************************************
    ! Create umac 
    !******************************************************************

    do i=is,ie+1
       ! extrapolate velocity to left face
       umacl(i) = utilde(i-1,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,1)))*slopex(i-1,1) &
            + dt2*force(i-1)
       ! extrapolate velocity to right face
       umacr(i) = utilde(i  ,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,1)))*slopex(i  ,1) &
            + dt2*force(i  )
    end do

    do i=is,ie+1
       ! solve Riemann problem using full velocity
       uavg = HALF*(umacl(i)+umacr(i))
       test = ((umacl(i)+w0(lev,i) .le. ZERO .and. umacr(i)+w0(lev,i) .ge. ZERO) .or. &
           (abs(umacl(i)+umacr(i)+TWO*w0(lev,i)) .lt. rel_eps))
       umac(i) = merge(umacl(i),umacr(i),uavg+w0(lev,i) .gt. ZERO)
       umac(i) = merge(ZERO,umac(i),test)
    enddo

    ! impose lo side bc's
    if (lo(1) .eq. domlo(1)) then
    select case(phys_bc(1,1))
    case (Inflow)
       umac(is) = utilde(is-1,1)
    case (SlipWall, NoSlipWall, Symmetry)
       umac(is) = ZERO
    case (Outflow)
       umac(is) = min(umacr(is),ZERO)
    case (Interior)
    case  default
       call bl_error("velpred_1d: invalid boundary type phys_bc(1,1)")
    end select
    end if

    ! impose hi side bc's
    if (hi(1) .eq. domhi(1)) then
    select case(phys_bc(1,2))
    case (Inflow)
       umac(ie+1) = utilde(ie+1,1)
    case (SlipWall, NoSlipWall, Symmetry)
       umac(ie+1) = ZERO
    case (Outflow)
       umac(ie+1) = max(umacl(ie+1),ZERO)
    case (Interior)
    case  default
       call bl_error("velpred_1d: invalid boundary type phys_bc(1,2)")
    end select
    end if

    deallocate(umacl,umacr)

  end subroutine velpred_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine velpred_2d(lev, domlo, domhi, lo, hi, &
                        utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
                        ufull,  uf_lo, uf_hi, nc_uf, &
                        utrans, uu_lo, uu_hi, &
                        vtrans, uv_lo, uv_hi, &
                        umac  , mu_lo, mu_hi, &
                        vmac  , mv_lo, mv_hi, &
                        force,   f_lo,  f_hi, nc_f, &
                        w0,dx,dt,adv_bc,phys_bc) bind(C,name="velpred_2d")

    integer         , intent(in   ) :: lev, domlo(2), domhi(2), lo(2), hi(2)
    integer         , intent(in   ) :: ut_lo(2), ut_hi(2), nc_ut, ng_ut
    integer         , intent(in   ) :: uf_lo(2), uf_hi(2), nc_uf
    integer         , intent(in   ) :: uu_lo(2), uu_hi(2)
    integer         , intent(in   ) :: uv_lo(2), uv_hi(2)
    integer         , intent(in   ) :: mu_lo(2), mu_hi(2)
    integer         , intent(in   ) :: mv_lo(2), mv_hi(2)
    integer         , intent(in   ) ::  f_lo(2),  f_hi(2), nc_f
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
    double precision :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)
    double precision :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)

    ! these correspond to u_L^x, etc.
    double precision, allocatable :: ulx(:,:,:),urx(:,:,:),uimhx(:,:,:)
    double precision, allocatable :: uly(:,:,:),ury(:,:,:),uimhy(:,:,:)

    ! these correspond to umac_L, etc.
    double precision, allocatable :: umacl(:,:),umacr(:,:)
    double precision, allocatable :: vmacl(:,:),vmacr(:,:)

    double precision :: hx, hy, dt2, dt4, uavg, maxu, minu
    double precision :: fl, fr

    integer :: i,j,is,js,ie,je

    logical :: test

    allocate(  ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(  urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))

    allocate(  uly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(  ury(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))

    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2)))

    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1))

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    call slopex_2d(utilde,slopex,domlo,domhi,lo,hi,ng_ut,2,adv_bc)
    call slopey_2d(utilde,slopey,domlo,domhi,lo,hi,ng_ut,2,adv_bc)
       
    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

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
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,1)")
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
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,2)")
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
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,1)")
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
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,2)")
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
          fl = force(i-1,j,1)
          fr = force(i,j  ,1)

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
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,1)")
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
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,2)")
    end select
    end if


    do j=js,je+1
       do i=is,ie
          fl = force(i,j-1,2)
          fr = force(i,j  ,2)
          
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
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,1)")
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
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,2)")
    end select
    end if

    deallocate(ulx,urx,uimhx,uly,ury,uimhy,umacl,umacr,vmacl,vmacr)

  end subroutine velpred_2d
#endif

#if (AMREX_SPACEDIM == 3)
  subroutine velpred_3d(lev, domlo, domhi, lo, hi, &
                        utilde, ut_lo, ut_hi, nc_ut, ng_ut, &
                        ufull,  uf_lo, uf_hi, nc_uf, &
                        utrans, uu_lo, uu_hi, &
                        vtrans, uv_lo, uv_hi, &
                        wtrans, uw_lo, uw_hi, &
                        umac,   mu_lo, mu_hi, &
                        vmac,   mv_lo, mv_hi, &
                        wmac,   mw_lo, mw_hi, &
                        force,   f_lo,  f_hi, nc_f, &
                        w0,dx,dt,adv_bc,phys_bc) bind(C,name="velpred_3d")

    integer         , intent(in   ) :: lev, domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: ut_lo(3), ut_hi(3), nc_ut, ng_ut
    integer         , intent(in   ) :: uf_lo(3), uf_hi(3), nc_uf
    integer         , intent(in   ) :: uu_lo(3), uu_hi(3)
    integer         , intent(in   ) :: uv_lo(3), uv_hi(3)
    integer         , intent(in   ) :: uw_lo(3), uw_hi(3)
    integer         , intent(in   ) :: mu_lo(3), mu_hi(3)
    integer         , intent(in   ) :: mv_lo(3), mv_hi(3)
    integer         , intent(in   ) :: mw_lo(3), mw_hi(3)
    integer         , intent(in   ) ::  f_lo(3),  f_hi(3), nc_f
    double precision, intent(in   ) :: utilde(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),nc_ut)
    double precision, intent(in   ) :: ufull (uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),uf_lo(3):uf_hi(3),nc_uf)
    double precision, intent(inout) :: utrans(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3))
    double precision, intent(inout) :: vtrans(uv_lo(1):uv_hi(1),uv_lo(2):uv_hi(2),uv_lo(3):uv_hi(3))
    double precision, intent(inout) :: wtrans(uw_lo(1):uw_hi(1),uw_lo(2):uw_hi(2),uw_lo(3):uw_hi(3))
    double precision, intent(inout) :: umac  (mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))
    double precision, intent(inout) :: vmac  (mv_lo(1):mv_hi(1),mv_lo(2):mv_hi(2),mv_lo(3):mv_hi(3))
    double precision, intent(inout) :: wmac  (mw_lo(1):mw_hi(1),mw_lo(2):mw_hi(2),mw_lo(3):mw_hi(3))
    double precision, intent(in   ) :: force ( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3),nc_f)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3), dt
    integer         , intent(in   ) :: adv_bc(3,2,3), phys_bc(3,2) ! dim, lohi, (comp)

    ! local variables
    double precision, allocatable :: slopex(:,:,:,:)
    double precision, allocatable :: slopey(:,:,:,:)
    double precision, allocatable :: slopez(:,:,:,:)

    ! these correspond to u_L^x, etc.
    double precision, allocatable:: ulx(:,:,:,:),urx(:,:,:,:),uimhx(:,:,:,:)
    double precision, allocatable:: uly(:,:,:,:),ury(:,:,:,:),uimhy(:,:,:,:)
    double precision, allocatable:: ulz(:,:,:,:),urz(:,:,:,:),uimhz(:,:,:,:)

    ! these correspond to u_L^{y|z}, etc.
    double precision, allocatable:: ulyz(:,:,:)
    double precision, allocatable:: uryz(:,:,:)
    double precision, allocatable:: uimhyz(:,:,:)

    double precision, allocatable:: ulzy(:,:,:)
    double precision, allocatable:: urzy(:,:,:)
    double precision, allocatable:: uimhzy(:,:,:)

    double precision, allocatable:: vlxz(:,:,:)
    double precision, allocatable:: vrxz(:,:,:)
    double precision, allocatable:: vimhxz(:,:,:)

    double precision, allocatable:: vlzx(:,:,:)
    double precision, allocatable:: vrzx(:,:,:)
    double precision, allocatable:: vimhzx(:,:,:)

    double precision, allocatable:: wlxy(:,:,:)
    double precision, allocatable:: wrxy(:,:,:)
    double precision, allocatable:: wimhxy(:,:,:)

    double precision, allocatable:: wlyx(:,:,:)
    double precision, allocatable:: wryx(:,:,:)
    double precision, allocatable:: wimhyx(:,:,:)

    ! these correspond to umac_L, etc.
    double precision, allocatable:: umacl(:,:,:),umacr(:,:,:)
    double precision, allocatable:: vmacl(:,:,:),vmacr(:,:,:)
    double precision, allocatable:: wmacl(:,:,:),wmacr(:,:,:)

    double precision :: hx, hy, hz, dt2, dt4, dt6, uavg, maxu, minu
    double precision :: fl, fr

    integer :: i,j,k,is,js,ie,je,ks,ke

    logical :: test

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

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

    do k = lo(3)-1,hi(3)+1
       call slopex_2d(utilde(:,:,k,:),slopex(:,:,k,:),domlo,domhi,lo,hi,ng_ut,3,adv_bc)
       call slopey_2d(utilde(:,:,k,:),slopey(:,:,k,:),domlo,domhi,lo,hi,ng_ut,3,adv_bc)
    end do
    call slopez_3d(utilde,slopez,domlo,domhi,lo,hi,ng_ut,3,adv_bc)

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

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

    deallocate(slopex)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select
    end if

    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

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

    ! normal predictor states

    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(uly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(ury(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))

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

    deallocate(slopey)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select
    end if

    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))

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

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(ulz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))
    allocate(urz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))

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

    deallocate(slopez)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select
    end if

    allocate(uimhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))

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

    !******************************************************************
    ! Create u_{\i-\half\e_y}^{y|z}, etc.
    !******************************************************************

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(ulyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uryz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uimhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

    ! uimhyz loop
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select
    end if

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

    deallocate(ulyz,uryz)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(ulzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(urzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(uimhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

    ! uimhzy loop
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select
    end if

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

    deallocate(ulzy,urzy)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(vlxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(vrxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    ! vimhxz loop
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

    deallocate(uimhz)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select
    end if

    allocate(vimhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

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

    deallocate(vlxz,vrxz)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(vlzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vrzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vimhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    ! vimhzx loop
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select
    end if

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

    deallocate(vlzx,vrzx)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(wlxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(wrxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    ! wimhxy loop
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

    deallocate(uimhy)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select
    end if

    allocate(wimhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

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

    deallocate(wlxy,wrxy)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(wlyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wryx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    ! wimhyx loop
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

    deallocate(uimhx)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select
    end if

    allocate(wimhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

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

    deallocate(wlyx,wryx)

    !******************************************************************
    ! Create umac, etc.
    !******************************************************************

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))

    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             fl = force(i-1,j,k,1)
             fr = force(i,j  ,k,1)
             
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

    deallocate(ulx,urx,uimhyz,uimhzy)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select
    end if

    deallocate(umacl,umacr)

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))

    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             fl = force(i,j-1,k,2)
             fr = force(i,j  ,k,2)

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

    deallocate(uly,ury,vimhxz,vimhzx)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select
    end if

    deallocate(vmacl,vmacr)

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(wmacl(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
    allocate(wmacr(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             fl = force(i,j,k-1,3)
             fr = force(i,j,k  ,3)

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

    deallocate(ulz,urz,wimhxy,wimhyx)

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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
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
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select
    end if

    deallocate(wmacl,wmacr)

  end subroutine velpred_3d
#endif

end module velpred_module
