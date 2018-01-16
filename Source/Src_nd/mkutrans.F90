
#include "AMReX_BC_TYPES.H"

module mkutrans_module

  use amrex_constants_module
  use slope_module
  use meth_params_module, only: rel_eps

  implicit none

  private

contains

  subroutine mkutrans_1d(utilde,ng_u,ufull,ng_uf,utrans,ng_ut,w0, &
                         lo,hi,dx,dt,adv_bc,phys_bc)

    integer         , intent(in   ) :: lo(:),hi(:),ng_u,ng_uf,ng_ut
    double precision, intent(in   ) :: utilde(lo(1)-ng_u :,:)
    double precision, intent(in   ) ::  ufull(lo(1)-ng_uf:,:)
    double precision, intent(inout) :: utrans(lo(1)-ng_ut:)
    double precision, intent(in   ) :: w0(0:)    
    double precision, intent(in   ) :: dt,dx(:)
    integer         , intent(in   ) :: adv_bc(:,:,:)
    integer         , intent(in   ) :: phys_bc(:,:)
    
    double precision :: slopex(lo(1)-1:hi(1)+1,1)

    double precision, allocatable :: ulx(:),urx(:)

    double precision hx,dt2,uavg

    integer :: i,is,ie

    logical :: test
    
    allocate(ulx(lo(1):hi(1)+1))
    allocate(urx(lo(1):hi(1)+1))

    is = lo(1)
    ie = hi(1)
    
    dt2 = HALF*dt
    
    hx = dx(1)
    
    call slopex_1d(utilde(:,1:),slopex,lo,hi,ng_u,1,adv_bc(:,:,1:))

    !******************************************************************
    ! create utrans
    !******************************************************************

    do i=is,ie+1
       ! extrapolate to edges
       ulx(i) = utilde(i-1,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,1)))*slopex(i-1,1)
       urx(i) = utilde(i  ,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,1)))*slopex(i  ,1)
    end do

    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (Inflow)
       ulx(is) = utilde(is-1,1)
       urx(is) = utilde(is-1,1)
    case (SlipWall, NoSlipWall, Symmetry)
       ulx(is) = ZERO
       urx(is) = ZERO
    case (Outflow)
       ulx(is) = min(urx(is),ZERO)
       urx(is) = ulx(is)
    case (Interior) 
    case  default
       call bl_error("mkutrans_1d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's    
    select case(phys_bc(1,2))
    case (Inflow)
       ulx(ie+1) = utilde(ie+1,1)
       urx(ie+1) = utilde(ie+1,1)
    case (SlipWall, NoSlipWall, Symmetry)
       ulx(ie+1) = ZERO
       urx(ie+1) = ZERO
    case (Outflow)
       ulx(ie+1) = max(ulx(ie+1),ZERO)
       urx(ie+1) = ulx(ie+1)
    case (Interior) 
    case  default
       call bl_error("mkutrans_1d: invalid boundary type phys_bc(1,2)")
    end select

    do i=is,ie+1
       ! solve Riemann problem using full velocity
       uavg = HALF*(ulx(i)+urx(i))
       test = ((ulx(i)+w0(i) .le. ZERO .and. urx(i)+w0(i) .ge. ZERO) .or. &
               (abs(ulx(i)+urx(i)+TWO*w0(i)) .lt. rel_eps))
       utrans(i) = merge(ulx(i),urx(i),uavg+w0(i) .gt. ZERO)
       utrans(i) = merge(ZERO,utrans(i),test)
    end do

    deallocate(ulx,urx)

  end subroutine mkutrans_1d

  subroutine mkutrans_2d(utilde,ng_u,ufull,ng_uf,utrans,vtrans,ng_ut,w0, &
                         lo,hi,dx,dt,adv_bc,phys_bc)

    integer         , intent(in   ) :: lo(:),hi(:),ng_u,ng_uf,ng_ut
    double precision, intent(in   ) :: utilde(lo(1)-ng_u :,lo(2)-ng_u :,:)
    double precision, intent(in   ) ::  ufull(lo(1)-ng_uf:,lo(2)-ng_uf:,:)
    double precision, intent(inout) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    double precision, intent(inout) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    double precision, intent(in   ) :: w0(0:)    
    double precision, intent(in   ) :: dt,dx(:)
    integer         , intent(in   ) :: adv_bc(:,:,:)
    integer         , intent(in   ) :: phys_bc(:,:)
    
    double precision :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    double precision :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)

    double precision, allocatable :: ulx(:,:),urx(:,:)
    double precision, allocatable :: vly(:,:),vry(:,:)

    double precision hx,hy,dt2,uavg

    integer :: i,j,is,js,ie,je

    logical :: test
    
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    dt2 = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    call slopex_2d(utilde(:,:,1:),slopex,lo,hi,ng_u,1,adv_bc(:,:,1:))
    call slopey_2d(utilde(:,:,2:),slopey,lo,hi,ng_u,1,adv_bc(:,:,2:))

    !******************************************************************
    ! create utrans
    !******************************************************************

    do j=js,je
       do i=is,ie+1
          ! extrapolate to edges
          ulx(i,j) = utilde(i-1,j,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,1)))*slopex(i-1,j,1)
          urx(i,j) = utilde(i  ,j,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,1)))*slopex(i  ,j,1)
       end do
    end do

    ! impose lo i side bc's
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
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's    
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
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(1,2)")
    end select

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

    do j=js,je+1
       do i=is,ie
          ! extrapolate to edges
          vly(i,j) = utilde(i,j-1,2) + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,2)))*slopey(i,j-1,1)
          vry(i,j) = utilde(i,j  ,2) - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,2)))*slopey(i,j  ,1)
       end do
    end do

    ! impose lo side bc's
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
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
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
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(2,2)")
    end select

    do j=js,je+1
       do i=is,ie
          ! solve Riemann problem using full velocity
          uavg = HALF*(vly(i,j)+vry(i,j))
          test = ((vly(i,j)+w0(j) .le. ZERO .and. vry(i,j)+w0(j) .ge. ZERO) .or. &
               (abs(vly(i,j)+vry(i,j)+TWO*w0(j)) .lt. rel_eps))
          vtrans(i,j) = merge(vly(i,j),vry(i,j),uavg+w0(j) .gt. ZERO)
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
       enddo
    enddo

    deallocate(ulx,urx,vly,vry)

  end subroutine mkutrans_2d
  
  subroutine mkutrans_3d(utilde,ng_u,ufull,ng_uf,utrans,vtrans,wtrans,ng_ut, &
                         w0,lo,hi,dx,dt,adv_bc,phys_bc)

    integer         , intent(in)    :: lo(:),hi(:),ng_u,ng_uf,ng_ut
    double precision, intent(in   ) :: utilde(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    double precision, intent(in   ) ::  ufull(lo(1)-ng_uf:,lo(2)-ng_uf:,lo(3)-ng_uf:,:)
    double precision, intent(inout) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    double precision, intent(inout) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    double precision, intent(inout) :: wtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    double precision, intent(in   ) :: w0(0:)
    double precision, intent(in   ) :: dt,dx(:)
    integer         , intent(in   ) :: adv_bc(:,:,:)
    integer         , intent(in   ) :: phys_bc(:,:)
    
    double precision, allocatable :: slopex(:,:,:,:)
    double precision, allocatable :: slopey(:,:,:,:)
    double precision, allocatable :: slopez(:,:,:,:)
    
    double precision hx,hy,hz,dt2,uavg
    
    logical :: test

    integer :: i,j,k,is,js,ks,ie,je,ke
    
    double precision, allocatable:: ulx(:,:,:),urx(:,:,:)
    double precision, allocatable:: vly(:,:,:),vry(:,:,:)
    double precision, allocatable:: wlz(:,:,:),wrz(:,:,:)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))

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
    
    do k = lo(3)-1,hi(3)+1
       call slopex_2d(utilde(:,:,k,1:),slopex(:,:,k,:),lo,hi,ng_u,1,adv_bc(:,:,1:))
       call slopey_2d(utilde(:,:,k,2:),slopey(:,:,k,:),lo,hi,ng_u,1,adv_bc(:,:,2:))
    end do
    call slopez_3d(utilde(:,:,:,3:),slopez,lo,hi,ng_u,1,adv_bc(:,:,3:))
    
    !******************************************************************
    ! create utrans
    !******************************************************************

    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

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
    deallocate(slopex)

    ! impose lo side bc's
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
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
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
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(1,2)")
    end select

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

    deallocate(ulx,urx)

    !******************************************************************
    ! create vtrans
    !******************************************************************

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

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

    deallocate(slopey)

    ! impose lo side bc's
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
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
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
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(2,2)")
    end select
    
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

    deallocate(vly,vry)

    !******************************************************************
    ! create wtrans
    !******************************************************************

    allocate(wlz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(wrz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

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

    deallocate(slopez)
    
    ! impose lo side bc's
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
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi side bc's
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
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(3,2)")
    end select
    
    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! solve Riemann problem using full velocity
             uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
             test = ((wlz(i,j,k)+w0(k).le.ZERO .and. wrz(i,j,k)+w0(k).ge.ZERO) .or. &
                  (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0(k)) .lt. rel_eps))
             wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0(k) .gt. ZERO)
             wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
          enddo
       enddo
    enddo

    deallocate(wlz,wrz)

  end subroutine mkutrans_3d
  
end module mkutrans_module
