
#include "AMReX_BC_TYPES.H"

module slope_module

  implicit none

  integer, private, parameter :: cen = 1, lim = 2, flag = 3, fromm = 4

contains

  subroutine slopex_2d(s,slx,domlo,domhi,lo,hi,ng,nvar,adv_bc)

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in   ) :: domlo(:),domhi(:),lo(:),hi(:),ng,nvar
    double precision, intent(in   ) ::   s(lo(1)-ng:, lo(2)-ng:,:)
    double precision, intent(  out) :: slx(lo(1)- 1:, lo(2)- 1:,:)
    integer         , intent(in)    :: adv_bc(:,:,:)

    ! Local variables
    integer :: is,js,ie,je
    integer :: i,j,comp

    double precision :: del,slim,sflag,dpls,dmin,ds
    double precision :: dxscr(lo(1)-2:hi(1)+2,4)

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       slx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,:) = zero

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do j = js-1,je+1
             do i = is-1,ie+1
                del = half*(s(i+1,j,comp) - s(i-1,j,comp))
                dpls = two*(s(i+1,j,comp) - s(i  ,j,comp))
                dmin = two*(s(i  ,j,comp) - s(i-1,j,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(i,j,comp)= sflag*min(slim,abs(del))
             enddo

             if (lo(1) .eq. domlo(1)) then
                if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                   slx(is-1,j,comp) = zero
                   del = (s(is+1,j,comp)+three*s(is,j,comp)- &
                        four*s(is-1,j,comp) ) * third
                   dpls = two*(s(is+1,j,comp) - s(is  ,j,comp))
                   dmin = two*(s(is  ,j,comp) - s(is-1,j,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(is,j,comp)= sflag*min(slim,abs(del))
                endif
             endif

             if (hi(1) .eq. domhi(1)) then
                if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                   slx(ie+1,j,comp) = zero
                   del = -(s(ie-1,j,comp)+three*s(ie,j,comp)- &
                        four*s(ie+1,j,comp) ) * third
                   dpls = two*(s(ie  ,j,comp) - s(ie-1,j,comp))
                   dmin = two*(s(ie+1,j,comp) - s(ie  ,j,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(ie,j,comp)= sflag*min(slim,abs(del))
                endif
             endif

          enddo
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          do j = js-1,je+1

             do i = is-2,ie+2
                dxscr(i,cen) = half*(s(i+1,j,comp)-s(i-1,j,comp))
                dmin = two*(s(i  ,j,comp)-s(i-1,j,comp))
                dpls = two*(s(i+1,j,comp)-s(i  ,j,comp))
                dxscr(i,lim)= min(abs(dmin),abs(dpls))
                dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
                dxscr(i,flag) = sign(one,dxscr(i,cen))
                dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim), &
                     abs(dxscr(i,cen)))
             enddo

             do i = is-1,ie+1
                ds = two * two3rd * dxscr(i,cen) - &
                     sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
                slx(i,j,comp) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
             enddo

             if (lo(1) .eq. domlo(1)) then
                if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                   slx(is-1,j,comp) = zero
                   del = -sixteen/fifteen*s(is-1,j,comp) + half*s(is,j,comp) + &
                        two3rd*s(is+1,j,comp) - tenth*s(is+2,j,comp)
                   dmin = two*(s(is  ,j,comp)-s(is-1,j,comp))
                   dpls = two*(s(is+1,j,comp)-s(is  ,j,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(is,j,comp)= sflag*min(slim,abs(del))

                   ! Recalculate the slope at is+1 using the revised dxscr(is,fromm)
                   dxscr(is,fromm) = slx(is,j,comp)
                   ds = two * two3rd * dxscr(is+1,cen) - &
                        sixth * (dxscr(is+2,fromm) + dxscr(is,fromm))
                   slx(is+1,j,comp) = dxscr(is+1,flag)*min(abs(ds),dxscr(is+1,lim))
                endif
             endif

             if (hi(1) .eq. domhi(1)) then
                if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                   slx(ie+1,j,comp) = zero
                   del = -( -sixteen/fifteen*s(ie+1,j,comp) + half*s(ie,j,comp) +  &
                        two3rd*s(ie-1,j,comp) - tenth*s(ie-2,j,comp) )
                   dmin = two*(s(ie  ,j,comp)-s(ie-1,j,comp))
                   dpls = two*(s(ie+1,j,comp)-s(ie  ,j,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(ie,j,comp)= sflag*min(slim,abs(del))

                   ! Recalculate the slope at ie-1 using the revised dxscr(ie,fromm)
                   dxscr(ie,fromm) = slx(ie,j,comp)
                   ds = two * two3rd * dxscr(ie-1,cen) - &
                        sixth * (dxscr(ie-2,fromm) + dxscr(ie,fromm))
                   slx(ie-1,j,comp) = dxscr(ie-1,flag)*min(abs(ds),dxscr(ie-1,lim))
                endif
             endif

          enddo
       enddo

    endif

  end subroutine slopex_2d

  subroutine slopey_2d(s,sly,domlo,domhi,lo,hi,ng,nvar,adv_bc)

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in)  :: domlo(:),domhi(:),lo(:),hi(:),ng,nvar
    integer         , intent(in)  :: adv_bc(:,:,:)
    double precision, intent( in) ::   s(lo(1)-ng:,lo(2)-ng:,:)
    double precision, intent(out) :: sly(lo(1)- 1:,lo(2)- 1:,:)

    ! local
    double precision :: dyscr(lo(2)-2:hi(2)+2,4)
    double precision :: dpls,dmin,ds,del,slim,sflag

    integer :: is,js,ie,je,i,j,comp

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)


    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       sly(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,:) = zero

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do j = js-1,je+1
             do i = is-1,ie+1

                del  = half*(s(i,j+1,comp) - s(i,j-1,comp))
                dpls = two *(s(i,j+1,comp) - s(i,j  ,comp))
                dmin = two *(s(i,j  ,comp) - s(i,j-1,comp))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                sly(i,j,comp)= sflag*min(slim,abs(del))

             enddo
          enddo

          if (lo(2) .eq. domlo(2)) then
             if (adv_bc(2,1,comp) .eq. EXT_DIR .or. adv_bc(2,1,comp) .eq. HOEXTRAP) then
                do i = is-1,ie+1
                   sly(i,js-1,comp) = zero
                   del = (s(i,js+1,comp)+three*s(i,js,comp)- &
                        four*s(i,js-1,comp)) * third
                   dpls = two*(s(i,js+1,comp) - s(i,js  ,comp))
                   dmin = two*(s(i,js  ,comp) - s(i,js-1,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   sly(i,js,comp)= sflag*min(slim,abs(del))
                enddo
             endif
          endif

          if (hi(2) .eq. domhi(2)) then
             if (adv_bc(2,2,comp) .eq. EXT_DIR .or. adv_bc(2,2,comp) .eq. HOEXTRAP) then
                do i = is-1, ie+1
                   sly(i,je+1,comp) = zero
                   del = -(s(i,je-1,comp)+three*s(i,je,comp)- &
                        four*s(i,je+1,comp)) * third
                   dpls = two*(s(i,je+1,comp) - s(i,je ,comp))
                   dmin = two*(s(i,je  ,comp) - s(i,je-1,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   sly(i,je,comp)= sflag*min(slim,abs(del))
                enddo
             endif
          endif

       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          do i = is-1,ie+1
             do j = js-2,je+2
                dyscr(j,cen) = half*(s(i,j+1,comp)-s(i,j-1,comp))
                dmin = two*(s(i,j  ,comp)-s(i,j-1,comp))
                dpls = two*(s(i,j+1,comp)-s(i,j  ,comp))
                dyscr(j,lim)  = min(abs(dmin),abs(dpls))
                dyscr(j,lim)  = merge(dyscr(j,lim),zero,dpls*dmin.gt.ZERO)
                dyscr(j,flag) = sign(one,dyscr(j,cen))
                dyscr(j,fromm)= dyscr(j,flag)*min(dyscr(j,lim),abs(dyscr(j,cen)))
             enddo

             do j = js-1,je+1
                ds = two * two3rd * dyscr(j,cen) -  &
                     sixth * (dyscr(j+1,fromm) + dyscr(j-1,fromm))
                sly(i,j,comp) = dyscr(j,flag)*min(abs(ds),dyscr(j,lim))
             enddo

             if (lo(2) .eq. domlo(2)) then
                if (adv_bc(2,1,comp) .eq. EXT_DIR .or. adv_bc(2,1,comp) .eq. HOEXTRAP) then
                   sly(i,js-1,comp) = zero
                   del = -sixteen/fifteen*s(i,js-1,comp) +  half*s(i,js ,comp) +  &
                        two3rd*s(i,js+1,comp) - tenth*s(i,js+2,comp)
                   dmin = two*(s(i,js  ,comp)-s(i,js-1,comp))
                   dpls = two*(s(i,js+1,comp)-s(i,js  ,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   sly(i,js,comp)= sflag*min(slim,abs(del))

                   ! Recalculate the slope at js+1 using the revised dyscr(js,fromm)
                   dyscr(js,fromm) = sly(i,js,comp)
                   ds = two * two3rd * dyscr(js+1,cen) - &
                        sixth * (dyscr(js+2,fromm) + dyscr(js,fromm))
                   sly(i,js+1,comp) = dyscr(js+1,flag)*min(abs(ds),dyscr(js+1,lim))
                endif
             endif

             if (hi(2) .eq. domhi(2)) then
                if (adv_bc(2,2,comp) .eq. EXT_DIR .or. adv_bc(2,2,comp) .eq. HOEXTRAP) then
                   sly(i,je+1,comp) = zero
                   del = -( -sixteen/fifteen*s(i,je+1,comp) +  half*s(i,je  ,comp) + &
                        two3rd*s(i,je-1,comp) - tenth*s(i,je-2,comp) )
                   dmin = two*(s(i,je ,comp)-s(i,je-1,comp))
                   dpls = two*(s(i,je+1,comp)-s(i,je ,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   sly(i,je,comp)= sflag*min(slim,abs(del))

                   ! Recalculate the slope at js+1 using the revised dyscr(js,fromm)
                   dyscr(je,fromm) = sly(i,je,comp)
                   ds = two * two3rd * dyscr(je-1,cen) -  &
                        sixth * (dyscr(je-2,fromm) + dyscr(je,fromm))
                   sly(i,je-1,comp) = dyscr(je-1,flag)*min(abs(ds),dyscr(je-1,lim))
                endif
             endif

          enddo
       enddo

    endif

  end subroutine slopey_2d

  subroutine slopez_3d(s,slz,domlo,domhi,lo,hi,ng,nvar,adv_bc)

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in)  :: domlo(:),domhi(:),lo(:),hi(:),ng,nvar
    integer         , intent(in)  :: adv_bc(:,:,:)
    double precision, intent( in) ::   s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    double precision, intent(out) :: slz(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)

    double precision :: dzscr(lo(3)-2:hi(3)+2,4)
    double precision :: dpls,dmin,ds,del,slim,sflag

    integer :: is,js,ks,ie,je,ke,i,j,k,comp

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       slz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:) = zero

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do k = ks-1,ke+1
             do j = js-1,je+1
                do i = is-1,ie+1

                   del  = half*(s(i,j,k+1,comp) - s(i,j,k-1,comp))
                   dpls = two *(s(i,j,k+1,comp) - s(i,j,k  ,comp))
                   dmin = two *(s(i,j,k  ,comp) - s(i,j,k-1,comp))
                   slim = min(abs(dpls),abs(dmin))
                   slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slz(i,j,k,comp)= sflag*min(slim,abs(del))

                enddo
             enddo
          enddo

          if (lo(3) .eq. domlo(3)) then
             if (adv_bc(3,1,comp) .eq. EXT_DIR .or. adv_bc(3,1,comp) .eq. HOEXTRAP) then
                do j = js-1,je+1
                   do i = is-1,ie+1
                      slz(i,j,ks-1,comp) = zero
                      del = (s(i,j,ks+1,comp)+three*s(i,j,ks,comp)- &
                           four*s(i,j,ks-1,comp)) * third
                      dpls = two*(s(i,j,ks+1,comp) - s(i,j,ks  ,comp))
                      dmin = two*(s(i,j,ks  ,comp) - s(i,j,ks-1,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slz(i,j,ks,comp)= sflag*min(slim,abs(del))
                   enddo
                enddo
             endif
          endif

          if (hi(3) .eq. domhi(3)) then
             if (adv_bc(3,2,comp) .eq. EXT_DIR .or. adv_bc(3,2,comp) .eq. HOEXTRAP) then
                do j = js-1, je+1
                   do i = is-1, ie+1
                      slz(i,j,ke+1,comp) = zero
                      del = -(s(i,j,ke-1,comp)+three*s(i,j,ke,comp)- &
                           four*s(i,j,ke+1,comp)) * third
                      dpls = two*(s(i,j,ke+1,comp) - s(i,j,ke ,comp))
                      dmin = two*(s(i,j,ke  ,comp) - s(i,j,ke-1,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slz(i,j,ke,comp)= sflag*min(slim,abs(del))
                   enddo
                enddo
             endif
          endif

       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          do j = js-1,je+1
             do i = is-1,ie+1
                do k = ks-2,ke+2
                   dzscr(k,cen) = half*(s(i,j,k+1,comp)-s(i,j,k-1,comp))
                   dmin = two*(s(i,j,k  ,comp)-s(i,j,k-1,comp))
                   dpls = two*(s(i,j,k+1,comp)-s(i,j,k  ,comp))
                   dzscr(k,lim)  = min(abs(dmin),abs(dpls))
                   dzscr(k,lim)  = merge(dzscr(k,lim),zero,dpls*dmin.gt.ZERO)
                   dzscr(k,flag) = sign(one,dzscr(k,cen))
                   dzscr(k,fromm)= dzscr(k,flag)*min(dzscr(k,lim),abs(dzscr(k,cen)))
                enddo

                do k = ks-1,ke+1
                   ds = two * two3rd * dzscr(k,cen) -  &
                        sixth * (dzscr(k+1,fromm) + dzscr(k-1,fromm))
                   slz(i,j,k,comp) = dzscr(k,flag)*min(abs(ds),dzscr(k,lim))
                enddo

                if (lo(3) .eq. domlo(3)) then
                   if (adv_bc(3,1,comp) .eq. EXT_DIR .or. adv_bc(3,1,comp) .eq. HOEXTRAP) then
                      slz(i,j,ks-1,comp) = zero
                      del = -sixteen/fifteen*s(i,j,ks-1,comp) +  half*s(i,j,ks ,comp) +  &
                           two3rd*s(i,j,ks+1,comp) - tenth*s(i,j,ks+2,comp)
                      dmin = two*(s(i,j,ks  ,comp)-s(i,j,ks-1,comp))
                      dpls = two*(s(i,j,ks+1,comp)-s(i,j,ks  ,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slz(i,j,ks,comp)= sflag*min(slim,abs(del))

                      ! Recalculate the slope at js+1 using the revised dzscr(js,fromm)
                      dzscr(ks,fromm) = slz(i,j,ks,comp)
                      ds = two * two3rd * dzscr(ks+1,cen) - &
                           sixth * (dzscr(ks+2,fromm) + dzscr(ks,fromm))
                      slz(i,j,ks+1,comp) = dzscr(ks+1,flag)*min(abs(ds),dzscr(ks+1,lim))
                   endif
                endif

                if (hi(3) .eq. domhi(3)) then
                   if (adv_bc(3,2,comp) .eq. EXT_DIR .or. adv_bc(3,2,comp) .eq. HOEXTRAP) then
                      slz(i,j,ke+1,comp) = zero
                      del = -( -sixteen/fifteen*s(i,j,ke+1,comp) +  half*s(i,j,ke  ,comp) + &
                           two3rd*s(i,j,ke-1,comp) - tenth*s(i,j,ke-2,comp) )
                      dmin = two*(s(i,j,ke ,comp)-s(i,j,ke-1,comp))
                      dpls = two*(s(i,j,ke+1,comp)-s(i,j,ke ,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slz(i,j,ke,comp)= sflag*min(slim,abs(del))

                      ! Recalculate the slope at ks+1 using the revised dzscr(ks,fromm)
                      dzscr(ke,fromm) = slz(i,j,ke,comp)
                      ds = two * two3rd * dzscr(ke-1,cen) -  &
                           sixth * (dzscr(ke-2,fromm) + dzscr(ke,fromm))
                      slz(i,j,ke-1,comp) = dzscr(ke-1,flag)*min(abs(ds),dzscr(ke-1,lim))
                   endif
                endif

             enddo
          enddo
       enddo

    endif

  end subroutine slopez_3d

end module slope_module
