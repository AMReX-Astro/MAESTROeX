
#include "AMReX_BC_TYPES.H"

module slope_module

  implicit none

  integer, private, parameter :: cen = 1, lim = 2, flag = 3, fromm = 4

contains

  subroutine slopex_1d(s,slx,domlo,domhi,lo,hi,ng,nvar,adv_bc)

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in   ) :: domlo(:),domhi(:),lo(:),hi(:),ng,nvar
    double precision, intent(in   ) ::   s(lo(1)-ng:,:)
    double precision, intent(  out) :: slx(lo(1)- 1:,:)
    integer         , intent(in)    :: adv_bc(:,:,:)

    ! Local variables
    integer :: is,ie
    integer :: i,comp

    double precision :: del,slim,sflag,dpls,dmin,ds
    double precision :: dxscr(lo(1)-2:hi(1)+2,4)

    is = lo(1)
    ie = hi(1)
    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       slx(lo(1)-1:hi(1)+1, :) = zero


    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do i = is-1,ie+1
             del = half*(s(i+1,comp) - s(i-1,comp))
             dpls = two*(s(i+1,comp) - s(i  ,comp))
             dmin = two*(s(i  ,comp) - s(i-1,comp))
             slim = min(abs(dpls), abs(dmin))
             slim = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(one,del)
             slx(i,comp)= sflag*min(slim,abs(del))
          enddo

          if (lo(1) .eq. domlo(1)) then
             if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                slx(is-1,comp) = zero
                del = (s(is+1,comp)+three*s(is,comp)- &
                     four*s(is-1,comp) ) * third
                dpls = two*(s(is+1,comp) - s(is  ,comp))
                dmin = two*(s(is  ,comp) - s(is-1,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(is,comp)= sflag*min(slim,abs(del))
             endif
          endif

          if (hi(1) .eq. domhi(1)) then
             if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                slx(ie+1,comp) = zero
                del = -(s(ie-1,comp)+three*s(ie,comp)- &
                     four*s(ie+1,comp) ) * third
                dpls = two*(s(ie  ,comp) - s(ie-1,comp))
                dmin = two*(s(ie+1,comp) - s(ie  ,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(ie,comp)= sflag*min(slim,abs(del))
             endif
          endif

       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar

          do i = is-2,ie+2
             dxscr(i,cen) = half*(s(i+1,comp)-s(i-1,comp))
             dmin = two*(s(i  ,comp)-s(i-1,comp))
             dpls = two*(s(i+1,comp)-s(i  ,comp))
             dxscr(i,lim)= min(abs(dmin),abs(dpls))
             dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
             dxscr(i,flag) = sign(one,dxscr(i,cen))
             dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim), &
                  abs(dxscr(i,cen)))
          enddo

          do i = is-1,ie+1
             ds = two * two3rd * dxscr(i,cen) - &
                  sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
             slx(i,comp) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
          enddo

          if (lo(1) .eq. domlo(1)) then
             if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                slx(is-1,comp) = zero
                del = -sixteen/fifteen*s(is-1,comp) + half*s(is,comp) + &
                     two3rd*s(is+1,comp) - tenth*s(is+2,comp)
                dmin = two*(s(is  ,comp)-s(is-1,comp))
                dpls = two*(s(is+1,comp)-s(is  ,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(is,comp)= sflag*min(slim,abs(del))

                ! Recalculate the slope at is+1 using the revised dxscr(is,fromm)
                dxscr(is,fromm) = slx(is,comp)
                ds = two * two3rd * dxscr(is+1,cen) - &
                     sixth * (dxscr(is+2,fromm) + dxscr(is,fromm))
                slx(is+1,comp) = dxscr(is+1,flag)*min(abs(ds),dxscr(is+1,lim))
             endif
          endif

          if (hi(1) .eq. domhi(1)) then
             if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                slx(ie+1,comp) = zero
                del = -( -sixteen/fifteen*s(ie+1,comp) + half*s(ie,comp) +  &
                     two3rd*s(ie-1,comp) - tenth*s(ie-2,comp) )
                dmin = two*(s(ie  ,comp)-s(ie-1,comp))
                dpls = two*(s(ie+1,comp)-s(ie  ,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(ie,comp)= sflag*min(slim,abs(del))

                ! Recalculate the slope at ie-1 using the revised dxscr(ie,fromm)
                dxscr(ie,fromm) = slx(ie,comp)
                ds = two * two3rd * dxscr(ie-1,cen) - &
                     sixth * (dxscr(ie-2,fromm) + dxscr(ie,fromm))
                slx(ie-1,comp) = dxscr(ie-1,flag)*min(abs(ds),dxscr(ie-1,lim))
             endif
          endif

       enddo

    endif

  end subroutine slopex_1d


  subroutine slopex_2d(s,slx,domlo,domhi,lo,hi,ng,nvar,adv_bc)

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in   ) :: domlo(:),domhi(:),lo(:),hi(:),ng,nvar
    double precision, intent(in   ) ::   s(lo(1)-ng+1:, lo(2)-ng+1:,:)
    double precision, intent(  out) :: slx(lo(1):, lo(2):,:)
    integer         , intent(in)    :: adv_bc(:,:,:)

    ! Local variables
    integer :: i,j,comp

    double precision :: del,slim,sflag,dpls,dmin,ds
    double precision :: dcen, dlim, dflag, dxl, dxr

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                slx(i,j,comp) = ZERO
             enddo
          enddo
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                del = half*(s(i+1,j,comp) - s(i-1,j,comp))
                dpls = two*(s(i+1,j,comp) - s(i  ,j,comp))
                dmin = two*(s(i  ,j,comp) - s(i-1,j,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(i,j,comp)= sflag*min(slim,abs(del))
             enddo


             if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                if (i .eq. lo(1) .and. lo(1)+1 .eq. domlo(1)) then
                   slx(i,j,comp) = ZERO
               elseif (i .eq. lo(1)+1 .and. lo(1)+1 .eq. domlo(1)) then
                   del = (s(i+1,j,comp)+three*s(i,j,comp)- &
                        four*s(i-1,j,comp) ) * third
                   dpls = two*(s(i+1,j,comp) - s(i,j,comp))
                   dmin = two*(s(i,j,comp) - s(i-1,j,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(i,j,comp)= sflag*min(slim,abs(del))
                endif
             endif


             if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                   slx(i,j,comp) = ZERO
               elseif (i .eq. hi(1)-1 .and. hi(1)-1 .eq. domhi(1)) then
                   del = -(s(i-1,j,comp)+three*s(i,j,comp)- &
                        four*s(i+1,j,comp) ) * third
                   dpls = two*(s(i,j,comp) - s(i-1,j,comp))
                   dmin = two*(s(i+1,j,comp) - s(i,j,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(i,j,comp)= sflag*min(slim,abs(del))
                endif
             endif

          enddo
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ! left
                dcen = half*(s(i,j,comp)-s(i-2,j,comp))
                dmin = two*(s(i-1,j,comp)-s(i-2,j,comp))
                dpls = two*(s(i,j,comp)-s(i-1,j,comp))
                dlim = min(abs(dmin),abs(dpls))
                dlim = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                dflag = sign(one,dcen)
                dxl = dflag*min(dlim, abs(dcen))

                ! right
                dcen = half*(s(i+2,j,comp)-s(i,j,comp))
                dmin = two*(s(i+1,j,comp)-s(i,j,comp))
                dpls = two*(s(i+2,j,comp)-s(i+1,j,comp))
                dlim = min(abs(dmin),abs(dpls))
                dlim = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                dflag = sign(one,dcen)
                dxr = dflag*min(dlim, abs(dcen))

                ! center
                dcen = half*(s(i+1,j,comp)-s(i-1,j,comp))
                dmin = two*(s(i  ,j,comp)-s(i-1,j,comp))
                dpls = two*(s(i+1,j,comp)-s(i  ,j,comp))
                dlim = min(abs(dmin),abs(dpls))
                dlim = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                dflag = sign(one,dcen)

                ds = two * two3rd * dcen - sixth * (dxr + dxl)
                slx(i,j,comp) = dflag*min(abs(ds),dlim)


                if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                   if (i .eq. lo(1) .and. lo(1)+1 .eq. domlo(1)) then
                      slx(i,j,comp) = ZERO
                  elseif (i .eq. lo(1)+1 .and. lo(1)+1 .eq. domlo(1)) then
                      del = -sixteen/fifteen*s(i-1,j,comp) + half*s(i,j,comp) + &
                           two3rd*s(i+1,j,comp) - tenth*s(i+2,j,comp)
                      dmin = two*(s(i,j,comp)-s(i-1,j,comp))
                      dpls = two*(s(i+1,j,comp)-s(i,j,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slx(i,j,comp)= sflag*min(slim,abs(del))
                  elseif (i .eq. lo(1)+2 .and. lo(1)+1 .eq. domlo(1)) then
                      ! Recalculate the slope at lo(1)+1 using the revised dxl
                      del = -sixteen/fifteen*s(i-2,j,comp) + half*s(i-1,j,comp) + &
                           two3rd*s(i,j,comp) - tenth*s(i+1,j,comp)
                      dmin = two*(s(i-1,j,comp)-s(i-2,j,comp))
                      dpls = two*(s(i,j,comp)-s(i-1,j,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      dxl = sflag*min(slim,abs(del))

                      ds = two * two3rd * dcen - sixth * (dxr + dxl)
                      slx(i,j,comp) = dflag*min(abs(ds),dlim)
                   endif
                endif


                if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                   if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                      slx(i,j,comp) = ZERO
                  elseif (i .eq. hi(1)-1 .and. hi(1)-1 .eq. domhi(1)) then
                      del = -( -sixteen/fifteen*s(i+1,j,comp) + half*s(i,j,comp) +  &
                           two3rd*s(i-1,j,comp) - tenth*s(i-2,j,comp) )
                      dmin = two*(s(i,j,comp)-s(i-1,j,comp))
                      dpls = two*(s(i+1,j,comp)-s(i,j,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slx(i,j,comp)= sflag*min(slim,abs(del))
                  elseif (i .eq. hi(1)-2 .and. hi(1)-1 .eq. domhi(1)) then
                      ! Recalculate the slope at hi(1)-1 using the revised dxr
                      del = -( -sixteen/fifteen*s(i+2,j,comp) + half*s(i+1,j,comp) +  &
                           two3rd*s(i,j,comp) - tenth*s(i-1,j,comp) )
                      dmin = two*(s(i+1,j,comp)-s(i,j,comp))
                      dpls = two*(s(i+2,j,comp)-s(i+1,j,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      dxr = sflag*min(slim,abs(del))

                      ds = two * two3rd * dcen - sixth * (dxl + dxr)
                      slx(i,j,comp) = dflag*min(abs(ds),dlim)
                   endif
                endif
             enddo

          enddo
       enddo

    endif

  end subroutine slopex_2d

  subroutine slopey_2d(s,sly,domlo,domhi,lo,hi,ng,nvar,adv_bc)

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in)  :: domlo(:),domhi(:),lo(:),hi(:),ng,nvar
    integer         , intent(in)  :: adv_bc(:,:,:)
    double precision, intent( in) ::   s(lo(1)-ng+1:,lo(2)-ng+1:,:)
    double precision, intent(out) :: sly(lo(1):,lo(2):,:)

    ! local
    double precision :: dpls,dmin,ds,del,slim,sflag
    double precision :: dcen,dlim,dflag,dyl,dyr

    integer :: i,j,comp

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                sly(i,j,comp) = ZERO
             enddo
          enddo
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                del  = half*(s(i,j+1,comp) - s(i,j-1,comp))
                dpls = two *(s(i,j+1,comp) - s(i,j  ,comp))
                dmin = two *(s(i,j  ,comp) - s(i,j-1,comp))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                sly(i,j,comp)= sflag*min(slim,abs(del))


                if (adv_bc(2,1,comp) .eq. EXT_DIR .or. adv_bc(2,1,comp) .eq. HOEXTRAP) then
                   if (j .eq. lo(2) .and. lo(2)+1 .eq. domlo(2)) then
                      sly(i,j,comp) = ZERO
                  elseif (j .eq. lo(2)+1 .and. lo(2)+1 .eq. domlo(2)) then
                      del = (s(i,j+1,comp)+three*s(i,j,comp)- &
                           four*s(i,j-1,comp)) * third
                      dpls = two*(s(i,j+1,comp) - s(i,j,comp))
                      dmin = two*(s(i,j,comp) - s(i,j-1,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      sly(i,j,comp)= sflag*min(slim,abs(del))
                   endif
                endif


                if (adv_bc(2,2,comp) .eq. EXT_DIR .or. adv_bc(2,2,comp) .eq. HOEXTRAP) then
                   if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                      sly(i,j,comp) = ZERO
                  elseif (j .eq. hi(2)-1 .and. hi(2)-1 .eq. domhi(2)) then
                      del = -(s(i,j-1,comp)+three*s(i,j,comp)- &
                           four*s(i,j+1,comp)) * third
                      dpls = two*(s(i,j+1,comp) - s(i,j,comp))
                      dmin = two*(s(i,j,comp) - s(i,j-1,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      sly(i,j,comp)= sflag*min(slim,abs(del))

                   endif
                endif
             enddo
          enddo
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          do i = lo(1),hi(1)
             do j = lo(2),hi(2)
                 ! left
                dcen = half*(s(i,j,comp)-s(i,j-2,comp))
                dmin = two*(s(i,j-1,comp)-s(i,j-2,comp))
                dpls = two*(s(i,j,comp)-s(i,j-1,comp))
                dlim  = min(abs(dmin),abs(dpls))
                dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                dflag = sign(one,dcen)
                dyl = dflag*min(dlim,abs(dcen))

                ! right
                dcen = half*(s(i,j+2,comp)-s(i,j,comp))
                dmin = two*(s(i,j+1,comp)-s(i,j,comp))
                dpls = two*(s(i,j+2,comp)-s(i,j+1,comp))
                dlim  = min(abs(dmin),abs(dpls))
                dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                dflag = sign(one,dcen)
                dyr = dflag*min(dlim,abs(dcen))

                ! center
                dcen = half*(s(i,j+1,comp)-s(i,j-1,comp))
                dmin = two*(s(i,j  ,comp)-s(i,j-1,comp))
                dpls = two*(s(i,j+1,comp)-s(i,j  ,comp))
                dlim  = min(abs(dmin),abs(dpls))
                dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                dflag = sign(one,dcen)

                ds = two * two3rd * dcen - sixth * (dyr + dyl)
                sly(i,j,comp) = dflag*min(abs(ds),dlim)


                if (adv_bc(2,1,comp) .eq. EXT_DIR .or. adv_bc(2,1,comp) .eq. HOEXTRAP) then
                   if (j .eq. lo(2) .and. lo(2)+1 .eq. domlo(2)) then
                      sly(i,j,comp) = ZERO
                  elseif (j .eq. lo(2)+1 .and. lo(2)+1 .eq. domlo(2)) then
                      del = -sixteen/fifteen*s(i,j-1,comp) +  half*s(i,j,comp) +  &
                           two3rd*s(i,j+1,comp) - tenth*s(i,j+2,comp)
                      dmin = two*(s(i,j,comp)-s(i,j-1,comp))
                      dpls = two*(s(i,j+1,comp)-s(i,j,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      sly(i,j,comp)= sflag*min(slim,abs(del))
                  elseif (j .eq. lo(2)+2 .and. lo(2)+1 .eq. domlo(2)) then

                      ! Recalculate the slope at lo(2)+1 using the revised dyl
                      del = -sixteen/fifteen*s(i,j-2,comp) +  half*s(i,j-1,comp) +  &
                           two3rd*s(i,j,comp) - tenth*s(i,j+1,comp)
                      dmin = two*(s(i,j-1,comp)-s(i,j-2,comp))
                      dpls = two*(s(i,j,comp)-s(i,j-1,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      dyl = sflag*min(slim,abs(del))
                      ds = two * two3rd * dcen - sixth * (dyr + dyl)
                      sly(i,j,comp) = dflag*min(abs(ds),dlim)
                   endif
                endif


                if (adv_bc(2,2,comp) .eq. EXT_DIR .or. adv_bc(2,2,comp) .eq. HOEXTRAP) then
                   if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                      sly(i,j,comp) = ZERO
                  elseif (j .eq. hi(2)-1 .and. hi(2)-1 .eq. domhi(2)) then
                      del = -( -sixteen/fifteen*s(i,j+1,comp) +  half*s(i,j,comp) + &
                           two3rd*s(i,j-1,comp) - tenth*s(i,j-2,comp) )
                      dmin = two*(s(i,j,comp)-s(i,j-1,comp))
                      dpls = two*(s(i,j+1,comp)-s(i,j,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      sly(i,j,comp)= sflag*min(slim,abs(del))
                  elseif (j .eq. hi(2)-2 .and. hi(2)-1 .eq. domhi(2)) then
                      ! Recalculate the slope at lo(2)+1 using the revised dyr
                      del = -( -sixteen/fifteen*s(i,j+2,comp) +  half*s(i,j+1,comp) + &
                           two3rd*s(i,j,comp) - tenth*s(i,j-1,comp) )
                      dmin = two*(s(i,j+1,comp)-s(i,j,comp))
                      dpls = two*(s(i,j+2,comp)-s(i,j+1,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      dyr = sflag*min(slim,abs(del))
                      ds = two * two3rd * dcen - sixth * (dyl + dyr)
                      sly(i,j,comp) = dflag*min(abs(ds),dlim)
                   endif
                endif

             enddo
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
