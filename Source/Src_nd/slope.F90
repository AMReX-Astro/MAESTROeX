
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
       slx(lo(1)-1:hi(1)+1, :) = ZERO


    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do i = is-1,ie+1
             del = half*(s(i+1,comp) - s(i-1,comp))
             dpls = two*(s(i+1,comp) - s(i  ,comp))
             dmin = two*(s(i  ,comp) - s(i-1,comp))
             slim = min(abs(dpls), abs(dmin))
             slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
             sflag = sign(one,del)
             slx(i,comp)= sflag*min(slim,abs(del))
          enddo

          if (lo(1) .eq. domlo(1)) then
             if (adv_bc(1,1,comp) .eq. EXT_DIR  .or. adv_bc(1,1,comp) .eq. HOEXTRAP) then
                slx(is-1,comp) = ZERO
                del = (s(is+1,comp)+three*s(is,comp)- &
                     four*s(is-1,comp) ) * third
                dpls = two*(s(is+1,comp) - s(is  ,comp))
                dmin = two*(s(is  ,comp) - s(is-1,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                sflag = sign(one,del)
                slx(is,comp)= sflag*min(slim,abs(del))
             endif
          endif

          if (hi(1) .eq. domhi(1)) then
             if (adv_bc(1,2,comp) .eq. EXT_DIR  .or. adv_bc(1,2,comp) .eq. HOEXTRAP) then
                slx(ie+1,comp) = ZERO
                del = -(s(ie-1,comp)+three*s(ie,comp)- &
                     four*s(ie+1,comp) ) * third
                dpls = two*(s(ie  ,comp) - s(ie-1,comp))
                dmin = two*(s(ie+1,comp) - s(ie  ,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
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
             dxscr(i,lim) = merge(dxscr(i,lim),ZERO,dpls*dmin.gt.ZERO)
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
                slx(is-1,comp) = ZERO
                del = -sixteen/fifteen*s(is-1,comp) + half*s(is,comp) + &
                     two3rd*s(is+1,comp) - tenth*s(is+2,comp)
                dmin = two*(s(is  ,comp)-s(is-1,comp))
                dpls = two*(s(is+1,comp)-s(is  ,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
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
                slx(ie+1,comp) = ZERO
                del = -( -sixteen/fifteen*s(ie+1,comp) + half*s(ie,comp) +  &
                     two3rd*s(ie-1,comp) - tenth*s(ie-2,comp) )
                dmin = two*(s(ie  ,comp)-s(ie-1,comp))
                dpls = two*(s(ie+1,comp)-s(ie  ,comp))
                slim = min(abs(dpls), abs(dmin))
                slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
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


  subroutine slopex_2d(lo,hi,s,s_lo,s_hi,nc_s,slx,sl_lo,sl_hi,nc_sl, &
       domlo,domhi,nvar,adv_bc,nbccomp,bccomp) bind(C,name="slopex_2d")

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in   ) :: domlo(3),domhi(3),lo(3),hi(3)
    integer         , intent(in   ) :: s_lo(3),s_hi(3),sl_lo(3),sl_hi(3)
    integer  , value, intent(in   ) :: nc_s,nc_sl,nvar,nbccomp,bccomp
    double precision, intent(in   ) ::   s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(  out) :: slx(sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3),nc_sl)
    integer         , intent(in)    :: adv_bc(AMREX_SPACEDIM,2,nbccomp)

    ! Local variables
    integer :: i,j,k,comp,bc_comp

    double precision :: del,slim,sflag,dpls,dmin,ds
    double precision :: dcen, dlim, dflag, dxl, dxr

    !$gpu

    k = s_lo(3)

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
          do k = lo(3), hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   slx(i,j,k,comp) = ZERO
                enddo
             enddo
          enddo
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
          do k = lo(3), hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   del = half*(s(i+1,j,k,comp) - s(i-1,j,k,comp))
                   dpls = two*(s(i+1,j,k,comp) - s(i  ,j,k,comp))
                   dmin = two*(s(i  ,j,k,comp) - s(i-1,j,k,comp))
                   slim = min(abs(dpls), abs(dmin))
                   slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slx(i,j,k,comp)= sflag*min(slim,abs(del))
                enddo


                if (adv_bc(1,1,bc_comp) .eq. EXT_DIR  .or. adv_bc(1,1,bc_comp) .eq. HOEXTRAP) then
                   if (i .eq. lo(1) .and. lo(1)+1 .eq. domlo(1)) then
                      slx(i,j,k,comp) = ZERO
                   elseif (i .eq. lo(1)+1 .and. lo(1)+1 .eq. domlo(1)) then
                      del = (s(i+1,j,k,comp)+three*s(i,j,k,comp)- &
                           four*s(i-1,j,k,comp) ) * third
                      dpls = two*(s(i+1,j,k,comp) - s(i,j,k,comp))
                      dmin = two*(s(i,j,k,comp) - s(i-1,j,k,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slx(i,j,k,comp)= sflag*min(slim,abs(del))
                   endif
                endif


                if (adv_bc(1,2,bc_comp) .eq. EXT_DIR  .or. adv_bc(1,2,bc_comp) .eq. HOEXTRAP) then
                   if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                      slx(i,j,k,comp) = ZERO
                   elseif (i .eq. hi(1)-1 .and. hi(1)-1 .eq. domhi(1)) then
                      del = -(s(i-1,j,k,comp)+three*s(i,j,k,comp)- &
                           four*s(i+1,j,k,comp) ) * third
                      dpls = two*(s(i,j,k,comp) - s(i-1,j,k,comp))
                      dmin = two*(s(i+1,j,k,comp) - s(i,j,k,comp))
                      slim = min(abs(dpls), abs(dmin))
                      slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                      sflag = sign(one,del)
                      slx(i,j,k,comp)= sflag*min(slim,abs(del))
                   endif
                endif

             enddo
          enddo
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
          do k = lo(3), hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   ! left
                   dcen = half*(s(i,j,k,comp)-s(i-2,j,k,comp))
                   dmin = two*(s(i-1,j,k,comp)-s(i-2,j,k,comp))
                   dpls = two*(s(i,j,k,comp)-s(i-1,j,k,comp))
                   dlim = min(abs(dmin),abs(dpls))
                   dlim = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)
                   dxl = dflag*min(dlim, abs(dcen))

                   ! right
                   dcen = half*(s(i+2,j,k,comp)-s(i,j,k,comp))
                   dmin = two*(s(i+1,j,k,comp)-s(i,j,k,comp))
                   dpls = two*(s(i+2,j,k,comp)-s(i+1,j,k,comp))
                   dlim = min(abs(dmin),abs(dpls))
                   dlim = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)
                   dxr = dflag*min(dlim, abs(dcen))

                   ! center
                   dcen = half*(s(i+1,j,k,comp)-s(i-1,j,k,comp))
                   dmin = two*(s(i  ,j,k,comp)-s(i-1,j,k,comp))
                   dpls = two*(s(i+1,j,k,comp)-s(i  ,j,k,comp))
                   dlim = min(abs(dmin),abs(dpls))
                   dlim = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)

                   ds = two * two3rd * dcen - sixth * (dxr + dxl)
                   slx(i,j,k,comp) = dflag*min(abs(ds),dlim)


                   if (adv_bc(1,1,bc_comp) .eq. EXT_DIR  .or. adv_bc(1,1,bc_comp) .eq. HOEXTRAP) then
                      if (i .eq. lo(1) .and. lo(1)+1 .eq. domlo(1)) then
                         slx(i,j,k,comp) = ZERO
                      elseif (i .eq. lo(1)+1 .and. lo(1)+1 .eq. domlo(1)) then
                         del = -sixteen/fifteen*s(i-1,j,k,comp) + half*s(i,j,k,comp) + &
                              two3rd*s(i+1,j,k,comp) - tenth*s(i+2,j,k,comp)
                         dmin = two*(s(i,j,k,comp)-s(i-1,j,k,comp))
                         dpls = two*(s(i+1,j,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slx(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (i .eq. lo(1)+2 .and. lo(1)+1 .eq. domlo(1)) then
                         ! Recalculate the slope at lo(1)+1 using the revised dxl
                         del = -sixteen/fifteen*s(i-2,j,k,comp) + half*s(i-1,j,k,comp) + &
                              two3rd*s(i,j,k,comp) - tenth*s(i+1,j,k,comp)
                         dmin = two*(s(i-1,j,k,comp)-s(i-2,j,k,comp))
                         dpls = two*(s(i,j,k,comp)-s(i-1,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         dxl = sflag*min(slim,abs(del))

                         ds = two * two3rd * dcen - sixth * (dxr + dxl)
                         slx(i,j,k,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif


                   if (adv_bc(1,2,bc_comp) .eq. EXT_DIR  .or. adv_bc(1,2,bc_comp) .eq. HOEXTRAP) then
                      if (i .eq. hi(1) .and. hi(1)-1 .eq. domhi(1)) then
                         slx(i,j,k,comp) = ZERO
                      elseif (i .eq. hi(1)-1 .and. hi(1)-1 .eq. domhi(1)) then
                         del = -( -sixteen/fifteen*s(i+1,j,k,comp) + half*s(i,j,k,comp) +  &
                              two3rd*s(i-1,j,k,comp) - tenth*s(i-2,j,k,comp) )
                         dmin = two*(s(i,j,k,comp)-s(i-1,j,k,comp))
                         dpls = two*(s(i+1,j,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slx(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (i .eq. hi(1)-2 .and. hi(1)-1 .eq. domhi(1)) then
                         ! Recalculate the slope at hi(1)-1 using the revised dxr
                         del = -( -sixteen/fifteen*s(i+2,j,k,comp) + half*s(i+1,j,k,comp) +  &
                              two3rd*s(i,j,k,comp) - tenth*s(i-1,j,k,comp) )
                         dmin = two*(s(i+1,j,k,comp)-s(i,j,k,comp))
                         dpls = two*(s(i+2,j,k,comp)-s(i+1,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         dxr = sflag*min(slim,abs(del))

                         ds = two * two3rd * dcen - sixth * (dxl + dxr)
                         slx(i,j,k,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo

    endif

  end subroutine slopex_2d

  subroutine slopey_2d(lo,hi,s,s_lo,s_hi,nc_s,sly,sl_lo,sl_hi,nc_sl, &
       domlo,domhi,nvar,adv_bc,nbccomp,bccomp)  bind(C,name="slopey_2d")

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in)  :: domlo(3),domhi(3),lo(3),hi(3)
    integer         , intent(in   ) :: s_lo(3),s_hi(3),sl_lo(3),sl_hi(3)
    integer  , value, intent(in   ) :: nc_s,nc_sl,nvar,nbccomp,bccomp
    double precision, intent(in   ) ::   s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(  out) :: sly(sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3),nc_sl)
    integer         , intent(in)  :: adv_bc(AMREX_SPACEDIM,2,nbccomp)

    ! local
    double precision :: dpls,dmin,ds,del,slim,sflag
    double precision :: dcen,dlim,dflag,dyl,dyr

    integer :: i,j,k,comp,bc_comp

    !$gpu

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   sly(i,j,k,comp) = ZERO
                enddo
             enddo
          enddo
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)

                   del  = half*(s(i,j+1,k,comp) - s(i,j-1,k,comp))
                   dpls = two *(s(i,j+1,k,comp) - s(i,j  ,k,comp))
                   dmin = two *(s(i,j  ,k,comp) - s(i,j-1,k,comp))
                   slim = min(abs(dpls),abs(dmin))
                   slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   sly(i,j,k,comp)= sflag*min(slim,abs(del))


                   if (adv_bc(2,1,bc_comp) .eq. EXT_DIR .or. adv_bc(2,1,bc_comp) .eq. HOEXTRAP) then
                      if (j .eq. lo(2) .and. lo(2)+1 .eq. domlo(2)) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. lo(2)+1 .and. lo(2)+1 .eq. domlo(2)) then
                         del = (s(i,j+1,k,comp)+three*s(i,j,k,comp)- &
                              four*s(i,j-1,k,comp)) * third
                         dpls = two*(s(i,j+1,k,comp) - s(i,j,k,comp))
                         dmin = two*(s(i,j,k,comp) - s(i,j-1,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         sly(i,j,k,comp)= sflag*min(slim,abs(del))
                      endif
                   endif


                   if (adv_bc(2,2,bc_comp) .eq. EXT_DIR .or. adv_bc(2,2,bc_comp) .eq. HOEXTRAP) then
                      if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. hi(2)-1 .and. hi(2)-1 .eq. domhi(2)) then
                         del = -(s(i,j-1,k,comp)+three*s(i,j,k,comp)- &
                              four*s(i,j+1,k,comp)) * third
                         dpls = two*(s(i,j+1,k,comp) - s(i,j,k,comp))
                         dmin = two*(s(i,j,k,comp) - s(i,j-1,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         sly(i,j,k,comp)= sflag*min(slim,abs(del))

                      endif
                   endif
                enddo
             enddo
          enddo
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   ! left
                   dcen = half*(s(i,j,k,comp)-s(i,j-2,k,comp))
                   dmin = two*(s(i,j-1,k,comp)-s(i,j-2,k,comp))
                   dpls = two*(s(i,j,k,comp)-s(i,j-1,k,comp))
                   dlim  = min(abs(dmin),abs(dpls))
                   dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)
                   dyl = dflag*min(dlim,abs(dcen))

                   ! right
                   dcen = half*(s(i,j+2,k,comp)-s(i,j,k,comp))
                   dmin = two*(s(i,j+1,k,comp)-s(i,j,k,comp))
                   dpls = two*(s(i,j+2,k,comp)-s(i,j+1,k,comp))
                   dlim  = min(abs(dmin),abs(dpls))
                   dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)
                   dyr = dflag*min(dlim,abs(dcen))

                   ! center
                   dcen = half*(s(i,j+1,k,comp)-s(i,j-1,k,comp))
                   dmin = two*(s(i,j  ,k,comp)-s(i,j-1,k,comp))
                   dpls = two*(s(i,j+1,k,comp)-s(i,j  ,k,comp))
                   dlim  = min(abs(dmin),abs(dpls))
                   dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)

                   ds = two * two3rd * dcen - sixth * (dyr + dyl)
                   sly(i,j,k,comp) = dflag*min(abs(ds),dlim)


                   if (adv_bc(2,1,bc_comp) .eq. EXT_DIR .or. adv_bc(2,1,bc_comp) .eq. HOEXTRAP) then
                      if (j .eq. lo(2) .and. lo(2)+1 .eq. domlo(2)) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. lo(2)+1 .and. lo(2)+1 .eq. domlo(2)) then
                         del = -sixteen/fifteen*s(i,j-1,k,comp) +  half*s(i,j,k,comp) +  &
                              two3rd*s(i,j+1,k,comp) - tenth*s(i,j+2,k,comp)
                         dmin = two*(s(i,j,k,comp)-s(i,j-1,k,comp))
                         dpls = two*(s(i,j+1,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         sly(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (j .eq. lo(2)+2 .and. lo(2)+1 .eq. domlo(2)) then

                         ! Recalculate the slope at lo(2)+1 using the revised dyl
                         del = -sixteen/fifteen*s(i,j-2,k,comp) +  half*s(i,j-1,k,comp) +  &
                              two3rd*s(i,j,k,comp) - tenth*s(i,j+1,k,comp)
                         dmin = two*(s(i,j-1,k,comp)-s(i,j-2,k,comp))
                         dpls = two*(s(i,j,k,comp)-s(i,j-1,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         dyl = sflag*min(slim,abs(del))
                         ds = two * two3rd * dcen - sixth * (dyr + dyl)
                         sly(i,j,k,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif


                   if (adv_bc(2,2,bc_comp) .eq. EXT_DIR .or. adv_bc(2,2,bc_comp) .eq. HOEXTRAP) then
                      if (j .eq. hi(2) .and. hi(2)-1 .eq. domhi(2)) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. hi(2)-1 .and. hi(2)-1 .eq. domhi(2)) then
                         del = -( -sixteen/fifteen*s(i,j+1,k,comp) +  half*s(i,j,k,comp) + &
                              two3rd*s(i,j-1,k,comp) - tenth*s(i,j-2,k,comp) )
                         dmin = two*(s(i,j,k,comp)-s(i,j-1,k,comp))
                         dpls = two*(s(i,j+1,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         sly(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (j .eq. hi(2)-2 .and. hi(2)-1 .eq. domhi(2)) then
                         ! Recalculate the slope at lo(2)+1 using the revised dyr
                         del = -( -sixteen/fifteen*s(i,j+2,k,comp) +  half*s(i,j+1,k,comp) + &
                              two3rd*s(i,j,k,comp) - tenth*s(i,j-1,k,comp) )
                         dmin = two*(s(i,j+1,k,comp)-s(i,j,k,comp))
                         dpls = two*(s(i,j+2,k,comp)-s(i,j+1,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         dyr = sflag*min(slim,abs(del))
                         ds = two * two3rd * dcen - sixth * (dyl + dyr)
                         sly(i,j,k,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif

                enddo
             enddo
          enddo
       enddo

    endif

  end subroutine slopey_2d


  subroutine slopez_3d(lo,hi,s,s_lo,s_hi,nc_s,slz,sl_lo,sl_hi,nc_sl, &
                       domlo,domhi,nvar,adv_bc,nbccomp,bccomp) &
                       bind(C, name="slopez_3d")

    use amrex_constants_module
    use meth_params_module, only : slope_order

    integer         , intent(in)  :: domlo(3),domhi(3),lo(3),hi(3)
    integer         , intent(in)  :: s_lo(3),s_hi(3),sl_lo(3),sl_hi(3)
    integer  , value, intent(in)  :: nvar,nbccomp,bccomp,nc_s,nc_sl
    double precision, intent(in   ) ::   s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(  out) :: slz(sl_lo(1):sl_hi(1),sl_lo(2):sl_hi(2),sl_lo(3):sl_hi(3),nc_sl)
    integer         , intent(in)  :: adv_bc(AMREX_SPACEDIM,2,nbccomp)

    double precision :: dpls,dmin,ds,del,slim,sflag
    double precision :: dcen,dlim,dflag,dzl,dzr

    integer :: i,j,k,comp

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
          do k = lo(3)-1,hi(3)+1
             do j = lo(2)-1,hi(2)+1
                do i = lo(1)-1,hi(1)+1
                   slz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,comp) = ZERO
                enddo
             enddo
          enddo
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do k = lo(3)-1,hi(3)+1
             do j = lo(2)-1,hi(2)+1
                do i = lo(1)-1,hi(1)+1

                   del  = half*(s(i,j,k+1,comp) - s(i,j,k-1,comp))
                   dpls = two *(s(i,j,k+1,comp) - s(i,j,k  ,comp))
                   dmin = two *(s(i,j,k  ,comp) - s(i,j,k-1,comp))
                   slim = min(abs(dpls),abs(dmin))
                   slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slz(i,j,k,comp)= sflag*min(slim,abs(del))


                   if (adv_bc(3,1,comp) .eq. EXT_DIR .or. adv_bc(3,1,comp) .eq. HOEXTRAP) then
                      if (k .eq. lo(3)-1 .and. lo(3) .eq. domlo(3)) then
                         slz(i,j,k,comp) = ZERO
                      elseif (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                         del = (s(i,j,k+1,comp)+three*s(i,j,k,comp)- &
                              four*s(i,j,k-1,comp)) * third
                         dpls = two*(s(i,j,k+1,comp) - s(i,j,k,comp))
                         dmin = two*(s(i,j,k,comp) - s(i,j,k-1,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slz(i,j,k,comp)= sflag*min(slim,abs(del))
                      endif
                   endif


                   if (adv_bc(3,2,comp) .eq. EXT_DIR .or. adv_bc(3,2,comp) .eq. HOEXTRAP) then
                      if (k .eq. hi(3)+1 .and. hi(3) .eq. domhi(3)) then
                         slz(i,j,k,comp) = ZERO
                      elseif (k .eq. hi(3)+1 .and. hi(3) .eq. domhi(3)) then
                         del = -(s(i,j,k-1,comp)+three*s(i,j,k,comp)- &
                              four*s(i,j,k+1,comp)) * third
                         dpls = two*(s(i,j,k+1,comp) - s(i,j,k,comp))
                         dmin = two*(s(i,j,k,comp) - s(i,j,k-1,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slz(i,j,k,comp)= sflag*min(slim,abs(del))
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                do k = lo(3)-1,hi(3)+1
                    ! left
                   dcen = half*(s(i,j,k,comp)-s(i,j,k-2,comp))
                   dmin = two*(s(i,j,k-1,comp)-s(i,j,k-2,comp))
                   dpls = two*(s(i,j,k,comp)-s(i,j,k-1,comp))
                   dlim  = min(abs(dmin),abs(dpls))
                   dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)
                   dzl = dflag*min(dlim,abs(dcen))

                   ! right
                   dcen = half*(s(i,j,k+2,comp)-s(i,j,k,comp))
                   dmin = two*(s(i,j,k+1,comp)-s(i,j,k,comp))
                   dpls = two*(s(i,j,k+2,comp)-s(i,j,k+1,comp))
                   dlim  = min(abs(dmin),abs(dpls))
                   dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)
                   dzr = dflag*min(dlim,abs(dcen))

                   ! center
                   dcen = half*(s(i,j,k+1,comp)-s(i,j,k-1,comp))
                   dmin = two*(s(i,j,k  ,comp)-s(i,j,k-1,comp))
                   dpls = two*(s(i,j,k+1,comp)-s(i,j,k  ,comp))
                   dlim  = min(abs(dmin),abs(dpls))
                   dlim  = merge(dlim,ZERO,dpls*dmin.gt.ZERO)
                   dflag = sign(one,dcen)

                   ds = two * two3rd * dcen -  &
                        sixth * (dzr + dzl)
                   slz(i,j,k,comp) = dflag*min(abs(ds),dlim)

                   if (adv_bc(3,1,comp) .eq. EXT_DIR .or. adv_bc(3,1,comp) .eq. HOEXTRAP) then
                      if (k .eq. lo(3)-1 .and. lo(3) .eq. domlo(3)) then
                         slz(i,j,k,comp) = ZERO
                      elseif (k .eq. lo(3) .and. lo(3) .eq. domlo(3)) then
                         del = -sixteen/fifteen*s(i,j,k-1,comp) +  half*s(i,j,k,comp) +  &
                              two3rd*s(i,j,k+1,comp) - tenth*s(i,j,k+2,comp)
                         dmin = two*(s(i,j,k,comp)-s(i,j,k-1,comp))
                         dpls = two*(s(i,j,k+1,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slz(i,j,k,comp)= sflag*min(slim,abs(del))

                      elseif (k .eq. lo(3)+1 .and. lo(3) .eq. domlo(3)) then
                         ! Recalculate the slope at lo(2)+1 using the revised dzl
                         del = -sixteen/fifteen*s(i,j,k-2,comp) +  half*s(i,j,k-1,comp) +  &
                              two3rd*s(i,j,k,comp) - tenth*s(i,j,k+1,comp)
                         dmin = two*(s(i,j,k-1,comp)-s(i,j,k-2,comp))
                         dpls = two*(s(i,j,k,comp)-s(i,j,k-1,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         dzl = slz(i,j,k-1,comp)
                         ds = two * two3rd * dcen - &
                              sixth * (dzr + dzl)
                         slz(i,j,lo(3)+1,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif


                   if (adv_bc(3,2,comp) .eq. EXT_DIR .or. adv_bc(3,2,comp) .eq. HOEXTRAP) then
                      if (k .eq. hi(3)+1 .and. hi(3) .eq. domhi(3)) then
                         slz(i,j,k,comp) = ZERO
                      elseif (k .eq. hi(3) .and. hi(3) .eq. domhi(3)) then
                         del = -( -sixteen/fifteen*s(i,j,hi(3)+1,comp) +  half*s(i,j,hi(3)  ,comp) + &
                              two3rd*s(i,j,k-1,comp) - tenth*s(i,j,k-2,comp) )
                         dmin = two*(s(i,j,k,comp)-s(i,j,k-1,comp))
                         dpls = two*(s(i,j,k+1,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slz(i,j,k,comp)= sflag*min(slim,abs(del))

                      elseif (k .eq. hi(3)-1 .and. hi(3) .eq. domhi(3)) then
                         ! Recalculate the slope at lo(3)+1 using the revised dzr
                         del = -( -sixteen/fifteen*s(i,j,k+2,comp) +  half*s(i,j,k+1,comp) + &
                              two3rd*s(i,j,k,comp) - tenth*s(i,j,k-1,comp) )
                         dmin = two*(s(i,j,k+1,comp)-s(i,j,k,comp))
                         dpls = two*(s(i,j,k+2,comp)-s(i,j,k+1,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         dzr = slz(i,j,k+1,comp)
                         ds = two * two3rd * dcen -  &
                              sixth * (dzl + dzr)
                         slz(i,j,k,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif

                enddo
             enddo
          enddo
       enddo

    endif

  end subroutine slopez_3d



end module slope_module
