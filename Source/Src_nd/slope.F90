
#include "AMReX_BC_TYPES.H"

module slope_module

  implicit none

  integer, private, parameter :: cen = 1, lim = 2, flag = 3, fromm = 4

contains

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
#if (AMREX_SPACEDIM == 3)
           do k = lo(3), hi(3)
#endif
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   slx(i,j,k,comp) = ZERO
                enddo
             enddo
#if (AMREX_SPACEDIM == 3)
         enddo
#endif
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
#if (AMREX_SPACEDIM == 3)
          do k = lo(3), hi(3)
#endif
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
                   if (i .eq. domlo(1)-1) then
                      slx(i,j,k,comp) = ZERO
                   elseif (i .eq. domlo(1)) then
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
                   if (i .eq. domhi(1)+1) then
                      slx(i,j,k,comp) = ZERO
                   elseif (i .eq. domhi(1)) then
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
#if (AMREX_SPACEDIM == 3)
         enddo
#endif
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
#if (AMREX_SPACEDIM == 3)
          do k = lo(3), hi(3)
#endif
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
                      if (i .eq. domlo(1)-1) then
                         slx(i,j,k,comp) = ZERO
                      elseif (i .eq. domlo(1)) then
                         del = -sixteen/fifteen*s(i-1,j,k,comp) + half*s(i,j,k,comp) + &
                              two3rd*s(i+1,j,k,comp) - tenth*s(i+2,j,k,comp)
                         dmin = two*(s(i,j,k,comp)-s(i-1,j,k,comp))
                         dpls = two*(s(i+1,j,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slx(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (i .eq. domlo(1)+1) then
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
                      if (i .eq. domhi(1)+1) then
                         slx(i,j,k,comp) = ZERO
                      elseif (i .eq. domhi(1)) then
                         del = -( -sixteen/fifteen*s(i+1,j,k,comp) + half*s(i,j,k,comp) +  &
                              two3rd*s(i-1,j,k,comp) - tenth*s(i-2,j,k,comp) )
                         dmin = two*(s(i,j,k,comp)-s(i-1,j,k,comp))
                         dpls = two*(s(i+1,j,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slx(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (i .eq. domhi(1)-1) then
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
#if (AMREX_SPACEDIM == 3)
         enddo
#endif
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

    k = s_lo(3)

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
#if (AMREX_SPACEDIM == 3)
           do k = lo(3), hi(3)
#endif
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   sly(i,j,k,comp) = ZERO
                enddo
             enddo
#if (AMREX_SPACEDIM == 3)
         enddo
#endif
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
#if (AMREX_SPACEDIM == 3)
          do k = lo(3), hi(3)
#endif
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
                      if (j .eq. domlo(2)-1) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. domlo(2)) then
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
                      if (j .eq. domhi(2)+1) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. domhi(2)) then
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
#if (AMREX_SPACEDIM == 3)
         enddo
#endif
       enddo

    else

       ! HERE DOING 4TH ORDER
       do comp=1,nvar
          bc_comp = bccomp + comp-1
#if (AMREX_SPACEDIM == 3)
          do k = lo(3), hi(3)
#endif
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
                      if (j .eq. domlo(2)-1) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. domlo(2)) then
                         del = -sixteen/fifteen*s(i,j-1,k,comp) +  half*s(i,j,k,comp) +  &
                              two3rd*s(i,j+1,k,comp) - tenth*s(i,j+2,k,comp)
                         dmin = two*(s(i,j,k,comp)-s(i,j-1,k,comp))
                         dpls = two*(s(i,j+1,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         sly(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (j .eq. domlo(2)+1) then

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
                      if (j .eq. domhi(2)+1) then
                         sly(i,j,k,comp) = ZERO
                      elseif (j .eq. domhi(2)) then
                         del = -( -sixteen/fifteen*s(i,j+1,k,comp) +  half*s(i,j,k,comp) + &
                              two3rd*s(i,j-1,k,comp) - tenth*s(i,j-2,k,comp) )
                         dmin = two*(s(i,j,k,comp)-s(i,j-1,k,comp))
                         dpls = two*(s(i,j+1,k,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         sly(i,j,k,comp)= sflag*min(slim,abs(del))
                      elseif (j .eq. domhi(2)-1) then
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
#if (AMREX_SPACEDIM == 3)
         enddo
#endif
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

    !$gpu

    if (slope_order .eq. 0) then

       ! HERE DOING 1ST ORDER
       do comp=1,nvar
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   slz(i,j,k,comp) = ZERO
                enddo
             enddo
          enddo
       enddo

    else if (slope_order .eq. 2) then

       ! HERE DOING 2ND ORDER
       do comp=1,nvar
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)

                   del  = half*(s(i,j,k+1,comp) - s(i,j,k-1,comp))
                   dpls = two *(s(i,j,k+1,comp) - s(i,j,k  ,comp))
                   dmin = two *(s(i,j,k  ,comp) - s(i,j,k-1,comp))
                   slim = min(abs(dpls),abs(dmin))
                   slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                   sflag = sign(one,del)
                   slz(i,j,k,comp)= sflag*min(slim,abs(del))


                   if (adv_bc(3,1,comp) .eq. EXT_DIR .or. adv_bc(3,1,comp) .eq. HOEXTRAP) then
                      if (k .eq. domlo(3)-1) then
                         slz(i,j,k,comp) = ZERO
                     elseif (k .eq. domlo(3)) then
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
                      if (k .eq. domhi(3)+1) then
                         slz(i,j,k,comp) = ZERO
                     elseif (k .eq. domhi(3)) then
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
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                do k = lo(3),hi(3)
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
                      if (k .eq. domlo(3)-1) then
                         slz(i,j,k,comp) = ZERO
                     elseif (k .eq. domlo(3)) then
                         del = -sixteen/fifteen*s(i,j,k-1,comp) +  half*s(i,j,k,comp) +  &
                              two3rd*s(i,j,k+1,comp) - tenth*s(i,j,k+2,comp)
                         dmin = two*(s(i,j,k,comp)-s(i,j,k-1,comp))
                         dpls = two*(s(i,j,k+1,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slz(i,j,k,comp)= sflag*min(slim,abs(del))

                     elseif (k .eq. domlo(3)+1) then
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
                         slz(i,j,k,comp) = dflag*min(abs(ds),dlim)
                      endif
                   endif


                   if (adv_bc(3,2,comp) .eq. EXT_DIR .or. adv_bc(3,2,comp) .eq. HOEXTRAP) then
                      if (k .eq. domhi(3)+1) then
                         slz(i,j,k,comp) = ZERO
                     elseif (k .eq. domhi(3)) then
                         del = -( -sixteen/fifteen*s(i,j,k+1,comp) +  half*s(i,j,k,comp) + &
                              two3rd*s(i,j,k-1,comp) - tenth*s(i,j,k-2,comp) )
                         dmin = two*(s(i,j,k,comp)-s(i,j,k-1,comp))
                         dpls = two*(s(i,j,k+1,comp)-s(i,j,k,comp))
                         slim = min(abs(dpls), abs(dmin))
                         slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                         sflag = sign(one,del)
                         slz(i,j,k,comp)= sflag*min(slim,abs(del))

                     elseif (k .eq. domhi(3)-1) then
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
