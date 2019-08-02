! compute the PPM integrals, Ip and Im.  These are the integrals under
! the parabolic profile of the reconstructed quantity over the domain
! that can reach the interface over the timestep dt.
!
! Ip captures the amount of the state that can reach the right
! interface of the cell and Im captures what can reach the left
! interface of the cell over the step.
!
! There are cases here: one (originally called the 'fpu version') uses
! the MAC velocity for the tracing while the non-fpu versions use the
! cell-centered velocity.

#include "AMReX_BC_TYPES.H"

module ppm_module

  use amrex_error_module
  use amrex_constants_module
  use meth_params_module, only: ppm_type, rel_eps

  implicit none

contains

#if (AMREX_SPACEDIM == 1)
  !===========================================================================
  ! 1-d version
  !===========================================================================
  subroutine ppm_1d(s,ng_s,u,ng_u,Ip,Im,domlo,domhi,lo,hi,adv_bc,dx,dt,is_umac)

    integer         , intent(in   ) :: domlo(:),domhi(:),lo(:),hi(:),ng_s,ng_u
    double precision, intent(in   ) ::  s(lo(1)-ng_s:)
    double precision, intent(in   ) ::  u(lo(1)-ng_u:)
    double precision, intent(inout) :: Ip(lo(1)-1   :)
    double precision, intent(inout) :: Im(lo(1)-1   :)
    integer         , intent(in   ) :: adv_bc(:,:)
    double precision, intent(in   ) :: dx(:),dt
    logical         , intent(in   ) :: is_umac

    ! local
    integer :: i
    logical :: extremum, bigp, bigm

    double precision :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision :: sgn, sigma, s6, amax, delam, delap, D2ABS
    double precision :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
    double precision :: dslv_l, dslv_r

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! s_{\ib,+}, s_{\ib,-}
    double precision, allocatable :: sp(:)
    double precision, allocatable :: sm(:)

    ! \delta s_{\ib}^{vL}
    double precision, allocatable :: dsvl(:)

    ! s_{i+\half}^{H.O.}
    double precision, allocatable :: sedge(:)

    ! cell-centered indexing
    allocate(sp(lo(1)-1:hi(1)+1))
    allocate(sm(lo(1)-1:hi(1)+1))

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(lo(1)-2:hi(1)+2))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(lo(1)-1:hi(1)+2))
    else if (ppm_type .eq. 2) then
       allocate(sedge(lo(1)-2:hi(1)+3))
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       !----------------------------------------------------------------------
       ! ppm_type = 1
       !----------------------------------------------------------------------

       ! compute van Leer slopes in x-direction
       dsvl = ZERO
       do i=lo(1)-2,hi(1)+2
          dsc = HALF * (s(i+1) - s(i-1))
          dsl = TWO  * (s(i  ) - s(i-1))
          dsr = TWO  * (s(i+1) - s(i  ))
          if (dsl*dsr .gt. ZERO) dsvl(i) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do

       ! interpolate s to x-edges
       do i=lo(1)-1,hi(1)+2
          sedge(i) = HALF*(s(i)+s(i-1)) - SIXTH*(dsvl(i)-dsvl(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i) = max(sedge(i),min(s(i),s(i-1)))
          sedge(i) = min(sedge(i),max(s(i),s(i-1)))
       end do

       ! copy sedge into sp and sm
       do i=lo(1)-1,hi(1)+1
          sp(i) = sedge(i+1)
          sm(i) = sedge(i  )
       end do

       ! modify using quadratic limiters
       do i=lo(1)-1,hi(1)+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(1) .eq. domlo(1)) then
          if (adv_bc(1,1) .eq. EXT_DIR  .or. adv_bc(1,1) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sm(lo(1)) = s(lo(1)-1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)+1) = &
                  -FIFTH        *s(lo(1)-1) &
                  + (THREE/FOUR)*s(lo(1)  ) &
                  + HALF        *s(lo(1)+1) &
                  - (ONE/20.0d0)*s(lo(1)+2)

             ! make sure sedge lies in between adjacent cell-centered values
             sedge(lo(1)+1) = max(sedge(lo(1)+1),min(s(lo(1)+1),s(lo(1))))
             sedge(lo(1)+1) = min(sedge(lo(1)+1),max(s(lo(1)+1),s(lo(1))))

             ! copy sedge into sp and sm
             sp(lo(1)  ) = sedge(lo(1)+1)
             sm(lo(1)+1) = sedge(lo(1)+1)

             ! reset sp on second interior edge
             sp(lo(1)+1) = sedge(lo(1)+2)

             ! modify using quadratic limiters
             i = lo(1)+1
             if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
                sp(i) = s(i)
                sm(i) = s(i)
             else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
                sp(i) = THREE*s(i) - TWO*sm(i)
             else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
                sm(i) = THREE*s(i) - TWO*sp(i)
             end if
          end if
       end if

       if (hi(1) .eq. domhi(1)) then
          if (adv_bc(1,2) .eq. EXT_DIR  .or. adv_bc(1,2) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sp(hi(1)) = s(hi(1)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(hi(1)) = &
                  -FIFTH        *s(hi(1)+1) &
                  + (THREE/FOUR)*s(hi(1)  ) &
                  + HALF        *s(hi(1)-1) &
                  - (ONE/20.0d0)*s(hi(1)-2)

             ! make sure sedge lies in between adjacent cell-centered values
             sedge(hi(1)) = max(sedge(hi(1)),min(s(hi(1)-1),s(hi(1))))
             sedge(hi(1)) = min(sedge(hi(1)),max(s(hi(1)-1),s(hi(1))))

             ! copy sedge into sp and sm
             sp(hi(1)-1) = sedge(hi(1))
             sm(hi(1)  ) = sedge(hi(1))

             ! reset sm on second interior edge
             sm(hi(1)-1) = sedge(hi(1)-1)

             ! modify using quadratic limiters
             i = hi(1)-1
             if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
                sp(i) = s(i)
                sm(i) = s(i)
             else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
                sp(i) = THREE*s(i) - TWO*sm(i)
             else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
                sm(i) = THREE*s(i) - TWO*sp(i)
             end if
          end if
       end if

    else if (ppm_type .eq. 2) then

       !----------------------------------------------------------------------
       ! ppm_type = 2
       !----------------------------------------------------------------------
#ifndef AMREX_USE_GPU
       if (ng_s .lt. 4) then
          call amrex_error("Need 4 ghost cells for ppm_type=2")
       end if
#endif

       ! interpolate s to x-edges
       do i=lo(1)-2,hi(1)+3
          sedge(i) = (7.d0/12.d0)*(s(i-1)+s(i)) - (1.d0/12.d0)*(s(i-2)+s(i+1))
          ! limit sedge
          if ((sedge(i)-s(i-1))*(s(i)-sedge(i)) .lt. ZERO) then
             D2  = THREE*(s(i-1)-TWO*sedge(i)+s(i))
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i) = HALF*(s(i-1)+s(i)) - SIXTH*D2LIM
          end if
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       do i=lo(1)-1,hi(1)+1

          alphap = sedge(i+1)-s(i)
          alpham = sedge(i  )-s(i)
          bigp = abs(alphap).gt.TWO*abs(alpham)
          bigm = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             dafacem = sedge(i) - sedge(i-1)
             dafacep = sedge(i+2) - sedge(i+1)
             dabarm = s(i) - s(i-1)
             dabarp = s(i+1) - s(i)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin= min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2  = SIX*(alpham + alphap)
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i)-TWO*s(i+1)+s(i+2)
             D2C = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             D2ABS = max(abs(D2),1.d-10)
             alpham = alpham*D2LIM/D2ABS
             alphap = alphap*D2LIM/D2ABS
          else
             if (bigp) then
                sgn = sign(ONE,alpham)
                amax = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1) - s(i)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn = sign(ONE,alphap)
                amax = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1) - s(i)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i) = s(i) + alpham
          sp(i) = s(i) + alphap

       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(1) .eq. domlo(1)) then
          if (adv_bc(1,1) .eq. EXT_DIR  .or. adv_bc(1,1) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sm(lo(1))    = s(lo(1)-1)
             sedge(lo(1)) = s(lo(1)-1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)+1) = &
                  -FIFTH        *s(lo(1)-1) &
                  + (THREE/FOUR)*s(lo(1)  ) &
                  + HALF        *s(lo(1)+1) &
                  - (ONE/20.0d0)*s(lo(1)+2)

             ! make sure sedge lies in between adjacent cell-centered values
             sedge(lo(1)+1) = max(sedge(lo(1)+1),min(s(lo(1)+1),s(lo(1))))
             sedge(lo(1)+1) = min(sedge(lo(1)+1),max(s(lo(1)+1),s(lo(1))))

             ! copy sedge into sp
             sp(lo(1)  ) = sedge(lo(1)+1)

             ! apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells
             do i=lo(1)+1,lo(1)+2

                alphap = sedge(i+1)-s(i)
                alpham = sedge(i  )-s(i)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i) - sedge(i-1)
                   dafacep = sedge(i+2) - sedge(i+1)
                   dabarm = s(i) - s(i-1)
                   dabarp = s(i+1) - s(i)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2)-TWO*s(i-1)+s(i)
                   D2R = s(i)-TWO*s(i+1)+s(i+2)
                   D2C = s(i-1)-TWO*s(i)+s(i+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   D2ABS = max(abs(D2),1.d-10)
                   alpham = alpham*D2LIM/D2ABS
                   alphap = alphap*D2LIM/D2ABS
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1) - s(i)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1) - s(i)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i) = s(i) + alpham
                sp(i) = s(i) + alphap

             end do
          end if
       end if

       if (hi(1) .eq. domhi(1)) then
          if (adv_bc(1,2) .eq. EXT_DIR  .or. adv_bc(1,2) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sp(hi(1)     ) = s(hi(1)+1)
             sedge(hi(1)+1) = s(hi(1)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(hi(1)) = &
                  -FIFTH        *s(hi(1)+1) &
                  + (THREE/FOUR)*s(hi(1)  ) &
                  + HALF        *s(hi(1)-1) &
                  - (ONE/20.0d0)*s(hi(1)-2)

             ! make sure sedge lies in between adjacent cell-centered values
             sedge(hi(1)) = max(sedge(hi(1)),min(s(hi(1)-1),s(hi(1))))
             sedge(hi(1)) = min(sedge(hi(1)),max(s(hi(1)-1),s(hi(1))))

             ! copy sedge into sm
             sm(hi(1)  ) = sedge(hi(1))

             ! reset sm on second interior edge
             sm(hi(1)-1) = sedge(hi(1)-1)

             ! apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells
             do i=hi(1)-2,hi(1)-1

                alphap = sedge(i+1)-s(i)
                alpham = sedge(i  )-s(i)
                bigp = abs(alphap).gt.TWO*abs(alpham)
                bigm = abs(alpham).gt.TWO*abs(alphap)
                extremum = .false.

                if (alpham*alphap .ge. ZERO) then
                   extremum = .true.
                else if (bigp .or. bigm) then
                   ! Possible extremum. We look at cell centered values and face
                   ! centered values for a change in sign in the differences adjacent to
                   ! the cell. We use the pair of differences whose minimum magnitude is the
                   ! largest, and thus least susceptible to sensitivity to roundoff.
                   dafacem = sedge(i) - sedge(i-1)
                   dafacep = sedge(i+2) - sedge(i+1)
                   dabarm = s(i) - s(i-1)
                   dabarp = s(i+1) - s(i)
                   dafacemin = min(abs(dafacem),abs(dafacep))
                   dabarmin= min(abs(dabarm),abs(dabarp))
                   if (dafacemin.ge.dabarmin) then
                      dachkm = dafacem
                      dachkp = dafacep
                   else
                      dachkm = dabarm
                      dachkp = dabarp
                   endif
                   extremum = (dachkm*dachkp .le. 0.d0)
                end if

                if (extremum) then
                   D2  = SIX*(alpham + alphap)
                   D2L = s(i-2)-TWO*s(i-1)+s(i)
                   D2R = s(i)-TWO*s(i+1)+s(i+2)
                   D2C = s(i-1)-TWO*s(i)+s(i+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   D2ABS = max(abs(D2),1.d-10)
                   alpham = alpham*D2LIM/D2ABS
                   alphap = alphap*D2LIM/D2ABS
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1) - s(i)
                      if (sgn*amax .ge. sgn*delam) then
                         if (sgn*(delam - alpham).ge.1.d-10) then
                            alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                         else
                            alphap = -TWO*alpham
                         endif
                      endif
                   end if
                   if (bigm) then
                      sgn = sign(ONE,alphap)
                      amax = -alpham**2 / (4*(alpham + alphap))
                      delap = s(i+1) - s(i)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i) = s(i) + alpham
                sp(i) = s(i) + alphap

             end do
          end if

       end if
    end if

    !-------------------------------------------------------------------------
    ! compute x-component of Ip and Im
    !-------------------------------------------------------------------------

    if (is_umac) then

       ! u is the MAC velocity, use edge-based indexing
       do i=lo(1)-1,hi(1)
          sigma = abs(u(i+1))*dt/dx(1)
          s6 = SIX*s(i) - THREE*(sm(i)+sp(i))
          if (u(i+1) .gt. rel_eps) then
             Ip(i) = sp(i) - (sigma/TWO)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma)*s6)
          else
             Ip(i) = s(i)
          end if
       end do

       do i=lo(1),hi(1)+1
          sigma = abs(u(i))*dt/dx(1)
          s6 = SIX*s(i) - THREE*(sm(i)+sp(i))
          if (u(i) .lt. -rel_eps) then
             Im(i) = sm(i) + (sigma/TWO)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i) = s(i)
          end if
       end do

    else
       do i=lo(1)-1,hi(1)
          sigma = abs(u(i))*dt/dx(1)
          s6 = SIX*s(i) - THREE*(sm(i)+sp(i))
          if (u(i) .gt. rel_eps) then
             Ip(i) = sp(i) - (sigma/TWO)*(sp(i)-sm(i)-(ONE-TWO3RD*sigma)*s6)
          else
             Ip(i) = s(i)
          end if
       end do
       do i=lo(1),hi(1)+1
          sigma = abs(u(i))*dt/dx(1)
          s6 = SIX*s(i) - THREE*(sm(i)+sp(i))
          if (u(i) .lt. -rel_eps) then
             Im(i) = sm(i) + (sigma/TWO)*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i) = s(i)
          end if
       end do
    endif

    deallocate(sp)
    deallocate(sm)
    deallocate(dsvl)
    deallocate(sedge)

  end subroutine ppm_1d

#elif (AMREX_SPACEDIM==2)
  !===========================================================================
  ! 2-d version
  !===========================================================================

  subroutine ppm_2d(lo,hi,s,s_lo,s_hi,nc_s,&
       u,u_lo,u_hi,v,v_lo,v_hi,&
       Ip,ip_lo,ip_hi,ip_dim,Im,im_lo,im_hi,im_dim,&
       domlo,domhi,adv_bc,&
       dx,dt,is_umac,start_comp,ncomp,start_bccomp,nbccomp) bind(C,name="ppm_2d")

    implicit none

    ! note that u,v here may be the normal cell-centered velocity,
    ! or the MAC velocity.  The is_umac argument tells us which it
    ! is.

    integer         , intent(in   ) :: domlo(3),domhi(3),lo(3),hi(3),s_lo(3),s_hi(3)
    integer         , intent(in   ) :: u_lo(3),u_hi(3),v_lo(3),v_hi(3)
    integer  , value, intent(in   ) :: nc_s,ip_dim, im_dim
    integer         , intent(in   ) :: im_lo(3),im_hi(3),ip_lo(3),ip_hi(3)
    double precision, intent(in   ) ::  s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) ::  u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) ::  v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),ip_dim,AMREX_SPACEDIM)
    double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),im_dim,AMREX_SPACEDIM)
    integer         , intent(in   ) :: adv_bc(AMREX_SPACEDIM,2,nbccomp)
    double precision, intent(in   ) :: dx(3)
    double precision, value, intent(in   ) :: dt
    integer, value, intent(in   ) :: is_umac,start_comp,ncomp,start_bccomp,nbccomp

    ! local
    integer :: i,j,k,n,m,bccomp

    logical :: extremum, bigp, bigm

    double precision :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision :: sgn, sigma, s6, amax, delam, delap, D2ABS
    double precision :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
    double precision :: dsvl_l, dsvl_r, sedge, sm, sp
    double precision :: sedgel, sedger, sedgerr

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    !$gpu

    do n = start_comp,start_comp+ncomp-1
       bccomp = start_bccomp + n - start_comp

       if (ip_dim > 1) then
          m = n
       else
          m = 1
       endif

       !-------------------------------------------------------------------------
       ! x-direction
       !-------------------------------------------------------------------------

       ! compute s at x-edges
       if (ppm_type .eq. 1) then

          ! compute van Leer slopes in x-direction
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)

                   ! sm
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j,k,n) - s(i-2,j,k,n))
                   dsl = TWO  * (s(i-1,j,k,n) - s(i-2,j,k,n))
                   dsr = TWO  * (s(i,j,k,n) - s(i-1,j,k,n))
                   if (dsl*dsr .gt. 0) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i-1,j,k,n))
                   dsr = TWO  * (s(i+1,j,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. 0) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! interpolate s to x-edges
                   sm = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)

                   ! make sure sedge lies in between adjacent cell-centered values
                   sm = max(sm,min(s(i,j,k,n),s(i-1,j,k,n)))
                   sm = min(sm,max(s(i,j,k,n),s(i-1,j,k,n)))

                   ! sp
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i-1,j,k,n))
                   dsr = TWO  * (s(i+1,j,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. 0) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i+2,j,k,n) - s(i,j,k,n))
                   dsl = TWO  * (s(i+1,j,k,n) - s(i,j,k,n))
                   dsr = TWO  * (s(i+2,j,k,n) - s(i+1,j,k,n))
                   if (dsl*dsr .gt. 0) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! interpolate s to x-edges
                   sp = HALF*(s(i+1,j,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)

                   ! make sure sedge lies in between adjacent cell-centered values
                   sp = max(sp,min(s(i+1,j,k,n),s(i,j,k,n)))
                   sp = min(sp,max(s(i+1,j,k,n),s(i,j,k,n)))

                   ! modify using quadratic limiters
                   if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                      sp = s(i,j,k,n)
                      sm = s(i,j,k,n)
                   else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                      sp = THREE*s(i,j,k,n) - TWO*sm
                   else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                      sm = THREE*s(i,j,k,n) - TWO*sp
                   end if

                   ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
                   if (i .eq. lo(1)+1 .and. lo(1)+1 .eq. domlo(1)) then
                      if (adv_bc(1,1,bccomp) .eq. EXT_DIR  .or. adv_bc(1,1,bccomp) .eq. HOEXTRAP) then

                         ! make sure sedge lies in between adjacent cell-centered values
                         ! the value in the first cc ghost cell represents the edge value
                         sm = s(i-1,j,k,n)

                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i-1,j,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i+1,j,k,n) &
                              - (ONE/20.0d0)*s(i+2,j,k,n)

                         sedge = max(sedge,min(s(i+1,j,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i+1,j,k,n),s(i,j,k,n)))

                         sp = sedge
                      end if
                   end if

                   if (i .eq. lo(1)+2 .and. lo(1)+1 .eq. domlo(1)) then
                      if (adv_bc(1,1,bccomp) .eq. EXT_DIR  .or. adv_bc(1,1,bccomp) .eq. HOEXTRAP) then

                         ! make sure sedge lies in between adjacent cell-centered values

                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i-2,j,k,n) &
                              + (THREE/FOUR)*s(i-1,j,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i+1,j,k,n)

                         sedge = max(sedge,min(s(i,j,k,n),s(i-1,j,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i-1,j,k,n)))

                         sm = sedge

                         ! modify using quadratic limiters
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if
                      end if
                   end if

                   if (i .eq. hi(1)-1 .and. hi(1)-1 .eq. domhi(1)) then
                      if (adv_bc(1,2,bccomp) .eq. EXT_DIR  .or. adv_bc(1,2,bccomp) .eq. HOEXTRAP) then
                         ! the value in the first cc ghost cell represents the edge value
                         sp = s(i+1,j,k,n)

                         ! make sure sedge lies in between adjacent cell-centered values
                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i+1,j,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i-1,j,k,n) &
                              - (ONE/20.0d0)*s(i-2,j,k,n)

                         sedge = max(sedge,min(s(i-1,j,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i-1,j,k,n),s(i,j,k,n)))

                         sm = sedge

                      end if
                   end if

                   if (i .eq. hi(1)-2 .and. hi(1)-1 .eq. domhi(1)) then
                      if (adv_bc(1,2,bccomp) .eq. EXT_DIR  .or. adv_bc(1,2,bccomp) .eq. HOEXTRAP) then

                         ! make sure sedge lies in between adjacent cell-centered values
                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i+2,j,k,n) &
                              + (THREE/FOUR)*s(i+1,j,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i-1,j,k,n)

                         sedge = max(sedge,min(s(i,j,k,n),s(i+1,j,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i+1,j,k,n)))

                         ! copy sedge into sp and sm
                         sp = sedge

                         ! modify using quadratic limiters
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if
                      end if
                   end if


                   !-------------------------------------------------------------------------
                   ! compute x-component of Ip and Im
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! u here is umac, so use edge-based indexing
                      sigma = abs(u(i+1,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i+1,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if

                   else
                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if
                   endif
                end do
             end do
          end do


       else if (ppm_type .eq. 2) then

          !----------------------------------------------------------------------
          ! ppm_type = 2
          !----------------------------------------------------------------------
          ! #ifndef AMREX_USE_GPU
          !        if (ng_s .lt. 4) then
          !           call amrex_error("Need 4 ghost cells for ppm_type=2")
          !        end if
          ! #endif

          ! interpolate s to x-edges
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ! -1
                   sedgel = (7.d0/12.d0)*(s(i-2,j,k,n)+s(i-1,j,k,n)) - (1.d0/12.d0)*(s(i-3,j,k,n)+s(i,j,k,n))
                   ! limit sedge
                   if ((sedgel-s(i-2,j,k,n))*(s(i-1,j,k,n)-sedgel) .lt. ZERO) then
                      D2  = THREE*(s(i-2,j,k,n)-TWO*sedgel+s(i-1,j,k,n))
                      D2L = s(i-3,j,k,n)-TWO*s(i-2,j,k,n)+s(i-1,j,k,n)
                      D2R = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgel = HALF*(s(i-2,j,k,n)+s(i-1,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! 0
                   sedge = (7.d0/12.d0)*(s(i-1,j,k,n)+s(i,j,k,n)) - (1.d0/12.d0)*(s(i-2,j,k,n)+s(i+1,j,k,n))
                   ! limit sedge
                   if ((sedge-s(i-1,j,k,n))*(s(i,j,k,n)-sedge) .lt. ZERO) then
                      D2  = THREE*(s(i-1,j,k,n)-TWO*sedge+s(i,j,k,n))
                      D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                      D2R = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedge = HALF*(s(i-1,j,k,n)+s(i,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +1
                   sedger = (7.d0/12.d0)*(s(i,j,k,n)+s(i+1,j,k,n)) - (1.d0/12.d0)*(s(i-1,j,k,n)+s(i+2,j,k,n))
                   ! limit sedge
                   if ((sedger-s(i,j,k,n))*(s(i+1,j,k,n)-sedger) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k,n)-TWO*sedger+s(i+1,j,k,n))
                      D2L = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedger = HALF*(s(i,j,k,n)+s(i+1,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +2
                   sedgerr = (7.d0/12.d0)*(s(i+1,j,k,n)+s(i+2,j,k,n)) - (1.d0/12.d0)*(s(i,j,k,n)+s(i+3,j,k,n))
                   ! limit sedge
                   if ((sedgerr-s(i+1,j,k,n))*(s(i+2,j,k,n)-sedgerr) .lt. ZERO) then
                      D2  = THREE*(s(i+1,j,k,n)-TWO*sedgerr+s(i+2,j,k,n))
                      D2L = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                      D2R = s(i+1,j,k,n)-TWO*s(i+2,j,k,n)+s(i+3,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgerr = HALF*(s(i+1,j,k,n)+s(i+2,j,k,n)) - SIXTH*D2LIM
                   end if
                   !
                   ! ! use Colella 2008 limiters
                   ! ! This is a new version of the algorithm
                   ! ! to eliminate sensitivity to roundoff.
                   !  do i=lo(1)-1,hi(1)+1

                   alphap = sedger-s(i,j,k,n)
                   alpham = sedge-s(i,j,k,n)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is the
                      ! largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge - sedgel
                      dafacep = sedgerr - sedger
                      dabarm = s(i,j,k,n) - s(i-1,j,k,n)
                      dabarp = s(i+1,j,k,n) - s(i,j,k,n)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                      D2C = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i+1,j,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm = s(i,j,k,n) + alpham
                   sp = s(i,j,k,n) + alphap


                   ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
                   ! if (i .eq. lo(1) .and. lo(1) .eq. domlo(1)) then
                   if (lo(1)+1 .eq. domlo(1)) then
                      if (adv_bc(1,1,bccomp) .eq. EXT_DIR  .or. adv_bc(1,1,bccomp) .eq. HOEXTRAP) then

                         if (i .eq. lo(1)+1) then

                            ! the value in the first cc ghost cell represents the edge value
                            sm    = s(i-1,j,k,n)
                            sedge = s(i-1,j,k,n)

                            ! use a modified stencil to get sedge on the first interior edge
                            sedger = &
                                 -FIFTH        *s(i-1,j,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i+1,j,k,n) &
                                 - (ONE/20.0d0)*s(i+2,j,k,n)

                            sedger = max(sedger,min(s(i+1,j,k,n),s(i,j,k,n)))
                            sedger = min(sedger,max(s(i+1,j,k,n),s(i,j,k,n)))

                            sp = sedger

                         elseif (i .eq. lo(1)+2) then

                            ! use a modified stencil to get sedge on the first interior edge
                            sedge = &
                                 -FIFTH        *s(i-2,j,k,n) &
                                 + (THREE/FOUR)*s(i-1,j,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i+1,j,k,n)

                            sedge = max(sedge,min(s(i,j,k,n),s(i-1,j,k,n)))
                            sedge = min(sedge,max(s(i,j,k,n),s(i-1,j,k,n)))

                         elseif (i .eq. lo(1)+3) then

                            ! use a modified stencil to get sedge on the first interior edge
                            sedgel = &
                                 -FIFTH        *s(i-3,j,k,n) &
                                 + (THREE/FOUR)*s(i-2,j,k,n) &
                                 + HALF        *s(i-1,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)

                            sedgel = max(sedgel,min(s(i-1,j,k,n),s(i-2,j,k,n)))
                            sedgel = min(sedgel,max(s(i-1,j,k,n),s(i-2,j,k,n)))

                         endif

                         ! apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells

                         if (i .eq. lo(1)+2 .or. i .eq. lo(1)+3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is the
                               ! largest, and thus least susceptible to sensitivity to roundoff.
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedge
                               dabarm = s(i,j,k,n) - s(i-1,j,k,n)
                               dabarp = s(i+1,j,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                               D2C = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i-1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i+1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if
                      end if
                   end if

                   if (hi(1)-1 .eq. domhi(1)) then
                      if (adv_bc(1,2,bccomp) .eq. EXT_DIR  .or. adv_bc(1,2,bccomp) .eq. HOEXTRAP) then

                         if (i .eq. hi(1)-1) then

                            ! the value in the first cc ghost cell represents the edge value
                            sp      = s(i+1,j,k,n)
                            sedge = s(i+1,j,k,n)

                            ! use a modified stencil to get sedge on the first interior edge
                            sedge = &
                                 -FIFTH        *s(i+1,j,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i-1,j,k,n) &
                                 - (ONE/20.0d0)*s(i-2,j,k,n)

                            sedge = max(sedge,min(s(i-1,j,k,n),s(i,j,k,n)))
                            sedge = min(sedge,max(s(i-1,j,k,n),s(i,j,k,n)))

                            sm = sedge

                         elseif (i .eq. hi(1)-2) then

                            sedgerr = s(i+2,j,k,n)

                            sedger = &
                                 -FIFTH        *s(i+2,j,k,n) &
                                 + (THREE/FOUR)*s(i+1,j,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i-1,j,k,n)

                            sedger = max(sedger,min(s(i,j,k,n),s(i-1,j,k,n)))
                            sedger = min(sedger,max(s(i,j,k,n),s(i-1,j,k,n)))

                         elseif (i .eq. hi(1)-3) then

                            sedgerr = &
                                 -FIFTH        *s(i+3,j,k,n) &
                                 + (THREE/FOUR)*s(i+2,j,k,n) &
                                 + HALF        *s(i+1,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)

                            sedgerr = max(sedgerr,min(s(i+1,j,k,n),s(i+2,j,k,n)))
                            sedgerr = min(sedgerr,max(s(i+1,j,k,n),s(i+2,j,k,n)))

                         endif

                         !
                         ! ! apply Colella 2008 limiters to compute sm and sp in the second
                         ! ! and third inner cells

                         if (i .eq. hi(1)-3 .or. i .eq. hi(1)-2) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is the
                               ! largest, and thus least susceptible to sensitivity to roundoff.
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i-1,j,k,n)
                               dabarp = s(i+1,j,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                               D2C = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i-1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i+1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if
                      end if
                   end if
                   !-------------------------------------------------------------------------
                   ! compute x-component of Ip and Im
                   !-------------------------------------------------------------------------
                   if (is_umac == 1) then

                      ! u here is umac, so use edge-based indexing
                      sigma = abs(u(i+1,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i+1,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if

                   else
                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if
                   endif
                end do
             end do
          end do

       end if


       !-------------------------------------------------------------------------
       ! y-direction
       !-------------------------------------------------------------------------

       ! compute s at y-edges
       if (ppm_type .eq. 1) then

          !----------------------------------------------------------------------
          ! ppm_type = 1
          !----------------------------------------------------------------------

          ! compute van Leer slopes in y-direction
          ! dsvl = ZERO
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ! sm
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j,k,n) - s(i,j-2,k,n))
                   dsl = TWO  * (s(i,j-1,k,n) - s(i,j-2,k,n))
                   dsr = TWO  * (s(i,j,k,n) - s(i,j-1,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i,j-1,k,n))
                   dsr = TWO  * (s(i,j+1,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   sm = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   ! make sure sedge lies in between adjacent cell-centered values
                   sm = max(sm,min(s(i,j,k,n),s(i,j-1,k,n)))
                   sm = min(sm,max(s(i,j,k,n),s(i,j-1,k,n)))

                   ! sp
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i,j-1,k,n))
                   dsr = TWO  * (s(i,j+1,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i,j+2,k,n) - s(i,j,k,n))
                   dsl = TWO  * (s(i,j+1,k,n) - s(i,j,k,n))
                   dsr = TWO  * (s(i,j+2,k,n) - s(i,j+1,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   sp = HALF*(s(i,j+1,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   ! make sure sedge lies in between adjacent cell-centered values
                   sp = max(sp,min(s(i,j+1,k,n),s(i,j,k,n)))
                   sp = min(sp,max(s(i,j+1,k,n),s(i,j,k,n)))

                   ! modify using quadratic limiters
                   if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                      sp = s(i,j,k,n)
                      sm = s(i,j,k,n)
                   else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                      sp = THREE*s(i,j,k,n) - TWO*sm
                   else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                      sm = THREE*s(i,j,k,n) - TWO*sp
                   end if

                   ! different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's
                   if (j .eq. lo(2)+1 .and. lo(2)+1 .eq. domlo(2)) then
                      if (adv_bc(2,1,bccomp) .eq. EXT_DIR  .or. adv_bc(2,1,bccomp) .eq. HOEXTRAP) then

                         ! the value in the first cc ghost cell represents the edge value
                         sm = s(i,j-1,k,n)

                         ! use a modified stencil to get sedge on the first interior edge
                         sp = &
                              -FIFTH        *s(i,j-1,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i,j+1,k,n) &
                              - (ONE/20.0d0)*s(i,j+2,k,n)

                         sp = max(sp,min(s(i,j+1,k,n),s(i,j,k,n)))
                         sp = min(sp,max(s(i,j+1,k,n),s(i,j,k,n)))
                      end if
                   end if

                   ! different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's
                   if (j .eq. lo(2)+2 .and. lo(2)+1 .eq. domlo(2)) then
                      if (adv_bc(2,1,bccomp) .eq. EXT_DIR  .or. adv_bc(2,1,bccomp) .eq. HOEXTRAP) then

                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i,j-2,k,n) &
                              + (THREE/FOUR)*s(i,j-1,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i,j+1,k,n)

                         sedge = max(sedge,min(s(i,j,k,n),s(i,j-1,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i,j-1,k,n)))

                         ! copy sedge into sp and sm
                         sm = sedge

                         ! modify using quadratic limiters
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if
                         ! end do
                      end if
                   end if

                   if (j .eq. hi(2)-1 .and. hi(2)-1 .eq. domhi(2)) then
                      if (adv_bc(2,2,bccomp) .eq. EXT_DIR  .or. adv_bc(2,2,bccomp) .eq. HOEXTRAP) then

                         ! the value in the first cc ghost cell represents the edge value
                         sp = s(i,j+1,k,n)

                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i,j+1,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i,j-1,k,n) &
                              - (ONE/20.0d0)*s(i,j-2,k,n)

                         sedge = max(sedge,min(s(i,j-1,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i,j-1,k,n),s(i,j,k,n)))
                         !
                         ! ! copy sedge into sp and sm
                         sm = sedge

                      end if
                   end if

                   if (j .eq. hi(2)-2 .and. hi(2)-1 .eq. domhi(2)) then
                      if (adv_bc(2,2,bccomp) .eq. EXT_DIR  .or. adv_bc(2,2,bccomp) .eq. HOEXTRAP) then

                         ! use a modified stencil to get sedge on the first interior edge
                         sedge = &
                              -FIFTH        *s(i,j+2,k,n) &
                              + (THREE/FOUR)*s(i,j+1,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i,j-1,k,n)

                         sedge = max(sedge,min(s(i,j,k,n),s(i,j+1,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i,j+1,k,n)))
                         !
                         ! ! copy sedge into sp and sm
                         sp = sedge
                         !
                         ! ! modify using quadratic limiters
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if
                         ! end do
                      end if
                   end if

                   !-------------------------------------------------------------------------
                   ! compute y-component of Ip and Im
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! v here is vmac, so use edge-based indexing
                      sigma = abs(v(i,j+1,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j+1,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if

                   else
                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if
                   endif

                end do
             end do
          end do

       else if (ppm_type .eq. 2) then

          !----------------------------------------------------------------------
          ! ppm_type = 2
          !---------------------------------------------------------

          ! interpolate s to y-edges
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ! -1
                   sedgel = (7.d0/12.d0)*(s(i,j-2,k,n)+s(i,j-1,k,n)) - (1.d0/12.d0)*(s(i,j-3,k,n)+s(i,j,k,n))
                   ! limit sedge
                   if ((sedgel-s(i,j-2,k,n))*(s(i,j-1,k,n)-sedgel) .lt. ZERO) then
                      D2  = THREE*(s(i,j-2,k,n)-TWO*sedgel+s(i,j-1,k,n))
                      D2L = s(i,j-3,k,n)-TWO*s(i,j-2,k,n)+s(i,j-1,k,n)
                      D2R = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgel = HALF*(s(i,j-2,k,n)+s(i,j-1,k,n)) - SIXTH*D2LIM
                   end if

                   ! 0
                   sedge = (7.d0/12.d0)*(s(i,j-1,k,n)+s(i,j,k,n)) - (1.d0/12.d0)*(s(i,j-2,k,n)+s(i,j+1,k,n))
                   ! limit sedge
                   if ((sedge-s(i,j-1,k,n))*(s(i,j,k,n)-sedge) .lt. ZERO) then
                      D2  = THREE*(s(i,j-1,k,n)-TWO*sedge+s(i,j,k,n))
                      D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                      D2R = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedge = HALF*(s(i,j-1,k,n)+s(i,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +1
                   sedger = (7.d0/12.d0)*(s(i,j,k,n)+s(i,j+1,k,n)) - (1.d0/12.d0)*(s(i,j-1,k,n)+s(i,j+2,k,n))
                   ! limit sedge
                   if ((sedger-s(i,j,k,n))*(s(i,j+1,k,n)-sedger) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k,n)-TWO*sedger+s(i,j+1,k,n))
                      D2L = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedger = HALF*(s(i,j,k,n)+s(i,j+1,k,n)) - SIXTH*D2LIM
                   end if

                   ! +2
                   sedgerr = (7.d0/12.d0)*(s(i,j+1,k,n)+s(i,j+2,k,n)) - (1.d0/12.d0)*(s(i,j,k,n)+s(i,j+3,k,n))
                   ! limit sedge
                   if ((sedgerr-s(i,j+1,k,n))*(s(i,j+2,k,n)-sedgerr) .lt. ZERO) then
                      D2  = THREE*(s(i,j+1,k,n)-TWO*sedgerr+s(i,j+2,k,n))
                      D2L = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                      D2R = s(i,j+1,k,n)-TWO*s(i,j+2,k,n)+s(i,j+3,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgerr = HALF*(s(i,j+1,k,n)+s(i,j+2,k,n)) - SIXTH*D2LIM
                   end if

                   ! use Colella 2008 limiters
                   ! This is a new version of the algorithm
                   ! to eliminate sensitivity to roundoff.

                   alphap = sedger-s(i,j,k,n)
                   alpham = sedge-s(i,j,k,n)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is the
                      ! largest, and thus least susceptible to sensitivity to roundoff.
                      dafacem = sedge - sedgel
                      dafacep = sedgerr - sedger
                      dabarm = s(i,j,k,n) - s(i,j-1,k,n)
                      dabarp = s(i,j+1,k,n) - s(i,j,k,n)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                      D2C = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j+1,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm = s(i,j,k,n) + alpham
                   sp = s(i,j,k,n) + alphap

                   ! end do
                   ! end do

                   ! different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's
                   if (lo(2)+1 .eq. domlo(2)) then
                      if (adv_bc(2,1,bccomp) .eq. EXT_DIR  .or. adv_bc(2,1,bccomp) .eq. HOEXTRAP) then

                         if (j .eq. lo(2)+1) then
                            ! the value in the first cc ghost cell represents the edge value
                            sm    = s(i,j-1,k,n)

                            ! use a modified stencil to get sedge on the first interior edge
                            sp = &
                                 -FIFTH        *s(i,j-1,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i,j+1,k,n) &
                                 - (ONE/20.0d0)*s(i,j+2,k,n)

                            sp = max(sp,min(s(i,j+1,k,n),s(i,j,k,n)))
                            sp = min(sp,max(s(i,j+1,k,n),s(i,j,k,n)))

                         elseif (j .eq. lo(2)+2) then

                            sedgel = s(i,j-2,k,n)

                            sedge = &
                                 -FIFTH        *s(i,j-2,k,n) &
                                 + (THREE/FOUR)*s(i,j-1,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j+1,k,n)

                            sedge = max(sedge,min(s(i,j,k,n),s(i,j-1,k,n)))
                            sedge = min(sedge,max(s(i,j,k,n),s(i,j-1,k,n)))

                         elseif (j .eq. lo(2)+3) then

                            ! use a modified stencil to get sedge on the first interior edge
                            sedgel = &
                                 -FIFTH        *s(i,j-3,k,n) &
                                 + (THREE/FOUR)*s(i,j-2,k,n) &
                                 + HALF        *s(i,j-1,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)

                            sedgel = max(sedgel,min(s(i,j-1,k,n),s(i,j-2,k,n)))
                            sedgel = min(sedgel,max(s(i,j-1,k,n),s(i,j-2,k,n)))

                         endif

                         ! apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells

                         if (j .eq. lo(2)+2 .or. j .eq. lo(2)+3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is the
                               ! largest, and thus least susceptible to sensitivity to roundoff.
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i,j-1,k,n)
                               dabarp = s(i,j+1,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                               D2C = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i,j-1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i,j+1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap
                         end if
                      end if
                   end if

                   if (hi(2)-1 .eq. domhi(2)) then
                      if (adv_bc(2,2,bccomp) .eq. EXT_DIR  .or. adv_bc(2,2,bccomp) .eq. HOEXTRAP) then

                         if (j .eq. hi(2)-1) then

                            sp     = s(i,j+1,k,n)

                            ! use a modified stencil to get sedge on the first interior edge
                            sm = &
                                 -FIFTH        *s(i,j+1,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i,j-1,k,n) &
                                 - (ONE/20.0d0)*s(i,j-2,k,n)

                            sm = max(sm,min(s(i,j-1,k,n),s(i,j,k,n)))
                            sm = min(sm,max(s(i,j-1,k,n),s(i,j,k,n)))

                         elseif (j .eq. hi(2)-2) then

                            sedgerr = s(i,j+2,k,n)

                            sedger = &
                                 -FIFTH        *s(i,j+2,k,n) &
                                 + (THREE/FOUR)*s(i,j+1,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j-1,k,n)

                            sedger = max(sedger,min(s(i,j,k,n),s(i,j+1,k,n)))
                            sedger = min(sedger,max(s(i,j,k,n),s(i,j+1,k,n)))

                         elseif (j .eq. hi(2)-3) then

                            sedgerr = &
                                 -FIFTH        *s(i,j+3,k,n) &
                                 + (THREE/FOUR)*s(i,j+2,k,n) &
                                 + HALF        *s(i,j+1,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)

                            sedgerr = max(sedgerr,min(s(i,j+1,k,n),s(i,j+2,k,n)))
                            sedgerr = min(sedgerr,max(s(i,j+1,k,n),s(i,j+2,k,n)))

                         endif

                         ! apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells
                         if (j .eq. hi(2)-3 .or. j .eq. hi(2)-2) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is the
                               ! largest, and thus least susceptible to sensitivity to roundoff.
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i,j-1,k,n)
                               dabarp = s(i,j+1,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                               D2C = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i,j-1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i,j+1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap
                         end if
                      end if
                   end if


                   !-------------------------------------------------------------------------
                   ! compute y-component of Ip and Im
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! v here is vmac, so use edge-based indexing

                      sigma = abs(v(i,j+1,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j+1,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if
                   else
                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if
                   endif

                end do
             end do
          end do

       end if
    end do

  end subroutine ppm_2d

#elif (AMREX_SPACEDIM == 3)
  !===========================================================================
  ! 3-d version
  !===========================================================================

  ! characteristics based on u
  subroutine ppm_3d(lo,hi,s,s_lo,s_hi,nc_s, &
       u,u_lo,u_hi,v,v_lo,v_hi,w,w_lo,w_hi, &
       Ip,ip_lo,ip_hi,ip_dim,Im,im_lo,im_hi,im_dim, &
       domlo,domhi, &
       adv_bc,dx,dt,is_umac,start_comp,ncomp,start_bccomp,nbccomp) bind(C,name="ppm_3d")

    integer         , intent(in   ) :: domlo(3),domhi(3),lo(3),hi(3),s_lo(3),s_hi(3)
    integer         , intent(in   ) :: u_lo(3),u_hi(3),v_lo(3),v_hi(3),w_lo(3),w_hi(3)
    integer         , intent(in   ) :: im_lo(3),im_hi(3),ip_lo(3),ip_hi(3)
    integer,   value, intent(in   ) :: nc_s,ip_dim,im_dim
    double precision, intent(in   ) ::  s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) ::  u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) ::  v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent(in   ) ::  w(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(inout) :: Ip(ip_lo(1):ip_hi(1),ip_lo(2):ip_hi(2),ip_lo(3):ip_hi(3),ip_dim,AMREX_SPACEDIM)
    double precision, intent(inout) :: Im(im_lo(1):im_hi(1),im_lo(2):im_hi(2),im_lo(3):im_hi(3),im_dim,AMREX_SPACEDIM)
    integer         , intent(in   ) :: adv_bc(AMREX_SPACEDIM,2,AMREX_SPACEDIM)
    double precision, intent(in   ) :: dx(3)
    double precision, value, intent(in   ) :: dt
    integer,   value, intent(in   ) :: is_umac,start_comp,ncomp,start_bccomp,nbccomp

    ! local
    integer :: i,j,k,n,m,bccomp

    logical :: extremum, bigp, bigm

    double precision :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision :: sgn, sigma, s6, D2ABS
    double precision :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
    double precision :: amax, delam, delap
    double precision :: dsvl_l, dsvl_r, sedge, sm, sp
    double precision :: sedgel, sedger, sedgerr

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    !$gpu

    do n = start_comp,start_comp+ncomp-1
       bccomp = start_bccomp + n - start_comp

       if (ip_dim > 1) then
          m = n
       else
          m = 1
       endif

       !-------------------------------------------------------------------------
       ! x-direction
       !-------------------------------------------------------------------------

       ! compute s at x-edges
       if (ppm_type .eq. 1) then

          !----------------------------------------------------------------------
          ! ppm_type = 1
          !----------------------------------------------------------------------

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   !
                   ! Compute van Leer slopes in x-direction.
                   !

                   ! sm
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j,k,n) - s(i-2,j,k,n))
                   dsl = TWO  * (s(i-1,j,k,n) - s(i-2,j,k,n))
                   dsr = TWO  * (s(i,j,k,n) - s(i-1,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i-1,j,k,n))
                   dsr = TWO  * (s(i+1,j,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   !
                   ! Interpolate s to x-edges.
                   !
                   sm = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sm = max(sm,min(s(i,j,k,n),s(i-1,j,k,n)))
                   sm = min(sm,max(s(i,j,k,n),s(i-1,j,k,n)))

                   ! sp
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i-1,j,k,n))
                   dsr = TWO  * (s(i+1,j,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i+2,j,k,n) - s(i,j,k,n))
                   dsl = TWO  * (s(i+1,j,k,n) - s(i,j,k,n))
                   dsr = TWO  * (s(i+2,j,k,n) - s(i+1,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   !
                   ! Interpolate s to x-edges.
                   !
                   sp = HALF*(s(i+1,j,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sp = max(sp,min(s(i+1,j,k,n),s(i,j,k,n)))
                   sp = min(sp,max(s(i+1,j,k,n),s(i,j,k,n)))

                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                      sp = s(i,j,k,n)
                      sm = s(i,j,k,n)
                   else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                      sp = THREE*s(i,j,k,n) - TWO*sm
                   else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                      sm = THREE*s(i,j,k,n) - TWO*sp
                   end if

                   ! Different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's.
                   !
                   if (i .eq. lo(1)+1 .and. lo(1)+1 .eq. domlo(1)) then
                      if (adv_bc(1,1,bccomp) .eq. EXT_DIR  .or. adv_bc(1,1,bccomp) .eq. HOEXTRAP) then

                         !
                         ! The value in the first cc ghost cell represents the edge value.
                         !
                         sm = s(i-1,j,k,n)
                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i-1,j,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i+1,j,k,n) &
                              - (ONE/20.0d0)*s(i+2,j,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i+1,j,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i+1,j,k,n),s(i,j,k,n)))
                         !
                         ! Copy sedge into sp and sm.
                         !
                         sp = sedge

                      end if
                   end if

                   if (i .eq. lo(1)+2 .and. lo(1)+1 .eq. domlo(1)) then
                      if (adv_bc(1,1,bccomp) .eq. EXT_DIR  .or. adv_bc(1,1,bccomp) .eq. HOEXTRAP) then

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i-2,j,k,n) &
                              + (THREE/FOUR)*s(i-1,j,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i+1,j,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k,n),s(i-1,j,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i-1,j,k,n)))

                         sm = sedge

                         !
                         ! Modify using quadratic limiters.
                         !
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if

                      end if
                   end if

                   if (i .eq. hi(1)-1 .and. hi(1)-1 .eq. domhi(1)) then
                      if (adv_bc(1,2,bccomp) .eq. EXT_DIR  .or. adv_bc(1,2,bccomp) .eq. HOEXTRAP) then

                         ! The value in the first cc ghost cell represents the edge value.
                         !
                         sp = s(i+1,j,k,n)

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i+1,j,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i-1,j,k,n) &
                              - (ONE/20.0d0)*s(i-2,j,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i-1,j,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i-1,j,k,n),s(i,j,k,n)))

                         sm = sedge

                      end if
                   end if

                   if (i .eq. hi(1)-2 .and. hi(1)-1 .eq. domhi(1)) then
                      if (adv_bc(1,2,bccomp) .eq. EXT_DIR  .or. adv_bc(1,2,bccomp) .eq. HOEXTRAP) then

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i+2,j,k,n) &
                              + (THREE/FOUR)*s(i+1,j,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i-1,j,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k,n),s(i+1,j,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i+1,j,k,n)))
                         !
                         ! Copy sedge into sp and sm.
                         !
                         sp = sedge

                         !
                         ! Modify using quadratic limiters.
                         !
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if
                      end if
                   end if

                   !-------------------------------------------------------------------------
                   ! Compute x-component of Ip and Im.
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! u is MAC velocity -- use edge-based indexing
                      sigma = abs(u(i+1,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i+1,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if

                   else

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if

                   endif
                end do
             end do
          end do


       else if (ppm_type .eq. 2) then

          !----------------------------------------------------------------------
          ! ppm_type = 2
          !----------------------------------------------------------------------
          !
          ! if (ng_s .lt. 4) then
          !    call amrex_error("Need 4 ghost cells for ppm_type=2")
          ! end if

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ! -1
                   ! Interpolate s to x-edges.
                   sedgel = (7.d0/12.d0)*(s(i-2,j,k,n)+s(i-1,j,k,n)) &
                        - (1.d0/12.d0)*(s(i-3,j,k,n)+s(i,j,k,n))

                   ! Limit sedge.
                   if ((sedgel-s(i-2,j,k,n))*(s(i-1,j,k,n)-sedgel) .lt. ZERO) then
                      D2  = THREE*(s(i-2,j,k,n)-TWO*sedgel+s(i-1,j,k,n))
                      D2L = s(i-3,j,k,n)-TWO*s(i-2,j,k,n)+s(i-1,j,k,n)
                      D2R = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgel = HALF*(s(i-2,j,k,n)+s(i-1,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! 0
                   ! Interpolate s to x-edges.
                   sedge = (7.d0/12.d0)*(s(i-1,j,k,n)+s(i,j,k,n)) &
                        - (1.d0/12.d0)*(s(i-2,j,k,n)+s(i+1,j,k,n))

                   ! Limit sedge.
                   if ((sedge-s(i-1,j,k,n))*(s(i,j,k,n)-sedge) .lt. ZERO) then
                      D2  = THREE*(s(i-1,j,k,n)-TWO*sedge+s(i,j,k,n))
                      D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                      D2R = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedge = HALF*(s(i-1,j,k,n)+s(i,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +1
                   ! Interpolate s to x-edges.
                   sedger = (7.d0/12.d0)*(s(i,j,k,n)+s(i+1,j,k,n)) &
                        - (1.d0/12.d0)*(s(i-1,j,k,n)+s(i+2,j,k,n))

                   ! Limit sedge.
                   if ((sedger-s(i,j,k,n))*(s(i+1,j,k,n)-sedger) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k,n)-TWO*sedger+s(i+1,j,k,n))
                      D2L = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedger = HALF*(s(i,j,k,n)+s(i+1,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +2
                   ! Interpolate s to x-edges.
                   sedgerr = (7.d0/12.d0)*(s(i+1,j,k,n)+s(i+2,j,k,n)) &
                        - (1.d0/12.d0)*(s(i,j,k,n)+s(i+3,j,k,n))

                   ! Limit sedge.
                   if ((sedgerr-s(i+1,j,k,n))*(s(i+2,j,k,n)-sedgerr) .lt. ZERO) then
                      D2  = THREE*(s(i+1,j,k,n)-TWO*sedgerr+s(i+2,j,k,n))
                      D2L = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                      D2R = s(i+1,j,k,n)-TWO*s(i+2,j,k,n)+s(i+3,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgerr = HALF*(s(i+1,j,k,n)+s(i+2,j,k,n)) - SIXTH*D2LIM
                   end if

                   alphap = sedger-s(i,j,k,n)
                   alpham = sedge-s(i,j,k,n)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      !
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is the
                      ! largest, and thus least susceptible to sensitivity to roundoff.
                      !
                      dafacem = sedge - sedgel
                      dafacep = sedgerr - sedger
                      dabarm = s(i,j,k,n) - s(i-1,j,k,n)
                      dabarp = s(i+1,j,k,n) - s(i,j,k,n)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                      D2C = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i+1,j,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm = s(i,j,k,n) + alpham
                   sp = s(i,j,k,n) + alphap

                   ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
                   if (lo(1)+1 .eq. domlo(1)) then
                      if (adv_bc(1,1,bccomp) .eq. EXT_DIR  .or. adv_bc(1,1,bccomp) .eq. HOEXTRAP) then

                         if (i .eq. lo(1)+1) then
                            !
                            ! The value in the first cc ghost cell represents the edge value.
                            !
                            sm    = s(i-1,j,k,n)

                            ! use a modified stencil to get sedge on the first interior edge
                            sp = -FIFTH        *s(i-1,j,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i+1,j,k,n) &
                                 - (ONE/20.0d0)*s(i+2,j,k,n)

                            ! make sure sedge lies in between adjacent cell-centered values
                            sp = max(sp,min(s(i+1,j,k,n),s(i,j,k,n)))
                            sp = min(sp,max(s(i+1,j,k,n),s(i,j,k,n)))

                         elseif (i .eq. lo(1)+2) then

                            sedgel = s(i-2,j,k,n)

                            ! use a modified stencil to get sedge on the first interior edge
                            sedge = -FIFTH     *s(i-2,j,k,n) &
                                 + (THREE/FOUR)*s(i-1,j,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i+1,j,k,n)

                            ! make sure sedge lies in between adjacent cell-centered values
                            sedge = max(sedge,min(s(i,j,k,n),s(i-1,j,k,n)))
                            sedge = min(sedge,max(s(i,j,k,n),s(i-1,j,k,n)))

                         elseif (i .eq. lo(1)+3) then

                            ! use a modified stencil to get sedge on the first interior edge
                            sedgel = -FIFTH    *s(i-3,j,k,n) &
                                 + (THREE/FOUR)*s(i-2,j,k,n) &
                                 + HALF        *s(i-1,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)

                            ! make sure sedge lies in between adjacent cell-centered values
                            sedgel = max(sedgel,min(s(i-1,j,k,n),s(i-2,j,k,n)))
                            sedgel = min(sedgel,max(s(i-1,j,k,n),s(i-2,j,k,n)))

                         endif

                         !
                         ! Apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells.

                         if (i .eq. lo(1)+2 .or. i .eq. lo(1)+3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is
                               ! the largest, and thus least susceptible to sensitivity to roundoff.
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i-1,j,k,n)
                               dabarp = s(i+1,j,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                               D2C = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i-1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i+1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if
                      end if
                   end if

                   if (hi(1)-1 .eq. domhi(1)) then
                      if (adv_bc(1,2,bccomp) .eq. EXT_DIR  .or. adv_bc(1,2,bccomp) .eq. HOEXTRAP) then

                         if (i .eq. hi(1)-1) then

                            !
                            ! The value in the first cc ghost cell represents the edge value.
                            !
                            sp = s(i+1,j,k,n)

                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sm = -FIFTH        *s(i+1,j,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i-1,j,k,n) &
                                 - (ONE/20.0d0)*s(i-2,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sm = max(sm,min(s(i-1,j,k,n),s(i,j,k,n)))
                            sm = min(sm,max(s(i-1,j,k,n),s(i,j,k,n)))


                         elseif (i .eq. hi(1)-2) then

                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedger = -FIFTH    *s(i+2,j,k,n) &
                                 + (THREE/FOUR)*s(i+1,j,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i-1,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedger = max(sedger,min(s(i,j,k,n),s(i+1,j,k,n)))
                            sedger = min(sedger,max(s(i,j,k,n),s(i+1,j,k,n)))

                         elseif (i .eq. hi(1)-3) then

                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedgerr = -FIFTH   *s(i+3,j,k,n) &
                                 + (THREE/FOUR)*s(i+2,j,k,n) &
                                 + HALF        *s(i+1,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedgerr = max(sedgerr,min(s(i+1,j,k,n),s(i+2,j,k,n)))
                            sedgerr = min(sedgerr,max(s(i+1,j,k,n),s(i+2,j,k,n)))

                         endif

                         !
                         ! Apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells.

                         if (i .eq. hi(1)-2 .or. i .eq. hi(1)-3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               !
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is
                               ! the largest, and thus least susceptible to sensitivity to roundoff.
                               !
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i-1,j,k,n)
                               dabarp = s(i+1,j,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i-2,j,k,n)-TWO*s(i-1,j,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i+1,j,k,n)+s(i+2,j,k,n)
                               D2C = s(i-1,j,k,n)-TWO*s(i,j,k,n)+s(i+1,j,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i-1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i+1,j,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if

                      end if

                   end if

                   !-------------------------------------------------------------------------
                   ! Compute x-component of Ip and Im.
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! u is MAC velocity -- use edge-based indexing
                      sigma = abs(u(i+1,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i+1,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if


                   else

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,1) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,1) = s(i,j,k,n)
                      end if

                      sigma = abs(u(i,j,k))*dt/dx(1)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (u(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,1) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,1) = s(i,j,k,n)
                      end if

                   endif
                end do
             end do
          end do

       end if


       !-------------------------------------------------------------------------
       ! y-direction
       !-------------------------------------------------------------------------

       !
       ! Compute s at y-edges.
       !
       if (ppm_type .eq. 1) then

          !----------------------------------------------------------------------
          ! ppm_type = 1
          !----------------------------------------------------------------------

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   !
                   ! Compute van Leer slopes in y-direction.
                   !

                   ! sm
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j,k,n) - s(i,j-2,k,n))
                   dsl = TWO  * (s(i,j-1,k,n) - s(i,j-2,k,n))
                   dsr = TWO  * (s(i,j,k,n) - s(i,j-1,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i,j-1,k,n))
                   dsr = TWO  * (s(i,j+1,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! Interpolate s to y-edges.
                   !
                   sm = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sm = max(sm,min(s(i,j,k,n),s(i,j-1,k,n)))
                   sm = min(sm,max(s(i,j,k,n),s(i,j-1,k,n)))

                   ! sp
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i,j-1,k,n))
                   dsr = TWO  * (s(i,j+1,k,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i,j+2,k,n) - s(i,j,k,n))
                   dsl = TWO  * (s(i,j+1,k,n) - s(i,j,k,n))
                   dsr = TWO  * (s(i,j+2,k,n) - s(i,j+1,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! Interpolate s to y-edges.
                   !
                   sp = HALF*(s(i,j+1,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sp = max(sp,min(s(i,j+1,k,n),s(i,j,k,n)))
                   sp = min(sp,max(s(i,j+1,k,n),s(i,j,k,n)))

                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                      sp = s(i,j,k,n)
                      sm = s(i,j,k,n)
                   else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                      sp = THREE*s(i,j,k,n) - TWO*sm
                   else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                      sm = THREE*s(i,j,k,n) - TWO*sp
                   end if
                   !
                   !
                   ! Different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's.
                   !
                   if (j .eq. lo(2)+1 .and. lo(2)+1 .eq. domlo(2)) then
                      if (adv_bc(2,1,bccomp) .eq. EXT_DIR  .or. adv_bc(2,1,bccomp) .eq. HOEXTRAP) then
                         !
                         ! The value in the first cc ghost cell represents the edge value.
                         !
                         sm = s(i,j-1,k,n)
                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j-1,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i,j+1,k,n) &
                              - (ONE/20.0d0)*s(i,j+2,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j+1,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i,j+1,k,n),s(i,j,k,n)))
                         !
                         ! Copy sedge into sp and sm.
                         !
                         sp = sedge

                      end if
                   end if

                   if (j .eq. lo(2)+2 .and. lo(2)+1 .eq. domlo(2)) then
                      if (adv_bc(2,1,bccomp) .eq. EXT_DIR  .or. adv_bc(2,1,bccomp) .eq. HOEXTRAP) then

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j-2,k,n) &
                              + (THREE/FOUR)*s(i,j-1,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i,j+1,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k,n),s(i,j-1,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i,j-1,k,n)))
                         !
                         sm = sedge
                         !
                         ! Modify using quadratic limiters.
                         !
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if

                      end if
                   end if

                   if (j .eq. hi(2)-1 .and. hi(2)-1 .eq. domhi(2)) then
                      if (adv_bc(2,2,bccomp) .eq. EXT_DIR  .or. adv_bc(2,2,bccomp) .eq. HOEXTRAP) then
                         !
                         ! The value in the first cc ghost cell represents the edge value.
                         !
                         sp = s(i,j+1,k,n)
                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j+1,k,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i,j-1,k,n) &
                              - (ONE/20.0d0)*s(i,j-2,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j-1,k,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i,j-1,k,n),s(i,j,k,n)))
                         !
                         sm = sedge

                      end if
                   end if

                   if (j .eq. hi(2)-2 .and. hi(2)-1 .eq. domhi(2)) then
                      if (adv_bc(2,2,bccomp) .eq. EXT_DIR  .or. adv_bc(2,2,bccomp) .eq. HOEXTRAP) then

                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j+2,k,n) &
                              + (THREE/FOUR)*s(i,j+1,k,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i,j-1,k,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k,n),s(i,j+1,k,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i,j+1,k,n)))
                         !
                         ! Copy sedge into sp and sm.
                         !
                         sp = sedge
                         !
                         ! Modify using quadratic limiters.
                         !
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if

                      end if
                   end if

                   !-------------------------------------------------------------------------
                   ! Compute y-component of Ip and Im.
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! v is MAC velocity -- use edge-based indexing

                      sigma = abs(v(i,j+1,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j+1,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if


                   else

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if


                   endif
                end do
             end do
          end do

       else if (ppm_type .eq. 2) then

          !----------------------------------------------------------------------
          ! ppm_type = 2
          !----------------------------------------------------------------------

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ! -1
                   ! Interpolate s to y-edges.
                   sedgel = (7.d0/12.d0)*(s(i,j-2,k,n)+s(i,j-1,k,n)) &
                        - (1.d0/12.d0)*(s(i,j-3,k,n)+s(i,j,k,n))
                   !
                   ! Limit sedge.
                   if ((sedgel-s(i,j-2,k,n))*(s(i,j-1,k,n)-sedgel) .lt. ZERO) then
                      D2  = THREE*(s(i,j-2,k,n)-TWO*sedgel+s(i,j-1,k,n))
                      D2L = s(i,j-3,k,n)-TWO*s(i,j-2,k,n)+s(i,j-1,k,n)
                      D2R = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgel = HALF*(s(i,j-2,k,n)+s(i,j-1,k,n)) - SIXTH*D2LIM
                   end if

                   ! 0
                   ! Interpolate s to y-edges.
                   sedge = (7.d0/12.d0)*(s(i,j-1,k,n)+s(i,j,k,n)) &
                        - (1.d0/12.d0)*(s(i,j-2,k,n)+s(i,j+1,k,n))
                   !
                   ! Limit sedge.
                   if ((sedge-s(i,j-1,k,n))*(s(i,j,k,n)-sedge) .lt. ZERO) then
                      D2  = THREE*(s(i,j-1,k,n)-TWO*sedge+s(i,j,k,n))
                      D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                      D2R = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedge = HALF*(s(i,j-1,k,n)+s(i,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +1
                   ! Interpolate s to y-edges.
                   sedger = (7.d0/12.d0)*(s(i,j,k,n)+s(i,j+1,k,n)) &
                        - (1.d0/12.d0)*(s(i,j-1,k,n)+s(i,j+2,k,n))
                   !
                   ! Limit sedge.
                   if ((sedger-s(i,j,k,n))*(s(i,j+1,k,n)-sedger) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k,n)-TWO*sedger+s(i,j+1,k,n))
                      D2L = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedger = HALF*(s(i,j,k,n)+s(i,j+1,k,n)) - SIXTH*D2LIM
                   end if

                   ! +2
                   ! Interpolate s to y-edges.
                   sedgerr = (7.d0/12.d0)*(s(i,j+1,k,n)+s(i,j+2,k,n)) &
                        - (1.d0/12.d0)*(s(i,j,k,n)+s(i,j+3,k,n))
                   !
                   ! Limit sedge.
                   if ((sedgerr-s(i,j+1,k,n))*(s(i,j+2,k,n)-sedgerr) .lt. ZERO) then
                      D2  = THREE*(s(i,j+1,k,n)-TWO*sedgerr+s(i,j+2,k,n))
                      D2L = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                      D2R = s(i,j+1,k,n)-TWO*s(i,j+2,k,n)+s(i,j+3,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgerr = HALF*(s(i,j+1,k,n)+s(i,j+2,k,n)) - SIXTH*D2LIM
                   end if

                   !
                   ! Use Colella 2008 limiters.
                   ! This is a new version of the algorithm
                   ! to eliminate sensitivity to roundoff.

                   alphap = sedger-s(i,j,k,n)
                   alpham = sedge-s(i,j,k,n)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      !
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is the
                      ! largest, and thus least susceptible to sensitivity to roundoff.
                      !
                      dafacem = sedge - sedgel
                      dafacep = sedgerr - sedger
                      dabarm = s(i,j,k,n) - s(i,j-1,k,n)
                      dabarp = s(i,j+1,k,n) - s(i,j,k,n)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                      D2C = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j+1,k,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm = s(i,j,k,n) + alpham
                   sp = s(i,j,k,n) + alphap

                   !
                   ! Different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's.
                   !
                   if (lo(2)+1 .eq. domlo(2)) then
                      if (adv_bc(2,1,bccomp) .eq. EXT_DIR  .or. adv_bc(2,1,bccomp) .eq. HOEXTRAP) then

                         if (j .eq. lo(2)+1) then

                            !
                            ! The value in the first cc ghost cell represents the edge value.
                            !
                            sm = s(i,j-1,k,n)
                            sedge = s(i,j-1,k,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sp = -FIFTH        *s(i,j-1,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i,j+1,k,n) &
                                 - (ONE/20.0d0)*s(i,j+2,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sp = max(sp,min(s(i,j+1,k,n),s(i,j,k,n)))
                            sp = min(sp,max(s(i,j+1,k,n),s(i,j,k,n)))

                         elseif (j .eq. lo(2)+2) then

                            sedgel = s(i,j-2,k,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedge = -FIFTH     *s(i,j-2,k,n) &
                                 + (THREE/FOUR)*s(i,j-1,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j+1,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedge = max(sedge,min(s(i,j,k,n),s(i,j-1,k,n)))
                            sedge = min(sedge,max(s(i,j,k,n),s(i,j-1,k,n)))

                         elseif (j .eq. lo(2)+3) then
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedgel = -FIFTH    *s(i,j-3,k,n) &
                                 + (THREE/FOUR)*s(i,j-2,k,n) &
                                 + HALF        *s(i,j-1,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedgel = max(sedgel,min(s(i,j-1,k,n),s(i,j-2,k,n)))
                            sedgel = min(sedgel,max(s(i,j-1,k,n),s(i,j-2,k,n)))

                         endif

                         !
                         ! Apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells.

                         if (j .eq. lo(2)+2 .or. j .eq. lo(2)+3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               !
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is
                               ! the largest, and thus least susceptible to sensitivity to roundoff.
                               !
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i,j-1,k,n)
                               dabarp = s(i,j+1,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                               D2C = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i,j-1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i,j+1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if

                      end if
                   end if

                   if (hi(2)-1 .eq. domhi(2)) then
                      if (adv_bc(2,2,bccomp) .eq. EXT_DIR  .or. adv_bc(2,2,bccomp) .eq. HOEXTRAP) then

                         if (j .eq. hi(2)-1) then

                            !
                            ! The value in the first cc ghost cell represents the edge value.
                            !
                            sp = s(i,j+1,k,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sm  = -FIFTH       *s(i,j+1,k,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i,j-1,k,n) &
                                 - (ONE/20.0d0)*s(i,j-2,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sm  = max(sm,min(s(i,j-1,k,n),s(i,j,k,n)))
                            sm = min(sm,max(s(i,j-1,k,n),s(i,j,k,n)))

                         elseif (j .eq. hi(2)-2) then

                            sedgerr = s(i,j+2,k,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedger = -FIFTH    *s(i,j+2,k,n) &
                                 + (THREE/FOUR)*s(i,j+1,k,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j-1,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedger = max(sedger,min(s(i,j,k,n),s(i,j+1,k,n)))
                            sedger = min(sedger,max(s(i,j,k,n),s(i,j+1,k,n)))

                         elseif (j .eq. hi(2)-3) then
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedgerr = -FIFTH   *s(i,j+3,k,n) &
                                 + (THREE/FOUR)*s(i,j+2,k,n) &
                                 + HALF        *s(i,j+1,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedgerr = max(sedgerr,min(s(i,j+1,k,n),s(i,j+2,k,n)))
                            sedgerr = min(sedgerr,max(s(i,j+1,k,n),s(i,j+2,k,n)))

                         endif

                         !
                         ! Apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells.

                         if (j .eq. hi(2)-2 .or. j .eq. hi(2)-3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               !
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is
                               ! the largest, and thus least susceptible to sensitivity to roundoff.
                               !
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i,j-1,k,n)
                               dabarp = s(i,j+1,k,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i,j-2,k,n)-TWO*s(i,j-1,k,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i,j+1,k,n)+s(i,j+2,k,n)
                               D2C = s(i,j-1,k,n)-TWO*s(i,j,k,n)+s(i,j+1,k,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i,j-1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i,j+1,k,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if

                      end if
                   end if

                   !-------------------------------------------------------------------------
                   ! Compute y-component of Ip and Im.
                   !-------------------------------------------------------------------------


                   if (is_umac == 1) then

                      ! v is MAC velocity -- use edge-based indexing

                      sigma = abs(v(i,j+1,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j+1,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if

                   else

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,2) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,2) = s(i,j,k,n)
                      end if

                      sigma = abs(v(i,j,k))*dt/dx(2)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (v(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,2) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,2) = s(i,j,k,n)
                      end if

                   endif
                end do
             end do
          end do

       end if

       !-------------------------------------------------------------------------
       ! z-direction
       !-------------------------------------------------------------------------

       !
       ! Compute s at z-edges.
       !
       if (ppm_type .eq. 1) then

          !----------------------------------------------------------------------
          ! ppm_type = 1
          !----------------------------------------------------------------------

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   !
                   ! Compute van Leer slopes in z-direction.
                   !

                   ! sm
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j,k,n) - s(i,j,k-2,n))
                   dsl = TWO  * (s(i,j,k-1,n) - s(i,j,k-2,n))
                   dsr = TWO  * (s(i,j,k,n) - s(i,j,k-1,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i,j,k-1,n))
                   dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   !
                   ! Interpolate s to z-edges.
                   !
                   sm = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) - SIXTH*(dsvl_r-dsvl_l)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sm = max(sm,min(s(i,j,k,n),s(i,j,k-1,n)))
                   sm = min(sm,max(s(i,j,k,n),s(i,j,k-1,n)))

                   ! sp
                   dsvl_l = ZERO
                   dsvl_r = ZERO

                   ! left side
                   dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                   dsl = TWO  * (s(i,j,k,n) - s(i,j,k-1,n))
                   dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k,n))
                   if (dsl*dsr .gt. ZERO) dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   ! right side
                   dsc = HALF * (s(i,j,k+2,n) - s(i,j,k,n))
                   dsl = TWO  * (s(i,j,k+1,n) - s(i,j,k,n))
                   dsr = TWO  * (s(i,j,k+2,n) - s(i,j,k+1,n))
                   if (dsl*dsr .gt. ZERO) dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

                   !
                   ! Interpolate s to z-edges.
                   !
                   sp = HALF*(s(i,j,k+1,n)+s(i,j,k,n)) - SIXTH*(dsvl_r-dsvl_l)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sp = max(sp,min(s(i,j,k+1,n),s(i,j,k,n)))
                   sp = min(sp,max(s(i,j,k+1,n),s(i,j,k,n)))

                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                      sp = s(i,j,k,n)
                      sm = s(i,j,k,n)
                   else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                      sp = THREE*s(i,j,k,n) - TWO*sm
                   else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                      sm = THREE*s(i,j,k,n) - TWO*sp
                   end if
                   !
                   ! Different stencil needed for z-component of EXT_DIR and HOEXTRAP adv_bc's.
                   !
                   if (k .eq. lo(3)+1 .and. lo(3)+1 .eq. domlo(3)) then
                      if (adv_bc(3,1,bccomp) .eq. EXT_DIR  .or. adv_bc(3,1,bccomp) .eq. HOEXTRAP) then

                         ! The value in the first cc ghost cell represents the edge value.
                         !
                         sm = s(i,j,k-1,n)

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j,k-1,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i,j,k+1,n) &
                              - (ONE/20.0d0)*s(i,j,k+2,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k+1,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i,j,k+1,n),s(i,j,k,n)))
                         !
                         ! Copy sedge into sp and sm.
                         !
                         sp = sedge

                      end if
                   end if

                   if (k .eq. lo(3)+2 .and. lo(3)+1 .eq. domlo(3)) then
                      if (adv_bc(3,1,bccomp) .eq. EXT_DIR  .or. adv_bc(3,1,bccomp) .eq. HOEXTRAP) then

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j,k-2,n) &
                              + (THREE/FOUR)*s(i,j,k-1,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i,j,k+1,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k,n),s(i,j,k-1,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i,j,k-1,n)))
                         !
                         sm = sedge
                         !
                         ! Modify using quadratic limiters.
                         !
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if

                      end if
                   end if

                   if (k .eq. hi(3)-1 .and. hi(3)-1 .eq. domhi(3)) then
                      if (adv_bc(3,2,bccomp) .eq. EXT_DIR  .or. adv_bc(3,2,bccomp) .eq. HOEXTRAP) then

                         !
                         ! The value in the first cc ghost cell represents the edge value.
                         !
                         sp = s(i,j,k+1,n)

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j,k+1,n) &
                              + (THREE/FOUR)*s(i,j,k,n) &
                              + HALF        *s(i,j,k-1,n) &
                              - (ONE/20.0d0)*s(i,j,k-2,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k-1,n),s(i,j,k,n)))
                         sedge = min(sedge,max(s(i,j,k-1,n),s(i,j,k,n)))
                         !
                         sm = sedge
                      end if
                   end if

                   if (k .eq. hi(3)-2 .and. hi(3)-1 .eq. domhi(3)) then
                      if (adv_bc(3,2,bccomp) .eq. EXT_DIR  .or. adv_bc(3,2,bccomp) .eq. HOEXTRAP) then

                         !
                         ! Use a modified stencil to get sedge on the first interior edge.
                         !
                         sedge = -FIFTH     *s(i,j,k+2,n) &
                              + (THREE/FOUR)*s(i,j,k+1,n) &
                              + HALF        *s(i,j,k,n) &
                              - (ONE/20.0d0)*s(i,j,k-1,n)
                         !
                         ! Make sure sedge lies in between adjacent cell-centered values.
                         !
                         sedge = max(sedge,min(s(i,j,k,n),s(i,j,k+1,n)))
                         sedge = min(sedge,max(s(i,j,k,n),s(i,j,k+1,n)))
                         !
                         ! Copy sedge into sp and sm.
                         !
                         sp = sedge

                         ! Modify using quadratic limiters.
                         !
                         if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                            sp = s(i,j,k,n)
                            sm = s(i,j,k,n)
                         else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                            sp = THREE*s(i,j,k,n) - TWO*sm
                         else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                            sm = THREE*s(i,j,k,n) - TWO*sp
                         end if
                      end if
                   end if

                   !-------------------------------------------------------------------------
                   ! Compute z-component of Ip and Im.
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! w is MAC velocity -- use edge-based indexing

                      sigma = abs(w(i,j,k+1))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k+1) .gt. rel_eps) then
                         Ip(i,j,k,m,3) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,3) = s(i,j,k,n)
                      end if

                      sigma = abs(w(i,j,k))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,3) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,3) = s(i,j,k,n)
                      end if

                   else

                      sigma = abs(w(i,j,k))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,3) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,3) = s(i,j,k,n)
                      end if

                      sigma = abs(w(i,j,k))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,3) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,3) = s(i,j,k,n)
                      end if

                   endif
                end do
             end do
          end do

       else if (ppm_type .eq. 2) then

          !----------------------------------------------------------------------
          ! ppm_type = 2
          !----------------------------------------------------------------------
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ! -1
                   ! Interpolate s to z-edges.
                   sedgel = (7.d0/12.d0)*(s(i,j,k-2,n)+s(i,j,k-1,n)) &
                        - (1.d0/12.d0)*(s(i,j,k-3,n)+s(i,j,k,n))
                   !
                   ! Limit sedge.
                   if ((sedgel-s(i,j,k-2,n))*(s(i,j,k-1,n)-sedgel) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k-2,n)-TWO*sedgel+s(i,j,k-1,n))
                      D2L = s(i,j,k-3,n)-TWO*s(i,j,k-2,n)+s(i,j,k-1,n)
                      D2R = s(i,j,k-2,n)-TWO*s(i,j,k-1,n)+s(i,j,k,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgel = HALF*(s(i,j,k-2,n)+s(i,j,k-1,n)) - SIXTH*D2LIM
                   end if

                   ! 0
                   ! Interpolate s to z-edges.
                   sedge = (7.d0/12.d0)*(s(i,j,k-1,n)+s(i,j,k,n)) &
                        - (1.d0/12.d0)*(s(i,j,k-2,n)+s(i,j,k+1,n))
                   !
                   ! Limit sedge.
                   if ((sedge-s(i,j,k-1,n))*(s(i,j,k,n)-sedge) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k-1,n)-TWO*sedge+s(i,j,k,n))
                      D2L = s(i,j,k-2,n)-TWO*s(i,j,k-1,n)+s(i,j,k,n)
                      D2R = s(i,j,k-1,n)-TWO*s(i,j,k,n)+s(i,j,k+1,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedge = HALF*(s(i,j,k-1,n)+s(i,j,k,n)) - SIXTH*D2LIM
                   end if

                   ! +1
                   ! Interpolate s to z-edges.
                   sedger = (7.d0/12.d0)*(s(i,j,k,n)+s(i,j,k+1,n)) &
                        - (1.d0/12.d0)*(s(i,j,k-1,n)+s(i,j,k+2,n))
                   !
                   ! Limit sedge.
                   if ((sedger-s(i,j,k,n))*(s(i,j,k+1,n)-sedger) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k,n)-TWO*sedger+s(i,j,k+1,n))
                      D2L = s(i,j,k-1,n)-TWO*s(i,j,k,n)+s(i,j,k+1,n)
                      D2R = s(i,j,k,n)-TWO*s(i,j,k+1,n)+s(i,j,k+2,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedger = HALF*(s(i,j,k,n)+s(i,j,k+1,n)) - SIXTH*D2LIM
                   end if

                   ! +2
                   ! Interpolate s to z-edges.
                   sedgerr = (7.d0/12.d0)*(s(i,j,k+1,n)+s(i,j,k+2,n)) &
                        - (1.d0/12.d0)*(s(i,j,k,n)+s(i,j,k+3,n))
                   !
                   ! Limit sedge.
                   if ((sedgerr-s(i,j,k+1,n))*(s(i,j,k+2,n)-sedgerr) .lt. ZERO) then
                      D2  = THREE*(s(i,j,k+1,n)-TWO*sedgerr+s(i,j,k+2,n))
                      D2L = s(i,j,k,n)-TWO*s(i,j,k+1,n)+s(i,j,k+2,n)
                      D2R = s(i,j,k+1,n)-TWO*s(i,j,k+2,n)+s(i,j,k+3,n)
                      sgn = sign(ONE,D2)
                      D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                      sedgerr = HALF*(s(i,j,k+1,n)+s(i,j,k+2,n)) - SIXTH*D2LIM
                   end if


                   alphap = sedger-s(i,j,k,n)
                   alpham = sedge-s(i,j,k,n)
                   bigp = abs(alphap).gt.TWO*abs(alpham)
                   bigm = abs(alpham).gt.TWO*abs(alphap)
                   extremum = .false.

                   if (alpham*alphap .ge. ZERO) then
                      extremum = .true.
                   else if (bigp .or. bigm) then
                      !
                      ! Possible extremum. We look at cell centered values and face
                      ! centered values for a change in sign in the differences adjacent to
                      ! the cell. We use the pair of differences whose minimum magnitude is the
                      ! largest, and thus least susceptible to sensitivity to roundoff.
                      !
                      dafacem = sedge - sedgel
                      dafacep = sedgerr - sedger
                      dabarm = s(i,j,k,n) - s(i,j,k-1,n)
                      dabarp = s(i,j,k+1,n) - s(i,j,k,n)
                      dafacemin = min(abs(dafacem),abs(dafacep))
                      dabarmin= min(abs(dabarm),abs(dabarp))
                      if (dafacemin.ge.dabarmin) then
                         dachkm = dafacem
                         dachkp = dafacep
                      else
                         dachkm = dabarm
                         dachkp = dabarp
                      endif
                      extremum = (dachkm*dachkp .le. 0.d0)
                   end if

                   if (extremum) then
                      D2  = SIX*(alpham + alphap)
                      D2L = s(i,j,k-2,n)-TWO*s(i,j,k-1,n)+s(i,j,k,n)
                      D2R = s(i,j,k,n)-TWO*s(i,j,k+1,n)+s(i,j,k+2,n)
                      D2C = s(i,j,k-1,n)-TWO*s(i,j,k,n)+s(i,j,k+1,n)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j,k-1,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delam) then
                            if (sgn*(delam - alpham).ge.1.d-10) then
                               alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                            else
                               alphap = -TWO*alpham
                            endif
                         endif
                      end if
                      if (bigm) then
                         sgn = sign(ONE,alphap)
                         amax = -alpham**2 / (4*(alpham + alphap))
                         delap = s(i,j,k+1,n) - s(i,j,k,n)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm = s(i,j,k,n) + alpham
                   sp = s(i,j,k,n) + alphap

                   !
                   ! Different stencil needed for z-component of EXT_DIR and HOEXTRAP adv_bc's.
                   !
                   if (lo(3)+1 .eq. domlo(3)) then
                      if (adv_bc(3,1,bccomp) .eq. EXT_DIR  .or. adv_bc(3,1,bccomp) .eq. HOEXTRAP) then

                         if (k .eq. lo(3)+1) then

                            !
                            ! The value in the first cc ghost cell represents the edge value.
                            !
                            sm = s(i,j,k-1,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sp = -FIFTH        *s(i,j,k-1,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i,j,k+1,n) &
                                 - (ONE/20.0d0)*s(i,j,k+2,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sp = max(sp,min(s(i,j,k+1,n),s(i,j,k,n)))
                            sp = min(sp,max(s(i,j,k+1,n),s(i,j,k,n)))

                         elseif (k .eq. lo(3)+2) then

                            sedgel = s(i,j,k-2,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedge = -FIFTH     *s(i,j,k-2,n) &
                                 + (THREE/FOUR)*s(i,j,k-1,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k+1,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedge = max(sedge,min(s(i,j,k,n),s(i,j,k-1,n)))
                            sedge = min(sedge,max(s(i,j,k,n),s(i,j,k-1,n)))

                         elseif (k .eq. lo(3)+3) then

                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedgel = -FIFTH    *s(i,j,k-3,n) &
                                 + (THREE/FOUR)*s(i,j,k-2,n) &
                                 + HALF        *s(i,j,k-1,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedgel = max(sedgel,min(s(i,j,k-1,n),s(i,j,k-2,n)))
                            sedgel = min(sedgel,max(s(i,j,k-1,n),s(i,j,k-2,n)))

                         endif

                         ! Apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells.

                         if (k .eq. lo(3)+2 .or. k .eq. lo(3)+3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               !
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is
                               ! the largest, and thus least susceptible to sensitivity to roundoff.
                               !
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i,j,k-1,n)
                               dabarp = s(i,j,k+1,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i,j,k-2,n)-TWO*s(i,j,k-1,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i,j,k+1,n)+s(i,j,k+2,n)
                               D2C = s(i,j,k-1,n)-TWO*s(i,j,k,n)+s(i,j,k+1,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i,j,k-1,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i,j,k+1,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if

                      end if
                   end if

                   if (hi(3)-1 .eq. domhi(3)) then
                      if (adv_bc(3,2,bccomp) .eq. EXT_DIR  .or. adv_bc(3,2,bccomp) .eq. HOEXTRAP) then

                         if (k .eq. hi(3)-1) then

                            !
                            ! The value in the first cc ghost cell represents the edge value.
                            !
                            sp = s(i,j,k+1,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sm = -FIFTH        *s(i,j,k+1,n) &
                                 + (THREE/FOUR)*s(i,j,k,n) &
                                 + HALF        *s(i,j,k-1,n) &
                                 - (ONE/20.0d0)*s(i,j,k-2,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sm= max(sm,min(s(i,j,k-1,n),s(i,j,k,n)))
                            sm = min(sm,max(s(i,j,k-1,n),s(i,j,k,n)))

                         elseif (k .eq. hi(3)-2) then

                            sedgerr = s(i,j,k+2,n)
                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedger = -FIFTH    *s(i,j,k+2,n) &
                                 + (THREE/FOUR)*s(i,j,k+1,n) &
                                 + HALF        *s(i,j,k,n) &
                                 - (ONE/20.0d0)*s(i,j,k-1,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedger = max(sedger,min(s(i,j,k,n),s(i,j,k+1,n)))
                            sedger = min(sedger,max(s(i,j,k,n),s(i,j,k+1,n)))

                         elseif (k .eq. hi(3)-3) then

                            !
                            ! Use a modified stencil to get sedge on the first interior edge.
                            !
                            sedgerr = -FIFTH   *s(i,j,k+3,n) &
                                 + (THREE/FOUR)*s(i,j,k+2,n) &
                                 + HALF        *s(i,j,k+1,n) &
                                 - (ONE/20.0d0)*s(i,j,k,n)
                            !
                            ! Make sure sedge lies in between adjacent cell-centered values.
                            !
                            sedgerr = max(sedgerr,min(s(i,j,k+1,n),s(i,j,k+2,n)))
                            sedgerr = min(sedgerr,max(s(i,j,k+1,n),s(i,j,k+2,n)))

                         endif

                         !
                         ! Apply Colella 2008 limiters to compute sm and sp in the second
                         ! and third inner cells.
                         !
                         if (k .eq. hi(3)-2 .or. k .eq. hi(3)-3) then

                            alphap = sedger-s(i,j,k,n)
                            alpham = sedge-s(i,j,k,n)
                            bigp = abs(alphap).gt.TWO*abs(alpham)
                            bigm = abs(alpham).gt.TWO*abs(alphap)
                            extremum = .false.

                            if (alpham*alphap .ge. ZERO) then
                               extremum = .true.
                            else if (bigp .or. bigm) then
                               !
                               ! Possible extremum. We look at cell centered values and face
                               ! centered values for a change in sign in the differences adjacent to
                               ! the cell. We use the pair of differences whose minimum magnitude is
                               ! the largest, and thus least susceptible to sensitivity to roundoff.
                               !
                               dafacem = sedge - sedgel
                               dafacep = sedgerr - sedger
                               dabarm = s(i,j,k,n) - s(i,j,k-1,n)
                               dabarp = s(i,j,k+1,n) - s(i,j,k,n)
                               dafacemin = min(abs(dafacem),abs(dafacep))
                               dabarmin= min(abs(dabarm),abs(dabarp))
                               if (dafacemin.ge.dabarmin) then
                                  dachkm = dafacem
                                  dachkp = dafacep
                               else
                                  dachkm = dabarm
                                  dachkp = dabarp
                               endif
                               extremum = (dachkm*dachkp .le. 0.d0)
                            end if

                            if (extremum) then
                               D2  = SIX*(alpham + alphap)
                               D2L = s(i,j,k-2,n)-TWO*s(i,j,k-1,n)+s(i,j,k,n)
                               D2R = s(i,j,k,n)-TWO*s(i,j,k+1,n)+s(i,j,k+2,n)
                               D2C = s(i,j,k-1,n)-TWO*s(i,j,k,n)+s(i,j,k+1,n)
                               sgn = sign(ONE,D2)
                               D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                               D2ABS = max(abs(D2),1.d-10)
                               alpham = alpham*D2LIM/D2ABS
                               alphap = alphap*D2LIM/D2ABS
                            else
                               if (bigp) then
                                  sgn = sign(ONE,alpham)
                                  amax = -alphap**2 / (4*(alpham + alphap))
                                  delam = s(i,j,k-1,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delam) then
                                     if (sgn*(delam - alpham).ge.1.d-10) then
                                        alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                                     else
                                        alphap = -TWO*alpham
                                     endif
                                  endif
                               end if
                               if (bigm) then
                                  sgn = sign(ONE,alphap)
                                  amax = -alpham**2 / (4*(alpham + alphap))
                                  delap = s(i,j,k+1,n) - s(i,j,k,n)
                                  if (sgn*amax .ge. sgn*delap) then
                                     if (sgn*(delap - alphap).ge.1.d-10) then
                                        alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                                     else
                                        alpham = -TWO*alphap
                                     endif
                                  endif
                               end if
                            end if

                            sm = s(i,j,k,n) + alpham
                            sp = s(i,j,k,n) + alphap

                         end if
                      end if
                   end if

                   !-------------------------------------------------------------------------
                   ! Compute z-component of Ip and Im.
                   !-------------------------------------------------------------------------

                   if (is_umac == 1) then

                      ! w is MAC velocity -- use edge-based indexing

                      sigma = abs(w(i,j,k+1))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k+1) .gt. rel_eps) then
                         Ip(i,j,k,m,3) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,3) = s(i,j,k,n)
                      end if

                      sigma = abs(w(i,j,k))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,3) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,3) = s(i,j,k,n)
                      end if

                   else
                      sigma = abs(w(i,j,k))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k) .gt. rel_eps) then
                         Ip(i,j,k,m,3) = sp - &
                              (sigma/TWO)*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                      else
                         Ip(i,j,k,m,3) = s(i,j,k,n)
                      end if

                      sigma = abs(w(i,j,k))*dt/dx(3)
                      s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)
                      if (w(i,j,k) .lt. -rel_eps) then
                         Im(i,j,k,m,3) = sm + &
                              (sigma/TWO)*(sp-sm+(ONE-TWO3RD*sigma)*s6)
                      else
                         Im(i,j,k,m,3) = s(i,j,k,n)
                      end if

                   endif
                end do
             end do
          end do

       end if
    end do

  end subroutine ppm_3d

#endif


end module ppm_module
