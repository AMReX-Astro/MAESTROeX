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
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use meth_params_module, only: ppm_type, rel_eps

  implicit none

  private

  public :: ppm_2d, ppm_3d

contains

  !===========================================================================
  ! 2-d version
  !===========================================================================

  subroutine ppm_2d(s,ng_s,u,v,ng_u,Ip,Im,domlo,domhi,lo,hi,adv_bc,dx,dt,is_umac)

    ! note that u,v here may be the normal cell-centered velocity,
    ! or the MAC velocity.  The is_umac argument tells us which it
    ! is.

    integer         , intent(in   ) :: domlo(:),domhi(:),lo(:),hi(:),ng_s,ng_u
    double precision, intent(in   ) ::  s(lo(1)-ng_s:,lo(2)-ng_s:)
    double precision, intent(in   ) ::  u(lo(1)-ng_u:,lo(2)-ng_u:)
    double precision, intent(in   ) ::  v(lo(1)-ng_u:,lo(2)-ng_u:)
    double precision, intent(inout) :: Ip(lo(1)-1   :,lo(2)-1   :,:)
    double precision, intent(inout) :: Im(lo(1)-1   :,lo(2)-1   :,:)
    integer         , intent(in   ) :: adv_bc(:,:)
    double precision, intent(in   ) :: dx(:),dt
    logical         , intent(in   ) :: is_umac

    ! local
    integer :: i,j

    logical :: extremum, bigp, bigm

    double precision :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision :: sgn, sigma, s6, amax, delam, delap, D2ABS
    double precision :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! s_{\ib,+}, s_{\ib,-}
    double precision, pointer :: sp(:,:)
    double precision, pointer :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    double precision, pointer :: dsvl(:,:)

    ! s_{i+\half}^{H.O.}
    double precision, pointer :: sedge(:,:)

    ! cell-centered indexing
    call bl_allocate(sp,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1)
    call bl_allocate(sm,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1)

    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! cell-centered indexing w/extra x-ghost cell
    call bl_allocate(dsvl,lo(1)-2,hi(1)+2,lo(2)-1,hi(2)+1)

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+2,lo(2)-1,hi(2)+1)
    else if (ppm_type .eq. 2) then
       call bl_allocate(sedge,lo(1)-2,hi(1)+3,lo(2)-1,hi(2)+1)
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       !----------------------------------------------------------------------
       ! ppm_type = 1
       !----------------------------------------------------------------------

       ! compute van Leer slopes in x-direction
       dsvl = ZERO
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+2
             dsc = HALF * (s(i+1,j) - s(i-1,j))
             dsl = TWO  * (s(i  ,j) - s(i-1,j))
             dsr = TWO  * (s(i+1,j) - s(i  ,j))
             if (dsl*dsr .gt. 0) dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do

       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+2
             sedge(i,j) = HALF*(s(i,j)+s(i-1,j)) - SIXTH*(dsvl(i,j)-dsvl(i-1,j))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i-1,j)))
             sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i-1,j)))
          end do
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j) = sedge(i+1,j)
             sm(i,j) = sedge(i  ,j)
          end do
       end do

       ! modify using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(1) .eq. domlo(1)) then
          if (adv_bc(1,1) .eq. EXT_DIR  .or. adv_bc(1,1) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sm(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)+1,lo(2)-1:hi(2)+1) = &
                  -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
                  + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1) &
                  + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1) &
                  - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

             ! make sure sedge lies in between adjacent cell-centered values
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j) = max(sedge(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
                sedge(lo(1)+1,j) = min(sedge(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
             end do

             ! copy sedge into sp and sm
             do j=lo(2)-1,hi(2)+1
                sp(lo(1)  ,j) = sedge(lo(1)+1,j)
                sm(lo(1)+1,j) = sedge(lo(1)+1,j)
             end do

             ! reset sp on second interior edge
             do j=lo(2)-1,hi(2)+1
                sp(lo(1)+1,j) = sedge(lo(1)+2,j)
             end do

             ! modify using quadratic limiters
             do j=lo(2)-1,hi(2)+1
                i = lo(1)+1
                if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                   sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
                else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                   sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
                end if
             end do
          end if
       end if

       if (hi(1) .eq. domhi(1)) then
          if (adv_bc(1,2) .eq. EXT_DIR  .or. adv_bc(1,2) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sp(hi(1),lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(hi(1),lo(2)-1:hi(2)+1) = &
                  -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
                  + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1) &
                  + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1) &
                  - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

             ! make sure sedge lies in between adjacent cell-centered values
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j) = max(sedge(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
                sedge(hi(1),j) = min(sedge(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
             end do

             ! copy sedge into sp and sm
             do j=lo(2)-1,hi(2)+1
                sp(hi(1)-1,j) = sedge(hi(1),j)
                sm(hi(1)  ,j) = sedge(hi(1),j)
             end do

             ! reset sm on second interior edge
             do j=lo(2)-1,hi(2)+1
                sm(hi(1)-1,j) = sedge(hi(1)-1,j)
             end do

             ! modify using quadratic limiters
             do j=lo(2)-1,hi(2)+1
                i = hi(1)-1
                if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                   sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
                else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                   sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
                end if
             end do
          end if
       end if

    else if (ppm_type .eq. 2) then

       !----------------------------------------------------------------------
       ! ppm_type = 2
       !----------------------------------------------------------------------

       if (ng_s .lt. 4) then
          call amrex_error("Need 4 ghost cells for ppm_type=2")
       end if

       ! interpolate s to x-edges
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-2,hi(1)+3
             sedge(i,j) = (7.d0/12.d0)*(s(i-1,j)+s(i,j)) - (1.d0/12.d0)*(s(i-2,j)+s(i+1,j))
             ! limit sedge
             if ((sedge(i,j)-s(i-1,j))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i-1,j)-TWO*sedge(i,j)+s(i,j))
                D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                D2R = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedge(i,j) = HALF*(s(i-1,j)+s(i,j)) - SIXTH*D2LIM
             end if
          end do
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             alphap = sedge(i+1,j)-s(i,j)
             alpham = sedge(i  ,j)-s(i,j)
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
                dafacem = sedge(i,j) - sedge(i-1,j)
                dafacep = sedge(i+2,j) - sedge(i+1,j)
                dabarm = s(i,j) - s(i-1,j)
                dabarp = s(i+1,j) - s(i,j)
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
                D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                D2ABS = max(abs(D2),1.d-10)
                alpham = alpham*D2LIM/D2ABS
                alphap = alphap*D2LIM/D2ABS
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1,j) - s(i,j)
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
                   delap = s(i+1,j) - s(i,j)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm(i,j) = s(i,j) + alpham
             sp(i,j) = s(i,j) + alphap

          end do
       end do

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(1) .eq. domlo(1)) then
          if (adv_bc(1,1) .eq. EXT_DIR  .or. adv_bc(1,1) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sm(lo(1),lo(2)-1:hi(2)+1)    = s(lo(1)-1,lo(2)-1:hi(2)+1)
             sedge(lo(1),lo(2)-1:hi(2)+1) = s(lo(1)-1,lo(2)-1:hi(2)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)+1,lo(2)-1:hi(2)+1) = &
                  -FIFTH        *s(lo(1)-1,lo(2)-1:hi(2)+1) &
                  + (THREE/FOUR)*s(lo(1)  ,lo(2)-1:hi(2)+1) &
                  + HALF        *s(lo(1)+1,lo(2)-1:hi(2)+1) &
                  - (ONE/20.0d0)*s(lo(1)+2,lo(2)-1:hi(2)+1)

             ! make sure sedge lies in between adjacent cell-centered values
             do j=lo(2)-1,hi(2)+1
                sedge(lo(1)+1,j) = max(sedge(lo(1)+1,j),min(s(lo(1)+1,j),s(lo(1),j)))
                sedge(lo(1)+1,j) = min(sedge(lo(1)+1,j),max(s(lo(1)+1,j),s(lo(1),j)))
             end do

             ! copy sedge into sp
             do j=lo(2)-1,hi(2)+1
                sp(lo(1)  ,j) = sedge(lo(1)+1,j)
             end do

             ! apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)+1,lo(1)+2

                   alphap = sedge(i+1,j)-s(i,j)
                   alpham = sedge(i  ,j)-s(i,j)
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
                      dafacem = sedge(i,j) - sedge(i-1,j)
                      dafacep = sedge(i+2,j) - sedge(i+1,j)
                      dabarm = s(i,j) - s(i-1,j)
                      dabarp = s(i+1,j) - s(i,j)
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
                      D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                      D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                      D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j) - s(i,j)
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
                         delap = s(i+1,j) - s(i,j)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j) = s(i,j) + alpham
                   sp(i,j) = s(i,j) + alphap

                end do
             end do
          end if
       end if

       if (hi(1) .eq. domhi(1)) then
          if (adv_bc(1,2) .eq. EXT_DIR  .or. adv_bc(1,2) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sp(hi(1),lo(2)-1:hi(2)+1)      = s(hi(1)+1,lo(2)-1:hi(2)+1)
             sedge(hi(1)+1,lo(2)-1:hi(2)+1) = s(hi(1)+1,lo(2)-1:hi(2)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(hi(1),lo(2)-1:hi(2)+1) = &
                  -FIFTH        *s(hi(1)+1,lo(2)-1:hi(2)+1) &
                  + (THREE/FOUR)*s(hi(1)  ,lo(2)-1:hi(2)+1) &
                  + HALF        *s(hi(1)-1,lo(2)-1:hi(2)+1) &
                  - (ONE/20.0d0)*s(hi(1)-2,lo(2)-1:hi(2)+1)

             ! make sure sedge lies in between adjacent cell-centered values
             do j=lo(2)-1,hi(2)+1
                sedge(hi(1),j) = max(sedge(hi(1),j),min(s(hi(1)-1,j),s(hi(1),j)))
                sedge(hi(1),j) = min(sedge(hi(1),j),max(s(hi(1)-1,j),s(hi(1),j)))
             end do

             ! copy sedge into sm
             do j=lo(2)-1,hi(2)+1
                sm(hi(1)  ,j) = sedge(hi(1),j)
             end do

             ! reset sm on second interior edge
             do j=lo(2)-1,hi(2)+1
                sm(hi(1)-1,j) = sedge(hi(1)-1,j)
             end do

             ! apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells
             do j=lo(2)-1,hi(2)+1
                do i=hi(1)-2,hi(1)-1

                   alphap = sedge(i+1,j)-s(i,j)
                   alpham = sedge(i  ,j)-s(i,j)
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
                      dafacem = sedge(i,j) - sedge(i-1,j)
                      dafacep = sedge(i+2,j) - sedge(i+1,j)
                      dabarm = s(i,j) - s(i-1,j)
                      dabarp = s(i+1,j) - s(i,j)
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
                      D2L = s(i-2,j)-TWO*s(i-1,j)+s(i,j)
                      D2R = s(i,j)-TWO*s(i+1,j)+s(i+2,j)
                      D2C = s(i-1,j)-TWO*s(i,j)+s(i+1,j)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i-1,j) - s(i,j)
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
                         delap = s(i+1,j) - s(i,j)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j) = s(i,j) + alpham
                   sp(i,j) = s(i,j) + alphap

                end do
             end do
          end if
       end if

    end if

    !-------------------------------------------------------------------------
    ! compute x-component of Ip and Im
    !-------------------------------------------------------------------------

    if (is_umac) then

       ! u here is umac, so use edge-based indexing
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)
             sigma = abs(u(i+1,j))*dt/dx(1)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (u(i+1,j) .gt. rel_eps) then
                Ip(i,j,1) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,1) = s(i,j)
             end if
          end do

          do i=lo(1),hi(1)+1
             sigma = abs(u(i,j))*dt/dx(1)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (u(i,j) .lt. -rel_eps) then
                Im(i,j,1) = sm(i,j) + (sigma/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
             else
                Im(i,j,1) = s(i,j)
             end if
          end do
       end do

    else
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)
             sigma = abs(u(i,j))*dt/dx(1)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (u(i,j) .gt. rel_eps) then
                Ip(i,j,1) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,1) = s(i,j)
             end if
          end do

          do i=lo(1),hi(1)+1
             sigma = abs(u(i,j))*dt/dx(1)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (u(i,j) .lt. -rel_eps) then
                Im(i,j,1) = sm(i,j) + (sigma/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
             else
                Im(i,j,1) = s(i,j)
             end if
          end do
       end do
    endif

    call bl_deallocate(sedge)
    call bl_deallocate(dsvl)


    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------

    ! cell-centered indexing w/extra y-ghost cell
    call bl_allocate(dsvl,lo(1)-1,hi(1)+1,lo(2)-2,hi(2)+2)

    ! edge-centered indexing for y-faces
    if (ppm_type .eq. 1) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+2)
    else if (ppm_type .eq. 2) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+1,lo(2)-2,hi(2)+3)
    end if

    ! compute s at y-edges
    if (ppm_type .eq. 1) then

       !----------------------------------------------------------------------
       ! ppm_type = 1
       !----------------------------------------------------------------------

       ! compute van Leer slopes in y-direction
       dsvl = ZERO
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             dsc = HALF * (s(i,j+1) - s(i,j-1))
             dsl = TWO  * (s(i,j  ) - s(i,j-1))
             dsr = TWO  * (s(i,j+1) - s(i,j  ))
             if (dsl*dsr .gt. ZERO) dsvl(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          end do
       end do

       ! interpolate s to y-edges
       do j=lo(2)-1,hi(2)+2
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = HALF*(s(i,j)+s(i,j-1)) - SIXTH*(dsvl(i,j)-dsvl(i,j-1))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j) = max(sedge(i,j),min(s(i,j),s(i,j-1)))
             sedge(i,j) = min(sedge(i,j),max(s(i,j),s(i,j-1)))
          end do
       end do

       ! copy sedge into sp and sm
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sp(i,j) = sedge(i,j+1)
             sm(i,j) = sedge(i,j  )
          end do
       end do

       ! modify using quadratic limiters
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                sp(i,j) = s(i,j)
                sm(i,j) = s(i,j)
             else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
             else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
             end if
          end do
       end do

       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(2) .eq. domlo(2)) then
          if (adv_bc(2,1) .eq. EXT_DIR  .or. adv_bc(2,1) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sm(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)-1:hi(1)+1,lo(2)+1) = &
                  -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
                  + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ) &
                  + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1) &
                  - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

             ! make sure sedge lies in between adjacent cell-centered values
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1) = max(sedge(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
                sedge(i,lo(2)+1) = min(sedge(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
             end do

             ! copy sedge into sp and sm
             do i=lo(1)-1,hi(1)+1
                sp(i,lo(2)  ) = sedge(i,lo(2)+1)
                sm(i,lo(2)+1) = sedge(i,lo(2)+1)
             end do

             ! reset sp on second interior edge
             do i=lo(1)-1,hi(1)+1
                sp(i,lo(2)+1) = sedge(i,lo(2)+2)
             end do

             ! modify using quadratic limiters
             do i=lo(1)-1,hi(1)+1
                j = lo(2)+1
                if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                   sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
                else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                   sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
                end if
             end do
          end if
       end if

       if (hi(2) .eq. domhi(2)) then
          if (adv_bc(2,2) .eq. EXT_DIR  .or. adv_bc(2,2) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sp(lo(1)-1:hi(1)+1,hi(2)) = s(lo(1)-1:hi(1)+1,hi(2)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)-1:hi(1)+1,hi(2)) = &
                  -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
                  + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ) &
                  + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1) &
                  - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

             ! make sure sedge lies in between adjacent cell-centered values
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2)) = max(sedge(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
                sedge(i,hi(2)) = min(sedge(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
             end do

             ! copy sedge into sp and sm
             do i=lo(1)-1,hi(1)+1
                sp(i,hi(2)-1) = sedge(i,hi(2))
                sm(i,hi(2)  ) = sedge(i,hi(2))
             end do

             ! reset sm on second interior edge
             do i=lo(1)-1,hi(1)+1
                sm(i,hi(2)-1) = sedge(i,hi(2)-1)
             end do

             ! modify using quadratic limiters
             do i=lo(1)-1,hi(1)+1
                j = hi(2)-1
                if ((sp(i,j)-s(i,j))*(s(i,j)-sm(i,j)) .le. ZERO) then
                   sp(i,j) = s(i,j)
                   sm(i,j) = s(i,j)
                else if (abs(sp(i,j)-s(i,j)) .ge. TWO*abs(sm(i,j)-s(i,j))) then
                   sp(i,j) = THREE*s(i,j) - TWO*sm(i,j)
                else if (abs(sm(i,j)-s(i,j)) .ge. TWO*abs(sp(i,j)-s(i,j))) then
                   sm(i,j) = THREE*s(i,j) - TWO*sp(i,j)
                end if
             end do
          end if
       end if

    else if (ppm_type .eq. 2) then

       !----------------------------------------------------------------------
       ! ppm_type = 2
       !----------------------------------------------------------------------

       ! interpolate s to y-edges
       do j=lo(2)-2,hi(2)+3
          do i=lo(1)-1,hi(1)+1
             sedge(i,j) = (7.d0/12.d0)*(s(i,j-1)+s(i,j)) - (1.d0/12.d0)*(s(i,j-2)+s(i,j+1))
             ! limit sedge
             if ((sedge(i,j)-s(i,j-1))*(s(i,j)-sedge(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i,j-1)-TWO*sedge(i,j)+s(i,j))
                D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                D2R = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedge(i,j) = HALF*(s(i,j-1)+s(i,j)) - SIXTH*D2LIM
             end if
          end do
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             alphap = sedge(i,j+1)-s(i,j)
             alpham = sedge(i,j  )-s(i,j)
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
                dafacem = sedge(i,j) - sedge(i,j-1)
                dafacep = sedge(i,j+2) - sedge(i,j+1)
                dabarm = s(i,j) - s(i,j-1)
                dabarp = s(i,j+1) - s(i,j)
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
                D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                sgn = sign(ONE,D2)
                D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                D2ABS = max(abs(D2),1.d-10)
                alpham = alpham*D2LIM/D2ABS
                alphap = alphap*D2LIM/D2ABS
             else
                if (bigp) then
                   sgn = sign(ONE,alpham)
                   amax = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j-1) - s(i,j)
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
                   delap = s(i,j+1) - s(i,j)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm(i,j) = s(i,j) + alpham
             sp(i,j) = s(i,j) + alphap

          end do
       end do

       ! different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(2) .eq. domlo(2)) then
          if (adv_bc(2,1) .eq. EXT_DIR  .or. adv_bc(2,1) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sm(lo(1)-1:hi(1)+1,lo(2))    = s(lo(1)-1:hi(1)+1,lo(2)-1)
             sedge(lo(1)-1:hi(1)+1,lo(2)) = s(lo(1)-1:hi(1)+1,lo(2)-1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)-1:hi(1)+1,lo(2)+1) = &
                  -FIFTH        *s(lo(1)-1:hi(1)+1,lo(2)-1) &
                  + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,lo(2)  ) &
                  + HALF        *s(lo(1)-1:hi(1)+1,lo(2)+1) &
                  - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,lo(2)+2)

             ! make sure sedge lies in between adjacent cell-centered values
             do i=lo(1)-1,hi(1)+1
                sedge(i,lo(2)+1) = max(sedge(i,lo(2)+1),min(s(i,lo(2)+1),s(i,lo(2))))
                sedge(i,lo(2)+1) = min(sedge(i,lo(2)+1),max(s(i,lo(2)+1),s(i,lo(2))))
             end do

             ! copy sedge into sp
             do i=lo(1)-1,hi(1)+1
                sp(i,lo(2)  ) = sedge(i,lo(2)+1)
             end do

             ! apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells
             do j=lo(2)+1,lo(2)+2
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j+1)-s(i,j)
                   alpham = sedge(i,j  )-s(i,j)
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
                      dafacem = sedge(i,j) - sedge(i,j-1)
                      dafacep = sedge(i,j+2) - sedge(i,j+1)
                      dabarm = s(i,j) - s(i,j-1)
                      dabarp = s(i,j+1) - s(i,j)
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
                      D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                      D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                      D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1) - s(i,j)
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
                         delap = s(i,j+1) - s(i,j)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j) = s(i,j) + alpham
                   sp(i,j) = s(i,j) + alphap

                end do
             end do
          end if
       end if

       if (hi(2) .eq. domhi(2)) then
          if (adv_bc(2,2) .eq. EXT_DIR  .or. adv_bc(2,2) .eq. HOEXTRAP) then
             ! the value in the first cc ghost cell represents the edge value
             sp(lo(1)-1:hi(1)+1,hi(2))      = s(lo(1)-1:hi(1)+1,hi(2)+1)
             sedge(lo(1)-1:hi(1)+1,hi(2)+1) = s(lo(1)-1:hi(1)+1,hi(2)+1)

             ! use a modified stencil to get sedge on the first interior edge
             sedge(lo(1)-1:hi(1)+1,hi(2)) = &
                  -FIFTH        *s(lo(1)-1:hi(1)+1,hi(2)+1) &
                  + (THREE/FOUR)*s(lo(1)-1:hi(1)+1,hi(2)  ) &
                  + HALF        *s(lo(1)-1:hi(1)+1,hi(2)-1) &
                  - (ONE/20.0d0)*s(lo(1)-1:hi(1)+1,hi(2)-2)

             ! make sure sedge lies in between adjacent cell-centered values
             do i=lo(1)-1,hi(1)+1
                sedge(i,hi(2)) = max(sedge(i,hi(2)),min(s(i,hi(2)-1),s(i,hi(2))))
                sedge(i,hi(2)) = min(sedge(i,hi(2)),max(s(i,hi(2)-1),s(i,hi(2))))
             end do

             ! copy sedge into sm
             do i=lo(1)-1,hi(1)+1
                sm(i,hi(2)  ) = sedge(i,hi(2))
             end do

             ! apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells
             do j=hi(2)-2,hi(2)-1
                do i=lo(1)-1,hi(1)+1

                   alphap = sedge(i,j+1)-s(i,j)
                   alpham = sedge(i,j  )-s(i,j)
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
                      dafacem = sedge(i,j) - sedge(i,j-1)
                      dafacep = sedge(i,j+2) - sedge(i,j+1)
                      dabarm = s(i,j) - s(i,j-1)
                      dabarp = s(i,j+1) - s(i,j)
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
                      D2L = s(i,j-2)-TWO*s(i,j-1)+s(i,j)
                      D2R = s(i,j)-TWO*s(i,j+1)+s(i,j+2)
                      D2C = s(i,j-1)-TWO*s(i,j)+s(i,j+1)
                      sgn = sign(ONE,D2)
                      D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                      D2ABS = max(abs(D2),1.d-10)
                      alpham = alpham*D2LIM/D2ABS
                      alphap = alphap*D2LIM/D2ABS
                   else
                      if (bigp) then
                         sgn = sign(ONE,alpham)
                         amax = -alphap**2 / (4*(alpham + alphap))
                         delam = s(i,j-1) - s(i,j)
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
                         delap = s(i,j+1) - s(i,j)
                         if (sgn*amax .ge. sgn*delap) then
                            if (sgn*(delap - alphap).ge.1.d-10) then
                               alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                            else
                               alpham = -TWO*alphap
                            endif
                         endif
                      end if
                   end if

                   sm(i,j) = s(i,j) + alpham
                   sp(i,j) = s(i,j) + alphap

                end do
             end do
          end if

       end if
    end if

    !-------------------------------------------------------------------------
    ! compute y-component of Ip and Im
    !-------------------------------------------------------------------------

    if (is_umac) then

       ! v here is vmac, so use edge-based indexing
       do j=lo(2)-1,hi(2)
          do i=lo(1)-1,hi(1)+1
             sigma = abs(v(i,j+1))*dt/dx(2)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (v(i,j+1) .gt. rel_eps) then
                Ip(i,j,2) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,2) = s(i,j)
             end if
          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sigma = abs(v(i,j))*dt/dx(2)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (v(i,j) .lt. -rel_eps) then
                Im(i,j,2) = sm(i,j) + (sigma/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
             else
                Im(i,j,2) = s(i,j)
             end if
          end do
       end do

    else
       do j=lo(2)-1,hi(2)
          do i=lo(1)-1,hi(1)+1
             sigma = abs(v(i,j))*dt/dx(2)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (v(i,j) .gt. rel_eps) then
                Ip(i,j,2) = sp(i,j) - (sigma/TWO)*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,2) = s(i,j)
             end if
          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1)-1,hi(1)+1
             sigma = abs(v(i,j))*dt/dx(2)
             s6 = SIX*s(i,j) - THREE*(sm(i,j)+sp(i,j))
             if (v(i,j) .lt. -rel_eps) then
                Im(i,j,2) = sm(i,j) + (sigma/TWO)*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
             else
                Im(i,j,2) = s(i,j)
             end if
          end do
       end do
    endif

    call bl_deallocate(sp)
    call bl_deallocate(sm)
    call bl_deallocate(dsvl)
    call bl_deallocate(sedge)

  end subroutine ppm_2d


  !===========================================================================
  ! 3-d version
  !===========================================================================

  ! characteristics based on u
  subroutine ppm_3d(s,ng_s,u,v,w,ng_u,Ip,Im,domlo,domhi,lo,hi,adv_bc,dx,dt,is_umac)

    integer         , intent(in   ) :: domlo(:),domhi(:),lo(:),hi(:),ng_s,ng_u
    double precision, intent(in   ) ::  s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
    double precision, intent(in   ) ::  u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    double precision, intent(in   ) ::  v(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    double precision, intent(in   ) ::  w(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
    double precision, intent(inout) :: Ip(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :,:)
    double precision, intent(inout) :: Im(lo(1)-1    :,lo(2)-1    :,lo(3)-1    :,:)
    integer         , intent(in   ) :: adv_bc(:,:)
    double precision, intent(in   ) :: dx(:),dt
    logical         , intent(in   ) :: is_umac

    ! local
    integer :: i,j,k

    logical :: extremum, bigp, bigm

    double precision :: dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision :: sgn, sigma, s6, D2ABS
    double precision :: dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp
    double precision :: amax, delam, delap

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! s_{\ib,+}, s_{\ib,-}
    double precision, pointer :: sp(:,:,:)
    double precision, pointer :: sm(:,:,:)

    ! \delta s_{\ib}^{vL}
    double precision, pointer :: dsvl(:,:,:)

    ! s_{i+\half}^{H.O.}
    double precision, pointer :: sedge(:,:,:)

    ! cell-centered indexing
    call bl_allocate(sp,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(sm,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! cell-centered indexing w/extra x-ghost cell
    call bl_allocate(dsvl,lo(1)-2,hi(1)+2,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+2,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    else if (ppm_type .eq. 2) then
       call bl_allocate(sedge,lo(1)-2,hi(1)+3,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       !----------------------------------------------------------------------
       ! ppm_type = 1
       !----------------------------------------------------------------------

       dsvl = ZERO

       !$OMP PARALLEL PRIVATE(i,j,k,dsc,dsl,dsr)
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+2
                !
                ! Compute van Leer slopes in x-direction.
                !
                dsc = HALF * (s(i+1,j,k) - s(i-1,j,k))
                dsl = TWO  * (s(i  ,j,k) - s(i-1,j,k))
                dsr = TWO  * (s(i+1,j,k) - s(i  ,j,k))
                if (dsl*dsr .gt. ZERO) dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                !
                ! Interpolate s to x-edges.
                !
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i-1,j,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i-1,j,k))
                !
                ! Make sure sedge lies in between adjacent cell-centered values.
                !
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i-1,j,k)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i-1,j,k)))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                !
                ! Copy sedge into sp and sm.
                !
                sp(i,j,k) = sedge(i+1,j,k)
                sm(i,j,k) = sedge(i  ,j,k)
                !
                ! Modify using quadratic limiters.
                !
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       !
       ! Different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's.
       !
       if (lo(1) .eq. domlo(1)) then
          if (adv_bc(1,1) .eq. EXT_DIR  .or. adv_bc(1,1) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

             i = lo(1)+1

             !$OMP PARALLEL DO PRIVATE(j,k)
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)-1,hi(2)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(lo(1)+1,j,k) = -FIFTH        *s(lo(1)-1,j,k) &
                        + (THREE/FOUR)*s(lo(1)  ,j,k) &
                        + HALF        *s(lo(1)+1,j,k) &
                        - (ONE/20.0d0)*s(lo(1)+2,j,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                   sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))
                   !
                   ! Copy sedge into sp and sm.
                   !
                   sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
                   sm(lo(1)+1,j,k) = sedge(lo(1)+1,j,k)
                   !
                   ! Reset sp on second interior edge.
                   !
                   sp(lo(1)+1,j,k) = sedge(lo(1)+2,j,k)
                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                      sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                   else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                      sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                   end if
                end do
             end do
             !$OMP END PARALLEL DO

          end if
       end if


       if (hi(1) .eq. domhi(1)) then
          if (adv_bc(1,2) .eq. EXT_DIR  .or. adv_bc(1,2) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
                  s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

             i = hi(1)-1

             !$OMP PARALLEL DO PRIVATE(j,k)
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)-1,hi(2)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(hi(1),j,k) = -FIFTH        *s(hi(1)+1,j,k) &
                        + (THREE/FOUR)*s(hi(1)  ,j,k) &
                        + HALF        *s(hi(1)-1,j,k) &
                        - (ONE/20.0d0)*s(hi(1)-2,j,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                   sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
                   !
                   ! Copy sedge into sp and sm.
                   !
                   sp(hi(1)-1,j,k) = sedge(hi(1),j,k)
                   sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
                   !
                   ! Reset sm on second interior edge.
                   !
                   sm(hi(1)-1,j,k) = sedge(hi(1)-1,j,k)
                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                      sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                   else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                      sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                   end if
                end do
             end do
             !$OMP END PARALLEL DO

          end if
       end if

    else if (ppm_type .eq. 2) then

       !----------------------------------------------------------------------
       ! ppm_type = 2
       !----------------------------------------------------------------------

       if (ng_s .lt. 4) then
          call amrex_error("Need 4 ghost cells for ppm_type=2")
       end if

       !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
       !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
       !$OMP PRIVATE(delam,delap,D2ABS)
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-2,hi(1)+3
                !
                ! Interpolate s to x-edges.
                !
                sedge(i,j,k) = (7.d0/12.d0)*(s(i-1,j,k)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i-2,j,k)+s(i+1,j,k))
                !
                ! Limit sedge.
                !
                if ((sedge(i,j,k)-s(i-1,j,k))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i-1,j,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                   D2R = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i-1,j,k)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
       !$OMP END DO
       !
       ! Use Colella 2008 limiters.
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       !
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i+1,j,k)-s(i,j,k)
                alpham = sedge(i  ,j,k)-s(i,j,k)
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
                   dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                   dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                   dabarm = s(i,j,k) - s(i-1,j,k)
                   dabarp = s(i+1,j,k) - s(i,j,k)
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
                   D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                   D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   D2ABS = max(abs(D2),1.d-10)
                   alpham = alpham*D2LIM/D2ABS
                   alphap = alphap*D2LIM/D2ABS
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i-1,j,k) - s(i,j,k)
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
                      delap = s(i+1,j,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

       ! different stencil needed for x-component of EXT_DIR and HOEXTRAP adv_bc's
       if (lo(1) .eq. domlo(1)) then
          if (adv_bc(1,1) .eq. EXT_DIR  .or. adv_bc(1,1) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sm(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)    = &
                  s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
             sedge(lo(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

             !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
             !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
             !$OMP PRIVATE(delam,delap,D2ABS)

             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)-1,hi(2)+1
                   ! use a modified stencil to get sedge on the first interior edge
                   sedge(lo(1)+1,j,k) = -FIFTH        *s(lo(1)-1,j,k) &
                        + (THREE/FOUR)*s(lo(1)  ,j,k) &
                        + HALF        *s(lo(1)+1,j,k) &
                        - (ONE/20.0d0)*s(lo(1)+2,j,k)

                   ! make sure sedge lies in between adjacent cell-centered values
                   sedge(lo(1)+1,j,k) = max(sedge(lo(1)+1,j,k),min(s(lo(1)+1,j,k),s(lo(1),j,k)))
                   sedge(lo(1)+1,j,k) = min(sedge(lo(1)+1,j,k),max(s(lo(1)+1,j,k),s(lo(1),j,k)))

                   ! copy sedge into sp
                   sp(lo(1)  ,j,k) = sedge(lo(1)+1,j,k)
                end do
             end do
             !$OMP END DO
             !
             ! Apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells.
             !
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)-1,hi(2)+1
                   do i=lo(1)+1,lo(1)+2

                      alphap = sedge(i+1,j,k)-s(i,j,k)
                      alpham = sedge(i  ,j,k)-s(i,j,k)
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
                         dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                         dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                         dabarm = s(i,j,k) - s(i-1,j,k)
                         dabarp = s(i+1,j,k) - s(i,j,k)
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
                         D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                         D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                         D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                         sgn = sign(ONE,D2)
                         D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                         D2ABS = max(abs(D2),1.d-10)
                         alpham = alpham*D2LIM/D2ABS
                         alphap = alphap*D2LIM/D2ABS
                      else
                         if (bigp) then
                            sgn = sign(ONE,alpham)
                            amax = -alphap**2 / (4*(alpham + alphap))
                            delam = s(i-1,j,k) - s(i,j,k)
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
                            delap = s(i+1,j,k) - s(i,j,k)
                            if (sgn*amax .ge. sgn*delap) then
                               if (sgn*(delap - alphap).ge.1.d-10) then
                                  alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                               else
                                  alpham = -TWO*alphap
                               endif
                            endif
                         end if
                      end if

                      sm(i,j,k) = s(i,j,k) + alpham
                      sp(i,j,k) = s(i,j,k) + alphap

                   end do
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
          end if
       end if

       if (hi(1) .eq. domhi(1)) then
          if (adv_bc(1,2) .eq. EXT_DIR  .or. adv_bc(1,2) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sp(hi(1),lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) = &
                  s(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

             !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
             !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
             !$OMP PRIVATE(delam,delap,D2ABS)
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)-1,hi(2)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(hi(1),j,k) = -FIFTH        *s(hi(1)+1,j,k) &
                        + (THREE/FOUR)*s(hi(1)  ,j,k) &
                        + HALF        *s(hi(1)-1,j,k) &
                        - (ONE/20.0d0)*s(hi(1)-2,j,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(hi(1),j,k) = max(sedge(hi(1),j,k),min(s(hi(1)-1,j,k),s(hi(1),j,k)))
                   sedge(hi(1),j,k) = min(sedge(hi(1),j,k),max(s(hi(1)-1,j,k),s(hi(1),j,k)))
                   !
                   ! Copy sedge into sm.
                   !
                   sm(hi(1)  ,j,k) = sedge(hi(1),j,k)
                end do
             end do
             !$OMP END DO
             !
             ! Apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells.
             !
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)-1,hi(2)+1
                   do i=hi(1)-2,hi(1)-1

                      alphap = sedge(i+1,j,k)-s(i,j,k)
                      alpham = sedge(i  ,j,k)-s(i,j,k)
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
                         dafacem = sedge(i,j,k) - sedge(i-1,j,k)
                         dafacep = sedge(i+2,j,k) - sedge(i+1,j,k)
                         dabarm = s(i,j,k) - s(i-1,j,k)
                         dabarp = s(i+1,j,k) - s(i,j,k)
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
                         D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                         D2R = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                         D2C = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                         sgn = sign(ONE,D2)
                         D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                         D2ABS = max(abs(D2),1.d-10)
                         alpham = alpham*D2LIM/D2ABS
                         alphap = alphap*D2LIM/D2ABS
                      else
                         if (bigp) then
                            sgn = sign(ONE,alpham)
                            amax = -alphap**2 / (4*(alpham + alphap))
                            delam = s(i-1,j,k) - s(i,j,k)
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
                            delap = s(i+1,j,k) - s(i,j,k)
                            if (sgn*amax .ge. sgn*delap) then
                               if (sgn*(delap - alphap).ge.1.d-10) then
                                  alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                               else
                                  alpham = -TWO*alphap
                               endif
                            endif
                         end if
                      end if

                      sm(i,j,k) = s(i,j,k) + alpham
                      sp(i,j,k) = s(i,j,k) + alphap

                   end do
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
          end if

       end if
    end if

    !-------------------------------------------------------------------------
    ! Compute x-component of Ip and Im.
    !-------------------------------------------------------------------------

    if (is_umac) then

       ! u is MAC velocity -- use edge-based indexing

       !$OMP PARALLEL DO PRIVATE(i,j,k,sigma,s6)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)
                sigma = abs(u(i+1,j,k))*dt/dx(1)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (u(i+1,j,k) .gt. rel_eps) then
                   Ip(i,j,k,1) = sp(i,j,k) - &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                else
                   Ip(i,j,k,1) = s(i,j,k)
                end if
             end do

             do i=lo(1),hi(1)+1
                sigma = abs(u(i,j,k))*dt/dx(1)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (u(i,j,k) .lt. -rel_eps) then
                   Im(i,j,k,1) = sm(i,j,k) + &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
                else
                   Im(i,j,k,1) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,sigma,s6)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)
                sigma = abs(u(i,j,k))*dt/dx(1)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (u(i,j,k) .gt. rel_eps) then
                   Ip(i,j,k,1) = sp(i,j,k) - &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                else
                   Ip(i,j,k,1) = s(i,j,k)
                end if
             end do

             do i=lo(1),hi(1)+1
                sigma = abs(u(i,j,k))*dt/dx(1)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (u(i,j,k) .lt. -rel_eps) then
                   Im(i,j,k,1) = sm(i,j,k) + &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
                else
                   Im(i,j,k,1) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    call bl_deallocate(sedge)
    call bl_deallocate(dsvl)


    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------

    ! cell-centered indexing w/extra y-ghost cell
    call bl_allocate( dsvl,lo(1)-1,hi(1)+1,lo(2)-2,hi(2)+2,lo(3)-1,hi(3)+1)

    ! edge-centered indexing for y-faces
    if (ppm_type .eq. 1) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+2,lo(3)-1,hi(3)+1)
    else if (ppm_type .eq. 2) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+1,lo(2)-2,hi(2)+3,lo(3)-1,hi(3)+1)
    end if
    !
    ! Compute s at y-edges.
    !
    if (ppm_type .eq. 1) then

       !----------------------------------------------------------------------
       ! ppm_type = 1
       !----------------------------------------------------------------------

       dsvl = ZERO

       !$OMP PARALLEL PRIVATE(i,j,k,dsc,dsl,dsr)
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                !
                ! Compute van Leer slopes in y-direction.
                !
                dsc = HALF * (s(i,j+1,k) - s(i,j-1,k))
                dsl = TWO  * (s(i,j  ,k) - s(i,j-1,k))
                dsr = TWO  * (s(i,j+1,k) - s(i,j  ,k))
                if (dsl*dsr .gt. ZERO) dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                !
                ! Interpolate s to y-edges.
                !
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j-1,k)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j-1,k))
                !
                ! Make sure sedge lies in between adjacent cell-centered values.
                !
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i,j-1,k)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i,j-1,k)))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                !
                ! Copy sedge into sp and sm.
                !
                sp(i,j,k) = sedge(i,j+1,k)
                sm(i,j,k) = sedge(i,j  ,k)
                !
                ! Modify using quadratic limiters.
                !
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       !
       !
       ! Different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's.
       !
       if (lo(2) .eq. domlo(2)) then
          if (adv_bc(2,1) .eq. EXT_DIR  .or. adv_bc(2,1) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

             j = lo(2)+1

             !$OMP PARALLEL DO PRIVATE(i,k)
             do k=lo(3)-1,hi(3)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,lo(2)+1,k) = -FIFTH        *s(i,lo(2)-1,k) &
                        + (THREE/FOUR)*s(i,lo(2)  ,k) &
                        + HALF        *s(i,lo(2)+1,k) &
                        - (ONE/20.0d0)*s(i,lo(2)+2,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                   sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
                   !
                   ! Copy sedge into sp and sm.
                   !
                   sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
                   sm(i,lo(2)+1,k) = sedge(i,lo(2)+1,k)
                   !
                   ! Reset sp on second interior edge.
                   !
                   sp(i,lo(2)+1,k) = sedge(i,lo(2)+2,k)
                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                      sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                   else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                      sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                   end if
                end do
             end do
             !$OMP END PARALLEL DO

          end if
       end if

       if (hi(2) .eq. domhi(2)) then
          if (adv_bc(2,2) .eq. EXT_DIR  .or. adv_bc(2,2) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

             j = hi(2)-1

             !$OMP PARALLEL DO PRIVATE(i,k)
             do k=lo(3)-1,hi(3)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,hi(2),k) = -FIFTH        *s(i,hi(2)+1,k) &
                        + (THREE/FOUR)*s(i,hi(2)  ,k) &
                        + HALF        *s(i,hi(2)-1,k) &
                        - (ONE/20.0d0)*s(i,hi(2)-2,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                   sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
                   !
                   ! Copy sedge into sp and sm.
                   !
                   sp(i,hi(2)-1,k) = sedge(i,hi(2),k)
                   sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
                   !
                   ! Reset sm on second interior edge.
                   !
                   sm(i,hi(2)-1,k) = sedge(i,hi(2)-1,k)
                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                      sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                   else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                      sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                   end if
                end do
             end do
             !$OMP END PARALLEL DO

          end if
       end if

    else if (ppm_type .eq. 2) then

       !----------------------------------------------------------------------
       ! ppm_type = 2
       !----------------------------------------------------------------------

       !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
       !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
       !$OMP PRIVATE(delam,delap,D2ABS)
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-2,hi(2)+3
             do i=lo(1)-1,hi(1)+1
                !
                ! Interpolate s to y-edges.
                !
                sedge(i,j,k) = (7.d0/12.d0)*(s(i,j-1,k)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i,j-2,k)+s(i,j+1,k))
                !
                ! Limit sedge.
                !
                if ((sedge(i,j,k)-s(i,j-1,k))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i,j-1,k)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                   D2R = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i,j-1,k)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
       !$OMP END DO
       !
       ! Use Colella 2008 limiters.
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       !
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j+1,k)-s(i,j,k)
                alpham = sedge(i,j  ,k)-s(i,j,k)
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
                   dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                   dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                   dabarm = s(i,j,k) - s(i,j-1,k)
                   dabarp = s(i,j+1,k) - s(i,j,k)
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
                   D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                   D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   D2ABS = max(abs(D2),1.d-10)
                   alpham = alpham*D2LIM/D2ABS
                   alphap = alphap*D2LIM/D2ABS
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j-1,k) - s(i,j,k)
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
                      delap = s(i,j+1,k) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       !
       ! Different stencil needed for y-component of EXT_DIR and HOEXTRAP adv_bc's.
       !
       if (lo(2) .eq. domlo(2)) then
          if (adv_bc(2,1) .eq. EXT_DIR  .or. adv_bc(2,1) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sm(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)
             sedge(lo(1)-1:hi(1)+1,lo(2),lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1)

             !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
             !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
             !$OMP PRIVATE(delam,delap,D2ABS)
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,lo(2)+1,k) = -FIFTH        *s(i,lo(2)-1,k) &
                        + (THREE/FOUR)*s(i,lo(2)  ,k) &
                        + HALF        *s(i,lo(2)+1,k) &
                        - (ONE/20.0d0)*s(i,lo(2)+2,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,lo(2)+1,k) = max(sedge(i,lo(2)+1,k),min(s(i,lo(2)+1,k),s(i,lo(2),k)))
                   sedge(i,lo(2)+1,k) = min(sedge(i,lo(2)+1,k),max(s(i,lo(2)+1,k),s(i,lo(2),k)))
                   !
                   ! Copy sedge into sp.
                   !
                   sp(i,lo(2)  ,k) = sedge(i,lo(2)+1,k)
                end do
             end do
             !$OMP END DO
             !
             ! Apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells.
             !
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do j=lo(2)+1,lo(2)+2
                   do i=lo(1)-1,hi(1)+1

                      alphap = sedge(i,j+1,k)-s(i,j,k)
                      alpham = sedge(i,j  ,k)-s(i,j,k)
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
                         dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                         dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                         dabarm = s(i,j,k) - s(i,j-1,k)
                         dabarp = s(i,j+1,k) - s(i,j,k)
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
                         D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                         D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                         D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                         sgn = sign(ONE,D2)
                         D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                         D2ABS = max(abs(D2),1.d-10)
                         alpham = alpham*D2LIM/D2ABS
                         alphap = alphap*D2LIM/D2ABS
                      else
                         if (bigp) then
                            sgn = sign(ONE,alpham)
                            amax = -alphap**2 / (4*(alpham + alphap))
                            delam = s(i,j-1,k) - s(i,j,k)
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
                            delap = s(i,j+1,k) - s(i,j,k)
                            if (sgn*amax .ge. sgn*delap) then
                               if (sgn*(delap - alphap).ge.1.d-10) then
                                  alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                               else
                                  alpham = -TWO*alphap
                               endif
                            endif
                         end if
                      end if

                      sm(i,j,k) = s(i,j,k) + alpham
                      sp(i,j,k) = s(i,j,k) + alphap

                   end do
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
          end if
       end if

       if (hi(2) .eq. domhi(2)) then
          if (adv_bc(2,2) .eq. EXT_DIR  .or. adv_bc(2,2) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sp(lo(1)-1:hi(1)+1,hi(2),lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)
             sedge(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1)

             !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
             !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
             !$OMP PRIVATE(delam,delap,D2ABS)
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,hi(2),k) = -FIFTH        *s(i,hi(2)+1,k) &
                        + (THREE/FOUR)*s(i,hi(2)  ,k) &
                        + HALF        *s(i,hi(2)-1,k) &
                        - (ONE/20.0d0)*s(i,hi(2)-2,k)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,hi(2),k) = max(sedge(i,hi(2),k),min(s(i,hi(2)-1,k),s(i,hi(2),k)))
                   sedge(i,hi(2),k) = min(sedge(i,hi(2),k),max(s(i,hi(2)-1,k),s(i,hi(2),k)))
                   !
                   ! Copy sedge into sm.
                   !
                   sm(i,hi(2)  ,k) = sedge(i,hi(2),k)
                end do
             end do
             !$OMP END DO
             !
             ! Apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells.
             !
             !$OMP DO
             do k=lo(3)-1,hi(3)+1
                do j=hi(2)-2,hi(2)-1
                   do i=lo(1)-1,hi(1)+1

                      alphap = sedge(i,j+1,k)-s(i,j,k)
                      alpham = sedge(i,j  ,k)-s(i,j,k)
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
                         dafacem = sedge(i,j,k) - sedge(i,j-1,k)
                         dafacep = sedge(i,j+2,k) - sedge(i,j+1,k)
                         dabarm = s(i,j,k) - s(i,j-1,k)
                         dabarp = s(i,j+1,k) - s(i,j,k)
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
                         D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                         D2R = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                         D2C = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                         sgn = sign(ONE,D2)
                         D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                         D2ABS = max(abs(D2),1.d-10)
                         alpham = alpham*D2LIM/D2ABS
                         alphap = alphap*D2LIM/D2ABS
                      else
                         if (bigp) then
                            sgn = sign(ONE,alpham)
                            amax = -alphap**2 / (4*(alpham + alphap))
                            delam = s(i,j-1,k) - s(i,j,k)
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
                            delap = s(i,j+1,k) - s(i,j,k)
                            if (sgn*amax .ge. sgn*delap) then
                               if (sgn*(delap - alphap).ge.1.d-10) then
                                  alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                               else
                                  alpham = -TWO*alphap
                               endif
                            endif
                         end if
                      end if

                      sm(i,j,k) = s(i,j,k) + alpham
                      sp(i,j,k) = s(i,j,k) + alphap

                   end do
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
          end if
       end if

    end if

    !-------------------------------------------------------------------------
    ! Compute y-component of Ip and Im.
    !-------------------------------------------------------------------------

    if (is_umac) then

       ! v is MAC velocity -- use edge-based indexing

       !$OMP PARALLEL DO PRIVATE(i,j,k,sigma,s6)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)
             do i=lo(1)-1,hi(1)+1
                sigma = abs(v(i,j+1,k))*dt/dx(2)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (v(i,j+1,k) .gt. rel_eps) then
                   Ip(i,j,k,2) = sp(i,j,k) - &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                else
                   Ip(i,j,k,2) = s(i,j,k)
                end if
             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sigma = abs(v(i,j,k))*dt/dx(2)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (v(i,j,k) .lt. -rel_eps) then
                   Im(i,j,k,2) = sm(i,j,k) + &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
                else
                   Im(i,j,k,2) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else

       !$OMP PARALLEL DO PRIVATE(i,j,k,sigma,s6)
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)
             do i=lo(1)-1,hi(1)+1
                sigma = abs(v(i,j,k))*dt/dx(2)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (v(i,j,k) .gt. rel_eps) then
                   Ip(i,j,k,2) = sp(i,j,k) - &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                else
                   Ip(i,j,k,2) = s(i,j,k)
                end if
             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sigma = abs(v(i,j,k))*dt/dx(2)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (v(i,j,k) .lt. -rel_eps) then
                   Im(i,j,k,2) = sm(i,j,k) + &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
                else
                   Im(i,j,k,2) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    call bl_deallocate(sedge)
    call bl_deallocate(dsvl)


    !-------------------------------------------------------------------------
    ! z-direction
    !-------------------------------------------------------------------------

    ! cell-centered indexing w/extra z-ghost cell
    call bl_allocate(dsvl,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-2,hi(3)+2)

    ! edge-centered indexing for z-faces
    if (ppm_type .eq. 1) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+2)
    else if (ppm_type .eq. 2) then
       call bl_allocate(sedge,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-2,hi(3)+3)
    end if
    !
    ! Compute s at z-edges.
    !
    if (ppm_type .eq. 1) then

       !----------------------------------------------------------------------
       ! ppm_type = 1
       !----------------------------------------------------------------------

       dsvl = ZERO

       !$OMP PARALLEL PRIVATE(i,j,k,dsc,dsl,dsr)
       !$OMP DO
       do k=lo(3)-2,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                !
                ! Compute van Leer slopes in z-direction.
                !
                dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
                dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
                dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
                if (dsl*dsr .gt. ZERO) dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                !
                ! Interpolate s to z-edges.
                !
                sedge(i,j,k) = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvl(i,j,k)-dsvl(i,j,k-1))
                !
                ! Make sure sedge lies in between adjacent cell-centered values.
                !
                sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k),s(i,j,k-1)))
                sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k),s(i,j,k-1)))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                !
                ! Copy sedge into sp and sm.
                !
                sp(i,j,k) = sedge(i,j,k+1)
                sm(i,j,k) = sedge(i,j,k  )
                !
                ! Modify using quadratic limiters.
                !
                if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                   sp(i,j,k) = s(i,j,k)
                   sm(i,j,k) = s(i,j,k)
                else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                   sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                   sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       !
       ! Different stencil needed for z-component of EXT_DIR and HOEXTRAP adv_bc's.
       !
       if (lo(3) .eq. domlo(3)) then
          if (adv_bc(3,1) .eq. EXT_DIR  .or. adv_bc(3,1) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

             k = lo(3)+1

             !$OMP PARALLEL DO PRIVATE(i,j)
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,j,lo(3)+1) = -FIFTH        *s(i,j,lo(3)-1) &
                        + (THREE/FOUR)*s(i,j,lo(3)  ) &
                        + HALF        *s(i,j,lo(3)+1) &
                        - (ONE/20.0d0)*s(i,j,lo(3)+2)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                   sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
                   !
                   ! Copy sedge into sp and sm.
                   !
                   sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
                   sm(i,j,lo(3)+1) = sedge(i,j,lo(3)+1)
                   !
                   ! Reset sp on second interior edge.
                   !
                   sp(i,j,lo(3)+1) = sedge(i,j,lo(3)+2)
                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                      sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                   else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                      sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                   end if
                end do
             end do
             !$OMP END PARALLEL DO

          end if
       end if

       if (hi(3) .eq. domhi(3)) then
          if (adv_bc(3,2) .eq. EXT_DIR  .or. adv_bc(3,2) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

             k = hi(3)-1

             !$OMP PARALLEL DO PRIVATE(i,j)
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,j,hi(3)) = -FIFTH        *s(i,j,hi(3)+1) &
                        + (THREE/FOUR)*s(i,j,hi(3)  ) &
                        + HALF        *s(i,j,hi(3)-1) &
                        - (ONE/20.0d0)*s(i,j,hi(3)-2)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                   sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
                   !
                   ! Copy sedge into sp and sm.
                   !
                   sp(i,j,hi(3)-1) = sedge(i,j,hi(3))
                   sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
                   !
                   ! Reset sm on second interior edge.
                   !
                   sm(i,j,hi(3)-1) = sedge(i,j,hi(3)-1)
                   !
                   ! Modify using quadratic limiters.
                   !
                   if ((sp(i,j,k)-s(i,j,k))*(s(i,j,k)-sm(i,j,k)) .le. ZERO) then
                      sp(i,j,k) = s(i,j,k)
                      sm(i,j,k) = s(i,j,k)
                   else if (abs(sp(i,j,k)-s(i,j,k)) .ge. TWO*abs(sm(i,j,k)-s(i,j,k))) then
                      sp(i,j,k) = THREE*s(i,j,k) - TWO*sm(i,j,k)
                   else if (abs(sm(i,j,k)-s(i,j,k)) .ge. TWO*abs(sp(i,j,k)-s(i,j,k))) then
                      sm(i,j,k) = THREE*s(i,j,k) - TWO*sp(i,j,k)
                   end if
                end do
             end do
             !$OMP END PARALLEL DO

          end if
       end if

    else if (ppm_type .eq. 2) then

       !----------------------------------------------------------------------
       ! ppm_type = 2
       !----------------------------------------------------------------------

       !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
       !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
       !$OMP PRIVATE(delam,delap,D2ABS)
       !$OMP DO
       do k=lo(3)-2,hi(3)+3
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                !
                ! Interpolate s to z-edges.
                !
                sedge(i,j,k) = (7.d0/12.d0)*(s(i,j,k-1)+s(i,j,k)) &
                     - (1.d0/12.d0)*(s(i,j,k-2)+s(i,j,k+1))
                !
                ! Limit sedge.
                !
                if ((sedge(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedge(i,j,k)) .lt. ZERO) then
                   D2  = THREE*(s(i,j,k-1)-TWO*sedge(i,j,k)+s(i,j,k))
                   D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                   D2R = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                   sgn = sign(ONE,D2)
                   D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                   sedge(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - SIXTH*D2LIM
                end if
             end do
          end do
       end do
       !$OMP END DO
       !
       ! Use Colella 2008 limiters.
       ! This is a new version of the algorithm
       ! to eliminate sensitivity to roundoff.
       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1

                alphap = sedge(i,j,k+1)-s(i,j,k)
                alpham = sedge(i,j,k  )-s(i,j,k)
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
                   dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                   dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                   dabarm = s(i,j,k) - s(i,j,k-1)
                   dabarp = s(i,j,k+1) - s(i,j,k)
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
                   D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                   D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                   D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                   sgn = sign(ONE,D2)
                   D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                   D2ABS = max(abs(D2),1.d-10)
                   alpham = alpham*D2LIM/D2ABS
                   alphap = alphap*D2LIM/D2ABS
                else
                   if (bigp) then
                      sgn = sign(ONE,alpham)
                      amax = -alphap**2 / (4*(alpham + alphap))
                      delam = s(i,j,k-1) - s(i,j,k)
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
                      delap = s(i,j,k+1) - s(i,j,k)
                      if (sgn*amax .ge. sgn*delap) then
                         if (sgn*(delap - alphap).ge.1.d-10) then
                            alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                         else
                            alpham = -TWO*alphap
                         endif
                      endif
                   end if
                end if

                sm(i,j,k) = s(i,j,k) + alpham
                sp(i,j,k) = s(i,j,k) + alphap

             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       !
       ! Different stencil needed for z-component of EXT_DIR and HOEXTRAP adv_bc's.
       !
       if (lo(3) .eq. domlo(3)) then
          if (adv_bc(3,1) .eq. EXT_DIR  .or. adv_bc(3,1) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sm(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)
             sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1)

             !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
             !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
             !$OMP PRIVATE(delam,delap,D2ABS)
             !$OMP DO
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,j,lo(3)+1) = -FIFTH        *s(i,j,lo(3)-1) &
                        + (THREE/FOUR)*s(i,j,lo(3)  ) &
                        + HALF        *s(i,j,lo(3)+1) &
                        - (ONE/20.0d0)*s(i,j,lo(3)+2)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,j,lo(3)+1) = max(sedge(i,j,lo(3)+1),min(s(i,j,lo(3)+1),s(i,j,lo(3))))
                   sedge(i,j,lo(3)+1) = min(sedge(i,j,lo(3)+1),max(s(i,j,lo(3)+1),s(i,j,lo(3))))
                   !
                   ! Copy sedge into sp.
                   !
                   sp(i,j,lo(3)  ) = sedge(i,j,lo(3)+1)
                end do
             end do
             !$OMP END DO
             !
             ! Apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells.
             !
             !$OMP DO
             do j=lo(2)-1,hi(2)+1
                do k=lo(3)+1,lo(3)+2
                   do i=lo(1)-1,hi(1)+1

                      alphap = sedge(i,j,k+1)-s(i,j,k)
                      alpham = sedge(i,j,k  )-s(i,j,k)
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
                         dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                         dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                         dabarm = s(i,j,k) - s(i,j,k-1)
                         dabarp = s(i,j,k+1) - s(i,j,k)
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
                         D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                         D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                         D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                         sgn = sign(ONE,D2)
                         D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                         D2ABS = max(abs(D2),1.d-10)
                         alpham = alpham*D2LIM/D2ABS
                         alphap = alphap*D2LIM/D2ABS
                      else
                         if (bigp) then
                            sgn = sign(ONE,alpham)
                            amax = -alphap**2 / (4*(alpham + alphap))
                            delam = s(i,j,k-1) - s(i,j,k)
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
                            delap = s(i,j,k+1) - s(i,j,k)
                            if (sgn*amax .ge. sgn*delap) then
                               if (sgn*(delap - alphap).ge.1.d-10) then
                                  alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                               else
                                  alpham = -TWO*alphap
                               endif
                            endif
                         end if
                      end if

                      sm(i,j,k) = s(i,j,k) + alpham
                      sp(i,j,k) = s(i,j,k) + alphap

                   end do
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
          end if
       end if

       if (hi(3) .eq. domhi(3)) then
          if (adv_bc(3,2) .eq. EXT_DIR  .or. adv_bc(3,2) .eq. HOEXTRAP) then
             !
             ! The value in the first cc ghost cell represents the edge value.
             !
             sp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)
             sedge(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1) = &
                  s(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1)

             !$OMP PARALLEL PRIVATE(i,j,k,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
             !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax) &
             !$OMP PRIVATE(delam,delap,D2ABS)
             !$OMP DO
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+1
                   !
                   ! Use a modified stencil to get sedge on the first interior edge.
                   !
                   sedge(i,j,hi(3)) = -FIFTH        *s(i,j,hi(3)+1) &
                        + (THREE/FOUR)*s(i,j,hi(3)  ) &
                        + HALF        *s(i,j,hi(3)-1) &
                        - (ONE/20.0d0)*s(i,j,hi(3)-2)
                   !
                   ! Make sure sedge lies in between adjacent cell-centered values.
                   !
                   sedge(i,j,hi(3)) = max(sedge(i,j,hi(3)),min(s(i,j,hi(3)-1),s(i,j,hi(3))))
                   sedge(i,j,hi(3)) = min(sedge(i,j,hi(3)),max(s(i,j,hi(3)-1),s(i,j,hi(3))))
                   !
                   ! Copy sedge into sm.
                   !
                   sm(i,j,hi(3)  ) = sedge(i,j,hi(3))
                end do
             end do
             !$OMP END DO
             !
             ! Apply Colella 2008 limiters to compute sm and sp in the second
             ! and third inner cells.
             !
             !$OMP DO
             do j=lo(2)-1,hi(2)+1
                do k=hi(3)-2,hi(3)-1
                   do i=lo(1)-1,hi(1)+1

                      alphap = sedge(i,j,k+1)-s(i,j,k)
                      alpham = sedge(i,j,k  )-s(i,j,k)
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
                         dafacem = sedge(i,j,k) - sedge(i,j,k-1)
                         dafacep = sedge(i,j,k+2) - sedge(i,j,k+1)
                         dabarm = s(i,j,k) - s(i,j,k-1)
                         dabarp = s(i,j,k+1) - s(i,j,k)
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
                         D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                         D2R = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                         D2C = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                         sgn = sign(ONE,D2)
                         D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                         D2ABS = max(abs(D2),1.d-10)
                         alpham = alpham*D2LIM/D2ABS
                         alphap = alphap*D2LIM/D2ABS
                      else
                         if (bigp) then
                            sgn = sign(ONE,alpham)
                            amax = -alphap**2 / (4*(alpham + alphap))
                            delam = s(i,j,k-1) - s(i,j,k)
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
                            delap = s(i,j,k+1) - s(i,j,k)
                            if (sgn*amax .ge. sgn*delap) then
                               if (sgn*(delap - alphap).ge.1.d-10) then
                                  alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                               else
                                  alpham = -TWO*alphap
                               endif
                            endif
                         end if
                      end if

                      sm(i,j,k) = s(i,j,k) + alpham
                      sp(i,j,k) = s(i,j,k) + alphap

                   end do
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
          end if
       end if

    end if

    !-------------------------------------------------------------------------
    ! Compute z-component of Ip and Im.
    !-------------------------------------------------------------------------

    if (is_umac) then

       ! w is MAC velocity -- use edge-based indexing

       !$OMP PARALLEL PRIVATE(i,j,k,sigma,s6)
       !$OMP DO
       do k=lo(3)-1,hi(3)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sigma = abs(w(i,j,k+1))*dt/dx(3)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (w(i,j,k+1) .gt. rel_eps) then
                   Ip(i,j,k,3) = sp(i,j,k) - &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                else
                   Ip(i,j,k,3) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k=lo(3),hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sigma = abs(w(i,j,k))*dt/dx(3)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (w(i,j,k) .lt. -rel_eps) then
                   Im(i,j,k,3) = sm(i,j,k) + &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
                else
                   Im(i,j,k,3) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    else

       !$OMP PARALLEL PRIVATE(i,j,k,sigma,s6)
       !$OMP DO
       do k=lo(3)-1,hi(3)
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sigma = abs(w(i,j,k))*dt/dx(3)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (w(i,j,k) .gt. rel_eps) then
                   Ip(i,j,k,3) = sp(i,j,k) - &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)-(ONE-TWO3RD*sigma)*s6)
                else
                   Ip(i,j,k,3) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k=lo(3),hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                sigma = abs(w(i,j,k))*dt/dx(3)
                s6 = SIX*s(i,j,k) - THREE*(sm(i,j,k)+sp(i,j,k))
                if (w(i,j,k) .lt. -rel_eps) then
                   Im(i,j,k,3) = sm(i,j,k) + &
                        (sigma/TWO)*(sp(i,j,k)-sm(i,j,k)+(ONE-TWO3RD*sigma)*s6)
                else
                   Im(i,j,k,3) = s(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    endif

    call bl_deallocate(sp)
    call bl_deallocate(sm)
    call bl_deallocate(dsvl)
    call bl_deallocate(sedge)

  end subroutine ppm_3d

end module ppm_module
