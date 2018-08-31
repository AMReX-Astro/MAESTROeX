
#include "AMReX_BC_TYPES.H"

module fill_umac_ghost_module

  implicit none

  private

contains

  subroutine fill_umac_ghost(domlo, domhi, lo, hi, &
                             umac, umac_lo, umac_hi, &
                             vmac, vmac_lo, vmac_hi, &
#if (AMREX_SPACEDIM == 3)
                             wmac, wmac_lo, wmac_hi, &
#endif
                             phys_bc) bind(C, name="fill_umac_ghost")


    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: umac_lo(3), umac_hi(3)
    integer         , intent(in   ) :: vmac_lo(3), vmac_hi(3)
    double precision, intent(inout) :: umac(umac_lo(1):umac_hi(1), &
                                            umac_lo(2):umac_hi(2), &
                                            umac_lo(3):umac_hi(3))
    double precision, intent(inout) :: vmac(vmac_lo(1):vmac_hi(1), &
                                            vmac_lo(2):vmac_hi(2), &
                                            vmac_lo(3):vmac_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: wmac_lo(3), wmac_hi(3)
    double precision, intent(inout) :: wmac(wmac_lo(1):wmac_hi(1), &
                                            wmac_lo(2):wmac_hi(2), &
                                            wmac_lo(3):wmac_hi(3))
#endif
    integer         , intent(in   ) :: phys_bc(AMREX_SPACEDIM,2)

    ! lo x-faces
    if (lo(1) .eq. domlo(1)) then
       select case (phys_bc(1,1))
       case (Inflow)
          umac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = umac(lo(1),lo(2):hi(2),lo(3):hi(3))
          vmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
          vmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#endif
       case (Outflow)
          umac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = umac(lo(1),lo(2):hi(2),lo(3):hi(3))
          vmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = vmac(lo(1),lo(2):hi(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = wmac(lo(1),lo(2):hi(2),lo(3):hi(3))
#endif
       case (Symmetry)
          umac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = -umac(lo(1)+1,lo(2):hi(2),lo(3):hi(3))
          vmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = vmac(lo(1),lo(2):hi(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1)-1,lo(2):hi(2),lo(3):hi(3)) = wmac(lo(1),lo(2):hi(2),lo(3):hi(3))
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! hi x-faces
    if (hi(1) .eq. domhi(1)) then
       select case (phys_bc(1,2))
       case (Inflow)
          umac(hi(1)+2,lo(2):hi(2),lo(3):hi(3)) = umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3))
          vmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(hi(1)+2,lo(2):hi(2),lo(3):hi(3)) = 0.d0
          vmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
#endif
       case (Outflow)
          umac(hi(1)+2,lo(2):hi(2),lo(3):hi(3)) = umac(hi(1)+1,lo(2):hi(2),lo(3):hi(3))
          vmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = vmac(hi(1),lo(2):hi(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = wmac(hi(1),lo(2):hi(2),lo(3):hi(3))
#endif
       case (Symmetry)
          umac(hi(1)+2,lo(2):hi(2),lo(3):hi(3)) = -umac(hi(1),lo(2):hi(2),lo(3):hi(3))
          vmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = vmac(hi(1),lo(2):hi(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = wmac(hi(1),lo(2):hi(2),lo(3):hi(3))
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! lo y-faces
    if (lo(2) .eq. domlo(2)) then
       select case (phys_bc(2,1))
       case (Inflow)
          umac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = 0.d0
          vmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = vmac(lo(1):hi(1),lo(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = 0.d0
          vmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = 0.d0
#endif
       case (Outflow)
          umac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = umac(lo(1):hi(1),lo(2),lo(3):hi(3))
          vmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = vmac(lo(1):hi(1),lo(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = wmac(lo(1):hi(1),lo(2),lo(3):hi(3))
#endif
       case (Symmetry)
          umac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = umac(lo(1):hi(1),lo(2),lo(3):hi(3))
          vmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = -vmac(lo(1):hi(1),lo(2)+1,lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),lo(2)-1,lo(3):hi(3)) = wmac(lo(1):hi(1),lo(2),lo(3):hi(3))
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! hi y-faces
    if (hi(2) .eq. domhi(2)) then
       select case (phys_bc(2,2))
       case (Inflow)
          umac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
          vmac(lo(1):hi(1),hi(2)+2,lo(3):hi(3)) = vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
#endif
       case (SlipWall, NoSlipWall)
          umac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
          vmac(lo(1):hi(1),hi(2)+2,lo(3):hi(3)) = 0.d0
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
#endif
       case (Outflow)
          umac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = umac(lo(1):hi(1),hi(2),lo(3):hi(3))
          vmac(lo(1):hi(1),hi(2)+2,lo(3):hi(3)) = vmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = wmac(lo(1):hi(1),hi(2),lo(3):hi(3))
#endif
       case (Symmetry)
          umac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = umac(lo(1):hi(1),hi(2),lo(3):hi(3))
          vmac(lo(1):hi(1),hi(2)+2,lo(3):hi(3)) = -vmac(lo(1):hi(1),hi(2),lo(3):hi(3))
#if (AMREX_SPACEDIM == 3)
          wmac(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = wmac(lo(1):hi(1),hi(2),lo(3):hi(3))
#endif
       case (Interior)
          ! do nothing
       case default
       end select
    end if

#if (AMREX_SPACEDIM == 3)

    ! lo z-faces
    if (lo(3) .eq. domlo(3)) then
       select case (phys_bc(3,1))
       case (Inflow)
          umac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = 0.d0
          vmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = 0.d0
          wmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = wmac(lo(1):hi(1),lo(2):hi(2),lo(3))
       case (SlipWall, NoSlipWall)
          umac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = 0.d0
          vmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = 0.d0
          wmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = 0.d0
       case (Outflow)
          umac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = umac(lo(1):hi(1),lo(2):hi(2),lo(3))
          vmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = vmac(lo(1):hi(1),lo(2):hi(2),lo(3))
          wmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = wmac(lo(1):hi(1),lo(2):hi(2),lo(3))
       case (Symmetry)
          umac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = umac(lo(1):hi(1),lo(2):hi(2),lo(3))
          vmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = vmac(lo(1):hi(1),lo(2):hi(2),lo(3))
          wmac(lo(1):hi(1),lo(2):hi(2),lo(3)-1) = -wmac(lo(1):hi(1),lo(2):hi(2),lo(3)+1)
       case (Interior)
          ! do nothing
       case default
       end select
    end if

    ! hi z-faces
    if (hi(3) .eq. domhi(3)) then
       select case (phys_bc(3,2))
       case (Inflow)
          umac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
          vmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
          wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+2) = wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1)
       case (SlipWall, NoSlipWall)
          umac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
          vmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
          wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+2) = 0.d0
       case (Outflow)
          umac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = umac(lo(1):hi(1),lo(2):hi(2),hi(3))
          vmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = vmac(lo(1):hi(1),lo(2):hi(2),hi(3))
          wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+2) = wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1)
       case (Symmetry)
          umac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = umac(lo(1):hi(1),lo(2):hi(2),hi(3))
          vmac(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = vmac(lo(1):hi(1),lo(2):hi(2),hi(3))
          wmac(lo(1):hi(1),lo(2):hi(2),hi(3)+2) = -wmac(lo(1):hi(1),lo(2):hi(2),hi(3))
       case (Interior)
          ! do nothing
       case default
       end select
    end if

#endif

  end subroutine fill_umac_ghost

  subroutine PC_EDGE_INTERP(lo, hi, nc, ratio, dir, &
                            crse, c_lo, c_hi, nc_c, &
                            fine, f_lo, f_hi, nc_f) bind(C,name="PC_EDGE_INTERP")

    implicit none

    integer         , intent(in   ) :: lo(3),hi(3), nc, ratio(0:AMREX_SPACEDIM-1), dir
    integer         , intent(in   ) :: c_lo(3), c_hi(3), nc_c
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    double precision, intent(in   ) :: crse(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nc_c)
    double precision, intent(inout) :: fine(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)

    integer i,j,k,ii,jj,n,L
#if (AMREX_SPACEDIM == 3)
    integer kk, P
#endif
    !
    !     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
    !     Operate only on faces that overlap--ie, only fill the fine faces that make up each
    !     coarse face, leave the in-between faces alone.
    !

#if (AMREX_SPACEDIM == 2)
    k = 0
    if (dir.eq.0) then
       do n=1,nc
          do j=lo(2),hi(2)
             jj = ratio(1)*j
             do i=lo(1),hi(1)
                ii = ratio(0)*i
                do L=0,ratio(1)-1
                   fine(ii,jj+L,k,n) = crse(i,j,k,n)
                enddo
             enddo
          enddo
       enddo
    else
       do n=1,nc
          do j=lo(2),hi(2)
             jj = ratio(1)*j
             do i=lo(1),hi(1)
                ii = ratio(0)*i
                do L=0,ratio(0)-1
                   fine(ii+L,jj,k,n) = crse(i,j,k,n)
                enddo
             enddo
          enddo
       enddo
    endif
#elif (AMREX_SPACEDIM == 3)
    if (dir.eq.0) then
       do n=1,nc
          do k=lo(3),hi(3)
             kk = ratio(2)*k
             do j=lo(2),hi(2)
                jj = ratio(1)*j
                do i=lo(1),hi(1)
                   ii = ratio(0)*i
                   do P=0,ratio(2)-1
                      do L=0,ratio(1)-1
                         fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else if (dir.eq.1) then
       do n=1,nc
          do k=lo(3),hi(3)
             kk = ratio(2)*k
             do j=lo(2),hi(2)
                jj = ratio(1)*j
                do i=lo(1),hi(1)
                   ii = ratio(0)*i
                   do P=0,ratio(2)-1
                      do L=0,ratio(0)-1
                         fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do n=1,nc
          do k=lo(3),hi(3)
             kk = ratio(2)*k
             do j=lo(2),hi(2)
                jj = ratio(1)*j
                do i=lo(1),hi(1)
                   ii = ratio(0)*i
                   do P=0,ratio(1)-1
                      do L=0,ratio(0)-1
                         fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
#endif

  end subroutine PC_EDGE_INTERP

  subroutine EDGE_INTERP(flo, fhi, nc, ratio, dir, &
                         fine, f_lo, f_hi, nc_f) bind(C,name="EDGE_INTERP")

    implicit none

    integer         , intent(in   ) :: flo(0:2), fhi(0:2), nc, ratio(0:AMREX_SPACEDIM-1), dir
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    double precision, intent(inout) :: fine(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)

    integer i,j,k,n,P,M
    DOUBLE PRECISION val, df
#if (AMREX_SPACEDIM == 3)
    integer L
#endif
    !
    !     Do linear in dir, pc transverse to dir, leave alone the fine values
    !     lining up with coarse edges--assume these have been set to hold the
    !     values you want to interpolate to the rest.
    !

#if (AMREX_SPACEDIM == 2)
    k = 0
      if (dir.eq.0) then
         do n=1,nc
            do j=flo(1),fhi(1),ratio(1)
               do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                  df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,k,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                        fine(i+M,P,k,n) = val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=flo(1),fhi(1)-ratio(dir),ratio(1)
               do i=flo(0),fhi(0)
                  df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,k,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                        fine(P,j+M,k,n) = val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
#elif (AMREX_SPACEDIM == 3)
    if (dir.eq.0) then
       do n=1,nc
          do k=flo(2),fhi(2),ratio(2)
             do j=flo(1),fhi(1),ratio(1)
                do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                   df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                   do M=1,ratio(dir)-1
                      val = fine(i,j,k,n) &
                           + df*dble(M)/dble(ratio(dir))
                      do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                         do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                            fine(i+M,P,L,n) = val
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else if (dir.eq.1) then
       do n=1,nc
          do k=flo(2),fhi(2),ratio(2)
             do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                do i=flo(0),fhi(0)
                   df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                   do M=1,ratio(dir)-1
                      val = fine(i,j,k,n) &
                           + df*dble(M)/dble(ratio(dir))
                      do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                         do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                            fine(P,j+M,L,n) = val
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do n=1,nc
          do k=flo(2),fhi(2)-ratio(dir),ratio(2)
             do j=flo(1),fhi(1),ratio(1)
                do i=flo(0),fhi(0),ratio(0)
                   df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                   do M=1,ratio(dir)-1
                      val = fine(i,j,k,n) &
                           + df*dble(M)/dble(ratio(dir))
                      do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                         do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                            fine(P,L,k+M,n) = val
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
#endif

  end subroutine EDGE_INTERP

end module fill_umac_ghost_module
