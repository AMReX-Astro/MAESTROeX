module debug_module


  implicit none

  private

contains

  subroutine print_cc(lev, lo, hi, &
                      mf, m_lo, m_hi, nc_m, ng_m) bind (C,name="print_cc")
    
    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3), nc_m, ng_m
    double precision, intent (in   ) :: mf(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nc_m)

    integer comp,i,j,k

    print*,"level",lev
    print*,"grid",lo(1:AMREX_SPACEDIM),hi(1:AMREX_SPACEDIM)
    
#if (AMREX_SPACEDIM == 1)

    k = 0
    j = 0

    ! loop over the data
    do comp = 1, nc_m
    do i = lo(1)-ng_m,hi(1)+ng_m
       print*,'i,comp',i,comp,mf(i,j,k,comp)
    enddo
    enddo

#elif (AMREX_SPACEDIM == 2)

    k = 0

    ! loop over the data
    do comp = 1, nc_m
    do j = lo(2)-ng_m,hi(2)+ng_m
    do i = lo(1)-ng_m,hi(1)+ng_m
       print*,'i,j,comp',i,j,comp,mf(i,j,k,comp)
    enddo
    enddo
    enddo

#elif (AMREX_SPACEDIM == 3)

    ! loop over the data
    do comp = 1, nc_m
    do k = lo(3)-ng_m,hi(3)+ng_m
    do j = lo(2)-ng_m,hi(2)+ng_m
    do i = lo(1)-ng_m,hi(1)+ng_m
       print*,'i,j,k,comp',i,j,k,comp,mf(i,j,k,comp)
    enddo
    enddo
    enddo
    enddo

#endif

    call flush(6)

  end subroutine print_cc

  subroutine print_edge(lev, dir, lo, hi, &
                        mf, m_lo, m_hi, nc_m, ng_m) bind (C,name="print_edge")
    
    integer         , intent (in   ) :: lev, dir, lo(3), hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3), nc_m, ng_m
    double precision, intent (in   ) :: mf(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nc_m)

    integer comp,i,j,k
    integer offset(3)

    offset(1:3) = 0
    offset(dir+1) = 1

    print*,"level",lev
    print*,"grid",lo(1:AMREX_SPACEDIM),hi(1:AMREX_SPACEDIM)
    
#if (AMREX_SPACEDIM == 1)

    k = 0
    j = 0

    ! loop over the data
    do comp = 1, nc_m
    do i = lo(1)-ng_m,hi(1)+ng_m+offset(1)
       print*,'i,comp',i,comp,mf(i,j,k,comp)
    enddo
    enddo

#elif (AMREX_SPACEDIM == 2)

    k = 0

    ! loop over the data
    do comp = 1, nc_m
    do j = lo(2)-ng_m,hi(2)+ng_m+offset(2)
    do i = lo(1)-ng_m,hi(1)+ng_m+offset(1)
       print*,'i,j,comp',i,j,comp,mf(i,j,k,comp)
    enddo
    enddo
    enddo

#elif (AMREX_SPACEDIM == 3)

    ! loop over the data
    do comp = 1, nc_m
    do k = lo(3)-ng_m,hi(3)+ng_m+offset(3)
    do j = lo(2)-ng_m,hi(2)+ng_m+offset(2)
    do i = lo(1)-ng_m,hi(1)+ng_m+offset(1)
       print*,'i,j,k,comp',i,j,k,comp,mf(i,j,k,comp)
    enddo
    enddo
    enddo
    enddo

#endif

    call flush(6)

  end subroutine print_edge

  subroutine print_nodal(lev, lo, hi, &
                         mf, m_lo, m_hi, nc_m, ng_m) bind (C,name="print_nodal")
    
    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3), nc_m, ng_m
    double precision, intent (in   ) :: mf(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nc_m)

    integer comp,i,j,k
 
    print*,"level",lev
    print*,"grid",lo(1:AMREX_SPACEDIM),hi(1:AMREX_SPACEDIM)
    
#if (AMREX_SPACEDIM == 1)

    k = 0
    j = 0

    ! loop over the data
    do comp = 1, nc_m
    do i = lo(1)-ng_m,hi(1)+ng_m+1
       print*,'i,comp',i,comp,mf(i,j,k,comp)
    enddo
    enddo

#elif (AMREX_SPACEDIM == 2)

    k = 0

    ! loop over the data
    do comp = 1, nc_m
    do j = lo(2)-ng_m,hi(2)+ng_m+1
    do i = lo(1)-ng_m,hi(1)+ng_m+1
       print*,'i,j,comp',i,j,comp,mf(i,j,k,comp)
    enddo
    enddo
    enddo

#elif (AMREX_SPACEDIM == 3)

    ! loop over the data
    do comp = 1, nc_m
    do k = lo(3)-ng_m,hi(3)+ng_m+1
    do j = lo(2)-ng_m,hi(2)+ng_m+1
    do i = lo(1)-ng_m,hi(1)+ng_m+1
       print*,'i,j,k,comp',i,j,k,comp,mf(i,j,k,comp)
    enddo
    enddo
    enddo
    enddo

#endif

    call flush(6)

  end subroutine print_nodal

end module debug_module
