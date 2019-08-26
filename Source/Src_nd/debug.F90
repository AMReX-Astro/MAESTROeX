module debug_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use base_state_geometry_module, only: max_radial_level, nr_fine, r_start_coord, &
                                        r_end_coord, numdisjointchunks

  implicit none

  private

contains

  subroutine print_base_cc(base) bind(C, name="print_base_cc")

    double precision, intent(in   ) :: base(0:max_radial_level,0:nr_fine-1)

    integer :: n,i,r

    if (parallel_IOProcessor()) then

       do n=0,max_radial_level
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                print*,'base lev,r',n,r,base(n,r)
             end do
          end do
       end do

    end if

  end subroutine print_base_cc

  subroutine print_base_edge(base) bind(C, name="print_base_edge")

    double precision, intent(in   ) :: base(0:max_radial_level,0:nr_fine)

    integer :: n,i,r

    if (parallel_IOProcessor()) then
       do n=0,max_radial_level
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)+1
                print*,'base lev,r',n,r,base(n,r)
             end do
          end do
       end do
    end if

  end subroutine print_base_edge

  subroutine print_mf(lev, lo, hi, mf, m_lo, m_hi, nc_m) bind (C,name="print_mf")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3), nc_m
    double precision, intent (inout) :: mf(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nc_m)

    integer comp,i,j,k

    print*,"level",lev
    print*,"valid box",lo(1:AMREX_SPACEDIM),hi(1:AMREX_SPACEDIM)

    ! loop over the data
    do comp = 1, nc_m
    do k = m_lo(3),m_hi(3)
    do j = m_lo(2),m_hi(2)
    do i = m_lo(1),m_hi(1)
#if (AMREX_SPACEDIM == 2)
       print*,'lev,i,j,comp',lev,i,j,comp,mf(i,j,k,comp)
#elif (AMREX_SPACEDIM == 3)
       print*,'lev,i,j,k,comp',lev,i,j,k,comp,mf(i,j,k,comp)
#endif
    enddo
    enddo
    enddo
    enddo

    call flush(6)

  end subroutine print_mf

  subroutine print_edge(lev, lo, hi, mf, m_lo, m_hi, nc_m) bind (C,name="print_edge")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3), nc_m
    double precision, intent (inout) :: mf(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nc_m)

    integer comp,i,j,k

    print*,"level",lev
    print*,"valid box",lo(1:AMREX_SPACEDIM),hi(1:AMREX_SPACEDIM)

    ! loop over the data
    do comp = 1, nc_m
    do k = m_lo(3),m_hi(3)
    do j = m_lo(2),m_hi(2)
    do i = m_lo(1),m_hi(1)
#if (AMREX_SPACEDIM == 2)
       print*,'lev,i,j,comp',lev,i,j,comp,mf(i,j,k,comp)
#elif (AMREX_SPACEDIM == 3)
       print*,'lev,i,j,k,comp',lev,i,j,k,comp,mf(i,j,k,comp)
#endif
    enddo
    enddo
    enddo
    enddo

    call flush(6)

  end subroutine print_edge

end module debug_module
