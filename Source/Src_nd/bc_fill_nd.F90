module bc_fill_module

  ! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  ! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  implicit none

contains

  subroutine phifill(phi,phi_lo,phi_hi,domlo,domhi,dx,gridlo,time,bc,icomp) &
       bind(C, name="phifill")

    integer, intent(in)      :: phi_lo(3),phi_hi(3)
    integer, intent(in)      :: bc(AMREX_SPACEDIM,2)
    integer, intent(in)      :: domlo(3), domhi(3)
    double precision, intent(in) :: dx(3), gridlo(3), time
    double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    integer, value, intent(in) :: icomp

#if (AMREX_SPACEDIM == 1)
    call filcc(phi,phi_lo(1),phi_hi(1),domlo,domhi,dx,gridlo,bc)
#elif (AMREX_SPACEDIM == 2)
    call filcc(phi,phi_lo(1),phi_lo(2),phi_hi(1),phi_hi(2),domlo,domhi,dx,gridlo,bc)
#else
    call filcc(phi,phi_lo(1),phi_lo(2),phi_lo(3),phi_hi(1),phi_hi(2),phi_hi(3),domlo,domhi,dx,gridlo,bc)
#endif

    call fill_scalar_ext_bc(phi_lo,phi_hi,phi,phi_lo,phi_hi,domlo,domhi,bc)

  end subroutine phifill

  subroutine velfill(vel,vel_lo,vel_hi,domlo,domhi,dx,gridlo,time,bc,icomp) &
       bind(C, name="velfill")

    integer, intent(in)      :: vel_lo(3),vel_hi(3)
    integer, intent(in)      :: bc(AMREX_SPACEDIM,2)
    integer, intent(in)      :: domlo(3), domhi(3)
    double precision, intent(in) :: dx(3), gridlo(3), time
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3))
    integer, value, intent(in) :: icomp

#if (AMREX_SPACEDIM == 1)
    call filcc(vel,vel_lo(1),vel_hi(1),domlo,domhi,dx,gridlo,bc)
#elif (AMREX_SPACEDIM == 2)
    call filcc(vel,vel_lo(1),vel_lo(2),vel_hi(1),vel_hi(2),domlo,domhi,dx,gridlo,bc)
#else
    call filcc(vel,vel_lo(1),vel_lo(2),vel_lo(3),vel_hi(1),vel_hi(2),vel_hi(3),domlo,domhi,dx,gridlo,bc)
#endif

    call fill_vel_ext_bc(vel_lo,vel_hi,vel,vel_lo,vel_hi,domlo,domhi,bc)

  end subroutine velfill

  subroutine fill_scalar_ext_bc(lo,hi,q,q_lo,q_hi,domlo,domhi,bc,icomp)

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(in   ) :: q_lo(3),q_hi(3)
    integer, intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer, intent(in   ) :: domlo(3), domhi(3)
    double precision, intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    integer, value, intent(in) :: icomp

    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, n
    integer :: imin, imax, jmin, jmax, kmin, kmax

    ! do nothing if there are no exterior boundaries
    if (.not. any(bc .eq. amrex_bc_ext_dir)) then
       return
    endif

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    ilo = domlo(1)
    ihi = domhi(1)

#if AMREX_SPACEDIM >= 2
    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)
#endif

#if AMREX_SPACEDIM == 3
    ks = max(q_lo(3), domlo(3))
    ke = min(q_hi(3), domhi(3))
    klo = domlo(3)
    khi = domhi(3)
#endif

    if (lo(1) < ilo) then
       imin = lo(1)
       imax = ilo-1

       if (bc(1,1) .eq. amrex_bc_ext_dir) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = imin, imax
                   q(i,j,k) = 0.d0
                end do
             end do
          end do

       end if
    end if

    if (hi(1) > ihi) then
       imin = ihi+1
       imax = hi(1)

       if (bc(1,2) .eq. amrex_bc_ext_dir) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = imin, imax
                   q(i,j,k) = 0.d0
                end do
             end do
          end do

       end if
    end if

#if AMREX_SPACEDIM >= 2

    if (lo(2) < jlo) then
       jmin = lo(2)
       jmax = jlo-1

       if (bc(2,1) .eq. amrex_bc_ext_dir) then

          do k = lo(3), hi(3)
             do j = jmin, jmax
                do i = lo(1), hi(1)
                   q(i,j,k) = 0.d0
                end do
             end do
          end do

       end if
    end if

    if (hi(2) > jhi) then
       jmin = jhi+1
       jmax = hi(2)

       if (bc(2,2) .eq. amrex_bc_ext_dir) then

          do k = lo(3), hi(3)
             do j = jmin, jmax
                do i = lo(1), hi(1)
                   q(i,j,k) = 0.d0
                end do
             end do
          end do

       end if
    end if
#endif

#if AMREX_SPACEDIM == 3

    if (lo(3) < klo) then
       kmin = lo(3)
       kmax = klo-1

       if (bc(3,1) .eq. amrex_bc_ext_dir) then

          do k = kmin, kmax
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   q(i,j,k) = 0.d0
                end do
             end do
          end do

       end if
    end if

    if (hi(3) > khi) then
       kmin = khi+1
       kmax = hi(3)

       if (bc(3,2) .eq. amrex_bc_ext_dir) then

          do k = kmin, kmax
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   q(i,j,k) = 0.d0
                end do
             end do
          end do

       end if
    end if
#endif

  end subroutine fill_scalar_ext_bc

  subroutine fill_vel_ext_bc(lo,hi,v,v_lo,v_hi,domlo,domhi,bc,icomp)

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(in   ) :: v_lo(3),v_hi(3)
    integer, intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer, intent(in   ) :: domlo(3), domhi(3)
    double precision, intent(inout) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    integer, value, intent(in) :: icomp

    call fill_scalar_ext_bc(lo,hi,v,v_lo,v_hi,domlo,domhi,bc,icomp)

  end subroutine fill_vel_ext_bc

end module bc_fill_module
