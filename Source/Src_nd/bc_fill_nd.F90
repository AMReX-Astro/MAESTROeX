module bc_fill_module

  ! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  ! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  use amrex_error_module
  use amrex_bc_types_module
  use meth_params_module, only: pi_comp

  implicit none

contains

  subroutine scalarfill(scal,scal_lo,scal_hi,domlo,domhi,dx,gridlo,time,bc,icomp) &
       bind(C, name="scalarfill")

    integer, intent(in)      :: scal_lo(3),scal_hi(3)
    integer, intent(in)      :: bc(AMREX_SPACEDIM,2)
    integer, intent(in)      :: domlo(3), domhi(3)
    double precision, intent(in) :: dx(3), gridlo(3), time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1),scal_lo(2):scal_hi(2),scal_lo(3):scal_hi(3))
    integer, value, intent(in) :: icomp

#if (AMREX_SPACEDIM == 2)
    call filcc(scal,scal_lo(1),scal_lo(2),scal_hi(1),scal_hi(2),domlo,domhi,dx,gridlo,bc)
#else
    call filcc(scal,scal_lo(1),scal_lo(2),scal_lo(3),scal_hi(1),scal_hi(2),scal_hi(3),domlo,domhi,dx,gridlo,bc)
#endif

    call fill_scalar_ext_bc(scal_lo,scal_hi,scal,scal_lo,scal_hi,domlo,domhi,bc,icomp)

  end subroutine scalarfill

  subroutine velfill(vel,vel_lo,vel_hi,domlo,domhi,dx,gridlo,time,bc,icomp) &
       bind(C, name="velfill")

    integer, intent(in)      :: vel_lo(3),vel_hi(3)
    integer, intent(in)      :: bc(AMREX_SPACEDIM,2)
    integer, intent(in)      :: domlo(3), domhi(3)
    double precision, intent(in) :: dx(3), gridlo(3), time
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3))
    integer, value, intent(in) :: icomp

#if (AMREX_SPACEDIM == 2)
    call filcc(vel,vel_lo(1),vel_lo(2),vel_hi(1),vel_hi(2),domlo,domhi,dx,gridlo,bc)
#else
    call filcc(vel,vel_lo(1),vel_lo(2),vel_lo(3),vel_hi(1),vel_hi(2),vel_hi(3),domlo,domhi,dx,gridlo,bc)
#endif

    call fill_vel_ext_bc(vel_lo,vel_hi,vel,vel_lo,vel_hi,domlo,domhi,bc,icomp)

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

    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)

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

          if (icomp .eq. pi_comp) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k) = 0.d0
                   end do
                end do
             end do
          else
             call amrex_error("bc_fill_nd.F90 bc(1,1) - must supply Dirichlet boundary conditions for scalar")
          end if

       end if
    end if

    if (hi(1) > ihi) then
       imin = ihi+1
       imax = hi(1)

       if (bc(1,2) .eq. amrex_bc_ext_dir) then

          if (icomp .eq. pi_comp) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k) = 0.d0
                   end do
                end do
             end do
          else
             call amrex_error("bc_fill_nd.F90 bc(1,2) - must supply Dirichlet boundary conditions for scalar")
          end if

       end if
    end if

    if (lo(2) < jlo) then
       jmin = lo(2)
       jmax = jlo-1

       if (bc(2,1) .eq. amrex_bc_ext_dir) then

          if (icomp .eq. pi_comp) then
             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k) = 0.d0
                   end do
                end do
             end do
          else
             call amrex_error("bc_fill_nd.F90 bc(2,1) - must supply Dirichlet boundary conditions for scalar")
          end if

       end if
    end if

    if (hi(2) > jhi) then
       jmin = jhi+1
       jmax = hi(2)

       if (bc(2,2) .eq. amrex_bc_ext_dir) then

          if (icomp .eq. pi_comp) then
             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k) = 0.d0
                   end do
                end do
             end do
          else
             call amrex_error("bc_fill_nd.F90 bc(2,2) - must supply Dirichlet boundary conditions for scalar")
          end if

       end if
    end if

#if AMREX_SPACEDIM == 3

    if (lo(3) < klo) then
       kmin = lo(3)
       kmax = klo-1

       if (bc(3,1) .eq. amrex_bc_ext_dir) then

          if (icomp .eq. pi_comp) then
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k) = 0.d0
                   end do
                end do
             end do
          else
             call amrex_error("bc_fill_nd.F90 bc(3,1) - must supply Dirichlet boundary conditions for scalar")
          end if

       end if
    end if

    if (hi(3) > khi) then
       kmin = khi+1
       kmax = hi(3)

       if (bc(3,2) .eq. amrex_bc_ext_dir) then

          if (icomp .eq. pi_comp) then
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k) = 0.d0
                   end do
                end do
             end do
          else
             call amrex_error("bc_fill_nd.F90 bc(3,2) - must supply Dirichlet boundary conditions for scalar")
          end if

       end if
    end if
#endif

  end subroutine fill_scalar_ext_bc

  subroutine fill_vel_ext_bc(lo,hi,q,q_lo,q_hi,domlo,domhi,bc,icomp)

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(in   ) :: q_lo(3),q_hi(3)
    integer, intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer, intent(in   ) :: domlo(3), domhi(3)
    double precision, intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    integer, value, intent(in) :: icomp

    !     Local variables
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

    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)

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

  end subroutine fill_vel_ext_bc

end module bc_fill_module
