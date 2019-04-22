module bc_fill_module

  ! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  ! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  use amrex_error_module
  use amrex_bc_types_module
  use amrex_constants_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, rhoh_comp, spec_comp, temp_comp
  use inlet_bc_module, only: INLET_RHO, INLET_RHOH, INLET_TEMP, INLET_VEL, INLET_RHOX
  use base_state_geometry_module, only: dr_fine
  ! use probin_module, only: inlet_mach

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

#if (AMREX_SPACEDIM == 1)
    call filcc(scal,scal_lo(1),scal_hi(1),domlo,domhi,dx,gridlo,bc)
#elif (AMREX_SPACEDIM == 2)
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

#if (AMREX_SPACEDIM == 1)
    call filcc(vel,vel_lo(1),vel_hi(1),domlo,domhi,dx,gridlo,bc)
#elif (AMREX_SPACEDIM == 2)
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
    double precision :: A, B, x, vy

    ! do nothing if there are no exterior boundaries
    if (.not. any(bc .eq. amrex_bc_ext_dir)) then
       return
    endif

    A = 4.5d-2
    B = 1.d2

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
                   q(i,j,k) = ONE
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
                   q(i,j,k) = ONE
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

           ! we only want to apply the BCs where the y-velocity is non-zero
           do i=lo(1),hi(1)
              if (INLET_VEL .ne. 0.0d0) then
                  ! rho
                  if (icomp .eq. rho_comp) then
                     q(i,jmin:jmax,lo(3):hi(3)) = INLET_RHO
                     ! rhoh
                  else if (icomp .eq. rhoh_comp) then
                     q(i,jmin:jmax,lo(3):hi(3)) = INLET_RHOH
                     ! species
                 else if (icomp >= spec_comp .and. icomp < spec_comp+nspec) then
                     q(i,jmin:jmax,lo(3):hi(3)) = INLET_RHOX(icomp-spec_comp+1)
                     ! temperature
                  else if (icomp .eq. temp_comp) then
                     q(i,jmin:jmax,lo(3):hi(3)) = INLET_TEMP
                  else
                     q(i,jmin:jmax,lo(3):hi(3)) = ZERO
                  endif
              else
                  q(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = ZERO
              endif
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
                   q(i,j,k) = ONE
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
                   q(i,j,k) = ONE
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
                   q(i,j,k) = ONE
                end do
             end do
          end do

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
    double precision :: A,B,x

    A = 4.5d-2
    B = 1.d2

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
                   q(i,j,k) = ONE
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
                   q(i,j,k) = ONE
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

          ! xvel
          if (icomp .eq. 1) q(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = ZERO
          ! yvel
          if (icomp .eq. 2) then
             do i=lo(1),hi(1)
                q(i,jmin:jmax,lo(3):hi(3)) = INLET_VEL
             end do
          end if

       end if
    end if

    if (hi(2) > jhi) then
       jmin = jhi+1
       jmax = hi(2)

       if (bc(2,2) .eq. amrex_bc_ext_dir) then

          do k = lo(3), hi(3)
             do j = jmin, jmax
                do i = lo(1), hi(1)
                   q(i,j,k) = ONE
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
                   q(i,j,k) = ONE
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
                   q(i,j,k) = ONE
                end do
             end do
          end do

       end if
    end if
#endif

  end subroutine fill_vel_ext_bc

end module bc_fill_module
