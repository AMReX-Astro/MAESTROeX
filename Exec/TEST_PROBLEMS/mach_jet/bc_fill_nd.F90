module bc_fill_module

  ! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  ! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  use amrex_fort_module, only : amrex_real, amrex_spacedim, amrex_get_loop_bounds
  use amrex_bc_types_module
  use amrex_constants_module
  use meth_params_module, only: rho_comp, rhoh_comp, spec_comp, temp_comp
  use inlet_bc_module
  use base_state_geometry_module, only: dr_fine, dr, nr, max_radial_level
  use probin_module, only: inlet_mach

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
    call amrex_error("physbc_1d not written")
#elif (AMREX_SPACEDIM == 2)
    call physbc_2d(phi,phi_lo,phi_hi,domlo,domhi,dx,&
         gridlo,bc,icomp)
#else
    call amrex_error("physbc_3d not written")
#endif

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
    call amrex_error("velphysbc_1d not written")
#elif (AMREX_SPACEDIM == 2)
    call velphysbc_2d(vel,vel_lo,vel_hi,domlo,domhi,dx,&
         gridlo,bc,icomp)
#else
    call amrex_error("velphysbc_3d not written")
#endif

  end subroutine velfill

  subroutine physbc_2d(phi,phi_lo,phi_hi,domlo,domhi,dx,gridlo,bc,icomp)

    ! use geometry, only: dr_fine
    ! use probin_module, only: inlet_mach

    integer, intent(in)      :: phi_lo(3),phi_hi(3)
    integer, intent(in)      :: bc(AMREX_SPACEDIM,2)
    integer, intent(in)      :: domlo(3), domhi(3)
    double precision, intent(in) :: dx(3), gridlo(3)
    double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    integer, value, intent(in) :: icomp

    !     Local variables
    integer :: i
    double precision :: A,B,x

    integer :: lo(3), hi(3)
    integer :: is, ie, js, je
    ! integer :: i, j, k, n
    integer :: jmin, jmax

    lo(:) = phi_lo(:)
    hi(:) = phi_hi(:)

    A = 4.5d-2
    B = 1.d2

    ! if (ng == 0) return

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) .eq. amrex_bc_int_dir) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("physbc_2d: bc(1,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) .eq. amrex_bc_int_dir) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("physbc_2d: bc(1,2) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------

    jmin = domlo(2)
    jmax = min(domhi(2),domlo(2)-1)

    if (bc(2,1) .eq. amrex_bc_ext_dir) then
       ! rho
       if (icomp .eq. rho_comp) phi(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = INLET_RHO
       ! rhoh
       if (icomp .eq. rhoh_comp) phi(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = INLET_RHOH
       ! species
       if (icomp .eq. spec_comp) phi(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = INLET_RHO
       ! temperature
       if (icomp .eq. temp_comp) phi(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = INLET_TEMP
       ! tracer
       ! if (icomp .eq. 7) phi(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = 0.d0
    else if (bc(2,1) .eq. amrex_bc_foextrap) then
       do i=lo(1),hi(1)
          phi(i,jmin:jmax,lo(3):hi(3)) = phi(i,lo(2),lo(3))
       end do
    else if (bc(2,1) .eq. amrex_bc_int_dir) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("physbc_2d: bc(2,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------

    jmin = max(domlo(2),domhi(2)+1)
    jmax = domhi(2)

    ! write(*,*) "bc(2,1) =, ", bc(2,1), "bc(2,2) = ", bc(2,2), "foextrap = ", amrex_bc_foextrap

    if (bc(2,2) .eq. amrex_bc_foextrap) then
       do i=lo(1),hi(1)
          phi(i,jmin:jmax,lo(3):hi(3)) = phi(i,hi(2),lo(3))
       end do
    else if (bc(2,2) .eq. amrex_bc_hoextrap) then
       do i=lo(1),hi(1)
          phi(i,jmin:jmax,lo(3):hi(3)) = &
               ( 15.d0 * phi(i,jmax  ,lo(3)) &
               -10.d0 * phi(i,jmax-1,lo(3)) &
               + 3.d0 * phi(i,jmax-2,lo(3)) ) * EIGHTH
       end do
    else if ((bc(2,2) .eq. amrex_bc_int_dir)) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("physbc_2d: bc(2,2) not yet supported")
    end if

  end subroutine physbc_2d

  subroutine velphysbc_2d(phi,phi_lo,phi_hi,domlo,domhi,dx,gridlo,bc,icomp)

    ! use geometry, only: dr_fine
    ! use probin_module, only: inlet_mach

    integer, intent(in)      :: phi_lo(3),phi_hi(3)
    integer, intent(in)      :: bc(AMREX_SPACEDIM,2)
    integer, intent(in)      :: domlo(3), domhi(3)
    double precision, intent(in) :: dx(3), gridlo(3)
    double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    integer, value, intent(in) :: icomp

    !     Local variables
    integer :: i
    double precision :: A,B,x

    integer :: lo(3), hi(3)
    integer :: is, ie, js, je
    ! integer :: i, j, k, n
    integer :: jmin, jmax

    lo(:) = phi_lo(:)
    hi(:) = phi_hi(:)

    A = 4.5d-2
    B = 1.d2

    ! if (ng == 0) return

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) .eq. amrex_bc_int_dir) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("velphysbc_2d: bc(1,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) .eq. amrex_bc_int_dir) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("velphysbc_2d: bc(1,2) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------

    jmin = domlo(2)
    jmax = min(domhi(2),domlo(2)-1)

    if (bc(2,1) .eq. amrex_bc_ext_dir) then
       ! xvel
       if (icomp .eq. 0) phi(lo(1):hi(1),jmin:jmax,lo(3):hi(3)) = 0.d0
       ! yvel
       if (icomp .eq. 1) then
          do i=lo(1),hi(1)
             x = (dble(i)+0.5d0)*dr_fine
             ! inflow is Mach number 0.01 front with a Mach number 0.1 bump in the middle
             phi(i,jmin:jmax,lo(3):hi(3)) = (inlet_mach/1.d-1)* &
                  INLET_CS*(1.d-2 + A*(tanh(B*(x-0.40d0)) + tanh(B*(0.6d0-x))))
          end do
       end if

    else if (bc(2,1) .eq. amrex_bc_foextrap) then
       do i=lo(1),hi(1)
          phi(i,jmin:jmax,lo(3):hi(3)) = phi(i,lo(2),lo(3))
       end do
    else if (bc(2,1) .eq. amrex_bc_int_dir) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("velphysbc_2d: bc(2,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------

    jmin = max(domlo(2),domhi(2)+1)
    jmax = domhi(2)

    ! write(*,*) "bc(2,1) =, ", bc(2,1), "bc(2,2) = ", bc(2,2), "foextrap = ", amrex_bc_foextrap

    if (bc(2,2) .eq. amrex_bc_foextrap) then
       do i=lo(1),hi(1)
          phi(i,jmin:jmax,lo(3):hi(3)) = phi(i,hi(2),lo(3))
       end do
    else if (bc(2,2) .eq. amrex_bc_hoextrap) then
       do i=lo(1),hi(1)
          phi(i,jmin:jmax,lo(3):hi(3)) = &
               ( 15.d0 * phi(i,jmax  ,lo(3)) &
               -10.d0 * phi(i,jmax-1,lo(3)) &
               + 3.d0 * phi(i,jmax-2,lo(3)) ) * EIGHTH
       end do
    else if ((bc(2,2) .eq. amrex_bc_int_dir)) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call amrex_error("velphysbc_2d: bc(2,2) not yet supported")
    end if

end subroutine velphysbc_2d

end module bc_fill_module
