module bc_fill_module

  ! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  ! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  use amrex_error_module
  use amrex_bc_types_module
  
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

  end subroutine velfill

end module bc_fill_module
