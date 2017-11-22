module bc_fill_module

! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  implicit none

contains

  subroutine phifill(phi,phi_lo,phi_hi,domlo,domhi,dx,gridlo,time,bc) &
       bind(C, name="phifill")

    integer      :: phi_lo(3),phi_hi(3)
    integer      :: bc(AMREX_SPACEDIM,2)
    integer      :: domlo(3), domhi(3)
    double precision :: dx(3), gridlo(3), time
    double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

#if (AMREX_SPACEDIM == 1)
       call filcc(phi,phi_lo(1),phi_hi(1),domlo,domhi,dx,gridlo,bc)
#elif (AMREX_SPACEDIM == 2)
       call filcc(phi,phi_lo(1),phi_lo(2),phi_hi(1),phi_hi(2),domlo,domhi,dx,gridlo,bc)
#else
       call filcc(phi,phi_lo(1),phi_lo(2),phi_lo(3),phi_hi(1),phi_hi(2),phi_hi(3),domlo,domhi,dx,gridlo,bc)
#endif

  end subroutine phifill
  
end module bc_fill_module
