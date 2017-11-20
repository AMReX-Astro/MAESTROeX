module bc_fill_module

! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
! for e.g., #if (AMREX_SPACEDIM == 1) statements.

  implicit none

contains

  subroutine phifill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="phifill")

    use amrex_fort_module, only: amrex_spacedim

    implicit none

    integer      :: phi_lo(3),phi_hi(3)
    integer      :: bc(amrex_spacedim,2)
    integer      :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

#if (AMREX_SPACEDIM == 1)
       call filcc(phi,phi_lo(1),phi_hi(1),domlo,domhi,delta,xlo,bc)
#elif (AMREX_SPACEDIM == 2)
       call filcc(phi,phi_lo(1),phi_lo(2),phi_hi(1),phi_hi(2),domlo,domhi,delta,xlo,bc)
#else
       call filcc(phi,phi_lo(1),phi_lo(2),phi_lo(3),phi_hi(1),phi_hi(2),phi_hi(3),domlo,domhi,delta,xlo,bc)
#endif

  end subroutine phifill
  
end module bc_fill_module
