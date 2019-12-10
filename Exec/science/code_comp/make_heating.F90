module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use meth_params_module, only: prob_lo, rho_comp
  use base_state_geometry_module, only: center
  use probin_module, only: heating_factor
  use amrex_constants_module, only: ZERO, HALF, ONE, M_PI

  implicit none

  private

contains

  subroutine make_heating(lo, hi, &
                          rho_Hext, r_lo, r_hi, &
                          scal,     s_lo, s_hi, nc_s, &
                          dx, time ) bind (C,name="make_heating")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3), nc_s
    double precision, intent (inout) :: rho_Hext(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)     )
    double precision, intent (in   ) ::     scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in   ) :: dx(3), time

    integer :: i, j, k, r
    double precision :: z, fheat

    rho_Hext(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
#if (AMREX_SPACEDIM == 2)
            r = j
#else
            r = k
#endif
            z = prob_lo(AMREX_SPACEDIM) + (dble(r) + HALF) * dx(AMREX_SPACEDIM)

            if (z < 1.125d0 * 4.d8) then 
                fheat = sin(8.d0 * M_PI * (z / 4.d8 - ONE))
    
                do i = lo(1), hi(1)
        
                   ! Source terms
                   rho_Hext(i,j,k) = heating_factor * fheat 
        
                end do
            endif
        enddo
    enddo

  end subroutine make_heating

end module make_heating_module
