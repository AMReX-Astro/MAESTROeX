module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use meth_params_module, only: prob_lo
  use base_state_geometry_module, only: center
  use probin_module, only: heating_factor

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

    integer :: i, j, k
    double precision :: x, y, z, r

    do k = lo(3), hi(3)
        z = prob_lo(3) + (dble(k) + 0.5d0) * dx(3) - center(3)
        do j = lo(2), hi(2)
            y = prob_lo(2) + (dble(j) + 0.5d0) * dx(2) - center(2)
            do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i) + 0.5d0) * dx(1) - center(1)

                r = sqrt(x**2 + y**2 + z**2)

                rho_Hext(i,j,k) = heating_factor * 6.7d5 * exp(-(r/4.525d10)**2)
            enddo
        enddo
    enddo

  end subroutine make_heating

end module make_heating_module
