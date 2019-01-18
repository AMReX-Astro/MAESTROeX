module make_analytic_solution_module

  use bl_error_module
  use bl_constants_module
  use base_state_geometry_module, only: center
  use probin_module, only: ambient_h, ambient_dens, &
       t0, peak_h, diffusion_coefficient
  use meth_params_module, only: prob_lo

  implicit none

  private

contains

  subroutine make_analytic_solution(lo, hi, solution, s_lo, s_hi, dx, time) bind(C, name="make_analytic_solution")

    integer, intent(in) :: lo(3), hi(3), s_lo(3), s_hi(3)
    double precision, intent(inout) :: solution(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ), value :: time

    integer :: n, i, j, k
    double precision :: xx, yy, dist

    solution(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3)) = ZERO

    do k = lo(3), hi(3)

       do j = lo(2), hi(2)

          yy = prob_lo(2) + (dble(j)+HALF) * dx(2)

          do i = lo(1), hi(1)

             xx = prob_lo(1) + (dble(i)+HALF) * dx(1)

             dist = sqrt((xx-center(1))**2 + (yy-center(2))**2)

             solution(i,j,k) = f(time,dist)

          enddo
       enddo
    enddo

  end subroutine make_analytic_solution

  function f(t,x) result(r)
    double precision :: t, x
    double precision :: r

    r = (peak_h-ambient_h)*(t0/(t+t0)) * &
         exp(-x*x/(FOUR*diffusion_coefficient*(t+t0))) + ambient_h

  end function f

end module make_analytic_solution_module
