
module initdata_module

  use meth_params_module, only: prob_lo
  use amrex_constants_module

  implicit none

  private

contains

  subroutine init_vel(lo, hi, &
       vel, vel_lo, vel_hi, nc_v, dx) bind(C, name="init_vel")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: dx(3)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3),1:nc_v)

    integer          :: i,j,k
    double precision :: x, y, z

    ! set velocity to zero
    vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3),1:nc_v) = 0.d0

    do k=lo(3),hi(3)
       z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)
       do j=lo(2),hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)
          do i=lo(1),hi(1)
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)

             vel(i,j,k,1) = TWO*M_PI*sin(FOUR*M_PI*x)*cos( TWO*M_PI*y) - &
                  FOUR*M_PI*sin( TWO*M_PI*x)*cos(FOUR*M_PI*z)
             vel(i,j,k,2) = TWO*M_PI*sin(FOUR*M_PI*y)*cos( TWO*M_PI*z) - &
                  FOUR*M_PI*cos(FOUR*M_PI*x)*sin( TWO*M_PI*y)
             vel(i,j,k,3) = TWO*M_PI*cos( TWO*M_PI*x)*sin(FOUR*M_PI*z) - &
                  FOUR*M_PI*cos(FOUR*M_PI*y)*sin( TWO*M_PI*z)

          end do
       end do
    end do

  end subroutine init_vel

end module initdata_module
