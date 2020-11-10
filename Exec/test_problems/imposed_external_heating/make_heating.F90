! This routine returns the externally imposed (i.e. not reactions)
! Heating source term to the enthalpy equation (rho * H )

module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use meth_params_module, only: prob_lo, rho_comp
  use amrex_constants_module

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
    double precision, parameter :: pi_ = 3.14159265358979323844d0
    ! local
    integer :: i,j,k
    double precision :: x,y,z,H,a,b,c,Ts, EE, BOT, FF
    H  = 10.d0
    a  = 2.d0
    b  =  2.d0
    Ts = 1.d0
    c  = H / pi_
    EE = exp(-time/Ts)
    do k = lo(3),hi(3)
      do j = lo(2),hi(2)
        y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)
        do i = lo(1),hi(1)
          x =  (dble(i)+0.5d0)*dx(1) + prob_lo(1)

          BOT = EE*cos(y/c) + b
          FF = exp(-2*time/Ts)*(sin(y/c)**2)/(BOT**2)
          FF = FF + (EE*cos(y/c) / BOT)
          rho_Hext(i,j,k) = scal(i,j,k,rho_comp) * FF/Ts
        end do
      end do
    end do

  end subroutine make_heating

end module make_heating_module

