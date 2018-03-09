module make_heating_module

  use parallel, only: parallel_IOProcessor
  use amrex_constants_module, only: M_PI
  use meth_params_module, only: prob_lo, rho_comp

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

    ! local
    integer :: i,j,k

    double precision :: x,y,z
    double precision :: y_layer,z_layer
    double precision :: ey,ez
    double precision :: L_x,L_y

#if (AMREX_SPACEDIM == 2)  
    
    L_x = 2.5d8

    if (time <= 200.0d0) then

       y_layer = 1.25d8

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)
          ey = exp(-(y-y_layer)*(y-y_layer)/1.d14)
          do i = lo(1),hi(1)
             x =  (dble(i)+0.5d0)*dx(1) + prob_lo(1)

             rho_Hext(i,j,k) = ey*(1.d0 + &
                  .00625d0 * sin(2*M_PI*x/L_x) &
                  + .01875d0 * sin((6*M_PI*x/L_x) + M_PI/3.d0) &
                  + .01250d0 * sin((8*M_PI*x/L_x) + M_PI/5.d0))*2.5d16

             rho_Hext(i,j,k) = rho_Hext(i,j,k) * scal(i,j,k,rho_comp)
          end do
       end do
       end do

    end if
#elif (AMREX_SPACEDIM == 3)

    L_x = 2.5d8
    L_y = 2.5d8

    if (time <= 200.0d0) then

       z_layer = 1.25d8

       do k = lo(3),hi(3)
          z = (dble(k)+0.5d0)*dx(3) + prob_lo(3)
          ez = exp(-(z-z_layer)*(z-z_layer)/1.d14)

          do j = lo(2),hi(2)
             y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)

             do i = lo(1),hi(1)
                x =  (dble(i)+0.5d0)*dx(1) + prob_lo(1)

                rho_Hext(i,j,k) = ez*(1.d0 + &
                     .00625d0 * sin(2*M_PI*x/L_x) * sin(2*M_PI*y/L_y) &
                     + .01875d0 * sin((6*M_PI*x/L_x) + M_PI/3.d0) * sin((6*M_PI*y/L_y) + M_PI/3.d0) &
                     + .01250d0 * sin((8*M_PI*x/L_x) + M_PI/5.d0) * sin((8*M_PI*y/L_y) + M_PI/5.d0))*2.5d16


                rho_Hext(i,j,k) = rho_Hext(i,j,k) * scal(i,j,k,rho_comp)
                
             end do
          end do
       enddo

    end if

#endif

  end subroutine make_heating

end module make_heating_module

