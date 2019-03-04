module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  
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

    if (parallel_IOProcessor()) then
       print*,"If make_heating=T, create a local copy of make_heating.F90 in your build directory"
    end if
    
    ! abort program
    call abort()
      
  end subroutine make_heating

end module make_heating_module

