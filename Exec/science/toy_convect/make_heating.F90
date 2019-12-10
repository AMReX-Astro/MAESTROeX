module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use amrex_constants_module
  use network
  ! use probin_module, only: drive_initial_convection
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, drive_initial_convection

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
    integer, save :: ih1, ic12, in14, io16

    double precision :: rho, temp, T6, T613, X_CNO, X_1, g14, eps_CNO

    logical, save :: firstCall = .true.

    if (firstCall) then
       ih1 =  network_species_index("hydrogen-1")
       ic12 = network_species_index("carbon-12")
       in14 = network_species_index("nitrogen-14")
       io16 = network_species_index("oxygen-16")

       firstCall = .false.
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rho = scal(i,j,k,rho_comp)

             ! if (drive_initial_convection) then
             !    temp = tempbar_init(j)
             ! else
                temp = scal(i,j,k,temp_comp)
             ! endif

             T6 = temp / 1.0d6
             T613 = T6**THIRD

             ! total CNO abundance
             X_CNO = (scal(i,j,k,spec_comp-1+ic12) + &
                  scal(i,j,k,spec_comp-1+in14) + &
                  scal(i,j,k,spec_comp-1+io16)) / rho

             ! H abundance
             X_1 = scal(i,j,k,spec_comp-1+ih1) / rho

             ! CNO heating from Kippenhahn & Weigert, Eq. 18.65
             g14 = 1.d0 + 2.7d-3*T613 - 7.78d-3*T613**2 - 1.49d-4*T6
             eps_CNO = 8.67d27 * g14 * X_CNO * X_1 * rho * exp(-152.28d0/T613) / T613**2

             rho_Hext(i,j,k) = rho * eps_CNO
          enddo
       enddo
    enddo

  end subroutine make_heating

end module make_heating_module
