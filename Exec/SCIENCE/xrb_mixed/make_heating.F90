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
    integer, save :: ilh1, ilc12, iln14, ilo16, ihe4

    double precision :: rho, temp, T6, T613, X_CNO, X_1, g14, eps_CNO

    logical, save :: firstCall = .true.

    if (firstCall) then
       ilh1 =  network_species_index("hydrogen-1")
       ihe4 =  network_species_index("helium-4")
       ilc12 = network_species_index("carbon-12")
       iln14 = network_species_index("nitrogen-14")
       ilo16 = network_species_index("oxygen-16")

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
             X_CNO = (scal(i,j,k,spec_comp-1+ilc12) + &
                  scal(i,j,k,spec_comp-1+iln14) + &
                  scal(i,j,k,spec_comp-1+ilo16)) / rho

             ! H abundance
             X_1 = scal(i,j,k,spec_comp-1+ilh1) / rho

             ! CNO heating from Kippenhahn & Weigert, Eq. 18.65
             g14 = 1.d0 + 2.7d-3*T613 - 7.78d-3*T613**2 - 1.49d-4*T6
             eps_CNO = HCNO_energy_generation(X_CNO)

             eps_CNO = eps_CNO + triple_alpha_energy_generation(scal(i,j,k,spec_comp+ihe4-1)/rho,rho,temp)

             rho_Hext(i,j,k) = rho * eps_CNO
          enddo
       enddo
    enddo

  end subroutine make_heating

  ! the temperature-insensitive Hot CNO cycle energy generation
  function HCNO_energy_generation(XCNO) result(r)
    double precision, intent(in   ) :: XCNO
    double precision :: r

    ! this factor is from Wallace & Woosley (1981) ApJS 45
    double precision, parameter :: factor = 5.86d15

    r = factor * XCNO
  end function HCNO_energy_generation

  ! the triple alpha reaction energy generation
  ! taken from Arnett's book pg 225
  ! needs more accurate screening factor; just setting it to unity now
  function triple_alpha_energy_generation(Y,rho,T) result(r)
    double precision, intent(in   ) :: Y, rho, T
    double precision :: r

    double precision, parameter :: eps_0 = 3.9d11, f3a = 1.0d0
    double precision :: t8i, t8i3

    t8i = 1.d8/T
    t8i3 = t8i*t8i*t8i

    r = eps_0*f3a*rho*rho*Y*Y*Y*exp(-42.94d0*t8i)*t8i3

  end function triple_alpha_energy_generation

end module make_heating_module
