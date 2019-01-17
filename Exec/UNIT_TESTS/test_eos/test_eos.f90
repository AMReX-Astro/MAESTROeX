module test_eos_module

  use amrex_constants_module
  use network
  use eos_module
  use eos_type_module
  use meth_params_module, only: nscal, spec_comp, rho_comp, rhoh_comp, temp_comp
  use probin_module, only: dens_min, dens_max, temp_min, temp_max, metalicity_max

  implicit none

  private

contains

  subroutine do_tests(lo, hi, &
      scal, s_lo, s_hi, &
      error, e_lo, e_hi, nc_e, &
      dom_lo, dom_hi) bind(C, name="do_tests")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: e_lo(3), e_hi(3), nc_e
    integer, intent(in) :: dom_lo(3), dom_hi(3)
    double precision, intent(inout) :: scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),1:nscal)
    double precision, intent(inout) :: error(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),0:nc_e-1)

    integer :: i, j, k
    integer :: ih1, ihe4
    double precision :: temp_zone, dens_zone, metalicity
    double precision :: dlogrho, dlogT, dmetal
    double precision :: xn_zone(nspec)

    double precision :: h, p, e, s
    double precision :: tfromrh, rfromtp, tfromrp, tfromre, tfromps

    type (eos_t) :: eos_state
    
    dlogrho   = (log10(dens_max) - log10(dens_min))/(dom_hi(1) - dom_lo(1))
    dlogT     = (log10(temp_max) - log10(temp_min))/(dom_hi(2) - dom_lo(2))
    dmetal    = (metalicity_max  - ZERO           )/(dom_hi(3) - dom_lo(3))

    ! initialize the thermodynamic cube and do the initial EOS call
    ih1 = network_species_index('hydrogen-1')
    ihe4 = network_species_index('helium-4')

    do k = lo(3), hi(3)

       ! set the composition -- approximately solar
       metalicity = ZERO + dble(k)*dmetal
       xn_zone(:) = metalicity/(nspec - 2)   ! all but H, He
       xn_zone(ih1)  = 0.75d0 - HALF*metalicity
       xn_zone(ihe4) = 0.25d0 - HALF*metalicity

       do j = lo(2), hi(2)

          ! set the temperature
          temp_zone = 10.0d0**(log10(temp_min) + dble(j)*dlogT)

          do i = lo(1), hi(1)

             ! set the density
             dens_zone = 10.0d0**(log10(dens_min) + dble(i)*dlogrho)

             !------------------------------------------------------------
             ! call the EOS -- rho, T, X directly
             ! input: eos_input_rt
             !------------------------------------------------------------
             eos_state%T     = temp_zone
             eos_state%rho   = dens_zone
             eos_state%xn(:) = xn_zone(:)

             call eos(eos_input_rt, eos_state)

             ! store the thermodynamic state
             scal(i,j,k,rho_comp) = dens_zone
             scal(i,j,k,rhoh_comp) = dens_zone*eos_state%h
             scal(i,j,k,temp_comp) = temp_zone
             scal(i,j,k,spec_comp:spec_comp-1+nspec) = xn_zone(:)

             h = eos_state%h
             p = eos_state%p
             e = eos_state%e
             s = eos_state%s

             !------------------------------------------------------------
             ! call the EOS with rho, h
             ! input: eos_input_rh
             !------------------------------------------------------------

             ! change initial T guess to make the root find do some work
             eos_state%T     = HALF*temp_zone
             eos_state%rho   = dens_zone
             eos_state%h     = h
             eos_state%xn(:) = xn_zone(:)

             call eos(eos_input_rh, eos_state)

             ! store the thermodynamic state
             tfromrh = eos_state%T
             error(i,j,k,0) = error(i,j,k,0) + (eos_state%T - temp_zone)/temp_zone

             !------------------------------------------------------------
             ! call the EOS with T, p
             ! input: eos_input_tp
             !------------------------------------------------------------

             ! change initial rho guess to make the root find do some work
             eos_state%T     = temp_zone
             eos_state%rho   = THIRD*dens_zone
             eos_state%p     = p
             eos_state%xn(:) = xn_zone(:)

             call eos(eos_input_tp, eos_state)

             ! store the thermodynamic state
             rfromtp = eos_state%rho
             error(i,j,k,1) = error(i,j,k,1) + (eos_state%rho - dens_zone)/dens_zone

             !------------------------------------------------------------
             ! call the EOS with rho, p
             ! input: eos_input_rp
             !------------------------------------------------------------

             ! change initial T guess to make the root find do some work
             eos_state%T     = HALF*temp_zone
             eos_state%rho   = dens_zone
             eos_state%p     = p
             eos_state%xn(:) = xn_zone(:)

             call eos(eos_input_rp, eos_state)

             ! store the thermodynamic state
             tfromrp = eos_state%T
             error(i,j,k,2) = error(i,j,k,2) + (eos_state%T - temp_zone)/temp_zone

             !------------------------------------------------------------
             ! call the EOS with rho, e
             ! input: eos_input_re
             !------------------------------------------------------------

             ! change initial T guess to make the root find do some work
             eos_state%T     = HALF*temp_zone
             eos_state%rho   = dens_zone
             eos_state%e     = e
             eos_state%xn(:) = xn_zone(:)

             call eos(eos_input_re, eos_state)

             ! store the thermodynamic state
             tfromre = eos_state%T
             error(i,j,k,3) = error(i,j,k,3) + (eos_state%T - temp_zone)/temp_zone

             !------------------------------------------------------------
             ! call the EOS with p, s
             ! input: eos_input_ps
             !------------------------------------------------------------

             ! change initial T and rho guess to make the root find do
             ! some work
             eos_state%T     = HALF*temp_zone
             eos_state%rho   = HALF*dens_zone
             eos_state%p     = p
             eos_state%s     = s
             eos_state%xn(:) = xn_zone(:)

             ! some EOSes don't have physically valid treatments
             ! of entropy throughout the entire rho-T plane
             if (eos_state%s > ZERO) then

                call eos(eos_input_ps, eos_state)

                ! store the thermodynamic state
                tfromps = eos_state%T
                error(i,j,k,4) = error(i,j,k,4) + (eos_state%T - temp_zone)/temp_zone

             else
                ! do nothing

             endif

          enddo
       enddo
    enddo

  end subroutine do_tests


end module test_eos_module
