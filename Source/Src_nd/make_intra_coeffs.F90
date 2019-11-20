! create the coefficients needed for I_T (intra for the SDC
! predict_T_* enthalpy updates).  We take an old and new state
! and return the time-centered quantities.


module make_intra_coeffs_module

  use meth_params_module, only: rho_comp, temp_comp, spec_comp, nscal
  use eos_module
  use eos_type_module
  use network, only: nspec
  use amrex_constants_module

  implicit none

  private

contains 
  
  subroutine make_intra_coeffs(lo, hi, &
       scalold, s1_lo, s1_hi, &
       scalnew, s2_lo, s2_hi, &
       cp, cp_lo, cp_hi, &
       xi, xi_lo, xi_hi) bind(C,name="make_intra_coeffs")
    ! input box is grown to include boundary cell information

    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module
    use eos_composition_module, only : eos_xderivs_t, composition_derivatives

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s1_lo(3), s1_hi(3)
    double precision, intent(in   ) :: scalold(s1_lo(1):s1_hi(1),s1_lo(2):s1_hi(2),s1_lo(3):s1_hi(3),nscal)
    integer         , intent(in   ) :: s2_lo(3), s2_hi(3)
    double precision, intent(in   ) :: scalnew(s2_lo(1):s2_hi(1),s2_lo(2):s2_hi(2),s2_lo(3):s2_hi(3),nscal)
    integer         , intent(in   ) :: cp_lo(3), cp_hi(3)
    double precision, intent(inout) :: cp(cp_lo(1):cp_hi(1),cp_lo(2):cp_hi(2),cp_lo(3):cp_hi(3))
    integer         , intent(in   ) :: xi_lo(3), xi_hi(3)
    double precision, intent(inout) :: xi(xi_lo(1):xi_hi(1),xi_lo(2):xi_hi(2),xi_lo(3):xi_hi(3),nspec)

    ! local
    integer :: i,j,k,comp
    type (eos_t) :: eos_state
    type (eos_xderivs_t) :: eos_xderivs

    !$gpu
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             ! old state first
             eos_state%rho   = scalold(i,j,k,rho_comp)
             eos_state%T     = scalold(i,j,k,temp_comp)
             eos_state%xn(:) = scalold(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state)

             call composition_derivatives(eos_state, eos_xderivs)

             cp(i,j,k) = eos_state%cp
             
             do comp=1,nspec
                xi(i,j,k,comp) = eos_xderivs % dhdX(comp)
             enddo

             ! new state now -- average results
             eos_state%rho   = scalnew(i,j,k,rho_comp)
             eos_state%T     = scalnew(i,j,k,temp_comp)
             eos_state%xn(:) = scalnew(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state)

             call composition_derivatives(eos_state, eos_xderivs)
             
             cp(i,j,k) = HALF*(eos_state%cp + cp(i,j,k))
             
             do comp=1,nspec
                xi(i,j,k,comp) = HALF*(eos_xderivs % dhdX(comp) + xi(i,j,k,comp))
             enddo

          enddo
       enddo
    enddo
    
  end subroutine make_intra_coeffs

end module make_intra_coeffs_module
