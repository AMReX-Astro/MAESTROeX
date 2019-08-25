module make_explicit_thermal_module

  use eos_type_module
  use eos_module
  use conductivity_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, nscal, &
       buoyancy_cutoff_factor, base_cutoff_density, &
       limit_conductivity
  use amrex_constants_module

  implicit none

  private

contains

  subroutine make_thermal_coeffs(lo, hi, &
       scal, s_lo, s_hi, &
       Tcoeff, t_lo, t_hi, &
       hcoeff, h_lo, h_hi, &
       Xkcoeff, xk_lo, xk_hi, &
       pcoeff, p_lo, p_hi) bind(C,name="make_thermal_coeffs")
    ! create the coefficients for grad{T}, grad{h}, grad{X_k}, and grad{p_0}
    ! for the thermal diffusion term in the enthalpy equation.
    !
    ! note: we explicitly fill the ghostcells by looping over them directly

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(inout) :: Tcoeff(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer         , intent(in   ) :: h_lo(3), h_hi(3)
    double precision, intent(inout) :: hcoeff(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    integer         , intent(in   ) :: xk_lo(3), xk_hi(3)
    double precision, intent(inout) :: Xkcoeff(xk_lo(1):xk_hi(1),xk_lo(2):xk_hi(2),xk_lo(3):xk_hi(3),nspec)
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(inout) :: pcoeff(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

    ! Local
    integer :: i,j,k,comp
    type (eos_t) :: eos_state

    !$gpu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! create Tcoeff = -kth,
    !        hcoeff = -kth/cp,
    !       Xkcoeff = xik*kth/cp,
    !        pcoeff = hp*kth/cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print *, "... Level ", lev, " create thermal coeffs:"

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! if we are outside the star, turn off the conductivity
             if (limit_conductivity .and. &
                  scal(i,j,k,rho_comp) < &
                  buoyancy_cutoff_factor*base_cutoff_density) then

                Tcoeff(i,j,k) = ZERO
                hcoeff(i,j,k) = ZERO
                pcoeff(i,j,k) = ZERO
                Xkcoeff(i,j,k,:) = ZERO

             else

                eos_state%rho   = scal(i,j,k,rho_comp)
                eos_state%T     = scal(i,j,k,temp_comp)
                eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! dens, temp, and xmass are inputs
                call conducteos(eos_input_rt, eos_state)

                Tcoeff(i,j,k) = -eos_state % conductivity
                hcoeff(i,j,k) = -eos_state % conductivity/eos_state%cp
                pcoeff(i,j,k) = (eos_state % conductivity/eos_state%cp)* &
                     ((1.0d0/eos_state%rho)* &
                     (1.0d0-eos_state%p/(eos_state%rho*eos_state%dpdr))+ &
                     eos_state%dedr/eos_state%dpdr)

                do comp=1,nspec
                   Xkcoeff(i,j,k,comp) = (eos_state % conductivity/eos_state%cp)*eos_state%dhdX(comp)
                enddo

             endif

          enddo
       enddo
    enddo

  end subroutine make_thermal_coeffs

end module make_explicit_thermal_module
