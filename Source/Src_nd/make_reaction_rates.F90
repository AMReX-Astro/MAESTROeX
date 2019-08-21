module make_reaction_rates_module

  use network, only: aion, nspec, nspec_evolve
  use actual_rhs_module, only: actual_rhs
  use meth_params_module, only: spec_comp, rhoh_comp, rho_comp, nscal
  use amrex_constants_module   , only: ZERO
  use burn_type_module, only: burn_t, burn_to_eos, eos_to_burn, net_ienuc
  use eos_type_module
  use eos_module
  
  implicit none

  private

contains

  subroutine instantaneous_reaction_rates(lo,hi, rho_omegadot,o_lo,o_hi, &
       rho_Hnuc,h_lo,h_hi, scal,s_lo,s_hi) bind(C,name="instantaneous_reaction_rates")

    integer,          intent(in   ) :: lo(3),hi(3)
    integer,          intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer,          intent(in   ) :: o_lo(3), o_hi(3)
    double precision, intent(inout) :: rho_omegadot(o_lo(1):o_hi(1),o_lo(2):o_hi(2),o_lo(3):o_hi(3),nspec)
    integer,          intent(in   ) :: h_lo(3), h_hi(3)
    double precision, intent(inout) :: rho_Hnuc(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))

    ! local
    integer :: i,j,k

    double precision :: temp_max, temp_min
    type (burn_t)    :: state
    type (eos_t)     :: eos_state
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! initialize state variables
             eos_state % rho = scal(i,j,k,rho_comp)
             eos_state % xn(1:nspec) = scal(i,j,k,spec_comp:spec_comp+nspec-1) / eos_state % rho
             eos_state % h   = scal(i,j,k,rhoh_comp) / eos_state % rho

             call eos_get_small_temp(temp_min)
             call eos_get_max_temp(temp_max)
             eos_state % T = sqrt(temp_min * temp_max)

             ! call the EOS with input rh to set T for rate evaluation
             call eos(eos_input_rh, eos_state)
             call eos_to_burn(eos_state, state)

             ! initialize arbitrary time
             state % time = ZERO

             call actual_rhs(state)

             rho_omegadot(i,j,k,1:nspec_evolve) = state % rho * aion(1:nspec_evolve) * &
                                                  state % ydot(1:nspec_evolve)
             rho_omegadot(i,j,k,nspec_evolve+1:nspec) = ZERO
             rho_Hnuc(i,j,k) = state % rho * state % ydot(net_ienuc)

          end do
       end do
    end do

  end subroutine instantaneous_reaction_rates

end module make_reaction_rates_module

