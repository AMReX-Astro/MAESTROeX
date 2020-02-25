module make_scal_force_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr, base_cutoff_density_coord
  use meth_params_module, only: enthalpy_pred_type, use_exact_base_state, base_cutoff_density, spherical, &
       nscal, rho_comp, temp_comp, spec_comp
  use amrex_constants_module

  implicit none

  private

contains

  subroutine mktempforce(lo, hi, lev, &
       temp_force, f_lo, f_hi, &
       scal, s_lo, s_hi, &
       umac, u_lo, u_hi, &
       vmac, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wmac, w_lo, w_hi, &
#endif
       thermal, t_lo, t_hi, &
       p0_cart, p_lo, p_hi, &
       psi_cart, ps_lo, ps_hi, &
       dx, domhi)  bind(C,name="mktempforce")

    use eos_module
    use eos_type_module
    use network, only: nspec
    
    integer         , intent(in   ) :: lo(3),hi(3)
    integer, value  , intent(in   ) :: lev
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: temp_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(in   ) :: thermal(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(in   ) :: p0_cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    integer         , intent(in   ) :: ps_lo(3), ps_hi(3)
    double precision, intent(in   ) :: psi_cart(ps_lo(1):ps_hi(1),ps_lo(2):ps_hi(2),ps_lo(3):ps_hi(3))
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: domhi(3)

    ! Local variable
    double precision :: gradp0,veladv,dhdp
    integer          :: i,j,k
    
    type (eos_t) :: eos_state

    ! spherical only variables
    double precision :: p0_lox,p0_hix,p0_loy,p0_hiy,p0_loz,p0_hiz
    double precision :: divup,p0divu,ugradp
    
    integer, parameter :: predict_rhoh             = 0;
    integer, parameter :: predict_rhohprime        = 1;
    integer, parameter :: predict_h                = 2;
    integer, parameter :: predict_T_then_rhohprime = 3;
    integer, parameter :: predict_T_then_h         = 4;
    integer, parameter :: predict_hprime           = 5;
    integer, parameter :: predict_Tprime_then_h    = 6;

    !$gpu

    ! For non-spherical, add wtilde d(p0)/dr
    ! For spherical, we make u grad p = div (u p) - p div (u)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            
#if (AMREX_SPACEDIM == 2)

             if (j .eq. 0) then
                gradp0 = ( p0_cart(i,j+1,k) - p0_cart(i,j,k) ) / dx(2)
             else if (j .eq. domhi(2)) then
                gradp0 = ( p0_cart(i,j,k) - p0_cart(i,j-1,k) ) / dx(2)
             else
                gradp0 = HALF*( p0_cart(i,j+1,k) - p0_cart(i,j-1,k) ) / dx(2)
             end if

             eos_state%T     = scal(i,j,k,temp_comp)
             eos_state%rho   = scal(i,j,k,rho_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1) / scal(i,j,k,rho_comp)

             ! dens, temp, xmass inputs
             call eos(eos_input_rt, eos_state)
             
             dhdp = ONE / scal(i,j,k,rho_comp) + ( scal(i,j,k,rho_comp) * eos_state%dedr - &
                                                   eos_state%p / scal(i,j,k,rho_comp) ) &
                                                 / ( scal(i,j,k,rho_comp) * eos_state%dpdr )

             veladv = HALF*(vmac(i,j,k)+vmac(i,j+1,k))
             
             temp_force(i,j,k) =  thermal(i,j,k) + (ONE - scal(i,j,k,rho_comp) * dhdp) * &
                                                   (veladv * gradp0 + psi_cart(i,j,k))
             temp_force(i,j,k) = temp_force(i,j,k) / (eos_state%cp * scal(i,j,k,rho_comp))

#else
             if (spherical .eq. 0) then 

                if (k .eq. 0) then
                   gradp0 = (p0_cart(i,j,k+1) - p0_cart(i,j,k) ) / dx(3)
                else if (k .eq. domhi(3)) then
                   gradp0 = (p0_cart(i,j,k) - p0_cart(i,j,k-1) ) / dx(3)
                else
                   gradp0 = HALF*(p0_cart(i,j,k+1) - p0_cart(i,j,k-1) ) / dx(3)
                end if
                
                eos_state%T     = scal(i,j,k,temp_comp)
                eos_state%rho   = scal(i,j,k,rho_comp)
                eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1) / scal(i,j,k,rho_comp)

                ! dens, temp, xmass inputs
                call eos(eos_input_rt, eos_state)
                
                dhdp = ONE / scal(i,j,k,rho_comp) + ( scal(i,j,k,rho_comp) * eos_state%dedr - &
                     eos_state%p / scal(i,j,k,rho_comp) ) / ( scal(i,j,k,rho_comp) * eos_state%dpdr )
                
                veladv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
                
                temp_force(i,j,k) =  thermal(i,j,k) + &
                     (ONE - scal(i,j,k,rho_comp) * dhdp) * (veladv * gradp0 + psi_cart(i,j,k))
                
                temp_force(i,j,k) = temp_force(i,j,k) / (eos_state%cp * scal(i,j,k,rho_comp))
             
             else 
                
                eos_state%T     = scal(i,j,k,temp_comp)
                eos_state%rho   = scal(i,j,k,rho_comp)
                eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1) / scal(i,j,k,rho_comp)
                                
                ! dens, temp, xmass inputs
                call eos(eos_input_rt, eos_state)
                
                dhdp = ONE / scal(i,j,k,rho_comp) + ( scal(i,j,k,rho_comp) * eos_state%dedr - &
                                                      eos_state%p / scal(i,j,k,rho_comp) ) &
                                                    / ( scal(i,j,k,rho_comp) * eos_state%dpdr )

                p0_lox = HALF * (p0_cart(i,j,k) + p0_cart(i-1,j,k)) 
                p0_hix = HALF * (p0_cart(i,j,k) + p0_cart(i+1,j,k)) 
                p0_loy = HALF * (p0_cart(i,j,k) + p0_cart(i,j-1,k)) 
                p0_hiy = HALF * (p0_cart(i,j,k) + p0_cart(i,j+1,k)) 
                p0_loz = HALF * (p0_cart(i,j,k) + p0_cart(i,j,k-1)) 
                p0_hiz = HALF * (p0_cart(i,j,k) + p0_cart(i,j,k+1)) 
                
                divup = (umac(i+1,j,k) * p0_hix - umac(i,j,k) * p0_lox) / dx(1) + &
                        (vmac(i,j+1,k) * p0_hiy - vmac(i,j,k) * p0_loy) / dx(2) + &
                        (wmac(i,j,k+1) * p0_hiz - wmac(i,j,k) * p0_loz) / dx(3)

                p0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                           (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                           (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * p0_cart(i,j,k)
             
                ugradp = divup - p0divu
                
                temp_force(i,j,k) =  thermal(i,j,k) + &
                     (ONE - scal(i,j,k,rho_comp) * dhdp) * (ugradp + psi_cart(i,j,k))
                
                temp_force(i,j,k) = temp_force(i,j,k) / (eos_state%cp * scal(i,j,k,rho_comp))
                
             endif
#endif
          end do
       end do
    end do

  end subroutine mktempforce
  
end module make_scal_force_module
