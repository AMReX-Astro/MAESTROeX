module rhoh_vs_t_module

  use amrex_constants_module
  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       nscal, use_eos_e_instead_of_h, use_pprime_in_tfromp
  use base_state_geometry_module, only:  max_radial_level, nr_fine

  implicit none

  private
  
  integer, parameter :: predict_rhoh             = 0;
  integer, parameter :: predict_rhohprime        = 1;
  integer, parameter :: predict_h                = 2;
  integer, parameter :: predict_T_then_rhohprime = 3;
  integer, parameter :: predict_T_then_h         = 4;
  integer, parameter :: predict_hprime           = 5;
  integer, parameter :: predict_Tprime_then_h    = 6;

  integer, parameter :: predict_rhoprime_and_X   = 1;
  integer, parameter :: predict_rhoX             = 2;
  integer, parameter :: predict_rho_and_X        = 3;

contains

  subroutine makeTfromRhoH(lo,hi,state,s_lo,s_hi,p0_cart,p0_lo,p0_hi) &
       bind(C,name="makeTfromRhoH")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: p0_lo(3), p0_hi(3)
    double precision, intent (inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in) :: p0_cart (p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! (rho, (h->e)) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0_cart(i,j,k) / state(i,j,k,rho_comp)

                pt_index(:) = (/i, j, k/)

                call eos(eos_input_re, eos_state, pt_index)

                state(i,j,k,temp_comp) = eos_state%T

             enddo
          enddo
       enddo

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! (rho, h) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)

                pt_index(:) = (/i, j, k/)

                call eos(eos_input_rh, eos_state, pt_index)

                state(i,j,k,temp_comp) = eos_state%T

             enddo
          enddo
       enddo

    endif

  end subroutine makeTfromRhoH


  subroutine makeTfromRhoP(lo,hi,state,s_lo,s_hi,p0_cart,p0_lo,p0_hi,updateRhoH) &
       bind(C,name="makeTfromRhoP")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: p0_lo(3), p0_hi(3)
    double precision, intent (inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in) :: p0_cart (p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))
    integer ,  value, intent (in   ) :: updateRhoH

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, p) --> T

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p     = p0_cart(i,j,k) + state(i,j,k,pi_comp)
             else
                eos_state%p     = p0_cart(i,j,k)
             endif

             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rp, eos_state, pt_index)

             state(i,j,k,temp_comp) = eos_state%T

             if (updateRhoH .eq. 1) then
                state(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h
             end if

          enddo
       enddo
    enddo

  end subroutine makeTfromRhoP

  !----------------------------------------------------------------------------
  ! makePfromRhoH
  !----------------------------------------------------------------------------
  subroutine makePfromRhoH(lo, hi, &
       state, s_lo, s_hi,  &
       temp_old, t_lo, t_hi, &
       peos, p_lo, p_hi) bind(C,name="makePfromRhoH")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(in   ) :: temp_old(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(inout) :: peos(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T, p
             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = temp_old(i,j,k)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rh, eos_state, pt_index)

             peos(i,j,k) = eos_state%p

          enddo
       enddo
    enddo

  end subroutine makePfromRhoH


  subroutine makeMachfromRhoH(lo,hi,state,s_lo,s_hi,u,u_lo,u_hi, &
       p0_cart,p0_lo,p0_hi,w0cart,w_lo,w_hi,mach,m_lo,m_hi) bind(C,name="makeMachfromRhoH")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    double precision, intent (in   ) ::  u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    integer         , intent (in   ) :: p0_lo(3), p0_hi(3)
    double precision, intent (in) :: p0_cart (p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
    double precision, intent (in   ) :: w0cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    double precision, intent (inout) :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    double precision :: vel
    type (eos_t) :: eos_state

    !$gpu

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! vel is the magnitude of the velocity, including w0
#if (AMREX_SPACEDIM == 2)
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + w0cart(i,j,k))**2 )
#elif (AMREX_SPACEDIM == 3)
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + w0cart(i,j,k))**2 )
#endif

                ! (rho, (h->e)) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0_cart(i,j,k) / state(i,j,k,rho_comp)

                pt_index(:) = (/i, j, k/)

                call eos(eos_input_re, eos_state, pt_index)

                mach(i,j,k) = vel / eos_state%cs

             enddo
          enddo
       enddo

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! vel is the magnitude of the velocity, including w0
#if (AMREX_SPACEDIM == 2)
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + w0cart(i,j,k) )**2 )
#elif (AMREX_SPACEDIM == 3)
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + w0cart(i,j,k))**2 )
#endif
                ! (rho, h) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)

                pt_index(:) = (/i, j, k/)

                call eos(eos_input_rh, eos_state, pt_index)

                mach(i,j,k) = vel / eos_state%cs

             enddo
          enddo
       enddo

    endif

  end subroutine makeMachfromRhoH


  subroutine makeCsfromRhoH(lo,hi,state,s_lo,s_hi,p0cart,p_lo,p_hi,cs,c_lo,c_hi) &
       bind(C,name="makeCsfromRhoH")


       use fill_3d_data_module, only: put_1d_array_on_cart_sphr

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    double precision, intent (in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in   ) :: p0cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),1)
    double precision, intent (inout) :: cs(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! (rho, (h->e)) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0cart(i,j,k,1) / state(i,j,k,rho_comp)

                pt_index(:) = (/i, j, k/)

                call eos(eos_input_re, eos_state, pt_index)

                cs(i,j,k) = eos_state%cs

             enddo
          enddo
       enddo

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ! (rho, h) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)

                pt_index(:) = (/i, j, k/)

                call eos(eos_input_rh, eos_state, pt_index)

                cs(i,j,k) = eos_state%cs

             enddo
          enddo
       enddo

    endif

  end subroutine makeCsfromRhoH


  subroutine makeHfromRhoT_edge(lo,hi,idir,sedge,x_lo,x_hi, &
       rho0_cart,r_lo,r_hi, &
       rhoh0_cart,rh_lo,rh_hi, &
       temp0_cart,t_lo,t_hi, &
       rho0_edge_cart,re_lo,re_hi, &
       rhoh0_edge_cart,rhe_lo,rhe_hi, &
       temp0_edge_cart,te_lo,te_hi) &
       bind(C,name="makeHfromRhoT_edge")
    
    use meth_params_module, only: species_pred_type, enthalpy_pred_type, small_temp
    
    integer         , intent (in   ) :: lo(3), hi(3)
    integer,   value, intent (in   ) :: idir
    integer         , intent (in   ) :: x_lo(3), x_hi(3)
    double precision, intent (inout) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nscal)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    double precision, intent (in   ) :: rho0_cart (r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    integer         , intent (in   ) :: rh_lo(3), rh_hi(3)
    double precision, intent (in   ) :: rhoh0_cart (rh_lo(1):rh_hi(1),rh_lo(2):rh_hi(2),rh_lo(3):rh_hi(3))
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    double precision, intent (in   ) :: temp0_cart (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer         , intent (in   ) :: re_lo(3), re_hi(3)
    double precision, intent (in   ) :: rho0_edge_cart (re_lo(1):re_hi(1),re_lo(2):re_hi(2),re_lo(3):re_hi(3))
    integer         , intent (in   ) :: rhe_lo(3), rhe_hi(3)
    double precision, intent (in   ) :: rhoh0_edge_cart (rhe_lo(1):rhe_hi(1),rhe_lo(2):rhe_hi(2),rhe_lo(3):rhe_hi(3))
    integer         , intent (in   ) :: te_lo(3), te_hi(3)
    double precision, intent (in   ) :: temp0_edge_cart (te_lo(1):te_hi(1),te_lo(2):te_hi(2),te_lo(3):te_hi(3))

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (idir == 1) then
        ! x-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    
                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_cart(i,j,k),small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then
                        
                        ! interface states are rho' and X
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_cart(i,j,k)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp) 

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp) 

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h

                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_cart(i,j,k)
                    end if

                enddo
            enddo
        enddo

    else if (idir == 2) then

#if (AMREX_SPACEDIM == 2)
        ! y-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    
                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_edge_cart(i,j,k),small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then
                        
                        ! interface states are rho' and X
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_edge_cart(i,j,k)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1) 

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h

                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge_cart(i,j,k)
                    end if
                    
                enddo
            enddo
        enddo
    
#elif (AMREX_SPACEDIM == 3)
        ! y-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    
                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_cart(i,j,k),small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then
                        
                        ! interface states are rho' and X
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_cart(i,j,k)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1) 

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h

                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_cart(i,j,k)
                    end if
                    
                enddo
            enddo
        enddo

    else if (idir == 3) then 

        ! z-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                        
                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_edge_cart(i,j,k),small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then

                        ! interface states are rho' and X
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_edge_cart(i,j,k)
                    
                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1) 

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h
                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge_cart(i,j,k)
                    end if
                    
                enddo
            enddo
        enddo
#endif

    endif

  end subroutine makeHfromRhoT_edge


  subroutine makeHfromRhoT_edge_sphr(lo,hi,idir,sedge,x_lo,x_hi, &
       rho0_cart,r_lo,r_hi, &
       rhoh0_cart,rh_lo,rh_hi, &
       temp0_cart,t_lo,t_hi) &
       bind(C,name="makeHfromRhoT_edge_sphr")
    
    use meth_params_module, only: species_pred_type, enthalpy_pred_type, small_temp
    
    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: idir
    integer         , intent (in   ) :: x_lo(3), x_hi(3)
    double precision, intent (inout) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nscal)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    double precision, intent (in   ) :: rho0_cart (r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    integer         , intent (in   ) :: rh_lo(3), rh_hi(3)
    double precision, intent (in   ) :: rhoh0_cart (rh_lo(1):rh_hi(1),rh_lo(2):rh_hi(2),rh_lo(3):rh_hi(3))
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    double precision, intent (in   ) :: temp0_cart (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! Local variables
    integer :: i, j, k
    double precision :: rho0_edge, rhoh0_edge, temp0_edge
    
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (idir == 1) then

        ! x-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        temp0_edge = HALF* (temp0_cart(i-1,j,k) + temp0_cart(i,j,k))
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_edge,small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then
                        
                        ! interface states are rho' and X
                        rho0_edge = HALF * (rho0_cart(i-1,j,k) + rho0_cart(i,j,k))
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_edge

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp) 

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp) 

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h

                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        rhoh0_edge = HALF * (rhoh0_cart(i-1,j,k) + rhoh0_cart(i,j,k))
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge
                    end if

                enddo
            enddo
        enddo

    else if (idir == 2) then
        ! y-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    
                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        temp0_edge = HALF * (temp0_cart(i,j-1,k) + temp0_cart(i,j,k))
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_edge,small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then
                        
                        ! interface states are rho' and X
                        rho0_edge = HALF * (rho0_cart(i,j-1,k) + rho0_cart(i,j,k))
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_edge

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1) 

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h

                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        rhoh0_edge = HALF * (rhoh0_cart(i,j-1,k) + rhoh0_cart(i,j,k))
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge
                    end if
                    
                enddo
            enddo
        enddo

    else if (idir == 3) then

        ! z-edge
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                        
                    ! get edge-centered temperature
                    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        temp0_edge = HALF * (temp0_cart(i,j,k-1) + temp0_cart(i,j,k))
                        eos_state%T = max(sedge(i,j,k,temp_comp)+temp0_edge,small_temp)
                    else
                        eos_state%T = max(sedge(i,j,k,temp_comp),small_temp)
                    end if

                    ! get edge-centered density and species
                    if (species_pred_type .eq. predict_rhoprime_and_X) then

                        ! interface states are rho' and X
                        rho0_edge = HALF * (rho0_cart(i,j,k-1) + rho0_cart(i,j,k))
                        eos_state%rho = sedge(i,j,k,rho_comp) + rho0_edge
                    
                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1) 

                    else if (species_pred_type .eq. predict_rhoX) then

                        ! interface states are rho and (rho X)
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                    else if (species_pred_type .eq. predict_rho_and_X) then

                        ! interface states are rho and X
                        eos_state%rho = sedge(i,j,k,rho_comp)

                        eos_state%xn(:) = sedge(i,j,k,spec_comp:spec_comp+nspec-1)

                    endif

                    pt_index(:) = (/i, j, k/)

                    call eos(eos_input_rt, eos_state, pt_index)
                    
                    if (enthalpy_pred_type .eq. predict_T_then_h .or. &
                        enthalpy_pred_type .eq. predict_Tprime_then_h) then
                        sedge(i,j,k,rhoh_comp) = eos_state%h
                    else if (enthalpy_pred_type .eq. predict_T_then_rhohprime) then
                        rhoh0_edge = HALF * (rhoh0_cart(i,j,k-1) + rhoh0_cart(i,j,k))
                        sedge(i,j,k,rhoh_comp) = eos_state%rho*eos_state%h - rhoh0_edge
                    end if
                    
                enddo
            enddo
        enddo

    endif
    
  end subroutine makeHfromRhoT_edge_sphr


end module rhoh_vs_t_module
