module rhoh_vs_t_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       nscal, use_eos_e_instead_of_h, use_pprime_in_tfromp
  use base_state_geometry_module, only:  max_radial_level, nr_fine

  implicit none

contains

  subroutine makeTfromRhoH(lo,hi,lev,state,s_lo,s_hi,p0) bind(C,name="makeTfromRhoH")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
                r = j
#elif (AMREX_SPACEDIM == 3)
                r = k
#endif

                ! (rho, (h->e)) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0(lev,r) / state(i,j,k,rho_comp)

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

#if (AMREX_SPACEDIM == 2)
                r = j
#elif (AMREX_SPACEDIM == 3)
                r = k
#endif

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

  subroutine makeTfromRhoH_sphr(lo,hi,state,s_lo,s_hi,p0_cart,p0_lo,p0_hi) &
       bind(C,name="makeTfromRhoH_sphr")

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

  end subroutine makeTfromRhoH_sphr

  subroutine makeTfromRhoP(lo,hi,lev,state,s_lo,s_hi,p0,updateRhoH) &
       bind(C,name="makeTfromRhoP")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    integer  , value, intent (in   ) :: updateRhoH

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
             r = j
#elif (AMREX_SPACEDIM == 3)
             r = k
#endif

             ! (rho, p) --> T

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p     = p0(lev,r) + state(i,j,k,pi_comp)
             else
                eos_state%p     = p0(lev,r)
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

  subroutine makeTfromRhoP_sphr(lo,hi,state,s_lo,s_hi,p0_cart,p0_lo,p0_hi,updateRhoH) &
       bind(C,name="makeTfromRhoP_sphr")

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

  end subroutine makeTfromRhoP_sphr

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

  subroutine makeMachfromRhoH(lo,hi,lev,state,s_lo,s_hi,u,u_lo,u_hi, &
       p0,w0,mach,m_lo,m_hi) bind(C,name="makeMachfromRhoH")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    double precision, intent (in   ) ::  u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    double precision, intent (inout) :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

    ! Local variables
    integer :: i, j, k, r
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
                r = j
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + 0.5d0*(w0(lev,r) + w0(lev,r+1)) )**2 )
#elif (AMREX_SPACEDIM == 3)
                r = k
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + 0.5d0*(w0(lev,r) + w0(lev,r+1)) )**2 )
#endif

                ! (rho, (h->e)) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0(lev,r) / state(i,j,k,rho_comp)

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
                r = j
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + 0.5d0*(w0(lev,r) + w0(lev,r+1)) )**2 )
#elif (AMREX_SPACEDIM == 3)
                r = k
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + 0.5d0*(w0(lev,r) + w0(lev,r+1)) )**2 )
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

  subroutine makeMachfromRhoH_sphr(lo,hi,state,s_lo,s_hi,u,u_lo,u_hi, &
       p0_cart,p0_lo,p0_hi,w0cart,w_lo,w_hi,mach,m_lo,m_hi) bind(C,name="makeMachfromRhoH_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    double precision, intent (in   ) ::  u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    integer         , intent (in   ) :: p0_lo(3), p0_hi(3)
    double precision, intent (in) :: p0_cart (p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
    double precision, intent (in   ) :: w0cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),1)
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    double precision, intent (inout) :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    ! Local variables
    integer :: i, j, k, r
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
                r = j
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + 0.5d0*(w0cart(i,j,k,1) + w0cart(i,j,k,1)) )**2 )
#elif (AMREX_SPACEDIM == 3)
                r = k
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + 0.5d0*(w0cart(i,j,k,1) + w0cart(i,j,k,1)) )**2 )
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
                r = j
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + w0cart(i,j,k,1) )**2 )
#elif (AMREX_SPACEDIM == 3)
                r = k
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + w0cart(i,j,k,1) )**2 )
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

  end subroutine makeMachfromRhoH_sphr

  subroutine makeCsfromRhoH(lo,hi,lev,state,s_lo,s_hi,p0,cs,c_lo,c_hi) &
       bind(C,name="makeCsfromRhoH")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    double precision, intent (inout) :: cs(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
                r = j
#elif (AMREX_SPACEDIM == 3)
                r = k
#endif

                ! (rho, (h->e)) --> T, p

                eos_state%rho   = state(i,j,k,rho_comp)
                eos_state%T     = state(i,j,k,temp_comp)
                eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                ! e = h - p/rho
                eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
                     p0(lev,r) / state(i,j,k,rho_comp)

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

#if (AMREX_SPACEDIM == 2)
                r = j
#elif (AMREX_SPACEDIM == 3)
                r = k
#endif
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

  subroutine makeCsfromRhoH_sphr(lo,hi,state,s_lo,s_hi,p0cart,p_lo,p_hi,cs,c_lo,c_hi) &
       bind(C,name="makeCsfromRhoH_sphr")


       use fill_3d_data_module, only: put_1d_array_on_cart_sphr

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    double precision, intent (in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal)
    double precision, intent (in   ) :: p0cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),1)
    double precision, intent (inout) :: cs(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
                r = j
#elif (AMREX_SPACEDIM == 3)
                r = k
#endif

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

#if (AMREX_SPACEDIM == 2)
                r = j
#elif (AMREX_SPACEDIM == 3)
                r = k
#endif
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

  end subroutine makeCsfromRhoH_sphr


end module rhoh_vs_t_module
