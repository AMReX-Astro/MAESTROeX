module make_S_module

  use amrex_error_module
  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, nscal, dpdt_factor, base_cutoff_density, use_delta_gamma1_term
  use base_state_geometry_module, only:  max_radial_level, nr_fine, base_cutoff_density_coord, anelastic_cutoff_density_coord, nr, dr

  implicit none

  private

contains

  subroutine make_S_cc(lo, hi, lev, &
       S_cc,  s_lo, s_hi, &
       delta_gamma1_term,  dg_lo, dg_hi, &
       delta_gamma1,  df_lo, df_hi, &
       scal,  c_lo, c_hi, &
       u,  u_lo, u_hi, &
       rodot, r_lo, r_hi, &
       rHnuc, n_lo, n_hi, &
       rHext, e_lo, e_hi, &
       therm, t_lo, t_hi, &
       p0, gamma1bar, dx) bind (C,name="make_S_cc")

    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    integer         , intent (in   ) :: df_lo(3), df_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: e_lo(3), e_hi(3)
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    double precision, intent (inout) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (inout) :: delta_gamma1_term(dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    double precision, intent (inout) :: delta_gamma1(df_lo(1):df_hi(1),df_lo(2):df_hi(2),df_lo(3):df_hi(3))
    double precision, intent (in   ) :: scal (c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nscal)
    double precision, intent (in   ) :: u (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),AMREX_SPACEDIM)
    double precision, intent (in   ) :: rodot(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec)
    double precision, intent (in   ) :: rHnuc(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3))
    double precision, intent (in   ) :: rHext(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    double precision, intent (in   ) :: therm(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer i,j,k,r
    integer pt_index(3)
    type(eos_t) :: eos_state

    integer comp
    double precision sigma, xi_term, pres_term, gradp0

    !$gpu

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state%rho   = scal(i,j,k,rho_comp)
             eos_state%T     = scal(i,j,k,temp_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, pt_index)

             sigma = eos_state%dpdt / &
                  (eos_state%rho * eos_state%cp * eos_state%dpdr)

             xi_term = 0.d0
             pres_term = 0.d0
             do comp = 1, nspec
                xi_term = xi_term - &
                     eos_state%dhdX(comp)*rodot(i,j,k,comp)/eos_state%rho

                pres_term = pres_term + &
                     eos_state%dpdX(comp)*rodot(i,j,k,comp)/eos_state%rho
             enddo

             S_cc(i,j,k) = (sigma/eos_state%rho) * &
                  ( rHext(i,j,k) + rHnuc(i,j,k) + therm(i,j,k) ) &
                  + sigma*xi_term &
                  + pres_term/(eos_state%rho*eos_state%dpdr)
#if (AMREX_SPACEDIM == 2)
             r = j
#else
             r = k
#endif
             if (use_delta_gamma1_term .and. r < anelastic_cutoff_density_coord(lev)) then
                if (r .eq. 0) then
                   gradp0 = (p0(lev,r+1) - p0(lev,r))/dx(AMREX_SPACEDIM)
                else if (r .eq. nr(lev)-1) then
                   gradp0 = (p0(lev,r) - p0(lev,r-1))/dx(AMREX_SPACEDIM)
                else
                   gradp0 = 0.5d0*(p0(lev,r+1) - p0(lev,r-1))/dx(AMREX_SPACEDIM)
                endif

                delta_gamma1(i,j,k) = eos_state%gam1 - gamma1bar(lev,r)

                delta_gamma1_term(i,j,k) = (eos_state%gam1 - gamma1bar(lev,r))*u(i,j,k,AMREX_SPACEDIM)* &
                     gradp0/(gamma1bar(lev,r)*gamma1bar(lev,r)*p0(lev,r))
             else
                delta_gamma1_term(i,j,k) = 0.0d0
                delta_gamma1(i,j,k) = 0.0d0
             endif

          enddo
       enddo
    enddo

  end subroutine make_S_cc

  subroutine make_S_cc_sphr(lo, hi, &
       S_cc,  s_lo, s_hi, &
       delta_gamma1_term,  dg_lo, dg_hi, &
       delta_gamma1,  df_lo, df_hi, &
       scal,  c_lo, c_hi, &
       u,  u_lo, u_hi, &
       rodot, r_lo, r_hi, &
       rHnuc, n_lo, n_hi, &
       rHext, e_lo, e_hi, &
       therm, t_lo, t_hi, &
       p0_cart, p0_lo, p0_hi, &
       gradp0_cart, gp0_lo, gp0_hi, &
       gamma1bar_cart, g1_lo, g1_hi, &
       normal, no_lo, no_hi) bind (C,name="make_S_cc_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    integer         , intent (in   ) :: df_lo(3), df_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: e_lo(3), e_hi(3)
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    integer         , intent (in   ) :: no_lo(3), no_hi(3)
    integer         , intent (in   ) :: p0_lo(3), p0_hi(3)
    integer         , intent (in   ) :: gp0_lo(3), gp0_hi(3)
    integer         , intent (in   ) :: g1_lo(3), g1_hi(3)
    double precision, intent (inout) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (inout) :: delta_gamma1_term(dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    double precision, intent (inout) :: delta_gamma1(df_lo(1):df_hi(1),df_lo(2):df_hi(2),df_lo(3):df_hi(3))
    double precision, intent (in   ) :: scal (c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nscal)
    double precision, intent (in   ) :: u (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),AMREX_SPACEDIM)
    double precision, intent (in   ) :: rodot(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec)
    double precision, intent (in   ) :: rHnuc(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3))
    double precision, intent (in   ) :: rHext(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    double precision, intent (in   ) :: therm(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent (in) :: p0_cart (p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))
    double precision, intent (in) :: gradp0_cart (gp0_lo(1):gp0_hi(1),gp0_lo(2):gp0_hi(2),gp0_lo(3):gp0_hi(3))
    double precision, intent (in) :: gamma1bar_cart (g1_lo(1):g1_hi(1),g1_lo(2):g1_hi(2),g1_lo(3):g1_hi(3))
    double precision, intent (in   ) :: normal(no_lo(1):no_hi(1),no_lo(2):no_hi(2),no_lo(3):no_hi(3),3)

    integer i,j,k,r
    integer pt_index(3)
    type(eos_t) :: eos_state

    integer comp
    double precision sigma, xi_term, pres_term, Ut_dot_er

    !$gpu

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state%rho   = scal(i,j,k,rho_comp)
             eos_state%T     = scal(i,j,k,temp_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, pt_index)

             sigma = eos_state%dpdt / &
                  (eos_state%rho * eos_state%cp * eos_state%dpdr)

             xi_term = 0.d0
             pres_term = 0.d0
             do comp = 1, nspec
                xi_term = xi_term - &
                     eos_state%dhdX(comp)*rodot(i,j,k,comp)/eos_state%rho

                pres_term = pres_term + &
                     eos_state%dpdX(comp)*rodot(i,j,k,comp)/eos_state%rho
             enddo

             S_cc(i,j,k) = (sigma/eos_state%rho) * &
                  ( rHext(i,j,k) + rHnuc(i,j,k) + therm(i,j,k) ) &
                  + sigma*xi_term &
                  + pres_term/(eos_state%rho*eos_state%dpdr)

             if (use_delta_gamma1_term) then
                delta_gamma1(i,j,k) = eos_state%gam1 - gamma1bar_cart(i,j,k)

                Ut_dot_er = &
                     u(i,j,k,1)*normal(i,j,k,1) + &
                     u(i,j,k,2)*normal(i,j,k,2) + &
                     u(i,j,k,3)*normal(i,j,k,3)

                delta_gamma1_term(i,j,k) = delta_gamma1(i,j,k)*Ut_dot_er* &
                     gradp0_cart(i,j,k)/ &
                     (gamma1bar_cart(i,j,k)**2*p0_cart(i,j,k))

             else
                delta_gamma1_term(i,j,k) = 0.0d0
                delta_gamma1(i,j,k) = 0.0d0
             end if

          enddo
       enddo
    enddo

  end subroutine make_S_cc_sphr

  subroutine make_rhcc_for_nodalproj(lo, hi, lev, &
       rhcc, c_lo, c_hi, &
       S_cc,  s_lo, s_hi, &
       Sbar, beta0, &
       delta_gamma1_term, dg_lo, dg_hi) bind (C,name="make_rhcc_for_nodalproj")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent (inout) :: rhcc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: Sbar (0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: delta_gamma1_term (dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))

    integer i,j,k,r

    !$gpu

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 2)
             r = j
#elif (AMREX_SPACEDIM == 3)
             r = k
#endif
             rhcc(i,j,k) = beta0(lev,r) * (S_cc(i,j,k) - Sbar(lev,r) + delta_gamma1_term(i,j,k))

          enddo
       enddo
    enddo

  end subroutine make_rhcc_for_nodalproj

  subroutine make_rhcc_for_nodalproj_sphr(lo, hi, &
       rhcc, c_lo, c_hi, &
       S_cc, s_lo, s_hi, &
       Sbar_cart, sb_lo, sb_hi, &
       beta0_cart, b_lo, b_hi, &
       delta_gamma1_term, dg_lo, dg_hi) &
       bind (C,name="make_rhcc_for_nodalproj_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: sb_lo(3), sb_hi(3)
    integer         , intent (in   ) :: b_lo(3), b_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent (inout) ::  rhcc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) ::  S_cc(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: Sbar_cart(sb_lo(1):sb_hi(1),sb_lo(2):sb_hi(2),sb_lo(3):sb_hi(3))
    double precision, intent (in   ) :: beta0_cart(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    double precision, intent (in   ) :: delta_gamma1_term (dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))

    ! Local variables
    integer :: i,j,k

    !$gpu

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhcc(i,j,k) = beta0_cart(i,j,k) * (S_cc(i,j,k) - Sbar_cart(i,j,k) + delta_gamma1_term(i,j,k))
          end do
       end do
    end do

  end subroutine make_rhcc_for_nodalproj_sphr

  subroutine create_correction_cc(lo, hi, lev, &
       correction_cc, c_lo, c_hi, &
       delta_p_term, dp_lo, dp_hi, &
       beta0, gamma1bar, &
       p0, dt) &
       bind (C,name="create_correction_cc")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    double precision, intent(  out) :: correction_cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    integer         , intent(in   ) :: dp_lo(3), dp_hi(3)
    double precision, intent(in   ) :: delta_p_term(dp_lo(1):dp_hi(1),dp_lo(2):dp_hi(2),dp_lo(3):dp_hi(3))
    double precision, intent(in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    double precision, value, intent(in   ) :: dt

    ! Local variables
    integer :: i, j, k, r
    double precision :: correction_factor

    !$gpu

#if (AMREX_SPACEDIM == 3)
    do k = lo(3),hi(3)
       if(k .lt. base_cutoff_density_coord(lev)) then
          correction_factor = beta0(lev,k)*(dpdt_factor/(gamma1bar(lev,k)*p0(lev,k))) / dt
       else
          correction_factor = 0.0d0
       end if
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             correction_cc(i,j,k) = correction_factor*delta_p_term(i,j,k)
          end do
       end do
    end do
#elif (AMREX_SPACEDIM == 2)
    k = lo(3)
    do j = lo(2),hi(2)
       if(j .lt. base_cutoff_density_coord(lev)) then
          correction_factor = beta0(lev,j)*(dpdt_factor/(gamma1bar(lev,j)*p0(lev,j))) / dt
       else
          correction_factor = 0.0d0
       end if
       do i = lo(1),hi(1)
          correction_cc(i,j,k) = correction_factor*delta_p_term(i,j,k)
       end do
    end do
#endif

  end subroutine create_correction_cc

  subroutine create_correction_cc_sphr(lo, hi, &
       correction_cc, c_lo, c_hi, &
       delta_p_term, dp_lo, dp_hi, &
       beta0_cart, b_lo, b_hi, &
       gamma1bar_cart, g_lo, g_hi, &
       p0_cart, p_lo, p_hi, &
       rho0_cart, r_lo, r_hi, &
       dt) bind (C,name="create_correction_cc_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    double precision, intent(  out) :: correction_cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    integer         , intent(in   ) :: dp_lo(3), dp_hi(3)
    double precision, intent(in   ) :: delta_p_term(dp_lo(1):dp_hi(1),dp_lo(2):dp_hi(2),dp_lo(3):dp_hi(3))
    integer         , intent(in   ) :: b_lo(3), b_hi(3)
    double precision, intent(in   ) :: beta0_cart(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    integer         , intent(in   ) :: g_lo(3), g_hi(3)
    double precision, intent(in   ) :: gamma1bar_cart(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3))
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(in   ) :: p0_cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    integer         , intent(in   ) :: r_lo(3), r_hi(3)
    double precision, intent(in   ) :: rho0_cart(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, value, intent(in   ) :: dt

    ! Local variables
    integer :: i, j, k
    double precision :: correction_factor

    !$gpu

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if(rho0_cart(i,j,k) .gt. base_cutoff_density) then
                correction_factor = beta0_cart(i,j,k) * &
                     (dpdt_factor/(gamma1bar_cart(i,j,k)*p0_cart(i,j,k))) / dt
             else
                correction_factor = 0.0d0
             end if
             correction_cc(i,j,k) = correction_factor*delta_p_term(i,j,k)
          end do
       end do
    end do

  end subroutine create_correction_cc_sphr

  subroutine make_rhcc_for_macproj(lo, hi, lev, &
       rhcc, c_lo, c_hi, &
       S_cc,  s_lo, s_hi, &
       Sbar, beta0, &
       delta_gamma1_term, dg_lo, dg_hi, &
       gamma1bar, p0, &
       delta_p_term, dp_lo, dp_hi, &
       delta_chi, dc_lo, dc_hi, &
       dt, is_predictor) bind (C,name="make_rhcc_for_macproj")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent (inout) :: rhcc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: Sbar (0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: delta_gamma1_term (dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    double precision, intent (in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: dp_lo(3), dp_hi(3)
    double precision, intent (in   ) :: delta_p_term(dp_lo(1):dp_hi(1),dp_lo(2):dp_hi(2),dp_lo(3):dp_hi(3))
    integer         , intent (in   ) :: dc_lo(3), dc_hi(3)
    double precision, intent (inout) :: delta_chi(dc_lo(1):dc_hi(1),dc_lo(2):dc_hi(2),dc_lo(3):dc_hi(3))
    double precision, value, intent (in   ) :: dt
    integer         , value, intent (in   ) :: is_predictor

    ! Local variables
    integer i,j,k,r

    !$gpu

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 2)
             r = j
#elif (AMREX_SPACEDIM == 3)
             r = k
#endif
             rhcc(i,j,k) = beta0(lev,r) * (S_cc(i,j,k) - Sbar(lev,r) + delta_gamma1_term(i,j,k))

          enddo
       enddo
    enddo

    if (dpdt_factor .gt. 0.0d0) then

       if (is_predictor .eq. 1) &
            delta_chi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0

#if (AMREX_SPACEDIM == 3)

       do k = lo(3),hi(3)
          if (k .lt. base_cutoff_density_coord(lev)) then
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                        (dt*gamma1bar(lev,k)*p0(lev,k))
                   rhcc(i,j,k) = rhcc(i,j,k) + beta0(lev,k) * delta_chi(i,j,k)
                end do
             end do
          end if
       end do

#elif (AMREX_SPACEDIM == 2)
       k = lo(3)
       do j = lo(2),hi(2)
          if (j .lt. base_cutoff_density_coord(lev)) then
             do i = lo(1),hi(1)
                delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                     (dt*gamma1bar(lev,j)*p0(lev,j))
                rhcc(i,j,k) = rhcc(i,j,k) + beta0(lev,j) * delta_chi(i,j,k)
             end do
          end if
       end do
#endif

    end if

  end subroutine make_rhcc_for_macproj

  subroutine make_rhcc_for_macproj_sphr(lo, hi, &
       rhcc, c_lo, c_hi, &
       S_cc,  s_lo, s_hi, &
       Sbar_cart, sb_lo, sb_hi, &
       beta0_cart, b_lo, b_hi, &
       rho0_cart, r_lo, r_hi, &
       delta_gamma1_term, dg_lo, dg_hi, &
       gamma1bar_cart, g1_lo, g1_hi, &
       p0_cart, p0_lo, p0_hi, &
       delta_p_term, dp_lo, dp_hi, &
       delta_chi, dc_lo, dc_hi, &
       dt, is_predictor) &
       bind (C,name="make_rhcc_for_macproj_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: sb_lo(3), sb_hi(3)
    integer         , intent (in   ) :: b_lo(3), b_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: g1_lo(3), g1_hi(3)
    integer         , intent (in   ) :: p0_lo(3), p0_hi(3)
    double precision, intent (inout) :: rhcc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) :: S_cc(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: Sbar_cart(sb_lo(1):sb_hi(1),sb_lo(2):sb_hi(2),sb_lo(3):sb_hi(3))
    double precision, intent (in   ) :: beta0_cart(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    double precision, intent (in   ) :: rho0_cart(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) :: delta_gamma1_term (dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    double precision, intent (in   ) :: gamma1bar_cart(g1_lo(1):g1_hi(1),g1_lo(2):g1_hi(2),g1_lo(3):g1_hi(3))
    double precision, intent (in   ) :: p0_cart(p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))
    integer         , intent (in   ) :: dp_lo(3), dp_hi(3)
    double precision, intent (in   ) :: delta_p_term(dp_lo(1):dp_hi(1),dp_lo(2):dp_hi(2),dp_lo(3):dp_hi(3))
    integer         , intent (in   ) :: dc_lo(3), dc_hi(3)
    double precision, intent (inout) :: delta_chi(dc_lo(1):dc_hi(1),dc_lo(2):dc_hi(2),dc_lo(3):dc_hi(3))
    double precision, value, intent (in   ) :: dt
    integer         , value, intent (in   ) :: is_predictor

    !     Local variables
    integer :: i, j, k

    !$gpu

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhcc(i,j,k) = beta0_cart(i,j,k) * (S_cc(i,j,k) - Sbar_cart(i,j,k) + delta_gamma1_term(i,j,k))
          end do
       end do
    end do

    if (dpdt_factor .gt. 0.0d0) then

       if (is_predictor .eq. 1) &
            delta_chi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (rho0_cart(i,j,k) .gt. base_cutoff_density) then
                   delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                        (dt*gamma1bar_cart(i,j,k)*p0_cart(i,j,k))
                   rhcc(i,j,k) = rhcc(i,j,k) + beta0_cart(i,j,k) * delta_chi(i,j,k)
                end if
             end do
          end do
       end do

    end if

  end subroutine make_rhcc_for_macproj_sphr

  subroutine create_correction_delta_gamma1_term(lo, hi, lev, &
       delta_gamma1_term, dg_lo, dg_hi, &
       delta_gamma1, df_lo, df_hi, &
       gamma1bar, psi, &
       p0) bind (C,name="create_correction_delta_gamma1_term")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent(inout) :: delta_gamma1_term(dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    integer         , intent(in   ) :: df_lo(3), df_hi(3)
    double precision, intent(in   ) :: delta_gamma1(df_lo(1):df_hi(1),df_lo(2):df_hi(2),dg_lo(3):df_hi(3))
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: psi(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer :: i, j, k, r

    !$gpu

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
#if (AMREX_SPACEDIM == 2)
             r = j
#else
             r = k
#endif
             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) &
                  + delta_gamma1(i,j,k)*psi(lev,r)/(gamma1bar(lev,r)**2*p0(lev,r))
          end do
       end do
    end do

  end subroutine create_correction_delta_gamma1_term

  subroutine create_correction_delta_gamma1_term_sphr(lo, hi, &
       delta_gamma1_term, dg_lo, dg_hi, &
       delta_gamma1, df_lo, df_hi, &
       gamma1bar_cart, g1_lo, g1_hi, &
       psi_cart, ps_lo, ps_hi, &
       p0_cart, p0_lo, p0_hi) bind (C,name="create_correction_delta_gamma1_term_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: p0_lo(3), p0_hi(3)
    integer         , intent(in   ) :: ps_lo(3), ps_hi(3)
    integer         , intent(in   ) :: g1_lo(3), g1_hi(3)
    integer         , intent(in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent(inout) :: delta_gamma1_term(dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    integer         , intent(in   ) :: df_lo(3), df_hi(3)
    double precision, intent(in   ) :: delta_gamma1(df_lo(1):df_hi(1),df_lo(2):df_hi(2),dg_lo(3):df_hi(3))
    double precision, intent(in   ) :: gamma1bar_cart (g1_lo(1):g1_hi(1),g1_lo(2):g1_hi(2),g1_lo(3):g1_hi(3))
    double precision, intent(in   ) :: psi_cart (ps_lo(1):ps_hi(1),ps_lo(2):ps_hi(2),ps_lo(3):ps_hi(3))
    double precision, intent(in   ) :: p0_cart (p0_lo(1):p0_hi(1),p0_lo(2):p0_hi(2),p0_lo(3):p0_hi(3))

    ! Local variables
    integer :: i, j, k

    !$gpu
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) &
                  + delta_gamma1(i,j,k)*psi_cart(i,j,k)/(gamma1bar_cart(i,j,k)**2*p0_cart(i,j,k))
          end do
       end do
    end do

  end subroutine create_correction_delta_gamma1_term_sphr

end module make_S_module
