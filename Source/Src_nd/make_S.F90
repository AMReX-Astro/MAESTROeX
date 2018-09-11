module make_S_module

  use amrex_error_module
  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, dpdt_factor, base_cutoff_density
  use base_state_geometry_module, only:  max_radial_level, nr_fine, base_cutoff_density_coord
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr
  use probin_module, only: use_delta_gamma1_term
  use geometry, only: anelastic_cutoff_coord, nr

  implicit none

  private

contains

  subroutine make_S_cc(lo, hi, &
       S_cc,  s_lo, s_hi, &
       delta_gamma1_term,  dg_lo, dg_hi, &
       delta_gamma1,  df_lo, df_hi, &
       scal,  c_lo, c_hi, nc_c, &
       rodot, r_lo, r_hi, nc_r, &
       rHnuc, n_lo, n_hi, &
       rHext, e_lo, e_hi, &
       therm, t_lo, t_hi, &
       p0, gamma1bar, dx) bind (C,name="make_S_cc")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    integer         , intent (in   ) :: df_lo(3), df_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3), nc_c
    integer         , intent (in   ) :: r_lo(3), r_hi(3), nc_r
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: e_lo(3), e_hi(3)
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    double precision, intent (inout) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (inout) :: delta_gamma1_term(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (inout) :: delta_gamma1(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: scal (c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nc_c)
    double precision, intent (in   ) :: rodot(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nc_r)
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
#if (AMREX_SPACEDIM == 1)
             r = i
#elif (AMREX_SPACEDIM == 2)
             r = j
#else
             r = k
#endif
             if (use_delta_gamma1_term .and. r < anelastic_cutoff_coord(n)) then
                if (r .eq. 0) then
                   gradp0 = (p0(lev,r+1) - p0(r))/dx(AMREX_SPACEDIM)
                else if (r .eq. nr(n)-1) then
                   gradp0 = (p0(lev,r) - p0(lev,r-1))/dx(AMREX_SPACEDIM)
                else
                   gradp0 = 0.5d0*(p0(lev,r+1) - p0(lev,r-1))/dx(AMREX_SPACEDIM)
                endif

                delta_gamma1(i,j,k) = eos_state%gam1 - gamma1bar(lev,r)

                delta_gamma1_term(i,j,k) = (eos_state%gam1 - gamma1bar(lev,r))*u(i,j,k,AMREX_SPACEDIM)* &
                     gradp0/(gamma1bar(lev,r)*gamma1bar(lev,r)*p0(lev,r))
             else
                delta_gamma1_term(i,j,k) = ZERO
                delta_gamma1(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine make_S_cc

  subroutine make_rhcc_for_nodalproj(lev, lo, hi, &
       rhcc, c_lo, c_hi, &
       S_cc,  s_lo, s_hi, &
       Sbar, beta0, &
       delta_gamma1_term, dg_lo, dg_hi) bind (C,name="make_rhcc_for_nodalproj")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent (inout) :: rhcc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: Sbar (0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: delta_gamma1_term (dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))

    integer i,j,k,r

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 1)
             r = i
#elif (AMREX_SPACEDIM == 2)
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

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhcc(i,j,k) = beta0_cart(i,j,k) * (S_cc(i,j,k) - Sbar_cart(i,j,k) + delta_gamma1_term(i,j,k))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_rhcc_for_nodalproj_sphr

  subroutine create_correction_cc(lev, lo, hi, &
       correction_cc, c_lo, c_hi, &
       delta_p_term, dp_lo, dp_hi, &
       beta0, gamma1bar, &
       p0, dt) &
       bind (C,name="create_correction_cc")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    double precision, intent(  out) :: correction_cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    integer         , intent(in   ) :: dp_lo(3), dp_hi(3)
    double precision, intent(in   ) :: delta_p_term(dp_lo(1):dp_hi(1),dp_lo(2):dp_hi(2),dp_lo(3):dp_hi(3))
    double precision, intent(in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dt

    ! Local variables
    integer :: i, j, k
    double precision :: correction_factor

#if (AMREX_SPACEDIM == 3)
    !$OMP PARALLEL DO PRIVATE(i,j,k,correction_factor)
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
    !$OMP END PARALLEL DO
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
#elif (AMREX_SPACEDIM == 1 )
    k = lo(3)
    j = lo(2)
    do i = lo(1),hi(1)
       if(i .lt. base_cutoff_density_coord(lev)) then
          correction_factor = beta0(lev,i)*(dpdt_factor/(gamma1bar(lev,i)*p0(lev,i))) / dt
       else
          correction_factor = 0.0d0
       end if
       correction_cc(i,j,k) = correction_factor*delta_p_term(i,j,k)
    end do
#endif

  end subroutine create_correction_cc

  subroutine create_correction_cc_sphr(lev, lo, hi, &
       correction_cc, c_lo, c_hi, &
       delta_p_term, dp_lo, dp_hi, &
       beta0_cart, b_lo, b_hi, &
       gamma1bar_cart, g_lo, g_hi, &
       p0_cart, p_lo, p_hi, &
       rho0_cart, r_lo, r_hi, &
       dt) bind (C,name="create_correction_cc_sphr")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
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
    double precision, intent(in   ) :: dt

    ! Local variables
    integer :: i, j, k
    double precision :: correction_factor

    !$OMP PARALLEL DO PRIVATE(i,j,k,correction_factor)
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
    !$OMP END PARALLEL DO

  end subroutine create_correction_cc_sphr

  subroutine make_rhcc_for_macproj(lev, lo, hi, &
       rhcc, c_lo, c_hi, &
       S_cc,  s_lo, s_hi, &
       Sbar, beta0, &
       delta_gamma1_term, dg_lo, dg_hi, &
       gamma1bar, p0, &
       delta_p_term, dp_lo, dp_hi, &
       delta_chi, dc_lo, dc_hi, &
       dt, is_predictor) bind (C,name="make_rhcc_for_macproj")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
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
    double precision, intent (in   ) :: dt
    integer         , intent (in   ) :: is_predictor

    ! Local variables
    integer i,j,k,r

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 1)
             r = i
#elif (AMREX_SPACEDIM == 2)
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
            delta_chi = 0.d0

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

#elif (AMREX_SPACEDIM == 1)
       k = lo(3)
       j = lo(2)
       do i = lo(1),hi(1)
          if (i .lt. base_cutoff_density_coord(n)) then
             delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                  (dt*gamma1bar(lev,i)*p0(lev,i))
             rhcc(i,j,k) = rhcc(i,j,k) + beta0(lev,i) * delta_chi(i,j,k)
          end if
       end do
#endif

    end if

  end subroutine make_rhcc_for_macproj

  subroutine make_rhcc_for_macproj_sphr(lo, hi, &
       rhcc, c_lo, c_hi, &
       S_cc,  s_lo, s_hi, &
       Sbar, beta0, &
       rho0, dx, &
       delta_gamma1_term, dg_lo, dg_hi, &
       gamma1bar, p0, &
       delta_p_term, dp_lo, dp_hi, &
       delta_chi, dc_lo, dc_hi, &
       dt, is_predictor, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) &
       bind (C,name="make_rhcc_for_macproj_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent (inout) :: rhcc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) :: S_cc(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) ::  Sbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::  rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: dx(3)
    double precision, intent (in   ) :: delta_gamma1_term (dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    double precision, intent (in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: dp_lo(3), dp_hi(3)
    double precision, intent (in   ) :: delta_p_term(dp_lo(1):dp_hi(1),dp_lo(2):dp_hi(2),dp_lo(3):dp_hi(3))
    integer         , intent (in   ) :: dc_lo(3), dc_hi(3)
    double precision, intent (inout) :: delta_chi(dc_lo(1):dc_hi(1),dc_lo(2):dc_hi(2),dc_lo(3):dc_hi(3))
    double precision, intent (in   ) :: dt
    integer         , intent (in   ) :: is_predictor
    double precision, intent (in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent (in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent (in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    !     Local variables
    integer :: i, j, k
    double precision, allocatable ::       div_cart(:,:,:,:)
    double precision, allocatable ::      Sbar_cart(:,:,:,:)
    double precision, allocatable :: gamma1bar_cart(:,:,:,:)
    double precision, allocatable ::        p0_cart(:,:,:,:)
    double precision, allocatable ::      rho0_cart(:,:,:,:)

    allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_sphr(lo,hi,div_cart,lo,hi,1,beta0,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)

    allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_sphr(lo,hi,Sbar_cart,lo,hi,1,Sbar,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r, ccr_lo, ccr_hi)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhcc(i,j,k) = div_cart(i,j,k,1) * (S_cc(i,j,k) - Sbar_cart(i,j,k,1) + delta_gamma1_term(i,j,k))
          end do
       end do
    end do

    deallocate(Sbar_cart)

    if (dpdt_factor .gt. 0.0d0) then

       allocate(gamma1bar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_sphr(lo,hi,gamma1bar_cart,lo,hi,1,gamma1bar,dx,0,0, &
            r_cc_loc,r_edge_loc, cc_to_r,ccr_lo,ccr_hi)

       allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_sphr(lo,hi,p0_cart,lo,hi,1,p0,dx,0,0,r_cc_loc,r_edge_loc, &
            cc_to_r,ccr_lo,ccr_hi)

       allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
       call put_1d_array_on_cart_sphr(lo,hi,rho0_cart,lo,hi,1,rho0,dx,0,0,r_cc_loc,r_edge_loc, &
            cc_to_r,ccr_lo,ccr_hi)

       if (is_predictor .eq. 1) &
            delta_chi = 0.d0

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (rho0_cart(i,j,k,1) .gt. base_cutoff_density) then
                   delta_chi(i,j,k) = delta_chi(i,j,k) + dpdt_factor * delta_p_term(i,j,k) / &
                        (dt*gamma1bar_cart(i,j,k,1)*p0_cart(i,j,k,1))
                   rhcc(i,j,k) = rhcc(i,j,k) + div_cart(i,j,k,1) * delta_chi(i,j,k)
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       deallocate(gamma1bar_cart,p0_cart,rho0_cart)

    end if

    deallocate(div_cart)

  end subroutine make_rhcc_for_macproj_sphr

  subroutine create_correction_delta_gamma1_term(lev, lo, hi, &
       delta_gamma1_term, dg_lo, dg_hi, &
       delta_gamma1, df_lo, df_hi, &
       gamma1bar, psi, &
       p0) bind (C,name="create_correction_delta_gamma1_term")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent(in   ) :: delta_gamma1_term(dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    integer         , intent(in   ) :: df_lo(3), df_hi(3)
    double precision, intent(in   ) :: delta_gamma1(df_lo(1):df_hi(1),df_lo(2):df_hi(2),dg_lo(3):df_hi(3))
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: psi(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer :: i, j, k

#if (AMREX_SPACEDIM == 3)
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
#else
       k = lo(3)
#endif
#if (AMREX_SPACEDIM >= 2)
       do j = lo(2),hi(2)
#else
          j = lo(2)
#endif
          do i = lo(1),hi(1)
             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) &
                  + delta_gamma1(i,j,k)*psi(lev,k)/(gamma1bar(lev,k)**2*p0(lev,k))
          end do
#if (AMREX_SPACEDIM >=2)
       end do
#if (AMREX_SPACEDIM == 3)
    end do
#endif
#endif
    !$OMP END PARALLEL DO

  end subroutine create_correction_delta_gamma1_term

  subroutine create_correction_delta_gamma1_term_sphr(lev, lo, hi, &
       delta_gamma1_term, dg_lo, dg_hi, &
       delta_gamma1, df_lo, df_hi, &
       gamma1bar, psi, &
       p0, dx, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind (C,name="create_correction_delta_gamma1_term_sphr")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: dg_lo(3), dg_hi(3)
    double precision, intent(in   ) :: delta_gamma1_term(dg_lo(1):dg_hi(1),dg_lo(2):dg_hi(2),dg_lo(3):dg_hi(3))
    integer         , intent(in   ) :: df_lo(3), df_hi(3)
    double precision, intent(in   ) :: delta_gamma1(df_lo(1):df_hi(1),df_lo(2):df_hi(2),dg_lo(3):df_hi(3))
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: psi(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    ! Local variables
    integer :: i, j, k
    double precision, allocatable :: gamma1bar_cart(:,:,:,:)
    double precision, allocatable ::        p0_cart(:,:,:,:)
    double precision, allocatable ::      psi_cart(:,:,:,:)

    allocate(gamma1bar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_sphr(lo,hi,gamma1bar_cart,lo,hi,1,gamma1bar,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,p0_cart,lo,hi,1,p0,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,psi_cart,lo,hi,1,psi,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)


#if (AMREX_SPACEDIM == 3)
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
#else
       k = lo(3)
#endif
#if (AMREX_SPACEDIM >= 2)
       do j = lo(2),hi(2)
#else
          j = lo(2)
#endif
          do i = lo(1),hi(1)
             delta_gamma1_term(i,j,k) = delta_gamma1_term(i,j,k) &
                  + delta_gamma1(i,j,k)*psi_cart(i,j,k,1)/(gamma1bar_cart(i,j,k,1)**2*p0_cart(i,j,k,1))
          end do
#if (AMREX_SPACEDIM >=2)
       end do
#if (AMREX_SPACEDIM == 3)
    end do
#endif
#endif
    !$OMP END PARALLEL DO

    deallocate(gamma1bar_cart, p0_cart, psi_cart)

  end subroutine create_correction_delta_gamma1_term_sphr

end module make_S_module
