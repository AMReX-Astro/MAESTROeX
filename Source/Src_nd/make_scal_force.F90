module make_scal_force_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr, base_cutoff_density_coord
  use meth_params_module, only: enthalpy_pred_type, use_exact_base_state, base_cutoff_density, spherical

  implicit none

  private

contains

  subroutine mkrhohforce(lo, hi, lev, &
       rhoh_force, f_lo, f_hi, &
       umac, u_lo, u_hi, &
       vmac, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wmac, w_lo, w_hi, &
#endif
       thermal, t_lo, t_hi, &
       grav_cart, g_lo, g_hi, &
       rho0_cart, r_lo, r_hi, &
       p0_cart, p_lo, p_hi, &
       p0macx, x_lo, x_hi, &
       p0macy, y_lo, y_hi, &
#if (AMREX_SPACEDIM == 3)
       p0macz, z_lo, z_hi, &
#endif
       psi_cart, ps_lo, ps_hi, &
       dx, domhi, &
       is_prediction, add_thermal) &
       bind(C,name="mkrhohforce")

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer         , intent(in   ) :: lo(3),hi(3)
    integer, value  , intent(in   ) :: lev
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: rhoh_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
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
    integer         , intent(in   ) :: g_lo(3), g_hi(3)
    double precision, intent(in   ) :: grav_cart(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3))
    integer         , intent(in   ) :: r_lo(3), r_hi(3)
    double precision, intent(in   ) :: rho0_cart(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(in   ) :: p0_cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) ::  p0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) ::  p0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) ::  p0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
#endif
    integer         , intent(in   ) :: ps_lo(3), ps_hi(3)
    double precision, intent(in   ) :: psi_cart(ps_lo(1):ps_hi(1),ps_lo(2):ps_hi(2),ps_lo(3):ps_hi(3))
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: domhi(3)
    integer  , value, intent(in   ) :: is_prediction, add_thermal

    ! Local variable
    double precision :: divup, p0divu, gradp0,veladv
    integer          :: i,j,k

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

            if (j .lt. base_cutoff_density_coord(lev)) then
                gradp0 = rho0_cart(i,j,k) * grav_cart(i,j,k)
            else if (j .eq. domhi(2)) then
                ! NOTE: this should be zero since p0 is constant up here
                gradp0 = ( p0_cart(i,j,k) - p0_cart(i,j-1,k) ) / dx(2)
            else
                ! NOTE: this should be zero since p0 is constant up here
                gradp0 = ( p0_cart(i,j+1,k) - p0_cart(i,j,k) ) / dx(2)
            end if

            veladv = 0.5d0*(vmac(i,j,k)+vmac(i,j+1,k))
            rhoh_force(i,j,k) =  veladv * gradp0
#else
            if (spherical .eq. 0) then 

                if (k .lt. base_cutoff_density_coord(lev)) then
                    gradp0 = rho0_cart(i,j,k) * grav_cart(i,j,k)
                 else if (k .eq. domhi(3)) then
                    ! NOTE: this should be zero since p0 is constant up here
                    gradp0 = (p0_cart(i,j,k) - p0_cart(i,j,k-1) ) / dx(3)
                 else
                    ! NOTE: this should be zero since p0 is constant up here
                    gradp0 = (p0_cart(i,j,k+1) - p0_cart(i,j,k) ) / dx(3)
                 end if
          
                veladv = 0.5d0*(wmac(i,j,k)+wmac(i,j,k+1))
                rhoh_force(i,j,k) = veladv * gradp0
            else 

                divup = (umac(i+1,j,k) * p0macx(i+1,j,k) - umac(i,j,k) * p0macx(i,j,k)) / dx(1) + &
                    (vmac(i,j+1,k) * p0macy(i,j+1,k) - vmac(i,j,k) * p0macy(i,j,k)) / dx(2) + &
                    (wmac(i,j,k+1) * p0macz(i,j,k+1) - wmac(i,j,k) * p0macz(i,j,k)) / dx(3)

                p0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                    (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                    (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * p0_cart(i,j,k)

                rhoh_force(i,j,k) = divup - p0divu

            endif
#endif
          end do
       end do
    end do
    !
    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h or rhoh
    !
    ! If use_exact_base_state or average_base_state is on, psi is instead dpdt term
    !
    if ((is_prediction .eq. 1 .AND. enthalpy_pred_type == predict_h) .OR. &
         (is_prediction .eq. 1 .AND. enthalpy_pred_type == predict_rhoh) .OR. &
         (is_prediction .eq. 0)) then

        do k = lo(3),hi(3)
           do j = lo(2),hi(2)
              do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi_cart(i,j,k)
             enddo
          enddo
       enddo
    endif

    if (add_thermal .eq. 1) then
        do k = lo(3),hi(3)
           do j = lo(2),hi(2)
              do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
    end if

  end subroutine mkrhohforce

  subroutine modify_scal_force(lo, hi, lev, &
    force, f_lo, f_hi, &
    scal,  s_lo, s_hi, &
    umac,  u_lo, u_hi, &
    vmac,  v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
    wmac,  w_lo, w_hi, &
#endif
    s0_cart, s0_lo, s0_hi, &
    s0_edge_cart, se_lo, se_hi, &
    w0_cart, w0_lo, w0_hi, &
    dx, do_fullform) &
    bind(C,name="modify_scal_force")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
#endif
    integer         , intent(in   ) :: s0_lo(3), s0_hi(3)
    integer         , intent(in   ) :: se_lo(3), se_hi(3)
    integer         , intent(in   ) :: w0_lo(3), w0_hi(3)
    double precision, intent(inout) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent(in   ) :: s0_cart(s0_lo(1):s0_hi(1),s0_lo(2):s0_hi(2),s0_lo(3):s0_hi(3))
    double precision, intent(in   ) :: s0_edge_cart(se_lo(1):se_hi(1),se_lo(2):se_hi(2),se_lo(3):se_hi(3))
    double precision, intent(in   ) :: w0_cart(w0_lo(1):w0_hi(1),w0_lo(2):w0_hi(2),w0_lo(3):w0_hi(3),AMREX_SPACEDIM)
    double precision, intent(in   ) :: dx(3)
    integer  , value, intent(in   ) :: do_fullform

    ! local
    integer :: i,j,k
    double precision :: divu,divs0u

    !$gpu

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)

                ! umac does not contain w0
#if (AMREX_SPACEDIM == 2)
                divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2)

                ! add w0 contribution
                divu = divu + (w0_cart(i,j+1,k,AMREX_SPACEDIM)-w0_cart(i,j,k,AMREX_SPACEDIM))/dx(AMREX_SPACEDIM)
#elif (AMREX_SPACEDIM == 3)
                divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                    +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)

                ! add w0 contribution
                divu = divu + (w0_cart(i,j,k+1,AMREX_SPACEDIM)-w0_cart(i,j,k,AMREX_SPACEDIM))/dx(AMREX_SPACEDIM)
#endif

                if (do_fullform .eq. 1) then
                    force(i,j,k) = force(i,j,k) - scal(i,j,k)*divu
                else

#if (AMREX_SPACEDIM == 2)
                    divs0u = s0_cart(i,j,k)*(  umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                        +(vmac(i,j+1,k) * s0_edge_cart(i,j+1,k) - &
                        vmac(i,j  ,k) * s0_edge_cart(i,j,k) )/ dx(2)
#elif (AMREX_SPACEDIM == 3)
                    divs0u = s0_cart(i,j,k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                        +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                        +(wmac(i,j,k+1) * s0_edge_cart(i,j,k+1) &
                        - wmac(i,j,k  ) * s0_edge_cart(i,j,k))/ dx(3)
#endif

                    force(i,j,k) = force(i,j,k) - (scal(i,j,k)-s0_cart(i,j,k))*divu - divs0u
                endif

            end do
        end do
    end do

  end subroutine modify_scal_force

  subroutine modify_scal_force_sphr(lo, hi, domlo, domhi, &
       force, f_lo, f_hi, &
       scal,  s_lo, s_hi, &
       umac,  u_lo, u_hi, &
       vmac,  v_lo, v_hi, &
       wmac,  w_lo, w_hi, &
       s0_cart,s0_lo, s0_hi, &
       dx, do_fullform, &
       divu_cart, d_lo, d_hi) &
       bind(C,name="modify_scal_force_sphr")

    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    integer         , intent(in   ) :: s0_lo(3), s0_hi(3)
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    double precision, intent(inout) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(in   ) :: s0_cart(s0_lo(1):s0_hi(1),s0_lo(2):s0_hi(2),s0_lo(3):s0_hi(3))
    double precision, intent(in   ) :: dx(3)
    integer  , value, intent(in   ) :: do_fullform
    double precision, intent(in   ) :: divu_cart(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i,j,k,r
    double precision :: divumac,divs0u
    double precision :: s0_xlo,s0_xhi
    double precision :: s0_ylo,s0_yhi
    double precision :: s0_zlo,s0_zhi

    !$gpu

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             ! umac does not contain w0
             divumac = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                  +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                  +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)


             if (do_fullform .eq. 1) then

                force(i,j,k) = force(i,j,k) - scal(i,j,k)*(divumac+divu_cart(i,j,k))
             else

                if (i.lt.domhi(1)) then
                   s0_xhi = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i+1,j,k))
                else
                   s0_xhi = s0_cart(i,j,k)
                end if
                if (i.gt.domlo(1)) then
                   s0_xlo = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i-1,j,k))
                else
                   s0_xlo = s0_cart(i,j,k)
                end if

                if (j.lt.domhi(2)) then
                   s0_yhi = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j+1,k))
                else
                   s0_yhi = s0_cart(i,j,k)
                end if
                if (j.gt.domlo(2)) then
                   s0_ylo = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j-1,k))
                else
                   s0_ylo = s0_cart(i,j,k)
                end if

                if (k.lt.domhi(3)) then
                   s0_zhi = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j,k+1))
                else
                   s0_zhi = s0_cart(i,j,k)
                end if
                if (k.gt.domlo(3)) then
                   s0_zlo = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j,k-1))
                else
                   s0_zlo = s0_cart(i,j,k)
                end if

                divs0u = (umac(i+1,j,k)*s0_xhi - umac(i,j,k)*s0_xlo)/dx(1) + &
                     (vmac(i,j+1,k)*s0_yhi - vmac(i,j,k)*s0_ylo)/dx(2) + &
                     (wmac(i,j,k+1)*s0_zhi - wmac(i,j,k)*s0_zlo)/dx(3)

                force(i,j,k) = force(i,j,k) - divs0u &
                     -(scal(i,j,k)-s0_cart(i,j,k))*(divumac+divu_cart(i,j,k))

             endif

          end do
       end do
    end do

  end subroutine modify_scal_force_sphr

end module make_scal_force_module
