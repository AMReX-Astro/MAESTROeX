module update_vel_module

  use amrex_constants_module
  use base_state_geometry_module, only:  max_radial_level, nr_fine
  use meth_params_module, only: do_sponge, spherical

  implicit none

  private

contains

  subroutine update_velocity(lo, hi, lev, &
       uold, uo_lo, uo_hi, &
       unew, un_lo, un_hi, &
       umac,  u_lo, u_hi, &
       vmac,  v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wmac,  w_lo, w_hi, &
#endif
       uedgex, x_lo, x_hi, &
       uedgey, y_lo, y_hi, &
#if (AMREX_SPACEDIM == 3)
       uedgez, z_lo, z_hi, &
#endif
       force,  f_lo, f_hi, &
       sponge, s_lo, s_hi, &
       w0, dx, dt) &
       bind(C,name="update_velocity")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: uo_lo(3), uo_hi(3)
    double precision, intent(in   ) :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: un_lo(3), un_hi(3)
    double precision, intent(inout) :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) :: uedgex (x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) :: uedgey (y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) :: uedgez (z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),AMREX_SPACEDIM)
#endif
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: sponge (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: w0     (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(AMREX_SPACEDIM)
    double precision, value, intent(in   ) :: dt

    integer :: i,j,k, dim
    double precision :: ubar, vbar, wbar, w0bar
    double precision :: ugradu, ugradv, ugradw

    !$gpu

    ! 1) Subtract (Utilde dot grad) Utilde term from old Utilde
    ! 2) Add forcing term to new Utilde
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             ! create cell-centered Utilde
             dim = i
             ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
             dim = j
             vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
#if (AMREX_SPACEDIM == 3)
             dim = k
             wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))
#endif

             ! create (Utilde dot grad) Utilde
             ugradu = ( ubar*(uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1) &
                  + vbar*(uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2) &
#if (AMREX_SPACEDIM == 3)
                  + wbar*(uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3) &
#endif
                  )

             ugradv = ( ubar*(uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1) &
                  + vbar*(uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2) &
#if (AMREX_SPACEDIM == 3)
                  + wbar*(uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3) &
#endif
                  )

#if (AMREX_SPACEDIM == 3)
             ugradw = ubar*(uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1) &
                  + vbar*(uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2) &
                  + wbar*(uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)
#endif

             ! update with (Utilde dot grad) Utilde and force
             unew(i,j,k,1) = uold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
             unew(i,j,k,2) = uold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
#if (AMREX_SPACEDIM == 3)
             unew(i,j,k,3) = uold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)
#endif

             ! subtract (w0 dot grad) Utilde term
             w0bar = HALF*(w0(lev,dim) + w0(lev,dim+1))
             unew(i,j,k,:) = unew(i,j,k,:) - dt * w0bar * &
#if (AMREX_SPACEDIM == 2)
             (uedgey(i,j+1,k,:) - uedgey(i,j,k,:))/dx(2)
#elif (AMREX_SPACEDIM == 3)
             (uedgez(i,j,k+1,:) - uedgez(i,j,k,:))/dx(3)
#endif

             ! Add the sponge
             if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

          end do
       end do
    end do

  end subroutine update_velocity

  subroutine update_velocity_sphr(lo, hi, &
       uold, uo_lo, uo_hi, &
       unew, un_lo, un_hi, &
       umac,  u_lo, u_hi, &
       vmac,  v_lo, v_hi, &
       wmac,  w_lo, w_hi, &
       uedgex, x_lo, x_hi, &
       uedgey, y_lo, y_hi, &
       uedgez, z_lo, z_hi, &
       force,  f_lo, f_hi, &
       sponge, s_lo, s_hi, &
       w0, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       dx, dt) bind(C,name="update_velocity_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: uo_lo(3), uo_hi(3)
    double precision, intent(in   ) :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: un_lo(3), un_hi(3)
    double precision, intent(inout) :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) :: uedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) :: uedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) :: uedgez(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: sponge(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: w0    (0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    double precision, intent(in   ) :: dx(3)
    double precision, value, intent(in   ) :: dt

    ! Local variables
    integer :: i, j, k
    double precision :: ubar, vbar, wbar
    double precision :: ugradu, ugradv, ugradw
    double precision :: gradux,graduy,graduz
    double precision :: gradvx,gradvy,gradvz
    double precision :: gradwx,gradwy,gradwz
    double precision :: w0_gradur,w0_gradvr,w0_gradwr

    !$gpu

    ! 1) Subtract (Utilde dot grad) Utilde term from old Utilde
    ! 2) Add forcing term to new Utilde
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! create cell-centered Utilde
             ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
             vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
             wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))

             ! create (Utilde dot grad) Utilde
             ugradu = ubar*(uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3)

             ugradv = ubar*(uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3)

             ugradw = ubar*(uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1) + &
                  vbar*(uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2) + &
                  wbar*(uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)

             ! update with (Utilde dot grad) Utilde and force
             unew(i,j,k,1) = uold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
             unew(i,j,k,2) = uold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
             unew(i,j,k,3) = uold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Subtract (w0 dot grad) Utilde term from new Utilde
             gradux = (uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1)
             gradvx = (uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1)
             gradwx = (uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1)

             graduy = (uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2)
             gradvy = (uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2)
             gradwy = (uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2)

             graduz = (uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3)
             gradvz = (uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3)
             gradwz = (uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)

             w0_gradur = gradux * HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) &
                  + graduy * HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) &
                  + graduz * HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))

             w0_gradvr = gradvx * HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) &
                  + gradvy * HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) &
                  + gradvz * HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))

             w0_gradwr = gradwx * HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) &
                  + gradwy * HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) &
                  + gradwz * HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))

             unew(i,j,k,1) = unew(i,j,k,1) - dt * w0_gradur
             unew(i,j,k,2) = unew(i,j,k,2) - dt * w0_gradvr
             unew(i,j,k,3) = unew(i,j,k,3) - dt * w0_gradwr

             ! Add the sponge
             if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

          enddo
       enddo
    enddo

  end subroutine update_velocity_sphr

end module update_vel_module
