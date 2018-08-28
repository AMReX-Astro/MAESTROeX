module make_vel_force_module

  use meth_params_module, only: base_cutoff_density,buoyancy_cutoff_factor, prob_lo, rotation_radius
  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr, center
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr, put_1d_array_on_cart
  use bl_constants_module
  use rotation_module, only: sin_theta, cos_theta, omega

  implicit none

  private

contains

  subroutine make_vel_force(lev, lo, hi, &
                            is_final_update, &
                            vel_force, f_lo, f_hi, nc_f, &
                            gpi, g_lo, g_hi, nc_g, &
                            rho, r_lo, r_hi, &
                            uedge, u_lo, u_hi, &
                            vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                            wedge, w_lo, w_hi, &
#endif
                            uold, uo_lo, uo_hi, nc_uo, &
                            w0,w0_force,rho0,grav, &
                            do_add_utilde_force) &
                            bind(C, name="make_vel_force")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: is_final_update
    integer         , intent (in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent (in   ) :: g_lo(3), g_hi(3), nc_g
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: v_lo(3), v_hi(3)
#if (AMREX_SPACEDIM == 3)
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
#endif
    integer         , intent (in   ) :: uo_lo(3), uo_hi(3), nc_uo
    double precision, intent (inout) :: vel_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent (in   ) ::       gpi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
    double precision, intent (in   ) ::       rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) ::     uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent (in   ) ::     vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent (in   ) ::     wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent (in   ) :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),nc_uo)
    double precision, intent (in   ) ::       w0(0:max_radial_level,0:nr_fine)
    double precision, intent (in   ) :: w0_force(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::     rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::     grav(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: do_add_utilde_force

    ! local
    integer :: i,j,k,r
    double precision :: rhopert

    real(kind=dp_t) :: coriolis_term(3), centrifugal_term(3)

    vel_force = 0.d0

    ! CURRENTLY for rotation in plane-parallel, we make the (bad) assueption
    ! that all points within the patch have the same centrifugal forcing terms.
    !
    ! We assume the centrifugal term applies at a constant radius,
    ! rotation_radius, for the patch.  In otherwords, the patch lives on the
    ! surface of a sphere of radius rotation_radius.
    !
    ! Furthermore, we assume the patch lives at longitude = 0.
    !
    ! Then the orientation of the patch is such that e_z is in the
    ! outward radial direction of the star, e_x is in the co_latitude (polar)
    ! angle direction and e_y is in the global y-direction.
    !
    ! centrifugal_term = omega x (omega x r) = (omega dot r) * omega
    !                                          - omega^2 * r
    ! where omega = (-|omega| sin_theta) e_x + (|omega| cos_theta) e_z
    !           r = rotation_radius e_z
    !
    ! See docs/rotation for derivation and figures.
    !

    centrifugal_term(1) = - omega**2 * rotation_radius * sin_theta * sin_theta
    centrifugal_term(2) = ZERO
    centrifugal_term(3) = omega**2 * rotation_radius * cos_theta * sin_theta &
                          - omega**2 * rotation_radius

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

       rhopert = rho(i,j,k) - rho0(lev,r)

       ! cutoff the buoyancy term if we are outside of the star
       if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
          rhopert = 0.d0
       end if

       ! the coriolis term is:
       !    TWO * omega x U
       ! where omega is given above and U = (u, v, w) is the velocity
#if (AMREX_SPACEDIM == 3)
       if (is_final_update) then

          ! use uedge so we are time-centered
          coriolis_term(1) = -TWO * omega * &
               HALF*(vedge(i,j,k) + vedge(i,j+1,k)) * cos_theta

          coriolis_term(2) =  TWO * omega * &
               (HALF*(wedge(i,j,k)   + w0(lev,k) + &
                      wedge(i,j,k+1) + w0(lev,k+1)) * sin_theta + &
                HALF*(uedge(i,j,k) + uedge(i+1,j,k)) * cos_theta)

          coriolis_term(3) = -TWO * omega * &
               HALF*(vedge(i,j,k) + vedge(i,j+1,k)) * sin_theta

       else
          coriolis_term(1) = -TWO * omega * uold(i,j,k,2) * cos_theta

          coriolis_term(2) =  TWO * omega * ((uold(i,j,k,3) + HALF*(w0(lev,k) + w0(lev,k+1))) * sin_theta + &
                                             uold(i,j,k,1) * cos_theta)

          coriolis_term(3) = -TWO * omega * uold(i,j,k,2) * sin_theta

       endif
#endif
       ! note: if use_alt_energy_fix = T, then gphi is already
       ! weighted by beta0
       vel_force(i,j,k,1:AMREX_SPACEDIM-1) =  -coriolis_term(1:AMREX_SPACEDIM-1) - centrifugal_term(1:AMREX_SPACEDIM-1) - &
            gpi(i,j,k,1:AMREX_SPACEDIM-1) / rho(i,j,k)

#if (AMREX_SPACEDIM == 3)
        vel_force(i,j,k,AMREX_SPACEDIM) = -coriolis_term(3) - centrifugal_term(3) + &
             ( rhopert * grav(lev,r) - gpi(i,j,k,AMREX_SPACEDIM) ) / rho(i,j,k) - w0_force(lev,r)
#else
        vel_force(i,j,k,AMREX_SPACEDIM) = &
             ( rhopert * grav(lev,r) - gpi(i,j,k,AMREX_SPACEDIM) ) / rho(i,j,k) - w0_force(lev,r)
#endif

    end do
    end do
    end do

    if (do_add_utilde_force .eq. 1) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

#if (AMREX_SPACEDIM == 1)
       r = i
#elif (AMREX_SPACEDIM == 2)
       r = j
#elif (AMREX_SPACEDIM == 3)
       r = k
#endif

             if (r .le. -1) then
                ! do not modify force since dw0/dr=0
             else if (r .ge. nr(lev)) then
                ! do not modify force since dw0/dr=0
             else

#if (AMREX_SPACEDIM == 2)
                vel_force(i,j,k,2) = vel_force(i,j,k,2) &
                     - (vedge(i,j+1,k)+vedge(i,j,k))*(w0(lev,r+1)-w0(lev,r)) / (2.d0*dr(lev))

#else
                vel_force(i,j,k,3) = vel_force(i,j,k,3) &
                     - (wedge(i,j,k+1)+wedge(i,j,k))*(w0(lev,r+1)-w0(lev,r)) / (2.d0*dr(lev))

#endif
             end if

       end do
       end do
       end do
    endif

  end subroutine make_vel_force

  subroutine make_vel_force_sphr(lo, hi, &
                                 is_final_update, &
                                 vel_force, f_lo, f_hi, nc_f, &
                                 gpi, g_lo, g_hi, nc_g, &
                                 rho, r_lo, r_hi, &
                                 uedge, u_lo, u_hi, &
                                 vedge, v_lo, v_hi, &
                                 wedge, w_lo, w_hi, &
                                 uold, uo_lo, uo_hi, nc_uo, &
                                 normal, n_lo, n_hi, nc_n, &
                                 w0x_cart, wcx_lo, wcx_hi, &
                                 w0y_cart, wcy_lo, wcy_hi, &
                                 gradw0_cart, gw_lo, gw_hi, &
                                 w0_force_cart, wf_lo, wf_hi, nc_wf, &
                                 w0macx, w0x_lo, w0x_hi, &
                                 w0macy, w0y_lo, w0y_hi, &
                                 rho0, grav, &
                                 dx, &
                                 r_cc_loc, r_edge_loc, &
                                 cc_to_r, ccr_lo, ccr_hi, &
                                 do_add_utilde_force) &
                                 bind(C, name="make_vel_force_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: is_final_update
    integer         , intent (in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent (in   ) :: g_lo(3), g_hi(3), nc_g
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: v_lo(3), v_hi(3)
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
    integer         , intent (in   ) :: uo_lo(3), uo_hi(3), nc_uo
    integer         , intent (in   ) :: n_lo(3), n_hi(3), nc_n
    integer         , intent (in   ) :: wcx_lo(3), wcx_hi(3)
    integer         , intent (in   ) :: wcy_lo(3), wcy_hi(3)
    integer         , intent (in   ) :: gw_lo(3), gw_hi(3)
    integer         , intent (in   ) :: wf_lo(3), wf_hi(3), nc_wf
    integer         , intent (in   ) :: w0x_lo(3), w0x_hi(3)
    integer         , intent (in   ) :: w0y_lo(3), w0y_hi(3)
    double precision, intent (inout) :: vel_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent (in   ) ::       gpi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
    double precision, intent (in   ) ::       rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) ::     uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent (in   ) ::     vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent (in   ) ::     wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent (in   ) ::    uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),nc_uo)
    double precision, intent (in   ) ::    normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),nc_n)
    double precision, intent (in   ) ::    w0x_cart(wcx_lo(1):wcx_hi(1),wcx_lo(2):wcx_hi(2),wcx_lo(3):wcx_hi(3))
    double precision, intent (in   ) ::    w0y_cart(wcy_lo(1):wcy_hi(1),wcy_lo(2):wcy_hi(2),wcy_lo(3):wcy_hi(3))
    double precision, intent (in   ) :: gradw0_cart(gw_lo(1):gw_hi(1),gw_lo(2):gw_hi(2),gw_lo(3):gw_hi(3))
    double precision, intent (in   ) :: w0_force_cart(wf_lo(1):wf_hi(1),wf_lo(2):wf_hi(2),wf_lo(3):wf_hi(3),nc_wf)
    double precision, intent (in   ) :: w0macx(w0x_lo(1):w0x_hi(1),w0x_lo(2):w0x_hi(2),w0x_lo(3):w0x_hi(3))
    double precision, intent (in   ) :: w0macy(w0y_lo(1):w0y_hi(1),w0y_lo(2):w0y_hi(2),w0y_lo(3):w0y_hi(3))
    double precision, intent (in   ) ::     rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::     grav(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: dx(3)
    double precision, intent (in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent (in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent (in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
                                               ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))
    integer         , intent (in   ) :: do_add_utilde_force


    integer         :: i,j,k

    double precision, allocatable :: rho0_cart(:,:,:,:)
    double precision, allocatable :: grav_cart(:,:,:,:)


    double precision :: rhopert
    double precision :: xx, yy, zz
    real(kind=dp_t) :: centrifugal_term(3), coriolis_term(3)

    real(kind=dp_t) :: Ut_dot_er

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(grav_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))

    vel_force = ZERO

    call put_1d_array_on_cart_sphr(lo,hi,rho0_cart,lo,hi,1,rho0,dx,0,0,r_cc_loc,r_edge_loc, &
                                      cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,grav_cart,lo,hi,3,grav,dx,0,1,r_cc_loc,r_edge_loc, &
                                      cc_to_r,ccr_lo,ccr_hi)

    !$OMP PARALLEL DO PRIVATE(i,j,k,xx,yy,zz,rhopert,centrifugal_term,coriolis_term)
    do k = lo(3),hi(3)
       zz = prob_lo(3) + (dble(k) + HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          yy = prob_lo(2) + (dble(j) + HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             xx = prob_lo(1) + (dble(i) + HALF)*dx(1) - center(1)

             rhopert = rho(i,j,k) - rho0_cart(i,j,k,1)

             ! cutoff the buoyancy term if we are outside of the star
             if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                rhopert = 0.d0
             end if

             ! Coriolis and centrifugal forces.  We assume that the
             ! rotation axis is the z direction, with angular velocity
             ! omega

             ! omega x (omega x r ) = - omega^2 x e_x  - omega^2 y e_y
             ! (with omega = omega e_z)
             centrifugal_term(1) = -omega * omega * xx
             centrifugal_term(2) = -omega * omega * yy
             centrifugal_term(3) = ZERO

             ! cutoff the centrifugal term if we are outside the star
             if (rho(i,j,k) .lt. buoyancy_cutoff_factor*base_cutoff_density) then
                centrifugal_term(:) = 0.d0
             end if

             ! 2 omega x U = - 2 omega v e_x  + 2 omega u e_y
             ! (with omega = omega e_z)
             if (is_final_update .eq. 1) then

                ! use uedge so we are time-centered
                coriolis_term(1) = -TWO * omega * &
                     HALF*(vedge(i,j,k)   + w0macy(i,j,k) + &
                           vedge(i,j+1,k) + w0macy(i,j+1,k))

                coriolis_term(2) =  TWO * omega * &
                     HALF*(uedge(i,j,k)   + w0macx(i,j,k) + &
                           uedge(i+1,j,k) + w0macx(i+1,j,k))

                coriolis_term(3) = ZERO

             else
                coriolis_term(1) = -TWO * omega * (uold(i,j,k,2) + w0y_cart(i,j,k))
                coriolis_term(2) =  TWO * omega * (uold(i,j,k,1) + w0x_cart(i,j,k))
                coriolis_term(3) = ZERO
             endif


             ! F_Coriolis = -2 omega x U
             ! F_centrifugal = - omega x (omega x r)

             ! we just computed the absolute value of the forces above, so use
             ! the right sign here


             ! note: if use_alt_energy_fix = T, then gphi is already weighted
             ! by beta0
             vel_force(i,j,k,1) = -coriolis_term(1) - centrifugal_term(1) + &
                  ( rhopert * grav_cart(i,j,k,1) - gpi(i,j,k,1) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,1)

             vel_force(i,j,k,2) = -coriolis_term(2) - centrifugal_term(2) + &
                  ( rhopert * grav_cart(i,j,k,2) - gpi(i,j,k,2) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,2)

             vel_force(i,j,k,3) = -coriolis_term(3) - centrifugal_term(3) + &
                  ( rhopert * grav_cart(i,j,k,3) - gpi(i,j,k,3) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,3)

          end do
       end do
    end do
    !$OMP END PARALLEL DO


    if (do_add_utilde_force .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,Ut_dot_er)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                Ut_dot_er = &
                     HALF*(uedge(i,j,k)+uedge(i+1,j  ,k  ))*normal(i,j,k,1) + &
                     HALF*(vedge(i,j,k)+vedge(i  ,j+1,k  ))*normal(i,j,k,2) + &
                     HALF*(wedge(i,j,k)+wedge(i  ,j,  k+1))*normal(i,j,k,3)

                vel_force(i,j,k,1) = vel_force(i,j,k,1) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                vel_force(i,j,k,2) = vel_force(i,j,k,2) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                vel_force(i,j,k,3) = vel_force(i,j,k,3) - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    deallocate(rho0_cart,grav_cart)

  end subroutine make_vel_force_sphr

end module make_vel_force_module
