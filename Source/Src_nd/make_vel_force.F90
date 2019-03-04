module make_vel_force_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only: base_cutoff_density,buoyancy_cutoff_factor, prob_lo
  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr, center
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr
  use amrex_constants_module

  implicit none

  private

contains

  subroutine make_vel_force(lev, lo, hi, &
       vel_force, f_lo, f_hi, nc_f, &
       gpi, g_lo, g_hi, nc_g, &
       rho, r_lo, r_hi, &
       uedge, u_lo, u_hi, &
       vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wedge, w_lo, w_hi, &
#endif
       w0,w0_force,rho0,grav, &
       do_add_utilde_force) &
       bind(C, name="make_vel_force")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent (in   ) :: g_lo(3), g_hi(3), nc_g
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: v_lo(3), v_hi(3)
#if (AMREX_SPACEDIM == 3)
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
#endif
    double precision, intent (inout) :: vel_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent (in   ) ::       gpi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
    double precision, intent (in   ) ::       rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) ::     uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent (in   ) ::     vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent (in   ) ::     wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent (in   ) ::       w0(0:max_radial_level,0:nr_fine)
    double precision, intent (in   ) :: w0_force(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::     rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::     grav(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: do_add_utilde_force

    ! local
    integer :: i,j,k,r
    double precision :: rhopert

    vel_force(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_f) = 0.d0

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

             ! note: if use_alt_energy_fix = T, then gphi is already
             ! weighted by beta0
             vel_force(i,j,k,1:AMREX_SPACEDIM-1) = - gpi(i,j,k,1:AMREX_SPACEDIM-1) / rho(i,j,k)

             vel_force(i,j,k,AMREX_SPACEDIM) = &
                  ( rhopert * grav(lev,r) - gpi(i,j,k,AMREX_SPACEDIM) ) / rho(i,j,k) - w0_force(lev,r)

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
#if (AMREX_SPACEDIM == 1)
                   vel_force(i,j,k,1) = vel_force(i,j,k,1) &
                        - (uedge(i+1,j,k)+uedge(i,j,k))*(w0(lev,r+1)-w0(lev,r)) / (2.d0*dr(lev))
#elif (AMREX_SPACEDIM == 2)
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
       vel_force, f_lo, f_hi, nc_f, &
       gpi, g_lo, g_hi, nc_g, &
       rho, r_lo, r_hi, &
       uedge, u_lo, u_hi, &
       vedge, v_lo, v_hi, &
       wedge, w_lo, w_hi, &
       normal, n_lo, n_hi, nc_n, &
       gradw0_cart, gw_lo, gw_hi, &
       w0_force_cart, wf_lo, wf_hi, nc_wf, &
       rho0, grav, &
       dx, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi, &
       do_add_utilde_force) &
       bind(C, name="make_vel_force_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent (in   ) :: g_lo(3), g_hi(3), nc_g
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: v_lo(3), v_hi(3)
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3), nc_n
    integer         , intent (in   ) :: gw_lo(3), gw_hi(3)
    integer         , intent (in   ) :: wf_lo(3), wf_hi(3), nc_wf
    double precision, intent (inout) :: vel_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent (in   ) ::       gpi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
    double precision, intent (in   ) ::       rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) ::     uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent (in   ) ::     vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent (in   ) ::     wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent (in   ) ::    normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),nc_n)
    double precision, intent (in   ) :: gradw0_cart(gw_lo(1):gw_hi(1),gw_lo(2):gw_hi(2),gw_lo(3):gw_hi(3))
    double precision, intent (in   ) :: w0_force_cart(wf_lo(1):wf_hi(1),wf_lo(2):wf_hi(2),wf_lo(3):wf_hi(3),nc_wf)
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

    double precision, pointer :: rho0_cart(:,:,:,:)
    double precision, pointer :: grav_cart(:,:,:,:)

    double precision :: rhopert
    double precision :: xx, yy, zz

    double precision :: Ut_dot_er

    call bl_allocate(rho0_cart,lo,hi,1)
    call bl_allocate(grav_cart,lo,hi,3)

    vel_force(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_f) = ZERO

    call put_1d_array_on_cart_sphr(lo,hi,rho0_cart,lo,hi,1,rho0,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,grav_cart,lo,hi,3,grav,dx,0,1,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)

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

             ! note: if use_alt_energy_fix = T, then gphi is already weighted
             ! by beta0
             vel_force(i,j,k,1) = ( rhopert * grav_cart(i,j,k,1) - gpi(i,j,k,1) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,1)

             vel_force(i,j,k,2) = ( rhopert * grav_cart(i,j,k,2) - gpi(i,j,k,2) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,2)

             vel_force(i,j,k,3) = ( rhopert * grav_cart(i,j,k,3) - gpi(i,j,k,3) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,3)

          end do
       end do
    end do


    if (do_add_utilde_force .eq. 1) then

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

    endif

    call bl_deallocate(rho0_cart)
    call bl_deallocate(grav_cart)

  end subroutine make_vel_force_sphr

  
  subroutine make_vel_force_noenergyfix_sphr(lo, hi, &
       vel_force, f_lo, f_hi, nc_f, &
       gpi, g_lo, g_hi, nc_g, &
       rho, r_lo, r_hi, &
       uedge, u_lo, u_hi, &
       vedge, v_lo, v_hi, &
       wedge, w_lo, w_hi, &
       normal, n_lo, n_hi, nc_n, &
       gradw0_cart, gw_lo, gw_hi, &
       w0_force_cart, wf_lo, wf_hi, nc_wf, &
       rho0, grav, beta0, &
       dx, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi, &
       do_add_utilde_force) &
       bind(C, name="make_vel_force_noenergyfix_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent (in   ) :: g_lo(3), g_hi(3), nc_g
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: v_lo(3), v_hi(3)
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3), nc_n
    integer         , intent (in   ) :: gw_lo(3), gw_hi(3)
    integer         , intent (in   ) :: wf_lo(3), wf_hi(3), nc_wf
    double precision, intent (inout) :: vel_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent (in   ) ::       gpi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
    double precision, intent (in   ) ::       rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) ::     uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent (in   ) ::     vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent (in   ) ::     wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent (in   ) ::    normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),nc_n)
    double precision, intent (in   ) :: gradw0_cart(gw_lo(1):gw_hi(1),gw_lo(2):gw_hi(2),gw_lo(3):gw_hi(3))
    double precision, intent (in   ) :: w0_force_cart(wf_lo(1):wf_hi(1),wf_lo(2):wf_hi(2),wf_lo(3):wf_hi(3),nc_wf)
    double precision, intent (in   ) ::     rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::     grav(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) ::    beta0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: dx(3)
    double precision, intent (in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent (in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent (in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))
    integer         , intent (in   ) :: do_add_utilde_force

    integer         :: i,j,k

    double precision, pointer :: rho0_cart(:,:,:,:)
    double precision, pointer :: grav_cart(:,:,:,:)
    double precision, pointer :: beta0_cart(:,:,:,:)

    double precision :: rhopert
    double precision :: xx, yy, zz

    double precision :: Ut_dot_er

    call bl_allocate(rho0_cart,lo,hi,1)
    call bl_allocate(grav_cart,lo,hi,3)
    call bl_allocate(beta0_cart,lo,hi,1)

    vel_force(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_f) = ZERO

    call put_1d_array_on_cart_sphr(lo,hi,rho0_cart,lo,hi,1,rho0,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,grav_cart,lo,hi,3,grav,dx,0,1,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,beta0_cart,lo,hi,1,beta0,dx,0,0,r_cc_loc,r_edge_loc, &
         cc_to_r,ccr_lo,ccr_hi)

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

             ! assume use_alt_energy_fix = T, but want use_alt_energy_fix = F
             ! hence, need to weigh gphi by beta0
             vel_force(i,j,k,1) = ( rhopert * grav_cart(i,j,k,1) - gpi(i,j,k,1)/beta0_cart(i,j,k,1) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,1)

             vel_force(i,j,k,2) = ( rhopert * grav_cart(i,j,k,2) - gpi(i,j,k,2)/beta0_cart(i,j,k,1) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,2)

             vel_force(i,j,k,3) = ( rhopert * grav_cart(i,j,k,3) - gpi(i,j,k,3)/beta0_cart(i,j,k,1) ) / rho(i,j,k) &
                  - w0_force_cart(i,j,k,3)

          end do
       end do
    end do


    if (do_add_utilde_force .eq. 1) then

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

    endif

    call bl_deallocate(rho0_cart)
    call bl_deallocate(grav_cart)
    call bl_deallocate(beta0_cart)

  end subroutine make_vel_force_noenergyfix_sphr
end module make_vel_force_module
