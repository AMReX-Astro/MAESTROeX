module make_vel_force_module

  use meth_params_module, only: base_cutoff_density,buoyancy_cutoff_factor
  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr

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

    vel_force = 0.d0

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

end module make_vel_force_module
