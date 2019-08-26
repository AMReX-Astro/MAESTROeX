module compute_dt_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, &
       cfl, use_soundspeed_firstdt, use_divu_firstdt, &
       use_exact_base_state
  use base_state_geometry_module, only:  max_radial_level, nr_fine, nr, dr

  implicit none

  private

contains

  subroutine estdt(lev, dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       dSdt,  t_lo, t_hi, &
       w0, p0, gamma1bar) bind (C,name="estdt")

    use amrex_fort_module, only: amrex_min, amrex_max

    integer  , value, intent(in   ) :: lev
    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u    (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: dSdt (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent(in   ) :: w0       (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)

    ! local variables
    double precision :: spdx, spdy, spdz, spdr, rho_min
    double precision :: fx, fy, fz, dt_temp
    double precision :: eps,denom,gradp0
    double precision :: a, b, c
    integer          :: i,j,k,r

    !$gpu

    rho_min = 1.d-20
    dt_temp = 1.e99

    eps = 1.d-8

    spdx    = 0.d0
    spdy    = 0.d0
    spdz    = 0.d0
    spdr    = 0.d0

    !
    ! Limit dt based on velocity terms.
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             spdx = max(spdx ,abs(u(i,j,k,1)))
#if (AMREX_SPACEDIM == 2)
             spdy = max(spdy ,abs(u(i,j,k,2)+0.5d0*(w0(lev,j)+w0(lev,j+1))))
#elif (AMREX_SPACEDIM == 3)
             spdy = max(spdy ,abs(u(i,j,k,2)))
             spdz = max(spdz ,abs(u(i,j,k,3)+0.5d0*(w0(lev,k)+w0(lev,k+1))))
#endif
          enddo
       enddo
    enddo

#if (AMREX_SPACEDIM == 2)
    do j = lo(2),hi(2)
       spdr = max(spdr ,abs(w0(lev,j)))
    enddo
#elif (AMREX_SPACEDIM == 3)
    do k = lo(3),hi(3)
       spdr = max(spdr ,abs(w0(lev,k)))
    enddo
#endif

    call amrex_max(umax, max(spdx,spdy,spdz,spdr))

    if (spdx > eps) dt_temp = min(dt_temp, dx(1)/spdx)
    if (spdy > eps) dt_temp = min(dt_temp, dx(2)/spdy)
    if (spdz > eps) dt_temp = min(dt_temp, dx(3)/spdz)
    if (spdr > eps) dt_temp = min(dt_temp, dx(AMREX_SPACEDIM)/spdr)

    dt_temp = dt_temp * cfl
    !
    ! Limit dt based on forcing terms
    !
    fx = 0.d0
    fy = 0.d0
    fz = 0.d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             fx = max(fx,abs(force(i,j,k,1)))
             fy = max(fy,abs(force(i,j,k,2)))
#if (AMREX_SPACEDIM == 3)
             fz = max(fz,abs(force(i,j,k,3)))
#endif
          enddo
       enddo
    enddo

    if (fx > eps) &
         dt_temp = min(dt_temp,sqrt(2.0d0*dx(1)/fx))

    if (fy > eps) &
         dt_temp = min(dt_temp,sqrt(2.0d0*dx(2)/fy))

#if (AMREX_SPACEDIM == 3)
    if (fz > eps) &
         dt_temp = min(dt_temp,sqrt(2.0d0*dx(3)/fz))
#endif

    !
    ! divU constraint
    !
    do k = lo(3), hi(3)

#if (AMREX_SPACEDIM == 3)
       if (k .eq. 0) then
          gradp0 = (p0(lev,k+1) - p0(lev,k))/dx(3)
       else if (k .eq. nr(lev)-1) then
          gradp0 = (p0(lev,k) - p0(lev,k-1))/dx(3)
       else
          gradp0 = 0.5d0*(p0(lev,k+1) - p0(lev,k-1))/dx(3)
       endif
       r = k
#endif

       do j = lo(2), hi(2)

#if (AMREX_SPACEDIM == 2)
          if (j .eq. 0) then
             gradp0 = (p0(lev,j+1) - p0(lev,j))/dx(2)
          else if (j .eq. nr(lev)-1) then
             gradp0 = (p0(lev,j) - p0(lev,j-1))/dx(2)
          else
             gradp0 = 0.5d0*(p0(lev,j+1) - p0(lev,j-1))/dx(2)
          endif
          r = j
#endif

          do i = lo(1), hi(1)

             denom = divU(i,j,k) - u(i,j,k,AMREX_SPACEDIM)*gradp0/(gamma1bar(lev,r)*p0(lev,r))

             if (denom > 0.d0) then
                dt_temp = min(dt_temp, 0.4d0*(1.d0 - rho_min/scal(i,j,k,rho_comp))/denom)
             endif

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             !
             ! An additional dS/dt timestep constraint originally
             ! used in nova
             ! solve the quadratic equation
             ! (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
             ! which is equivalent to
             ! (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
             ! which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
             !
             if (dSdt(i,j,k) .gt. 1.d-20) then
                a = 0.5d0*scal(i,j,k,rho_comp)*dSdt(i,j,k)
                b = scal(i,j,k,rho_comp)*divU(i,j,k)
                c = rho_min - scal(i,j,k,rho_comp)
                dt_temp = min(dt_temp,0.4d0*2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c)))
             endif

          enddo
       enddo
    enddo

    ! set dt to local min
    call amrex_min(dt, dt_temp)

  end subroutine estdt

  subroutine estdt_sphr(dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       dSdt,  t_lo, t_hi, &
       w0, &
       w0macx, x_lo, x_hi, &
       w0macy, y_lo, y_hi, &
       w0macz, z_lo, z_hi, &
       gp0_cart, g_lo, g_hi) bind (C,name="estdt_sphr")

    use amrex_fort_module, only: amrex_min, amrex_max

    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    integer         , intent(in   ) :: g_lo(3), g_hi(3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) :: scal  (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u     (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu  (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: dSdt  (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent(in   ) :: gp0_cart(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)
    double precision, intent(in   ) :: w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    double precision, intent(in   ) :: w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    double precision, intent(in   ) :: w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    double precision, intent(in   ) :: w0       (0:max_radial_level,0:nr_fine)

    double precision :: spdx, spdy, spdz, spdr, rho_min
    double precision :: gp_dot_u, dt_temp
    double precision :: fx, fy, fz, eps, denom, a, b, c
    integer          :: i,j,k,r

    !$gpu

    rho_min = 1.d-20
    dt_temp = 1.e99

    eps = 1.0d-8

    spdx = 0.d0
    spdy = 0.d0
    spdz = 0.d0
    spdr = 0.d0
    umax = 0.d0
    !
    ! Limit dt based on velocity terms
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             spdx = max(spdx ,abs(u(i,j,k,1)+0.5d0*(w0macx(i,j,k)+w0macx(i+1,j,k))))
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             spdy = max(spdy ,abs(u(i,j,k,2)+0.5d0*(w0macy(i,j,k)+w0macy(i,j+1,k))))
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             spdz = max(spdz ,abs(u(i,j,k,3)+0.5d0*(w0macz(i,j,k)+w0macz(i,j,k+1))))
          enddo
       enddo
    enddo

    do k=0,nr_fine
       spdr = max(spdr ,abs(w0(0,k)))
    enddo

    call amrex_max(umax,max(spdx,spdy,spdz,spdr))

    if (spdx > eps) dt_temp = min(dt_temp, dx(1)/spdx)
    if (spdy > eps) dt_temp = min(dt_temp, dx(2)/spdy)
    if (spdz > eps) dt_temp = min(dt_temp, dx(3)/spdz)
    if (spdr > eps) dt_temp = min(dt_temp, dr(0)/spdr)

    dt_temp = dt_temp*cfl

    ! Limit dt based on forcing terms
    fx = 0.d0
    fy = 0.d0
    fz = 0.d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             fx = max(fx,abs(force(i,j,k,1)))
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             fy = max(fy,abs(force(i,j,k,2)))
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             fz = max(fz,abs(force(i,j,k,3)))
          enddo
       enddo
    enddo

    if (fx > eps) &
         dt_temp = min(dt_temp,sqrt(2.0D0*dx(1)/fx))

    if (fy > eps) &
         dt_temp = min(dt_temp,sqrt(2.0D0*dx(2)/fy))

    if (fz > eps) &
         dt_temp = min(dt_temp,sqrt(2.0D0*dx(3)/fz))

    ! divU constraint
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             gp_dot_u = u(i,j,k,1) * gp0_cart(i,j,k,1) + &
                  u(i,j,k,2) * gp0_cart(i,j,k,2) + &
                  u(i,j,k,3) * gp0_cart(i,j,k,3)

             denom = divU(i,j,k) - gp_dot_u

             if (denom > 0.d0) then
                dt_temp = min(dt_temp,0.4d0*(1.d0 - rho_min/scal(i,j,k,rho_comp))/denom)
             endif
             !
             ! An additional dS/dt timestep constraint originally
             ! used in nova
             ! solve the quadratic equation
             ! (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
             ! which is equivalent to
             ! (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
             ! which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
             !
             if (dSdt(i,j,k) .gt. 1.d-20) then
                a = 0.5d0*scal(i,j,k,rho_comp)*dSdt(i,j,k)
                b = scal(i,j,k,rho_comp)*divU(i,j,k)
                c = rho_min - scal(i,j,k,rho_comp)
                dt_temp = min(dt_temp,0.4d0*2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c)))
             endif

          enddo
       enddo
    enddo

    ! set dt to local min
    call amrex_min(dt, dt_temp)

  end subroutine estdt_sphr

  subroutine firstdt(lev, dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       p0, gamma1bar) bind (C,name="firstdt")

    use amrex_fort_module, only: amrex_min, amrex_max

    integer, value  , intent(in   ) :: lev
    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u    (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)

    ! local variables
    double precision :: spdx,spdy,spdz,pforcex,pforcey,pforcez,ux,uy,uz
    double precision :: eps,dt_divu,dt_sound,gradp0,denom,rho_min
    integer          :: i,j,k

    integer :: pt_index(3)
    type(eos_t) :: eos_state

    !$gpu

    eps = 1.d-8

    rho_min = 1.d-20

    spdx    = 0.d0
    spdy    = 0.d0
    spdz    = 0.d0
    pforcex = 0.d0
    pforcey = 0.d0
    pforcez = 0.d0
    ux      = 0.d0
    uy      = 0.d0
    uz      = 0.d0

    dt = 1.d99
    umax = 0.d0

    ! loop over the data
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             ! compute the sound speed from rho and temp
             eos_state%rho = scal(i,j,k,rho_comp)
             eos_state%T = scal(i,j,k,temp_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, pt_index)

             spdx    = max(spdx,eos_state%cs)
             ux      = max(ux,abs(u(i,j,k,1)))
             pforcex = max(pforcex,abs(force(i,j,k,1)))
             spdy    = max(spdy,eos_state%cs)
             uy      = max(uy,abs(u(i,j,k,2)))
             pforcey = max(pforcey,abs(force(i,j,k,2)))
#if (AMREX_SPACEDIM == 3)
             spdz    = max(spdz,eos_state%cs)
             uz      = max(uz,abs(u(i,j,k,3)))
             pforcez = max(pforcez,abs(force(i,j,k,3)))
#endif

          enddo
       enddo
    enddo

    umax = max(umax,ux,uy,uz)

    ux = ux / dx(1)
    spdx = spdx / dx(1)
    uy = uy / dx(2)
    spdy = spdy / dx(2)
#if (AMREX_SPACEDIM == 3)
    uz = uz / dx(3)
    spdz = spdz / dx(3)
#endif

    ! use advective constraint unless velocities are zero everywhere
    ! in which case we use the sound speed
    if (ux .ne. 0.d0 .or. uy .ne. 0.d0 .or. uz .ne. 0.d0) then
       dt = cfl / max(ux,uy,uz)
    else if (spdx .ne. 0.d0 .or. spdy .ne. 0.d0 .or. spdz .ne. 0.d0) then
       dt = cfl / max(spdx,spdy,spdz)
    end if

    ! sound speed constraint
    if (use_soundspeed_firstdt) then
       if (spdx .eq. 0.d0 .and. spdy .eq. 0.d0 .and. spdz .eq. 0.d0) then
          dt_sound = 1.d99
       else
          dt_sound = cfl / max(spdx,spdy,spdz)
       end if
       dt = min(dt,dt_sound)
    end if

    ! force constraints
    if (pforcex > eps) dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
    if (pforcey > eps) dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
#if (AMREX_SPACEDIM == 3)
    if (pforcez > eps) dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))
#endif

    ! divU constraint
    if (use_divu_firstdt) then

       dt_divu = 1.d99

#if (AMREX_SPACEDIM == 2)

       do j = lo(2), hi(2)
          if (j .eq. 0) then
             gradp0 = (p0(lev,j+1) - p0(lev,j))/dx(2)
          else if (j .eq. nr(lev)-1) then
             gradp0 = (p0(lev,j) - p0(lev,j-1))/dx(2)
          else
             gradp0 = 0.5d0*(p0(lev,j+1) - p0(lev,j-1))/dx(2)
          endif

          do i = lo(1), hi(1)
             denom = divU(i,j,k) - u(i,j,k,2)*gradp0/(gamma1bar(lev,j)*p0(lev,j))
             if (denom > 0.d0) then
                dt_divu = min(dt_divu,0.4d0*(1.d0 - rho_min/scal(i,j,k,rho_comp))/denom)
             endif
          enddo
       enddo


#elif (AMREX_SPACEDIM == 3)

       do k = lo(3), hi(3)
          if (k .eq. 0) then
             gradp0 = (p0(lev,k+1) - p0(lev,k))/dx(3)
          else if (k .eq. nr(lev)-1) then
             gradp0 = (p0(lev,k) - p0(lev,k-1))/dx(3)
          else
             gradp0 = 0.5d0*(p0(lev,k+1) - p0(lev,k-1))/dx(3)
          endif

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                denom = divU(i,j,k) - u(i,j,k,3)*gradp0/(gamma1bar(lev,k)*p0(lev,k))
                if (denom > 0.d0) then
                   dt_divu = min(dt_divu,0.4d0*(1.d0 - rho_min/scal(i,j,k,rho_comp))/denom)
                endif
             enddo
          enddo
       enddo

#endif

       dt = min(dt,dt_divu)

    end if

  end subroutine firstdt

  subroutine firstdt_sphr(dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       gp0_cart, g_lo, g_hi) bind (C,name="firstdt_sphr")

    use amrex_fort_module, only: amrex_min, amrex_max

    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    integer         , intent(in   ) :: g_lo(3), g_hi(3)
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u    (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: gp0_cart(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)

    ! local variables
    double precision :: spdx,spdy,spdz,pforcex,pforcey,pforcez,ux,uy,uz
    double precision :: gp_dot_u,eps,dt_divu,dt_sound,denom,rho_min
    integer          :: i,j,k,r

    integer pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    eps = 1.0d-8

    rho_min = 1.d-20

    spdx    = 0.d0
    spdy    = 0.d0
    spdz    = 0.d0
    pforcex = 0.d0
    pforcey = 0.d0
    pforcez = 0.d0
    ux      = 0.d0
    uy      = 0.d0
    uz      = 0.d0

    dt = 1.d99
    umax = 0.d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! compute the sound speed from rho and temp
             eos_state%rho   = scal(i,j,k,rho_comp)
             eos_state%T     = scal(i,j,k,temp_comp)
             eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, eos_state, pt_index)

             spdx    = max(spdx,eos_state%cs)
             spdy    = max(spdy,eos_state%cs)
             spdz    = max(spdz,eos_state%cs)
             pforcex = max(pforcex,abs(force(i,j,k,1)))
             pforcey = max(pforcey,abs(force(i,j,k,2)))
             pforcez = max(pforcez,abs(force(i,j,k,3)))
             ux      = max(ux,abs(u(i,j,k,1)))
             uy      = max(uy,abs(u(i,j,k,2)))
             uz      = max(uz,abs(u(i,j,k,3)))

          enddo
       enddo
    enddo

    umax = max(umax,ux,uy,uz)

    ux = ux / dx(1)
    uy = uy / dx(2)
    uz = uz / dx(3)

    spdx = spdx / dx(1)
    spdy = spdy / dx(2)
    spdz = spdz / dx(3)

    ! advective constraint
    if (ux .ne. 0.d0 .or. uy .ne. 0.d0 .or. uz .ne. 0.d0) then
       dt = cfl / max(ux,uy,uz)
    else if (spdx .ne. 0.d0 .and. spdy .ne. 0.d0 .and. spdz .ne. 0.d0) then
       dt = cfl / max(spdx,spdy,spdz)
    end if

    ! sound speed constraint
    if (use_soundspeed_firstdt) then
       if (spdx .eq. 0.d0 .and. spdy .eq. 0.d0 .and. spdz .eq. 0.d0) then
          dt_sound = 1.d99
       else
          dt_sound = cfl / max(spdx,spdy,spdz)
       end if
       dt = min(dt,dt_sound)
    end if

    ! force constraints
    if (pforcex > eps) dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
    if (pforcey > eps) dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
    if (pforcez > eps) dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))

    ! divU constraint
    if (use_divu_firstdt) then

       dt_divu = 1.d99

       !REDUCTION(MIN : dt_divu)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                gp_dot_u = u(i,j,k,1) * gp0_cart(i,j,k,1) + &
                     u(i,j,k,2) * gp0_cart(i,j,k,2) + &
                     u(i,j,k,3) * gp0_cart(i,j,k,3)

                denom = divU(i,j,k) - gp_dot_u

                if (denom > 0.d0) then
                   dt_divu = min(dt_divu,0.4d0*(1.d0 - rho_min/scal(i,j,k,rho_comp))/denom)
                endif

             enddo
          enddo
       enddo

       dt = min(dt,dt_divu)

    end if

  end subroutine firstdt_sphr

  subroutine estdt_divu(gp0, &
       p0, gamma1bar, &
       r_cc_loc, r_edge_loc) bind (C,name="estdt_divu")

    double precision, intent(inout) :: gp0      (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)

    ! local variables
    double precision :: gamma1bar_p_avg
    integer          :: r


    !$gpu

    ! spherical divU constraint
    if (use_exact_base_state) then
       do r=1,nr_fine-1
          gamma1bar_p_avg = 0.5d0 * (gamma1bar(0,r)*p0(0,r) + gamma1bar(0,r-1)*p0(0,r-1))
          gp0(0,r) = ( (p0(0,r) - p0(0,r-1))/(r_cc_loc(0,r) - r_cc_loc(0,r-1)) ) / gamma1bar_p_avg
       end do
    else
       do r=1,nr_fine-1
          gamma1bar_p_avg = 0.5d0 * (gamma1bar(0,r)*p0(0,r) + gamma1bar(0,r-1)*p0(0,r-1))
          gp0(0,r) = ( (p0(0,r) - p0(0,r-1))/dr(0) ) / gamma1bar_p_avg
       end do
    end if

    gp0(0,nr_fine) = gp0(0,nr_fine-1)
    gp0(0,      0) = gp0(0,        1)

  end subroutine estdt_divu

end module compute_dt_module
