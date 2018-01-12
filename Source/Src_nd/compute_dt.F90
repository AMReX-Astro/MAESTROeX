module compute_dt_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, &
                                cfl, use_soundspeed_firstdt, use_divu_firstdt
  use base_state_geometry_module, only:  max_radial_level, nr_fine, nr

  implicit none

  private

contains

  subroutine firstdt(lev, dt, umax, lo, hi, dx, &
                     scal,  s_lo, s_hi, nc_s, &
                     u,     u_lo, u_hi, nc_u, &
                     force, f_lo, f_hi, nc_f, &
                     divu,  d_lo, d_hi, &
                     p0, gamma1bar) bind (C,name="firstdt")
    
    integer         , intent(in   ) :: lev
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
#if (AMREX_SPACEDIM >= 2)
       spdy    = max(spdy,eos_state%cs)
       uy      = max(uy,abs(u(i,j,k,2)))
       pforcey = max(pforcey,abs(force(i,j,k,2)))
#if (AMREX_SPACEDIM == 3)
       spdz    = max(spdz,eos_state%cs)
       uz      = max(uz,abs(u(i,j,k,3)))
       pforcez = max(pforcez,abs(force(i,j,k,3)))
#endif
#endif

    enddo
    enddo
    enddo

    umax = max(umax,ux,uy,uz)

    ux = ux / dx(1)
    spdx = spdx / dx(1)
#if (AMREX_SPACEDIM >= 2)
    uy = uy / dx(2)
    spdy = spdy / dx(2)
#if (AMREX_SPACEDIM == 3)
    uz = uz / dx(3)    
    spdz = spdz / dx(3)
#endif
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
                denom = divU(i,j,k) - u(i,j,k,3)*gradp0/(gamma1bar(k)*p0(lev,k))
                if (denom > 0.d0) then
                   dt_divu = min(dt_divu,0.4d0*(1.d0 - rho_min/s(i,j,k,rho_comp))/denom)
                endif
             enddo
          enddo
       enddo

#endif

       dt = min(dt,dt_divu)

    end if

  end subroutine firstdt

end module compute_dt_module
