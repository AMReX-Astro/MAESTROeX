module compute_dt_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp, &
                                cfl, use_soundspeed_firstdt

  implicit none

  private

contains

  subroutine firstdt(dt, umax, lo, hi, dx, &
                     scal, s_lo, s_hi, nc_s, &
                     u,    u_lo, u_hi, nc_u) bind (C,name="firstdt")
    
    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    double precision, intent(in   ) :: scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u   (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)

    ! local variables
    integer i,j,k
    integer pt_index(3)
    type(eos_t) :: eos_state

    double precision spdx, spdy, spdz, ux, uy, uz
    double precision dt_sound

    ! FIXME - need to add force and divu constraints

    spdx = 0.d0
    spdy = 0.d0
    spdz = 0.d0

    ux = 0.d0
    uy = 0.d0
    uz = 0.d0

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

       spdx    = max(spdx,eos_state%cs)
       ux      = max(ux,abs(u(i,j,k,1)))
#if (AMREX_SPACEDIM >= 2)
       spdy    = max(spdy,eos_state%cs)
       uy      = max(uy,abs(u(i,j,k,2)))
#if (AMREX_SPACEDIM == 3)
       spdz    = max(spdz,eos_state%cs)
       uz      = max(uz,abs(u(i,j,k,3)))
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

  end subroutine firstdt

end module compute_dt_module
