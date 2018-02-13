module update_vel_module

  use amrex_constants_module
  use base_state_geometry_module, only:  max_radial_level, nr_fine
  use meth_params_module, only: do_sponge

  implicit none

  private

contains

  subroutine update_velocity(lev, lo, hi, &
                               uold, uo_lo, uo_hi, &
                               unew, un_lo, un_hi, &
                               umac,  u_lo, u_hi, &
#if (AMREX_SPACEDIM >= 2)
                               vmac,  v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                               wmac,  w_lo, w_hi, &
#endif
#endif
                               uedgex, x_lo, x_hi, &
#if (AMREX_SPACEDIM >= 2)
                               uedgey, y_lo, y_hi, &
#if (AMREX_SPACEDIM == 3)
                               uedgez, z_lo, z_hi, &
#endif
#endif
                               force,  f_lo, f_hi, &
                               sponge, s_lo, s_hi, & 
                               w0, dx, dt) &
                               bind(C,name="update_velocity")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: uo_lo(3), uo_hi(3)
    double precision, intent(in   ) :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: un_lo(3), un_hi(3)
    double precision, intent(inout) :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
#if (AMREX_SPACEDIM >= 2)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
#endif
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) :: uedgex (x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),AMREX_SPACEDIM)
#if (AMREX_SPACEDIM >= 2)
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) :: uedgey (y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) :: uedgez (z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),AMREX_SPACEDIM)
#endif
#endif
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),AMREX_SPACEDIM)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    double precision, intent(in   ) :: sponge (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: w0     (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(AMREX_SPACEDIM), dt
    
    integer :: i,j,k, dim
    double precision :: ubar, vbar, wbar, w0bar
    double precision :: ugradu, ugradv, ugradw
    
    ! 1) Subtract (Utilde dot grad) Utilde term from old Utilde
    ! 2) Add forcing term to new Utilde
    !$OMP PARALLEL DO PRIVATE(i,j,k,ubar,vbar,wbar,ugradu,ugradv,ugradw)
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

       ! create cell-centered Utilde
       dim = i
       ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
#if (AMREX_SPACEDIM >= 2)
       dim = j
       vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
#if (AMREX_SPACEDIM == 3)
       dim = k
       wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))
#endif
#endif

       ! create (Utilde dot grad) Utilde
       ugradu = ( ubar*(uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1) &
#if (AMREX_SPACEDIM >= 2)
                + vbar*(uedgey(i,j+1,k,1) - uedgey(i,j,k,1))/dx(2) &
#if (AMREX_SPACEDIM == 3)
                + wbar*(uedgez(i,j,k+1,1) - uedgez(i,j,k,1))/dx(3) &
#endif
#endif
                  )

#if (AMREX_SPACEDIM >= 2)
       ugradv = ( ubar*(uedgex(i+1,j,k,2) - uedgex(i,j,k,2))/dx(1) &
                + vbar*(uedgey(i,j+1,k,2) - uedgey(i,j,k,2))/dx(2) &
#if (AMREX_SPACEDIM == 3)
                + wbar*(uedgez(i,j,k+1,2) - uedgez(i,j,k,2))/dx(3) &
#endif
#endif
                  ) 

#if (AMREX_SPACEDIM == 3)
       ugradw = ubar*(uedgex(i+1,j,k,3) - uedgex(i,j,k,3))/dx(1) &
                + vbar*(uedgey(i,j+1,k,3) - uedgey(i,j,k,3))/dx(2) &
                + wbar*(uedgez(i,j,k+1,3) - uedgez(i,j,k,3))/dx(3)
#endif 
                  
       ! update with (Utilde dot grad) Utilde and force
       unew(i,j,k,1) = uold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
#if (AMREX_SPACEDIM >= 2)
       unew(i,j,k,2) = uold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
#if (AMREX_SPACEDIM == 3)
       unew(i,j,k,3) = uold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)
#endif
#endif

       ! subtract (w0 dot grad) Utilde term
       w0bar = HALF*(w0(lev,dim) + w0(lev,dim+1))
       unew(i,j,k,:) = unew(i,j,k,:) - dt * w0bar * &
#if (AMREX_SPACEDIM == 1)
                          (uedgex(i+1,j,k,1) - uedgex(i,j,k,1))/dx(1)
#elif (AMREX_SPACEDIM == 2) 
                          (uedgey(i,j+1,k,:) - uedgey(i,j,k,:))/dx(2)
#elif (AMREX_SPACEDIM == 3)
                          (uedgez(i,j,k+1,:) - uedgez(i,j,k,:))/dx(3)
#endif

       ! Add the sponge
       if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

    end do
    end do
    end do
    !$OMP END PARALLEL DO

  end subroutine update_velocity

end module update_vel_module
