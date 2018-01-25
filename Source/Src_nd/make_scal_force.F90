module make_scal_force_module

  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr

  implicit none

  private

contains

  subroutine mkrhohforce()

  end subroutine mkrhohforce

  subroutine modify_scal_force(lev, lo, hi, &
                               force, f_lo, f_hi, &
                               scal,  s_lo, s_hi, &
                               umac,  u_lo, u_hi, &
#if (AMREX_SPACEDIM >= 2)
                               vmac,  v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                               wmac,  w_lo, w_hi, &
#endif
#endif
                               s0, s0_edge, w0, dx, do_fullform) &
                               bind(C,name="modify_scal_force")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
#if (AMREX_SPACEDIM >= 2)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
#endif
#endif
    double precision, intent(inout) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
#if (AMREX_SPACEDIM >= 2)
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
#endif
    double precision, intent(in   ) :: s0     (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: s0_edge(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0     (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: do_fullform
    
    integer :: i,j,k,r
    double precision :: divu,divs0u
    
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

       ! umac does not contain w0
#if (AMREX_SPACEDIM == 1)
       r = i
       divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1)
#elif (AMREX_SPACEDIM == 2)
       r = j
       divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
             +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2)
#elif (AMREX_SPACEDIM == 3)
       r = k
       divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
             +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
             +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
#endif

       ! add w0 contribution
       divu = divu + (w0(lev,r+1)-w0(lev,r))/dx(AMREX_SPACEDIM)

       if (do_fullform .eq. 1) then
          force(i,j,k) = force(i,j,k) - scal(i,j,k)*divu
       else

#if (AMREX_SPACEDIM == 1)
          divs0u = (umac(i+1,j,k) * s0_edge(lev,i+1) - &
                    umac(i  ,j,k) * s0_edge(lev,i  ) )/ dx(1)
#elif (AMREX_SPACEDIM == 2)
          divs0u = s0(lev,j)*(  umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                          +(vmac(i,j+1,k) * s0_edge(lev,j+1) - &
                            vmac(i,j  ,k) * s0_edge(lev,j  ) )/ dx(2)
#elif (AMREX_SPACEDIM == 3)
          divs0u = s0(lev,k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                          +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                          +(wmac(i,j,k+1) * s0_edge(lev,k+1) &
                          - wmac(i,j,k  ) * s0_edge(lev,k  ))/ dx(3)
#endif

          force(i,j,k) = force(i,j,k) - (scal(i,j,k)-s0(lev,r))*divu - divs0u 
       endif

    end do
    end do
    end do

  end subroutine modify_scal_force

end module make_scal_force_module
