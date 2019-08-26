
module compute_grad_phi_module

  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr

  implicit none

  private

contains

  subroutine compute_grad_phi(lo, hi, &
       phi,  p_lo, p_hi, &
       gphi, g_lo, g_hi, &
       dx) bind(C, name="compute_grad_phi")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    integer         , intent(in   ) :: g_lo(3), g_hi(3)
    double precision, intent(in   ) :: phi (p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision, intent(inout) :: gphi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),AMREX_SPACEDIM)
    double precision, intent(in   ) :: dx(3)

    integer :: i,j,k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
             gphi(i,j,k,1) = 0.5d0*(phi(i+1,j,k) + phi(i+1,j+1,k) - &
                  phi(i  ,j,k) - phi(i  ,j+1,k) ) /dx(1)
             gphi(i,j,k,2) = 0.5d0*(phi(i,j+1,k) + phi(i+1,j+1,k) - &
                  phi(i,j  ,k) - phi(i+1,j  ,k) ) /dx(2)
#elif (AMREX_SPACEDIM == 3)
             gphi(i,j,k,1) = 0.25d0*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                  +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) &
                  -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                  -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) /dx(1)
             gphi(i,j,k,2) = 0.25d0*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                  +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) &
                  -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                  -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) /dx(2)
             gphi(i,j,k,3) = 0.25d0*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                  +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) &
                  -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                  -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) /dx(3)
#endif

          end do
       end do
    end do

  end subroutine compute_grad_phi

  subroutine compute_grad_phi_rad(phi, gphi_rad) &
       bind(C, name="compute_grad_phi_rad")

    double precision, intent(in   ) ::      phi(0:max_radial_level,0:nr_fine)
    double precision, intent(inout) :: gphi_rad(0:max_radial_level,0:nr_fine-1)

    ! local
    integer :: r

    do r=0,nr_fine-1
       gphi_rad(0,r) = (phi(0,r+1) - phi(0,r)) / dr(0)
    end do

  end subroutine compute_grad_phi_rad

end module compute_grad_phi_module
