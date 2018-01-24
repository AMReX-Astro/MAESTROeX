
module compute_grad_phi_module

  implicit none

  private

contains

 subroutine compute_grad_phi(lo, hi, &
                             phi,  p_lo, p_hi, &
                             gphi, g_lo, g_hi, nc_g, &
                             dx) bind(C, name="compute_grad_phi")

   integer         , intent(in   ) :: lo(3), hi(3)
   integer         , intent(in   ) :: p_lo(3), p_hi(3)
   integer         , intent(in   ) :: g_lo(3), g_hi(3), nc_g
   double precision, intent(in   ) :: phi (p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
   double precision, intent(inout) :: gphi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
   double precision, intent(in   ) :: dx(3)

   integer :: i,j,k

   do k = lo(3), hi(3)
   do j = lo(2), hi(2)
   do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 1)
      gphi(i,j,k,1) = ( phi(i+1,j,k) - phi(i,j,k) ) / dx(1)
#elif (AMREX_SPACEDIM == 2)
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

end module compute_grad_phi_module
