
module make_pi_cc_module

  implicit none

  private

contains

 subroutine make_pi_cc(lo, hi, &
                       pi,    p_lo, p_hi, &
                       pi_cc, c_lo, c_hi) bind(C, name="make_pi_cc")

   integer         , intent(in   ) :: lo(3), hi(3)
   integer         , intent(in   ) :: p_lo(3), p_hi(3)
   integer         , intent(in   ) :: c_lo(3), c_hi(3)
   double precision, intent(in   ) :: pi   (p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
   double precision, intent(inout) :: pi_cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

   integer :: i,j,k

   do k = lo(3), hi(3)
   do j = lo(2), hi(2)
   do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 1)
      pi_cc(i,j,k) = (pi(i,j,k) + pi(i+1,j,k)) / 2.d0
#elif (AMREX_SPACEDIM == 2)
      pi_cc(i,j,k) = (pi(i,j,k) + pi(i+1,j,k) + pi(i,j+1,k) + pi(i+1,j+1,k)) / 4.d0

#elif (AMREX_SPACEDIM == 3)
      pi_cc(i,j,k) = (pi(i,j,k) + pi(i+1,j,k) + pi(i,j+1,k) + pi(i,j,k+1) &
                      + pi(i+1,j+1,k) + pi(i+1,j,k+1) + pi(i,j+1,k+1) &
                      + pi(i+1,j+1,k+1)) / 8.d0
#endif

   end do
   end do
   end do

 end subroutine make_pi_cc

end module make_pi_cc_module
