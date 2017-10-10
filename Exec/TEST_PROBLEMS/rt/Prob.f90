
subroutine initdata(level, time, lo, hi, &
                    scal, scal_lo, scal_hi, &
                    vel, vel_lo, vel_hi, &
                    dx, prob_lo, nscal) bind(C, name="initdata")

  use amrex_fort_module, only : amrex_spacedim

  implicit none
  integer, intent(in) :: level, lo(3), hi(3)
  integer, intent(in) :: scal_lo(3), scal_hi(3)
  integer, intent(in) :: vel_lo(3), vel_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
       &                                  scal_lo(2):scal_hi(2), &
       &                                  scal_lo(3):scal_hi(3), 1:nscal)
  double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
       &                                 vel_lo(2):vel_hi(2), &
       &                                 vel_lo(3):vel_hi(3), 1:amrex_spacedim)
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer, intent(in) :: nscal

  integer          :: i,j,k
  double precision :: x,y,z,r2

  vel = 0.d0
  
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
           
           if ( amrex_spacedim .eq. 2) then
              r2 = ((x-0.875d0)**2 + (y-0.5d0)**2) / 0.0025d0
              scal(i,j,k,1:nscal-1) = 1.d0 + exp(-r2)
              r2 = ((x-0.5d0)**2 + (y-0.875d0)**2) / 0.0025d0
              scal(i,j,k,nscal) = 1.d0 + exp(-r2)
           else
              r2 = ((x-0.875d0)**2 + (y-0.5d0)**2 + (z-0.5d0)**2) / 0.0025d0
              scal(i,j,k,1:nscal-1) = 1.d0 + exp(-r2)
              r2 = ((x-0.5d0)**2 + (y-0.875d0)**2 + (z-0.5d0)**2) / 0.0025d0
              scal(i,j,k,nscal) = 1.d0 + exp(-r2)
           end if
        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata
