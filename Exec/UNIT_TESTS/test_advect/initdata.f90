
module initdata_module

  use parallel, only: parallel_IOProcessor
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       prob_lo, prob_hi, base_cutoff_density

  implicit none

  private

  ! width of the gaussian
  double precision, parameter :: W = 0.05

contains

  subroutine initdata(lev, time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, &
       dx) bind(C, name="initdata")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), 1:nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), 1:nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k

    double precision :: x, y, z, dist

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.0d0
    scal(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_s) = 0.0d0

    do k=lo(3),hi(3)
       z = (dble(k)+0.5d0)*dx(3) + prob_lo(3) - center(3)
       do j=lo(2),hi(2)
          y = (dble(j)+0.5d0)*dx(2) + prob_lo(2) - center(2)
          do i=lo(1),hi(1)
             x = (dble(i)+0.5d0)*dx(1) + prob_lo(1) - center(1)

             dist = sqrt(x**2 + y**2 + z**2)

             scal(i,j,k,rho_comp)  = max(exp(-dist**2/W**2),base_cutoff_density)

             scal(i,j,k,spec_comp) = scal(i,j,k,rho_comp)

          end do
       end do
    end do

  end subroutine initdata

end module initdata_module
