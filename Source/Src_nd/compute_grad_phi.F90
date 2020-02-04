
module compute_grad_phi_module

  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr

  implicit none

  private

contains

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
