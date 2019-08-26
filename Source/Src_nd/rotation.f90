module rotation_module
  ! co_latitude, rotation_radius, theta_in_rad, sin_theta and cos_theta
  ! are only important when wanting to rotate a plane-parallel patch
  ! which lives at an angle co_latitude from the rotation axis and at a
  ! distance rotation_radius from center().  mk_vel_force should
  ! calculate the rotational forcing terms for the points within the
  ! patch.
  !
  ! for spherical problems, the only important variable from
  ! init_rotation() is omega, the angular frequency - mk_vel_force will
  ! calculate the appropriate terms for a given coordinate
  !
  ! NOTE: it is currently unclear how to handle BC's with a
  ! plane-parallel patch and it is not advisable to utilize
  ! rotation for such problems.

  use meth_params_module, only: rotational_frequency, co_latitude
  use amrex_constants_module

  implicit none

  private

  public :: init_rotation

  double precision, save :: sin_theta, cos_theta, omega


contains

  subroutine init_rotation()

    double precision :: theta_in_rad

    theta_in_rad = M_PI * co_latitude / 180.d0

    sin_theta = sin(theta_in_rad)
    cos_theta = cos(theta_in_rad)

    omega = 2 * M_PI * rotational_frequency

  end subroutine init_rotation

end module rotation_module
