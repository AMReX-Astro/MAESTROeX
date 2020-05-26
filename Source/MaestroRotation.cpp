
#include <Maestro.H>

using namespace amrex;

/*
 co_latitude, rotation_radius, theta_in_rad, sin_theta and cos_theta
 are only important when wanting to rotate a plane-parallel patch
 which lives at an angle co_latitude from the rotation axis and at a
 distance rotation_radius from center().  mk_vel_force should
 calculate the rotational forcing terms for the points within the
 patch.

 for spherical problems, the only important variable from
 init_rotation() is omega, the angular frequency - mk_vel_force will
 calculate the appropriate terms for a given coordinate

 NOTE: it is currently unclear how to handle BC's with a
 plane-parallel patch and it is not advisable to utilize
 rotation for such problems.
 */

void Maestro::RotationInit() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RotationInit()", RotationInit);

    Real theta_in_rad = M_PI * co_latitude / 180.0;

    sin_theta = std::sin(theta_in_rad);
    cos_theta = std::cos(theta_in_rad);

    omega = 2.0 * M_PI * rotational_frequency;
}
