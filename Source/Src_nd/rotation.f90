! rotation
!
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
!       plane-parallel patch and it is not advisable to utilize
!       rotation for such problems.
!

module rotation_module

  use meth_params_module, only: rotational_frequency, co_latitude, rho_comp
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  double precision, save :: sin_theta, cos_theta, omega

  public :: rotation_init, sin_theta, cos_theta, omega


contains

  subroutine rotation_init() bind(C, name="rotation_init")

    implicit none

    real(rt) :: theta_in_rad

    theta_in_rad = M_PI * co_latitude / 180.d0

    sin_theta = sin(theta_in_rad)
    cos_theta = cos(theta_in_rad)

    omega = 2 * M_PI * rotational_frequency

  end subroutine rotation_init


  subroutine derangmom(L,L_lo,L_hi,&
       state,s_lo,s_hi,ncomp_s, &
       u,u_lo,u_hi,&
       lo,hi,domlo,domhi, &
       dx,xlo) bind(C, name="derangmom")

    use math_module, only: cross_product
    use base_state_geometry_module, only: center

    implicit none

    integer, intent(in), value :: ncomp_s
    integer, intent(in) :: L_lo(3), L_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),3)
    real(rt), intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    real(rt), intent(in) :: dx(3), xlo(3)

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    do k = lo(3), hi(3)
       loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

             rho = state(i,j,k,rho_comp)
             mom = u(i,j,k,:) * rho
             ang_mom = cross_product(loc, mom)

             L(i,j,k,:) = ang_mom

          enddo
       enddo
    enddo

  end subroutine derangmom


end module rotation_module
