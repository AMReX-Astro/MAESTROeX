! Compute the timestep
module estdt_module

  use base_state_geometry_module, only:  max_radial_level, nr_fine, nr, dr
  use probin_module, only: diffusion_coefficient

  implicit none

  private

contains

  subroutine estdt(lev, dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       dSdt,  t_lo, t_hi, &
       w0, p0, gamma1bar) bind (C,name="estdt")

    integer         , intent(in   ) :: lev
    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u    (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: dSdt (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent(in   ) :: w0       (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)

    ! calculate the timestep
    dt = minval(dx*dx / diffusion_coefficient)

  end subroutine estdt

  subroutine firstdt(lev, dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       p0, gamma1bar) bind (C,name="firstdt")

    integer         , intent(in   ) :: lev
    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u    (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)

    ! local variables
    dt = minval(dx*dx / diffusion_coefficient)

  end subroutine firstdt

end module estdt_module
