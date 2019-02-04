module compute_dt_module

  use network, only: nspec
  use meth_params_module, only: cfl, prob_lo, prob_hi, drdxfac
  use probin_module, only: initial_dt, prob_type
  use base_state_geometry_module, only:  max_radial_level, nr_fine, nr, dr
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

  private

  double precision, parameter :: SMALL = 1.d-12

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

    double precision :: dr_fine

    dr_fine = (prob_hi(1) - prob_lo(1))/nr_fine

    dt = min(1.1d0 * dt, cfl * dr_fine / (maxval(abs(w0(1,0:nr_fine-2))) + SMALL))

  end subroutine estdt

  subroutine estdt_sphr(dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       dSdt,  t_lo, t_hi, &
       w0, &
       w0macx, x_lo, x_hi, &
       w0macy, y_lo, y_hi, &
       w0macz, z_lo, z_hi, &
       p0, gamma1bar, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind (C,name="estdt_sphr")

    double precision, intent(inout) :: dt, umax
    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    integer         , intent(in   ) :: u_lo(3), u_hi(3), nc_u
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) :: scal  (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(in   ) :: u     (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc_u)
    double precision, intent(in   ) :: force (f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(in   ) :: divu  (d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    double precision, intent(in   ) :: dSdt  (t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent(in   ) :: w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    double precision, intent(in   ) :: w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    double precision, intent(in   ) :: w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    double precision, intent(in   ) :: w0       (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    double precision :: dr_fine

    if (prob_type == 1 .or. prob_type == 3 &
         .or. prob_type == 4) then
       dr_fine = (prob_hi(1) - prob_lo(1))/nr_fine
    else
       ! need to compute it this way to agree with how the initial model was
       !  computed
       dr_fine = prob_hi(1) * dx(1) / drdxfac
    endif

    dt = min(1.1d0 * dt, cfl * dr_fine / (maxval(abs(w0(1,0:nr_fine-2))) + SMALL))

  end subroutine estdt_sphr

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

    dt = initial_dt

  end subroutine firstdt

  subroutine firstdt_sphr(dt, umax, lo, hi, dx, &
       scal,  s_lo, s_hi, nc_s, &
       u,     u_lo, u_hi, nc_u, &
       force, f_lo, f_hi, nc_f, &
       divu,  d_lo, d_hi, &
       p0, gamma1bar, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind (C,name="firstdt_sphr")

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
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    dt = initial_dt

  end subroutine firstdt_sphr

end module compute_dt_module
