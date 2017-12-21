module make_vel_force_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp
  use base_state_geometry_module, only:  max_radial_level, nr_fine

  implicit none

  private

contains

  subroutine make_vel_force(lev, lo, hi, &
                            vel_force, f_lo, f_hi, nc_f, &
                            gpi, g_lo, g_hi, nc_g, &
                            rho, r_lo, r_hi, &
                            uedge, u_lo, u_hi, &
                            vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                            wedge, w_lo, w_hi, &
#endif
                            w0,rho0,grav, &
                            is_final_update, do_add_utilde_force) &
                            bind(C, name="make_vel_force")
    
    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent (in   ) :: g_lo(3), g_hi(3), nc_g
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: u_lo(3), u_hi(3)
    integer         , intent (in   ) :: v_lo(3), v_hi(3)
#if (AMREX_SPACEDIM == 3)
    integer         , intent (in   ) :: w_lo(3), w_hi(3)
#endif
    double precision, intent (inout) :: vel_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent (in   ) ::       gpi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),nc_g)
    double precision, intent (in   ) ::       rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent (in   ) ::     uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent (in   ) ::     vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent (in   ) ::     wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent (in   ) ::   w0(0:max_radial_level,0:nr_fine)
    double precision, intent (in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: grav(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: is_final_update, do_add_utilde_force





  end subroutine make_vel_force

end module make_vel_force_module
