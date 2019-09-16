module make_grav_module

  use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
  use probin_module, only : g0
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level, finest_radial_level

  implicit none

  private

  public :: make_grav_cell, make_grav_edge, grav_zone

contains

  pure function grav_zone(y) result (g)

    double precision, intent(in) :: y
    double precision :: g
    double precision :: fg

    if (y < 1.0625d0 * 4.d8) then 
        fg = HALF * (ONE + sin(16.d0 * M_PI * (y/4.d8 - 1.03125d0)))
    else if (y > 2.9375d0 * 4.d8) then
        fg = HALF * (ONE - sin(16.d0 * M_PI * (y/4.d8 - 2.96875d0)))
    else
        fg = ONE
    endif

    g = fg * g0 / ((y / 4.d8)**1.25d0)

  end function grav_zone

  subroutine make_grav_cell(grav_cell,rho0,r_cc_loc,r_edge_loc) &
       bind(C, name="make_grav_cell")
    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell
    ! does too.

    double precision, intent(  out) ::  grav_cell(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::       rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

    ! Local variables
    integer :: r, n

    call bl_proffortfuncstart("Maestro::make_grav_cell")

    do n = 0, finest_radial_level
        do r = 0, nr(n)-1
           grav_cell(n,r) = grav_zone(r_cc_loc(n,r))
        enddo
     enddo

    call bl_proffortfuncstop("Maestro::make_grav_cell")

  end subroutine make_grav_cell

  subroutine make_grav_edge(grav_edge,rho0,r_edge_loc) &
       bind(C, name="make_grav_edge")
    ! compute the base state gravity at the cell edges
    ! grav_edge(0) is the gravitational acceleration at the left edge of zone 0).
    ! The base state uses 0-based indexing, so grav_edge does too.

    double precision, intent(  out) ::  grav_edge(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::       rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

    ! Local variables
    integer                      :: r, n

    call bl_proffortfuncstart("Maestro::make_grav_edge")

    do n = 0, finest_radial_level
        do r = 0, nr(n)-1
           grav_edge(n,r) = grav_zone(r_edge_loc(n,r))
        enddo
     enddo

    call bl_proffortfuncstop("Maestro::make_grav_edge")

  end subroutine make_grav_edge

end module make_grav_module
