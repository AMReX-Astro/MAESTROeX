module advect_base_module
  ! Create the psi term, where psi = D_0 p_0/Dt

  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, dr, &
       max_radial_level, numdisjointchunks, &
       r_start_coord, r_end_coord, &
       restrict_base, fill_ghost_base
  use meth_params_module, only: spherical, enthalpy_pred_type
  use make_edge_state_module

  implicit none

  private

contains

  subroutine advect_base_dens(w0,rho0_old,rho0_new, &
       rho0_predicted_edge,dt, &
       r_cc_loc, r_edge_loc) bind(C, name="advect_base_dens")

    double precision, intent(in   ) ::                  w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::            rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::            rho0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine  )
    double precision, value, intent(in   ) ::                  dt
    double precision, intent(in   ) ::            r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::          r_edge_loc(0:max_radial_level,0:nr_fine  )

    call bl_proffortfuncstart("Maestro::advect_base_dens")

    if (spherical .eq. 0) then
       call advect_base_dens_planar(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
       call restrict_base(rho0_new,1)
       call fill_ghost_base(rho0_new,1)
    else
       call advect_base_dens_spherical(w0,rho0_old,rho0_new,rho0_predicted_edge,dt, &
            r_cc_loc,r_edge_loc)
    end if

    call bl_proffortfuncstop("Maestro::advect_base_dens")

  end subroutine advect_base_dens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_dens_planar(w0,rho0_old,rho0_new, &
       rho0_predicted_edge,dt)

    double precision, intent(in   ) ::                  w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::            rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::            rho0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine  )
    double precision, value, intent(in   ) ::                  dt

    ! Local variables
    integer :: r, n, i

    double precision :: force(0:max_radial_level,0:nr_fine-1)
    double precision ::  edge(0:max_radial_level,0:nr_fine)

    rho0_predicted_edge = 0.

    ! zero the new density so we don't leave a non-zero density in fine radial
    ! regions that no longer have a corresponding full state
    rho0_new = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=0,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             force(n,r) = -rho0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dr(n)
          end do
       end do
    end do

    call make_edge_state_1d(rho0_old,edge,w0,force,dt)

    rho0_predicted_edge = edge

    ! update rho_0
    do n=0,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             rho0_new(n,r) = rho0_old(n,r) &
                  - dt / dr(n) * (edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r))
          end do
       end do
    end do

  end subroutine advect_base_dens_planar

  subroutine advect_base_dens_spherical(w0,rho0_old,rho0_new, &
       rho0_predicted_edge,dt, &
       r_cc_loc, r_edge_loc)

    double precision, intent(in   ) ::        w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::            rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::            rho0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine  )
    double precision, value, intent(in   ) ::                  dt
    double precision, intent(in   ) ::            r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::          r_edge_loc(0:max_radial_level,0:nr_fine  )

    ! local variables
    integer          :: r
    double precision :: dtdr
    double precision :: force(0:max_radial_level,0:nr_fine-1)
    double precision :: edge(0:max_radial_level,0:nr_fine)

    dtdr = dt / dr(0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Predict rho_0 to vertical edges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       force(0,r) = -rho0_old(0,r) * (w0(0,r+1) - w0(0,r)) / dr(0) - &
            2.0*rho0_old(0,r)*0.5*(w0(0,r) + w0(0,r+1))/r_cc_loc(0,r)
    end do
    !$OMP END PARALLEL DO

    call make_edge_state_1d(rho0_old,edge,w0,force,dt)

    rho0_predicted_edge = edge

    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       rho0_new(0,r) = rho0_old(0,r) - dtdr/r_cc_loc(0,r)**2 * &
            (r_edge_loc(0,r+1)**2 * edge(0,r+1) * w0(0,r+1) - &
            r_edge_loc(0,r  )**2 * edge(0,r  ) * w0(0,r  ))
    end do
    !$OMP END PARALLEL DO

  end subroutine advect_base_dens_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new,rho0_predicted_edge, &
       psi,dt,r_cc_loc,r_edge_loc) bind(C, name="advect_base_enthalpy")

    double precision, intent(in   ) ::                  w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::            rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           rhoh0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::           rhoh0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::                 psi(0:max_radial_level,0:nr_fine-1)
    double precision, value, intent(in   ) ::                  dt
    double precision, intent(in   ) ::            r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::          r_edge_loc(0:max_radial_level,0:nr_fine  )

    call bl_proffortfuncstart("Maestro::advect_base_enthalpy")

    if (spherical .eq. 0) then
       call advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
            rho0_predicted_edge,psi,dt)
       call restrict_base(rhoh0_new,1)
       call fill_ghost_base(rhoh0_new,1)
    else
       call advect_base_enthalpy_spherical(w0,rho0_old,rhoh0_old,rhoh0_new, &
            rho0_predicted_edge,psi,dt, &
            r_cc_loc, r_edge_loc)
    end if

    call bl_proffortfuncstop("Maestro::advect_base_enthalpy")

  end subroutine advect_base_enthalpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advect_base_enthalpy_planar(w0,rho0_old,rhoh0_old,rhoh0_new, &
       rho0_predicted_edge,psi,dt)

    double precision, intent(in   ) ::                  w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::            rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           rhoh0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::           rhoh0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::                 psi(0:max_radial_level,0:nr_fine-1)
    double precision, value, intent(in   ) ::                  dt

    ! Local variables
    integer :: r, i, n

    double precision :: force(0:max_radial_level,0:nr_fine-1)
    double precision ::  edge(0:max_radial_level,0:nr_fine)

    ! zero the new enthalpy so we don't leave a non-zero enthalpy in fine radial
    ! regions that no longer have a corresponding full state
    rhoh0_new = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Update (rho h)_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! here we predict (rho h)_0 on the edges
    do n=0,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             force(n,r) = -rhoh0_old(n,r) * (w0(n,r+1) - w0(n,r)) / dr(n) + psi(n,r)
          end do
       end do
    end do

    call make_edge_state_1d(rhoh0_old,edge,w0,force,dt)


    ! update (rho h)_0
    do n=0,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             rhoh0_new(n,r) = rhoh0_old(n,r) &
                  - dt/dr(n) * (edge(n,r+1) * w0(n,r+1) - edge(n,r) * w0(n,r)) + dt*psi(n,r)
          end do
       end do
    end do

  end subroutine advect_base_enthalpy_planar

  subroutine advect_base_enthalpy_spherical(w0,rho0_old,rhoh0_old,rhoh0_new, &
       rho0_predicted_edge,psi,dt, &
       r_cc_loc, r_edge_loc)

    double precision, intent(in   ) ::                  w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::            rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           rhoh0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::           rhoh0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::                 psi(0:max_radial_level,0:nr_fine-1)
    double precision, value, intent(in   ) ::                  dt
    double precision, intent(in   ) ::            r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::          r_edge_loc(0:max_radial_level,0:nr_fine  )

    ! Local variables
    integer :: r

    double precision :: dtdr
    double precision :: div_w0_cart

    double precision :: force(0:max_radial_level,0:nr_fine-1)
    double precision ::  edge(0:max_radial_level,0:nr_fine)

    dtdr = dt / dr(0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! UPDATE RHOH0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! here we predict (rho h)_0 on the edges
    !$OMP PARALLEL DO PRIVATE(r,div_w0_cart)
    do r=0,nr_fine-1

       div_w0_cart = (w0(0,r+1) - w0(0,r)) / dr(0)

       ! add psi at time-level n to the force for the prediction
       force(0,r) = -rhoh0_old(0,r) * div_w0_cart - &
            2.0*rhoh0_old(0,r)*0.5*(w0(0,r) + w0(0,r+1))/r_cc_loc(0,r) + psi(0,r)

    end do
    !$OMP END PARALLEL DO

    call make_edge_state_1d(rhoh0_old,edge,w0,force,dt)


    ! update (rho h)_0
    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       rhoh0_new(0,r) = rhoh0_old(0,r) - dtdr / r_cc_loc(0,r)**2 * &
            (r_edge_loc(0,r+1)**2 * edge(0,r+1) * w0(0,r+1) - &
            r_edge_loc(0,r  )**2 * edge(0,r  ) * w0(0,r  )) + dt * psi(0,r)
    end do
    !$OMP END PARALLEL DO

  end subroutine advect_base_enthalpy_spherical

end module advect_base_module
