module make_psi_module
  ! Create the psi term, where psi = D_0 p_0/Dt

  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, dr, &
                                        max_radial_level, numdisjointchunks, &
                                        r_start_coord, r_end_coord, finest_radial_level, &
                                        restrict_base, fill_ghost_base, &
                                        base_cutoff_density_coord, anelastic_cutoff_density_coord
  use meth_params_module, only: grav_const

  implicit none

  private

contains

  subroutine make_psi_planar(etarho_cc,psi) bind(C, name="make_psi_planar")
    ! Binds to C function ``make_psi_planar``

    double precision, intent(in   ) :: etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::       psi(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer         :: r,i,n

    psi = ZERO

    do n=0,finest_radial_level
       do i=1,numdisjointchunks(n)
          do r = r_start_coord(n,i), r_end_coord(n,i)
             if (r .lt. base_cutoff_density_coord(n)) then
                psi(n,r) = etarho_cc(n,r) * abs(grav_const)
             end if
          end do
       end do
    end do

    call restrict_base(psi,1)
    call fill_ghost_base(psi,1)

  end subroutine make_psi_planar

  subroutine make_psi_spherical(psi,w0,gamma1bar,p0_avg,Sbar_in,r_cc_loc,r_edge_loc) bind(C, name="make_psi_spherical")
    ! Binds to C function ``make_psi_spherical``

    double precision, intent(inout) ::       psi(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::        w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    p0_avg(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   Sbar_in(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::  r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine )

    ! local variables
    integer :: r
    double precision :: div_w0_sph

    call bl_proffortfuncstart("Maestro::make_psi_spherical")

    psi = ZERO

    !$OMP PARALLEL DO PRIVATE(r,div_w0_sph)
    do r=0,base_cutoff_density_coord(0)-1

       div_w0_sph = one/(r_cc_loc(0,r)**2)* &
            (r_edge_loc(0,r+1)**2 * w0(0,r+1) - &
             r_edge_loc(0,r  )**2 * w0(0,r  )) / dr(0)

       psi(0,r) = gamma1bar(0,r) * p0_avg(0,r) * (Sbar_in(0,r) - div_w0_sph)

    enddo
    !$OMP END PARALLEL DO

    call bl_proffortfuncstop("Maestro::make_psi_spherical")

  end subroutine make_psi_spherical

  subroutine make_psi_irreg(etarho_cc,grav_cell,psi) bind(C, name="make_psi_irreg")
    ! Binds to C function ``make_psi_irreg``

    double precision, intent(in   ) :: etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: grav_cell(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::       psi(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer         :: r

    call bl_proffortfuncstart("Maestro::make_psi_irreg")

    psi = ZERO

    do r=0,base_cutoff_density_coord(0)
       psi(0,r) = etarho_cc(0,r) * grav_cell(0,r)
    end do

    do r=base_cutoff_density_coord(0)+1,nr_fine-1
       psi(0,r) = psi(0,r-1)
    enddo

    call restrict_base(psi,1)
    call fill_ghost_base(psi,1)

    call bl_proffortfuncstop("Maestro::make_psi_irreg")

  end subroutine make_psi_irreg

!! OLD PSI SUBROUTINE
!!$  subroutine make_psi_irreg(psi,p0_old,p0_new,dt) bind(C, name="make_psi_irreg")
!!$
!!$    double precision, intent(inout) ::       psi(0:max_radial_level,0:nr_fine-1)
!!$    double precision, intent(in   ) ::    p0_old(0:max_radial_level,0:nr_fine-1)
!!$    double precision, intent(in   ) ::    p0_new(0:max_radial_level,0:nr_fine-1)
!!$    double precision, intent(in   ) ::        dt
!!$
!!$    ! local variables
!!$    integer :: r
!!$
!!$    psi = ZERO
!!$
!!$    !$OMP PARALLEL DO PRIVATE(r)
!!$    do r=0,base_cutoff_density_coord(0)
!!$       psi(0,r) = (p0_new(0,r) - p0_old(0,r))/dt
!!$    enddo
!!$    !$OMP END PARALLEL DO
!!$
!!$    !$OMP PARALLEL DO PRIVATE(r)
!!$    do r=base_cutoff_density_coord(0)+1,nr_fine-1
!!$       psi(0,r) = psi(0,r-1)
!!$    enddo
!!$    !$OMP END PARALLEL DO
!!$
!!$  end subroutine make_psi_irreg

end module make_psi_module
