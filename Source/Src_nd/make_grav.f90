module make_grav_module

  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, &
                                        max_radial_level, nr, numdisjointchunks, &
                                        r_start_coord, r_end_coord, finest_radial_level, &
                                        restrict_base, fill_ghost_base
  use meth_params_module, only: spherical, grav_const, base_cutoff_density, &
                                do_planar_invsq_grav, planar_invsq_mass, do_2d_planar_octant
  use fundamental_constants_module, only: Gconst

  implicit none

  private

  public :: make_grav_cell, make_grav_edge

contains

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
    integer :: r, n, i
    double precision, allocatable :: m(:,:)
    double precision              :: term1, term2

    call bl_proffortfuncstart("Maestro::make_grav_cell")

    if (spherical .eq. 0) then

       if (do_planar_invsq_grav)  then

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do n=0,finest_radial_level
             do r = 0, nr(n)-1
                grav_cell(n,r) = -Gconst*planar_invsq_mass / r_cc_loc(n,r)**2
             enddo
          enddo

       else if (do_2d_planar_octant .eq. 1) then

          ! compute gravity as in the spherical case

          allocate(m(0:finest_radial_level,0:nr_fine-1))

          n = 0
          m(n,0) = FOUR3RD*M_PI*rho0(n,0)*r_cc_loc(n,0)**3
          grav_cell(n,0) = -Gconst * m(n,0) / r_cc_loc(n,0)**2

          do r=1,nr(n)-1

             ! the mass is defined at the cell-centers, so to compute
             ! the mass at the current center, we need to add the
             ! contribution of the upper half of the zone below us and
             ! the lower half of the current zone.

             ! don't add any contributions from outside the star --
             ! i.e.  rho < base_cutoff_density
             if (rho0(n,r-1) > base_cutoff_density) then
                term1 = FOUR3RD*M_PI*rho0(n,r-1) * &
                     (r_edge_loc(n,r) - r_cc_loc(n,r-1)) * &
                     (r_edge_loc(n,r)**2 + &
                     r_edge_loc(n,r)*r_cc_loc(n,r-1) + &
                     r_cc_loc(n,r-1)**2)
             else
                term1 = ZERO
             endif

             if (rho0(n,r) > base_cutoff_density) then
                term2 = FOUR3RD*M_PI*rho0(n,r  )*&
                     (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
                     (r_cc_loc(n,r)**2 + &
                     r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                     r_edge_loc(n,r  )**2)
             else
                term2 = ZERO
             endif

             m(n,r) = m(n,r-1) + term1 + term2

             grav_cell(n,r) = -Gconst * m(n,r) / r_cc_loc(n,r)**2

          enddo

          do n=1,finest_radial_level
             do i=1,numdisjointchunks(n)

                if (r_start_coord(n,i) .eq. 0) then
                   m(n,0) = FOUR3RD*M_PI*rho0(n,0)*r_cc_loc(n,0)**3
                   grav_cell(n,0) = -Gconst * m(n,0) / r_cc_loc(n,0)**2
                else
                   r = r_start_coord(n,i)
                   m(n,r) = m(n-1,r/2-1)

                   ! the mass is defined at the cell-centers, so to compute
                   ! the mass at the current center, we need to add the
                   ! contribution of the upper half of the zone below us and
                   ! the lower half of the current zone.

                   ! don't add any contributions from outside the star --
                   ! i.e.  rho < base_cutoff_density
                   if (rho0(n-1,r/2-1) > base_cutoff_density) then
                      term1 = FOUR3RD*M_PI*rho0(n-1,r/2-1) * &
                           (r_edge_loc(n-1,r/2) - r_cc_loc(n-1,r/2-1)) * &
                           (r_edge_loc(n-1,r/2)**2 + &
                           r_edge_loc(n-1,r/2)*r_cc_loc(n-1,r/2-1) + &
                           r_cc_loc(n-1,r/2-1)**2)
                   else
                      term1 = ZERO
                   endif

                   if (rho0(n,r) > base_cutoff_density) then
                      term2 = FOUR3RD*M_PI*rho0(n,r  )*&
                           (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
                           (r_cc_loc(n,r)**2 + &
                           r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                           r_edge_loc(n,r  )**2)
                   else
                      term2 = ZERO
                   endif

                   m(n,r) = m(n,r) + term1 + term2

                   grav_cell(n,r) = -Gconst * m(n,r) / r_cc_loc(n,r)**2

                end if

                do r=r_start_coord(n,i)+1,r_end_coord(n,i)

                   ! the mass is defined at the cell-centers, so to compute
                   ! the mass at the current center, we need to add the
                   ! contribution of the upper half of the zone below us and
                   ! the lower half of the current zone.

                   ! don't add any contributions from outside the star --
                   ! i.e.  rho < base_cutoff_density
                   if (rho0(n,r-1) > base_cutoff_density) then
                      term1 = FOUR3RD*M_PI*rho0(n,r-1) * &
                           (r_edge_loc(n,r) - r_cc_loc(n,r-1)) * &
                           (r_edge_loc(n,r)**2 + &
                           r_edge_loc(n,r)*r_cc_loc(n,r-1) + &
                           r_cc_loc(n,r-1)**2)
                   else
                      term1 = ZERO
                   endif

                   if (rho0(n,r) > base_cutoff_density) then
                      term2 = FOUR3RD*M_PI*rho0(n,r  )*&
                           (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
                           (r_cc_loc(n,r)**2 + &
                           r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                           r_edge_loc(n,r  )**2)
                   else
                      term2 = ZERO
                   endif

                   m(n,r) = m(n,r-1) + term1 + term2

                   grav_cell(n,r) = -Gconst * m(n,r) / r_cc_loc(n,r)**2

                end do
             enddo
          end do

          call restrict_base(grav_cell,1)
          call fill_ghost_base(grav_cell,1)

       else

          ! constant gravity
          grav_cell = grav_const

       endif

    else  ! spherical = 1

       allocate(m(0:0,0:nr_fine-1))

       m(0,0) = FOUR3RD*M_PI*rho0(0,0)*r_cc_loc(0,0)**3
       grav_cell(0,0) = -Gconst * m(0,0) / r_cc_loc(0,0)**2

       do r=1,nr_fine-1

          ! the mass is defined at the cell-centers, so to compute
          ! the mass at the current center, we need to add the
          ! contribution of the upper half of the zone below us and
          ! the lower half of the current zone.

          ! don't add any contributions from outside the star --
          ! i.e.  rho < base_cutoff_density
          if (rho0(0,r-1) > base_cutoff_density) then
             term1 = FOUR3RD*M_PI*rho0(0,r-1) * &
                  (r_edge_loc(0,r) - r_cc_loc(0,r-1)) * &
                  (r_edge_loc(0,r)**2 + &
                   r_edge_loc(0,r)*r_cc_loc(0,r-1) + &
                   r_cc_loc(0,r-1)**2)
          else
             term1 = ZERO
          endif

          if (rho0(0,r) > base_cutoff_density) then
             term2 = FOUR3RD*M_PI*rho0(0,r  )*&
                  (r_cc_loc(0,r) - r_edge_loc(0,r  )) * &
                  (r_cc_loc(0,r)**2 + &
                   r_cc_loc(0,r)*r_edge_loc(0,r  ) + &
                   r_edge_loc(0,r  )**2)
          else
             term2 = ZERO
          endif

          m(0,r) = m(0,r-1) + term1 + term2

          grav_cell(0,r) = -Gconst * m(0,r) / r_cc_loc(0,r)**2

       enddo

       deallocate(m)

    end if

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
    integer                      :: r, n, i
    double precision              :: mencl
    double precision, allocatable :: m(:,:)

    call bl_proffortfuncstart("Maestro::make_grav_edge")

    if (spherical .eq. 0) then

       if (do_planar_invsq_grav)  then

          ! we are doing a plane-parallel geometry with a 1/r**2
          ! gravitational acceleration.  The mass is assumed to be
          ! at the origin.  The mass in the computational domain
          ! does not contribute to the gravitational acceleration.
          do n=0,finest_radial_level
             do r = 0, nr(n)-1
                grav_edge(n,r) = -Gconst*planar_invsq_mass / r_edge_loc(n,r)**2
             enddo
          enddo

       else if (do_2d_planar_octant .eq. 1) then

          ! compute gravity as in spherical geometry

          allocate(m(0:finest_radial_level,0:nr_fine))

          grav_edge(0,0) = zero
          m(0,0) = ZERO

          do r=1,nr(0)-1

             ! only add to the enclosed mass if the density is
             ! > base_cutoff_density
             if (rho0(0,r-1) > base_cutoff_density) then
                m(0,r) = m(0,r-1) + FOUR3RD*M_PI * &
                     (r_edge_loc(0,r) - r_edge_loc(0,r-1)) * &
                     (r_edge_loc(0,r)**2 + &
                     r_edge_loc(0,r)*r_edge_loc(0,r-1) + &
                     r_edge_loc(0,r-1)**2) * rho0(0,r-1)
             else
                m(0,r) = m(0,r-1)
             endif

             grav_edge(0,r) = -Gconst * m(0,r) / r_edge_loc(0,r)**2

          end do

          do n=1,finest_radial_level
             do i=1,numdisjointchunks(n)

                if (r_start_coord(n,i) .eq. 0) then

                   m(n,0) = ZERO

                else

                   m(n,r_start_coord(n,i)) = m(n-1,r_start_coord(n,i)/2)
                   grav_edge(n,r_start_coord(n,i)) = grav_edge(n-1,r_start_coord(n,i)/2)

                end if

                do r=r_start_coord(n,i)+1,r_end_coord(n,i)+1

                   ! only add to the enclosed mass if the density is
                   ! > base_cutoff_density
                   if (rho0(n,r-1) > base_cutoff_density) then
                      m(n,r) = m(n,r-1) + FOUR3RD*M_PI * &
                           (r_edge_loc(n,r) - r_edge_loc(n,r-1)) * &
                           (r_edge_loc(n,r)**2 + &
                           r_edge_loc(n,r)*r_edge_loc(n,r-1) + &
                           r_edge_loc(n,r-1)**2) * rho0(n,r-1)
                   else
                      m(n,r) = m(n,r-1)
                   endif

                   grav_edge(n,r) = -Gconst * m(n,r) / r_edge_loc(n,r)**2

                end do
             enddo
          end do

          deallocate(m)


          call restrict_base(grav_edge,0)
          call fill_ghost_base(grav_edge,0)


       else

          ! constant gravity
          grav_edge = grav_const

       endif

    else

       grav_edge(0,0) = zero
       mencl = ZERO

       do r=1,nr_fine

          ! only add to the enclosed mass if the density is
          ! > base_cutoff_density
          if (rho0(0,r-1) > base_cutoff_density) then
             mencl = mencl + FOUR3RD*M_PI * &
                  (r_edge_loc(0,r) - r_edge_loc(0,r-1)) * &
                  (r_edge_loc(0,r)**2 + &
                   r_edge_loc(0,r)*r_edge_loc(0,r-1) + &
                   r_edge_loc(0,r-1)**2) * rho0(0,r-1)
          endif

          grav_edge(0,r) = -Gconst * mencl / r_edge_loc(0,r)**2

       end do

    end if
    call bl_proffortfuncstop("Maestro::make_grav_edge")

  end subroutine make_grav_edge

end module make_grav_module
