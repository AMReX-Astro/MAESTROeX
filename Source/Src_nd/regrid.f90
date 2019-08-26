module regrid_base_module
  ! Regrid base state variables (ex. psi, etarho, rho0, etc.)
  !
  ! We copy the coarsest level only, interpolate to all
  ! the other levels and then copy the valid data from
  ! the old arrays onto the new.

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, dr, nr, &
       max_radial_level, numdisjointchunks, &
       r_start_coord, r_end_coord, finest_radial_level

  implicit none

  private

contains

  subroutine regrid_base_state_cc(state_cc) &
       bind(C, name="regrid_base_state_cc")
    ! Binds to C function ``regrid_base_state_cc``

    double precision, intent(inout) :: state_cc(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer :: r,n,i
    double precision, allocatable :: state_temp(:,:)

    allocate(state_temp(0:max_radial_level,0:nr_fine-1))

    ! copy the coarsest level of the real arrays into the
    ! temp arrays
    state_temp(0,:) = state_cc(0,:)

    ! piecewise linear interpolation to fill the cc temp arrays
    do n=1,max_radial_level
       do r=0,nr(n)-1
          if (r .eq. 0 .or. r .eq. nr(n)-1) then
             state_temp(n,r) = state_temp(n-1,r/2)
          else
             if (mod(r,2) .eq. 0) then
                state_temp(n,r) = 0.75d0*state_temp(n-1,r/2) &
                                  + 0.25d0*state_temp(n-1,r/2-1)
             else
                state_temp(n,r) = 0.75d0*state_temp(n-1,r/2) &
                                  + 0.25d0*state_temp(n-1,r/2+1)
             end if
          end if
       end do
    end do

    ! copy valid data into temp
    do n=1,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             state_temp(n,r) = state_cc(n,r)
          end do
       end do
    end do


    ! copy temp array back into the real thing
    state_cc = state_temp

    deallocate(state_temp)

  end subroutine regrid_base_state_cc

  subroutine regrid_base_state_edge(state_ec) &
       bind(C, name="regrid_base_state_edge")
    ! Binds to C function ``regrid_base_state_edge``

    double precision, intent(inout) :: state_ec(0:max_radial_level,0:nr_fine)

    ! Local variables
    integer :: r,n,i
    double precision, allocatable :: state_temp(:,:)

    allocate(state_temp(0:max_radial_level,0:nr_fine))

    ! copy the coarsest level of the real arrays into the
    ! temp arrays
    state_temp(0,:) = state_ec(0,:)

    ! piecewise linear interpolation to fill the edge-centered temp arrays
    do n=1,max_radial_level
       do r=0,nr(n)
          if (mod(r,2) .eq. 0) then
             state_temp(n,r) = state_temp(n-1,r/2)
          else
             state_temp(n,r) = HALF*(state_temp(n-1,r/2)+state_temp(n-1,r/2+1))
          end if
       end do
    end do

    ! copy valid data into temp
    do n=1,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)+1
             state_temp(n,r) = state_ec(n,r)
          end do
       end do
    end do


    ! copy temp array back into the real thing
    state_ec = state_temp

    deallocate(state_temp)

  end subroutine regrid_base_state_edge

end module regrid_base_module
