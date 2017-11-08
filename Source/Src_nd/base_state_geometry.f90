! a module for storing the geometric information so we don't have to pass it
!
! This module provides the coordinate value for the left edge of a base-state
! zone (r_edge_loc) and the zone center (r_cc_loc).  As always, it is assumed that 
! the base state arrays begin with index 0, not 1.

module base_state_geometry_module

  use bl_types
  use amrex_error_module
  use amrex_constants_module
  use parallel, only: parallel_IOProcessor
  use amrex_fort_module, only: amrex_spacedim
  use meth_params_module, only: prob_lo, prob_hi, spherical, octant, &
                                anelastic_cutoff, base_cutoff_density, burning_cutoff_density

  implicit none

  integer   , save :: nlevs_radial
  real(dp_t), save :: center(3)
  integer   , save :: nr_fine, nr_irreg
  real(dp_t), save :: dr_fine

  real(dp_t), allocatable, save :: dr(:), r_cc_loc(:,:), r_edge_loc(:,:)
  integer   , allocatable, save :: nr(:)

  integer   , allocatable, save :: numdisjointchunks(:)
  integer   , allocatable, save :: r_start_coord(:,:), r_end_coord(:,:)

  integer   , allocatable, save :: anelastic_cutoff_coord(:)
  integer   , allocatable, save :: base_cutoff_density_coord(:)
  integer   , allocatable, save :: burning_cutoff_density_coord(:)

  logical   , save :: radial_initialized = .false.
  logical   , save :: cutoff_initialized = .false.
  logical   , save :: multilevel_initialized = .false.

  public :: init_base_state_geometry, compute_cutoff_coords, destroy_base_state_geometry

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_base_state_geometry(max_levs,dx_fine,domhi_fine) &
       bind(C, name="init_base_state_geometry")

    integer          , intent(in   ) :: max_levs
    double precision , intent(in   ) ::    dx_fine(1:amrex_spacedim)
    integer          , intent(in   ) :: domhi_fine(1:amrex_spacedim)

    ! local
    integer :: n,i

    if ( parallel_IOProcessor() ) then
       print*,'Calling init_base_state_geometry()'
    end if

    ! compute center(:)
    if (octant) then
       if (.not. (spherical == 1 .and. amrex_spacedim == 3 .and. &
                  prob_lo(1) == ZERO .and. &
                  prob_lo(2) == ZERO .and. &
                  prob_lo(3) == ZERO)) then
          call amrex_error("ERROR: octant requires spherical with prob_lo = 0.0")
       endif
       center(1:amrex_spacedim) = ZERO
    else
       center(1:amrex_spacedim) = HALF * (prob_lo(1:amrex_spacedim) + prob_hi(1:amrex_spacedim))
    endif


    ! computes dr, nr, r_cc_loc, r_edge_loc

    if (spherical .eq. 0) then

       ! cartesian case

       ! allocate space for dr, nr, r_cc_loc, r_edge_loc
       allocate(dr(max_levs))
       allocate(nr(max_levs))
       allocate(  r_cc_loc(max_levs,0:nr_fine-1))
       allocate(r_edge_loc(max_levs,0:nr_fine))
       
       ! compute nr_fine and dr_fine
       dr_fine = dx_fine(amrex_spacedim)
       nr_fine = domhi_fine(amrex_spacedim) + 1

       ! compute nr(:) and dr(:) assuming refinement ratio = 2
       nr(max_levs) = nr_fine
       dr(max_levs) = dr_fine
       do n=max_levs-1,1,-1
          nr(n) = nr(n+1)/2
          dr(n) = dr(n+1)*2.d0
       enddo
       
       ! compute r_cc_loc, r_edge_loc
       do n=1,max_levs
          do i = 0,nr(n)-1
             r_cc_loc(n,i) = prob_lo(amrex_spacedim) + (dble(i)+HALF)*dr(n)
          end do
          do i = 0,nr(n)
             r_edge_loc(n,i) = prob_lo(amrex_spacedim) + (dble(i))*dr(n)
          end do
       enddo

    else

       ! spherical case

       ! allocate space for dr, nr, r_cc_loc, r_edge_loc
       allocate(dr(1))
       allocate(nr(1))
       allocate(  r_cc_loc(1,0:nr_fine-1))
       allocate(r_edge_loc(1,0:nr_fine))
       
       ! compute nr_fine and dr_fine
       ! FIXME
       call amrex_error("base_state_geometry.f90: FIXME")
       
       ! compute nr(:) and dr(:)
       nr(1) = nr_fine
       dr(1) = dr_fine

       do i=0,nr_fine-1
          r_cc_loc(1,i) = (dble(i)+HALF)*dr(1)
       end do
       
       do i=0,nr_fine
          r_edge_loc(1,i) = (dble(i))*dr(1)
       end do

    end if

    radial_initialized = .true.

    if (spherical .eq. 0) then
       allocate(      anelastic_cutoff_coord(max_levs))
       allocate(   base_cutoff_density_coord(max_levs))
       allocate(burning_cutoff_density_coord(max_levs))
    else
       allocate(      anelastic_cutoff_coord(1))
       allocate(   base_cutoff_density_coord(1))
       allocate(burning_cutoff_density_coord(1))
    end if

    cutoff_initialized = .true.

  end subroutine init_base_state_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_cutoff_coords(rho0)

    real(kind=dp_t), intent(in   ) :: rho0(:,0:)

    ! local
    integer :: i,n,r,which_lev
    logical :: found

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the coordinates of the anelastic cutoff
    found = .false.

    ! find the finest level containing the anelastic cutoff density,
    ! and set the anelastic cutoff coord for this level
    do n=nlevs_radial,1,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. anelastic_cutoff) then
                   anelastic_cutoff_coord(n) = r
                   which_lev = n
                   found = .true.
                   exit
                end if
             end do
          end if

       end do
    end do

    ! if the anelastic cutoff density was not found anywhere, then set
    ! it to above the top of the domain on the finest level
    if (.not. found) then
       which_lev = nlevs_radial
       anelastic_cutoff_coord(nlevs_radial) = nr(nlevs_radial)
    endif

    ! set the anelastic cutoff coordinate on the finer levels
    do n=which_lev+1,nlevs_radial
       anelastic_cutoff_coord(n) = 2*anelastic_cutoff_coord(n-1)+1
    end do

    ! set the anelastic cutoff coordinate on the coarser levels
    do n=which_lev-1,1,-1
       if (mod(anelastic_cutoff_coord(n+1),2) .eq. 0) then
          anelastic_cutoff_coord(n) = anelastic_cutoff_coord(n+1) / 2
       else
          anelastic_cutoff_coord(n) = anelastic_cutoff_coord(n+1) / 2 + 1
       end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the coordinates of the base cutoff density
    found = .false.

    ! find the finest level containing the base cutoff density,
    ! and set the base cutoff coord for this level
    do n=nlevs_radial,1,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. base_cutoff_density) then
                   base_cutoff_density_coord(n) = r
                   which_lev = n
                   found = .true.
                   exit
                end if
             end do
          end if

       end do
    end do

    ! if the base cutoff density was not found anywhere, then set
    ! it to above the top of the domain on the finest level
    if (.not. found) then
       which_lev = nlevs_radial
       base_cutoff_density_coord(nlevs_radial) = nr(nlevs_radial)
    endif

    ! set the base cutoff coordinate on the finer levels
    do n=which_lev+1,nlevs_radial
       base_cutoff_density_coord(n) = 2*base_cutoff_density_coord(n-1)+1
    end do

    ! set the base cutoff coordinate on the coarser levels
    do n=which_lev-1,1,-1
       if (mod(base_cutoff_density_coord(n+1),2) .eq. 0) then
          base_cutoff_density_coord(n) = base_cutoff_density_coord(n+1) / 2
       else
          base_cutoff_density_coord(n) = base_cutoff_density_coord(n+1) / 2 + 1
       end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the coordinates of the burning cutoff density
    found = .false.

    ! find the finest level containing the burning cutoff density,
    ! and set the burning cutoff coord for this level
    do n=nlevs_radial,1,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. burning_cutoff_density) then
                   burning_cutoff_density_coord(n) = r
                   which_lev = n
                   found = .true.
                   exit
                end if
             end do
          end if

       end do
    end do
 
    ! if the burning cutoff density was not found anywhere, then set
    ! it to above the top of the domain on the finest level
    if (.not. found) then
       which_lev = nlevs_radial
       burning_cutoff_density_coord(nlevs_radial) = nr(nlevs_radial)
    endif

    ! set the burning cutoff coordinate on the finer levels 
    do n=which_lev+1,nlevs_radial
       burning_cutoff_density_coord(n) = 2*burning_cutoff_density_coord(n-1)+1
    end do

    ! set the burning cutoff coordinate on the coarser levels
    do n=which_lev-1,1,-1
       if (mod(burning_cutoff_density_coord(n+1),2) .eq. 0) then
          burning_cutoff_density_coord(n) = burning_cutoff_density_coord(n+1) / 2
       else
          burning_cutoff_density_coord(n) = burning_cutoff_density_coord(n+1) / 2 + 1
       end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine compute_cutoff_coords

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! FIXME
  ! need a routine similar to init_multilevel that computes
  ! numdisjointchunks, r_start_coord, r_end_coord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine destroy_base_state_geometry() bind(C, name="destroy_base_state_geometry")

    if(radial_initialized) then
      deallocate(dr, nr, r_cc_loc, r_edge_loc)
    endif
    if(cutoff_initialized) then
      deallocate(anelastic_cutoff_coord)
      deallocate(base_cutoff_density_coord)
      deallocate(burning_cutoff_density_coord)
    endif
    if(multilevel_initialized) then
      deallocate(numdisjointchunks,r_start_coord,r_end_coord)
    endif

  end subroutine destroy_base_state_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module base_state_geometry_module
