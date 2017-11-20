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
  use meth_params_module, only: spherical, octant, &
                                anelastic_cutoff, base_cutoff_density, burning_cutoff_density

  implicit none

  integer         , save :: max_radial_level
  integer         , save :: finest_radial_level
  integer         , save :: nr_fine
  double precision, save :: dr_fine

  integer         , save :: nr_irreg
  double precision, save :: center(0:amrex_spacedim-1)

  real(dp_t), allocatable, save :: dr(:)
  integer   , allocatable, save :: nr(:)

  integer   , allocatable, save :: numdisjointchunks(:)
  integer   , allocatable, save :: r_start_coord(:,:)
  integer   , allocatable, save :: r_end_coord(:,:)

  integer   , allocatable, save :: anelastic_cutoff_coord(:)
  integer   , allocatable, save :: base_cutoff_density_coord(:)
  integer   , allocatable, save :: burning_cutoff_density_coord(:)

  public :: init_base_state_geometry, compute_cutoff_coords, destroy_base_state_geometry

contains

  subroutine init_base_state_geometry(max_radial_level_in,nr_fine_in,dr_fine_in, &
                                      r_cc_loc,r_edge_loc, &
                                      dx_fine,domhi_fine,prob_lo,prob_hi) &
       bind(C, name="init_base_state_geometry")

    integer          , intent(in   ) :: max_radial_level_in
    integer          , intent(in   ) :: nr_fine_in
    double precision , intent(in   ) :: dr_fine_in
    double precision , intent(inout) ::   r_cc_loc(0:max_radial_level_in,0:nr_fine_in-1)
    double precision , intent(inout) :: r_edge_loc(0:max_radial_level_in,0:nr_fine_in  )
    double precision , intent(in   ) ::    dx_fine(0:amrex_spacedim-1)
    integer          , intent(in   ) :: domhi_fine(0:amrex_spacedim-1)
    double precision , intent(in   ) ::    prob_lo(0:amrex_spacedim-1)
    double precision , intent(in   ) ::    prob_hi(0:amrex_spacedim-1)

    ! local
    integer :: n,i,domhi

    if ( parallel_IOProcessor() ) then
       print*,'Calling init_base_state_geometry()'
    end if

    max_radial_level = max_radial_level_in
    nr_fine = nr_fine_in
    dr_fine = dr_fine_in

    ! FIXME - we want to set this after regridding 
    finest_radial_level = max_radial_level

    ! compute center(:)
    if (octant) then
       if (.not. (spherical == 1 .and. amrex_spacedim == 3 .and. &
                  all(prob_lo(0:amrex_spacedim-1) == 0.d0) ) ) then
          call amrex_error("ERROR: octant requires spherical with prob_lo = 0.0")
       endif
       center = 0.d0
    else
       center = 0.5d0*(prob_lo + prob_hi)
    endif

    ! computes dr, nr, r_cc_loc, r_edge_loc
    if (spherical .eq. 0) then
       ! cartesian case
       
       ! allocate space for dr, nr
       allocate(dr(0:max_radial_level))
       allocate(nr(0:max_radial_level))

       ! compute nr(:) and dr(:) assuming refinement ratio = 2
       nr(max_radial_level) = nr_fine
       dr(max_radial_level) = dr_fine
       do n=max_radial_level-1,0,-1
          nr(n) = nr(n+1)/2
          dr(n) = dr(n+1)*2.d0
       enddo
       
       ! compute r_cc_loc, r_edge_loc
       do n = 0,max_radial_level
          do i = 0,nr(n)-1
             r_cc_loc(n,i) = prob_lo(amrex_spacedim-1) + (dble(i)+HALF)*dr(n)
          end do
          do i = 0,nr(n)
             r_edge_loc(n,i) = prob_lo(amrex_spacedim-1) + (dble(i))*dr(n)
          end do
       enddo

    else

       ! spherical case

       ! allocate space for dr, nr
       allocate(dr(0:0))
       allocate(nr(0:0))
       
       ! compute nr(:) and dr(:)
       nr(0) = nr_fine
       dr(0) = dr_fine

       ! compute r_cc_loc, r_edge_loc
       do i=0,nr_fine-1
          r_cc_loc(0,i) = (dble(i)+HALF)*dr(0)
       end do       
       do i=0,nr_fine
          r_edge_loc(0,i) = (dble(i))*dr(0)
       end do

       ! compute nr_irreg
       domhi = domhi_fine(0)+1
       if (.not. octant) then
          nr_irreg = (3*(domhi/2-0.5d0)**2-0.75d0)/2.d0
       else
          nr_irreg = (3*(domhi-0.5d0)**2-0.75d0)/2.d0
       endif
       
    end if

    allocate(      anelastic_cutoff_coord(0:max_radial_level))
    allocate(   base_cutoff_density_coord(0:max_radial_level))
    allocate(burning_cutoff_density_coord(0:max_radial_level))

  end subroutine init_base_state_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_cutoff_coords(rho0) bind(C, name="compute_cutoff_coords") 

     double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)

     ! local
     integer :: i,n,r,which_lev
     logical :: found


    ! compute the coordinates of the anelastic cutoff
    found = .false.

    ! find the finest level containing the anelastic cutoff density,
    ! and set the anelastic cutoff coord for this level
    do n=finest_radial_level,0,-1
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
       which_lev = finest_radial_level
       anelastic_cutoff_coord(finest_radial_level) = nr(finest_radial_level)
    endif

    ! set the anelastic cutoff coordinate on the finer levels
    do n=which_lev+1,finest_radial_level
       anelastic_cutoff_coord(n) = 2*anelastic_cutoff_coord(n-1)+1
    end do

    ! set the anelastic cutoff coordinate on the coarser levels
    do n=which_lev-1,0,-1
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
    do n=finest_radial_level,0,-1
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
       which_lev = finest_radial_level
       base_cutoff_density_coord(finest_radial_level) = nr(finest_radial_level)
    endif

    ! set the base cutoff coordinate on the finer levels
    do n=which_lev+1,finest_radial_level
       base_cutoff_density_coord(n) = 2*base_cutoff_density_coord(n-1)+1
    end do

    ! set the base cutoff coordinate on the coarser levels
    do n=which_lev-1,0,-1
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
    do n=finest_radial_level,0,-1
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
       which_lev = finest_radial_level
       burning_cutoff_density_coord(finest_radial_level) = nr(finest_radial_level)
    endif

    ! set the burning cutoff coordinate on the finer levels 
    do n=which_lev+1,finest_radial_level
       burning_cutoff_density_coord(n) = 2*burning_cutoff_density_coord(n-1)+1
    end do

    ! set the burning cutoff coordinate on the coarser levels
    do n=which_lev-1,0,-1
       if (mod(burning_cutoff_density_coord(n+1),2) .eq. 0) then
          burning_cutoff_density_coord(n) = burning_cutoff_density_coord(n+1) / 2
       else
          burning_cutoff_density_coord(n) = burning_cutoff_density_coord(n+1) / 2 + 1
       end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine compute_cutoff_coords

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_multilevel() bind(C, name="init_multilevel")
    
    ! compute numdisjointchunks, r_start_coord, r_end_coord
    ! FIXME - right now there is one chunk at each level that spans the domain
    
    integer :: n

    call amrex_abort("init_multilevel does not work with AMR")

    if (allocated(numdisjointchunks)) then
       deallocate(numdisjointchunks)
    end if
    allocate(numdisjointchunks(0:finest_radial_level))

    if (allocated(r_start_coord)) then
       deallocate(r_start_coord)
    end if
    allocate(r_start_coord(0:finest_radial_level,1)) ! FIXME - for > 1 chunk case

    if (allocated(r_end_coord)) then
       deallocate(r_end_coord)
    end if
    allocate(r_end_coord(0:finest_radial_level,1)) ! FIXME - for > 1 chunk case

    do n=0,finest_radial_level
       numdisjointchunks(n) = 1
       r_start_coord(n,1) = 0
       r_end_coord(n,1) = nr(n)
    end do

  end subroutine init_multilevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine destroy_base_state_geometry() bind(C, name="destroy_base_state_geometry")

      deallocate(dr, nr)
      deallocate(anelastic_cutoff_coord)
      deallocate(base_cutoff_density_coord)
      deallocate(burning_cutoff_density_coord)
      deallocate(numdisjointchunks,r_start_coord,r_end_coord)

  end subroutine destroy_base_state_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module base_state_geometry_module
