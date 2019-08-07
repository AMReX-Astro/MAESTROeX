! a module for storing the geometric information so we don't have to pass it
!
! This module provides the coordinate value for the left edge of a base-state
! zone (r_edge_loc) and the zone center (r_cc_loc).  As always, it is assumed that
! the base state arrays begin with index 0, not 1.

module base_state_geometry_module

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use amrex_constants_module
  use amrex_fort_module, only: amrex_spacedim
  use meth_params_module, only: spherical, octant, anelastic_cutoff_density, base_cutoff_density, &
       burning_cutoff_density, prob_lo, prob_hi, &
       use_exact_base_state

  implicit none

  private

  public :: restrict_base, fill_ghost_base

  ! note that for spherical problems the base state only has one level of refiment,
  ! for for spherical, max_radial_level = finest_radial_level = 0
  ! for planar, max_radial_level is the index of the finest possible radial level
  ! finest_radial_level is the current finest level index, which may be less than
  ! max_radial_level depending on the refinement criteria

  integer         , allocatable, save, public :: max_radial_level
  integer         , allocatable, save, public :: finest_radial_level
  integer         , allocatable, save, public :: nr_fine   ! number of zones associated with *max_radial_level*
  double precision, allocatable, save, public :: dr_fine   ! base state grid spacing associated with *max_radial_level*

  integer         , allocatable, save, public  :: nr_irreg
  double precision, allocatable, save, public  :: center(:)

  double precision, allocatable, save, public  :: dr(:)
  integer         , allocatable, save, public  :: nr(:)

  integer         , allocatable, save, public  :: numdisjointchunks(:)
  integer         , pointer, save, public  :: r_start_coord(:,:)
  integer         , pointer, save, public  :: r_end_coord(:,:)

  integer         , allocatable, save, public  :: anelastic_cutoff_density_coord(:)
  integer         , allocatable, save, public  :: base_cutoff_density_coord(:)
  integer         , pointer, save, public  :: burning_cutoff_density_coord(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: max_radial_level, finest_radial_level, nr_fine, dr_fine
  attributes(managed) ::  nr_irreg, center, dr, nr
  attributes(managed) :: base_cutoff_density_coord, anelastic_cutoff_density_coord
  attributes(managed) :: numdisjointchunks
#endif

contains

  subroutine init_base_state_geometry(max_radial_level_in,nr_fine_in,dr_fine_in, &
       r_cc_loc,r_edge_loc, &
       dx_fine, &
       nr_irreg_in) &
       bind(C, name="init_base_state_geometry")

    integer          , intent(in   ) :: max_radial_level_in
    integer          , intent(in   ) :: nr_fine_in
    double precision , intent(in   ) :: dr_fine_in
    integer          , intent(in   ) :: nr_irreg_in
    double precision , intent(inout) ::   r_cc_loc(0:max_radial_level_in,0:nr_fine_in-1)
    double precision , intent(inout) :: r_edge_loc(0:max_radial_level_in,0:nr_fine_in  )
    double precision , intent(in   ) :: dx_fine(0:amrex_spacedim-1)

    ! local
    integer :: n,i

    if ( parallel_IOProcessor() ) then
       print*,'Calling init_base_state_geometry()'
    end if

    allocate(max_radial_level)
    allocate(finest_radial_level)
    allocate(nr_fine)
    allocate(dr_fine)
    allocate(nr_irreg)

    max_radial_level = max_radial_level_in
    finest_radial_level = max_radial_level ! FIXME - we want to set this after regridding
    nr_fine = nr_fine_in
    dr_fine = dr_fine_in
    nr_irreg = nr_irreg_in

    allocate(center(3))
    allocate(dr(0:max_radial_level))
    allocate(nr(0:max_radial_level))
    allocate(base_cutoff_density_coord(0:max_radial_level))
    allocate(anelastic_cutoff_density_coord(0:max_radial_level))

    ! compute center(:)
    if (octant) then
       if (.not. (spherical == 1 .and. amrex_spacedim == 3 .and. &
            all(prob_lo(1:amrex_spacedim) == 0.d0) ) ) then
          call amrex_error("ERROR: octant requires spherical with prob_lo = 0.0")
       endif
       center = 0.d0
    else
       center = 0.5d0*(prob_lo + prob_hi)
    endif

    ! ! allocate space for dr, nr
    ! call bl_allocate(nr,0,max_radial_level)

    ! compute nr(:) and dr(:)
    nr(max_radial_level) = nr_fine
    dr(max_radial_level) = dr_fine

    ! computes dr, nr, r_cc_loc, r_edge_loc
    if (spherical .eq. 0) then
       ! cartesian case

       ! compute nr(:) and dr(:) assuming refinement ratio = 2
       do n=max_radial_level-1,0,-1
          nr(n) = nr(n+1)/2
          dr(n) = dr(n+1)*2.d0
       enddo

       ! compute r_cc_loc, r_edge_loc
       do n = 0,max_radial_level
          do i = 0,nr(n)-1
             r_cc_loc(n,i) = prob_lo(amrex_spacedim) + (dble(i)+HALF)*dr(n)
          end do
          do i = 0,nr(n)
             r_edge_loc(n,i) = prob_lo(amrex_spacedim) + (dble(i))*dr(n)
          end do
       enddo

    else

       ! spherical case
       ! compute r_cc_loc, r_edge_loc
       if (use_exact_base_state) then
          ! nr_fine = nr_irreg + 1
          do i=0,nr_fine-1
             r_cc_loc(0,i) = sqrt(0.75d0+2.d0*i)*dx_fine(0)
          end do
          r_edge_loc(0,0) = 0.d0
          do i=0,nr_fine-1
             r_edge_loc(0,i+1) = sqrt(0.75d0+2.d0*(i+0.5d0))*dx_fine(0)
          end do
       else
          do i=0,nr_fine-1
             r_cc_loc(0,i) = (dble(i)+HALF)*dr(0)
          end do
          do i=0,nr_fine
             r_edge_loc(0,i) = (dble(i))*dr(0)
          end do
       end if

    end if

    ! call bl_allocate(   base_cutoff_density_coord,0,max_radial_level)
    call bl_allocate(burning_cutoff_density_coord,0,max_radial_level)

  end subroutine init_base_state_geometry

  subroutine init_base_state_map_sphr(cc_to_r, lo, hi, &
       dx_fine, dx_lev) &
       bind(C, name="init_base_state_map_sphr")

    integer          , intent(in   ) :: lo(3), hi(3)
    double precision , intent(inout) :: cc_to_r(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision , intent(in   ) :: dx_fine(0:amrex_spacedim-1)
    double precision , intent(in   ) ::  dx_lev(3)

    ! local
    integer :: i,j,k
    double precision :: index,x,y,z

    if ( spherical .eq. 0 ) then
       print*,'init_base_state_map_sphr() does not work for planar'
       call abort()
    end if

    ! map cell centers to base state indices
    if (use_exact_base_state) then

       do k = lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx_lev(3) - center(3)
          do j = lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx_lev(2) - center(2)
             do i = lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx_lev(1) - center(1)

                index = (x**2 + y**2 + z**2)/(2.0d0*dx_fine(0)**2) - 0.375d0
                cc_to_r(i,j,k) = nint(index)
             end do
          end do
       end do

    end if

  end subroutine init_base_state_map_sphr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_cutoff_coords(rho0) bind(C, name="compute_cutoff_coords")

    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)

    ! local
    integer :: i,n,r,which_lev
    logical :: found

    call bl_proffortfuncstart("Maestro::compute_cutoff_coords")

    ! compute the coordinates of the anelastic cutoff
    found = .false.

    ! find the finest level containing the anelastic cutoff density,
    ! and set the anelastic cutoff coord for this level
    do n=finest_radial_level,0,-1
       do i=1,numdisjointchunks(n)

          if (.not. found) then
             do r=r_start_coord(n,i),r_end_coord(n,i)
                if (rho0(n,r) .le. anelastic_cutoff_density) then
                   anelastic_cutoff_density_coord(n) = r
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
       anelastic_cutoff_density_coord(finest_radial_level) = nr(finest_radial_level)
    endif

    ! set the anelastic cutoff coordinate on the finer levels
    do n=which_lev+1,finest_radial_level
       anelastic_cutoff_density_coord(n) = 2*anelastic_cutoff_density_coord(n-1)+1
    end do

    ! set the anelastic cutoff coordinate on the coarser levels
    do n=which_lev-1,0,-1
       if (mod(anelastic_cutoff_density_coord(n+1),2) .eq. 0) then
          anelastic_cutoff_density_coord(n) = anelastic_cutoff_density_coord(n+1) / 2
       else
          anelastic_cutoff_density_coord(n) = anelastic_cutoff_density_coord(n+1) / 2 + 1
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

    call bl_proffortfuncstop("Maestro::compute_cutoff_coords")

  end subroutine compute_cutoff_coords

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_multilevel(tag_array, finest_radial_level_in) bind(C, name="init_multilevel")

    ! compute numdisjointchunks, r_start_coord, r_end_coord
    ! FIXME - right now there is one chunk at each level that spans the domain

    integer, intent(in   ) :: tag_array(0:max_radial_level,0:nr_fine-1)
    integer, intent(in   ) :: finest_radial_level_in

    integer :: n, r
    integer :: nchunks, maxchunks
    logical :: chunk_start

    if (spherical .eq. 1) then
       finest_radial_level = 0
    else
       finest_radial_level = finest_radial_level_in
    end if

    if (allocated(numdisjointchunks)) then
       deallocate(numdisjointchunks)
    end if
    allocate(numdisjointchunks(0:finest_radial_level))

    ! loop through tag_array first to determine the maximum number of chunks
    ! to use for allocating r_start_coord and r_end_coord
    maxchunks = 1
    do n=1,finest_radial_level

       ! initialize variables
       chunk_start = .false.
       nchunks = 0

       ! increment nchunks at beginning of each chunk
       ! (ex. when the tagging index changes from 0 to 1)
       do r=0,nr(n-1)-1
          if (tag_array(n-1,r).gt.0 .AND. .not.chunk_start) then
             chunk_start = .true.
             nchunks = nchunks + 1
          elseif (tag_array(n-1,r).eq.0 .AND. chunk_start) then
             chunk_start = .false.
          end if
       end do

       maxchunks = max(nchunks,maxchunks)

    end do

    if (associated(r_start_coord)) then
       call bl_deallocate(r_start_coord)
    end if
    call bl_allocate(r_start_coord,0,finest_radial_level,1,maxchunks)

    if (associated(r_end_coord)) then
       call bl_deallocate(r_end_coord)
    end if
    call bl_allocate(r_end_coord,0,finest_radial_level,1,maxchunks)

    if (spherical .eq. 0) then

       ! coarsest grid always has 1 chunk of data
       numdisjointchunks(0) = 1
       r_start_coord(0,1) = 0
       r_end_coord(0,1) = nr(0)-1

       ! for > 1 chunks (multilevel)
       do n=1,finest_radial_level
          ! initialize variables
          chunk_start = .false.
          numdisjointchunks(n) = 0

          ! increment numdisjointchunks at beginning of each chunk
          ! (ex. when the tagging index changes from 0 to 1)
          do r=0,nr(n-1)-1
             if (tag_array(n-1,r).gt.0 .AND. .not.chunk_start) then
                chunk_start = .true.
                numdisjointchunks(n) = numdisjointchunks(n) + 1
                r_start_coord(n,numdisjointchunks(n)) = 2*r
             elseif (tag_array(n-1,r).eq.0 .AND. chunk_start) then
                r_end_coord(n,numdisjointchunks(n)) = 2*r-1
                chunk_start = .false.
             elseif (r.eq.nr(n-1)-1 .AND. chunk_start) then
                ! if last chunk is at the end of array
                r_end_coord(n,numdisjointchunks(n)) = 2*r-1
             end if
          end do
       end do

    else

       numdisjointchunks(0) = 1
       r_start_coord(0,1) = 0
       r_end_coord(0,1) = nr(0)-1

    end if

!!$    print *,"hack,",numdisjointchunks
!!$    print *,"hack,",r_start_coord(1,:),r_end_coord(1,:)

  end subroutine init_multilevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restrict_base(s0,is_cell_centered)

    double precision, intent(inout) :: s0(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) :: is_cell_centered

    ! local
    integer :: n, r, i

    if (is_cell_centered .eq. 1) then

       do n=finest_radial_level,1,-1
          ! for level n, make the coarser cells underneath simply the average of the fine
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)-1,2
                s0(n-1,r/2) = 0.5d0 * (s0(n,r) + s0(n,r+1))
             end do
          end do
       end do

    else

       do n=finest_radial_level,1,-1
          ! for level n, make the coarse edge underneath equal to the fine edge value
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,1),r_end_coord(n,1)+1,2
                s0(n-1,r/2) = s0(n,r)
             end do
          end do
       end do

    end if

  end subroutine restrict_base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fill_ghost_base(s0,is_cell_centered)

    double precision, intent(inout) :: s0(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) :: is_cell_centered

    ! local
    integer          :: n,i,r_crse
    double precision :: del,dpls,dmin,slim,slope

    if (is_cell_centered .eq. 1) then

       ! compute limited slopes at the coarse level and set 4 fine ghost cells
       ! maintaining conservation
       do n=finest_radial_level,1,-1
          do i=1,numdisjointchunks(n)

             ! lo side
             if (r_start_coord(n,i) .ne. 0) then
                r_crse = r_start_coord(n,i)/2-1
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_start_coord(n,i)-1) = s0(n-1,r_crse) + FOURTH*slope
                s0(n,r_start_coord(n,i)-2) = s0(n-1,r_crse) - FOURTH*slope

                r_crse = r_start_coord(n,i)/2-2
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_start_coord(n,i)-3) = s0(n-1,r_crse) + FOURTH*slope
                s0(n,r_start_coord(n,i)-4) = s0(n-1,r_crse) - FOURTH*slope
             end if

             ! hi side
             if (r_end_coord(n,i) .ne. nr(n)-1) then
                r_crse = (r_end_coord(n,i)+1)/2
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_end_coord(n,i)+1) = s0(n-1,r_crse) - FOURTH*slope
                s0(n,r_end_coord(n,i)+2) = s0(n-1,r_crse) + FOURTH*slope

                r_crse = (r_end_coord(n,i)+1)/2+1
                del = HALF*(s0(n-1,r_crse+1)-s0(n-1,r_crse-1))
                dpls = TWO*(s0(n-1,r_crse+1)-s0(n-1,r_crse  ))
                dmin = TWO*(s0(n-1,r_crse  )-s0(n-1,r_crse-1))
                slim = min(abs(dpls),abs(dmin))
                slim = merge(slim,ZERO,dpls*dmin.gt.ZERO)
                slope=sign(ONE,del)*min(slim,abs(del))
                s0(n,r_end_coord(n,i)+3) = s0(n-1,r_crse) - FOURTH*slope
                s0(n,r_end_coord(n,i)+4) = s0(n-1,r_crse) + FOURTH*slope
             end if

          end do
       end do

    else

       do n=finest_radial_level,1,-1
          do i=1,numdisjointchunks(n)

             if (r_start_coord(n,i) .ne. 0) then
                ! quadratic interpolation from the three closest points
                s0(n,r_start_coord(n,i)-1) = -THIRD*s0(n,r_start_coord(n,i)+1) &
                     + s0(n,r_start_coord(n,i)) + THIRD*s0(n-1,r_start_coord(n,i)/2-1)
                ! copy the next ghost cell value directly in from the coarser level
                s0(n,r_start_coord(n,i)-2) = s0(n-1,(r_start_coord(n,i)-2)/2)
             end if

             if (r_end_coord(n,i)+1 .ne. nr(n)) then
                ! quadratic interpolation from the three closest points
                s0(n,r_end_coord(n,i)+2) = -THIRD*s0(n,r_end_coord(n,i)) &
                     + s0(n,r_end_coord(n,i)+1) + THIRD*s0(n-1,(r_end_coord(n,i)+3)/2)
                ! copy the next ghost cell value directly in from the coarser level
                s0(n,r_end_coord(n,i)+3) = s0(n-1,(r_end_coord(n,i)+3)/2)
             end if

          end do
       end do

    end if

  end subroutine fill_ghost_base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine destroy_base_state_geometry() bind(C, name="destroy_base_state_geometry")

    deallocate(max_radial_level)
    deallocate(finest_radial_level)
    deallocate(nr_fine)
    deallocate(dr_fine)
    deallocate(nr_irreg)
    deallocate(center)
    deallocate(dr)

    deallocate(nr)
    deallocate(anelastic_cutoff_density_coord)
    deallocate(base_cutoff_density_coord)
    call bl_deallocate(burning_cutoff_density_coord)
    deallocate(numdisjointchunks)
    call bl_deallocate(r_start_coord)
    call bl_deallocate(r_end_coord)

  end subroutine destroy_base_state_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module base_state_geometry_module
