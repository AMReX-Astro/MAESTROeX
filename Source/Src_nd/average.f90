module average_module

  use base_state_geometry_module, only: max_radial_level, finest_radial_level, nr_fine, &
       restrict_base, fill_ghost_base, center, nr_irreg, dr, base_cutoff_density_coord, &
       numdisjointchunks, r_start_coord, r_end_coord
  use amrex_fort_module, only: amrex_spacedim
  use meth_params_module, only: spherical, prob_lo, drdxfac

  implicit none

  private

contains

  subroutine average(lev,lo,hi,phi,p_lo,p_hi,phisum) bind (C,name="average")


    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    double precision, intent (in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision, intent (inout) :: phisum(0:max_radial_level,0:nr_fine-1)

    ! local
    integer          :: j,k

    if (amrex_spacedim .eq. 2) then
       do j=lo(2),hi(2)
          phisum(lev,j) = phisum(lev,j) + sum(phi(lo(1):hi(1),j,0))
       end do
    else if (amrex_spacedim .eq. 3) then
       do k=lo(3),hi(3)
          phisum(lev,k) = phisum(lev,k) + sum(phi(lo(1):hi(1),lo(2):hi(2),k))
       end do
    end if

  end subroutine average

  subroutine divide_phisum_by_ncell(phisum,ncell) bind (C,name="divide_phisum_by_ncell")

    double precision, intent(inout) :: phisum(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) ::  ncell(0:max_radial_level)

    integer :: n,i,r

    ! compute phibar by normalizing phisum
    do n=0,max_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             phisum(n,r) = phisum(n,r) / dble(ncell(n))
          end do
       end do
    end do

    call restrict_base(phisum,1)
    call fill_ghost_base(phisum,1)

  end subroutine divide_phisum_by_ncell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! spherical subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_sphr_irreg(lev,lo,hi,phi,p_lo,p_hi,phisum,ncell, &
       cc_to_r,ccr_lo,ccr_hi) bind (C,name="average_sphr_irreg")


    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    double precision, intent (in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision, intent (inout) :: phisum(0:max_radial_level,0:nr_fine-1)
    integer         , intent (inout) ::  ncell(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent (in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    ! local
    integer          :: i,j,k,index

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             index = cc_to_r(i,j,k)
             phisum(lev,index) = phisum(lev,index) + phi(i,j,k)
             ncell(lev,index) = ncell(lev,index) + 1
          end do
       end do
    end do

  end subroutine average_sphr_irreg

  subroutine divide_phisum_by_ncell_irreg(phisum,ncell) &
       bind (C,name="divide_phisum_by_ncell_irreg")

    double precision, intent(inout) :: phisum(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) ::  ncell(0:max_radial_level,0:nr_fine-1)

    integer :: n,r

    do n=0,max_radial_level
       do r=0,nr_fine-1
          ! divide only if ncell>0
          if (ncell(n,r) > 0) then
             phisum(n,r) = phisum(n,r) / ncell(n,r)
          else
             ! keep value constant if it is outside the cutoff coords
             phisum(n,r) = phisum(n,r-1)
          end if
       end do
    end do

    call restrict_base(phisum,1)
    call fill_ghost_base(phisum,1)

  end subroutine divide_phisum_by_ncell_irreg


  subroutine average_sphr(phisum,phibar,ncell,radii,finest_level) &
       bind (C,name="average_sphr")

    integer         , intent(in   ) :: finest_level
    double precision, intent(inout) :: phisum(0:finest_level,-1:nr_irreg)
    double precision, intent(inout) :: phibar(0:max_radial_level,0:nr_fine-1)
    integer         , intent(inout) ::  ncell(0:finest_level,-1:nr_irreg)
    double precision, intent(inout) :: radii(0:finest_level,-1:nr_irreg+1)

    ! local
    integer          :: j,n,r
    integer          :: min_all, min_lev
    integer          :: max_rcoord(0:finest_level), rcoord(0:finest_level)
    integer          :: stencil_coord(0:finest_level)
    double precision :: radius

    integer          :: which_lev(0:nr_fine-1)

    logical          :: limit

    ! normalize phisum so it actually stores the average at a radius
    do n=0,finest_level
       do r=0,nr_irreg
          if (ncell(n,r) .ne. 0.d0) then
             phisum(n,r) = phisum(n,r) / dble(ncell(n,r))
          end if
       end do
    end do

    ! compute center point for the finest level
    phisum(finest_level,-1) =  (11.d0/8.d0)*phisum(finest_level,0) &
         - (3.d0/8.d0)*phisum(finest_level,1)
    ncell(finest_level,-1) = 1

    ! choose which level to interpolate from
    do n=0,finest_level
       rcoord(n) = 0
    end do

    !$OMP PARALLEL DO PRIVATE(r,radius,n,j,min_all,min_lev) FIRSTPRIVATE(rcoord)
    do r=0,nr_fine-1

       radius = (dble(r)+0.5d0)*dr(0)

       ! for each level, find the closest coordinate
       do n=0,finest_level
          do j=rcoord(n),nr_irreg
             if (abs(radius-radii(n,j)) .lt. abs(radius-radii(n,j+1))) then
                rcoord(n) = j
                exit
             end if
          end do
       end do

       ! make sure closest coordinate is in bounds
       do n=0,finest_level-1
          rcoord(n) = max(rcoord(n),1)
       end do
       do n=0,finest_level
          rcoord(n) = min(rcoord(n),nr_irreg-1)
       end do

       ! choose the level with the largest min over the ncell interpolation points
       which_lev(r) = 0
       min_all = min(ncell(0,rcoord(0)-1), &
            ncell(0,rcoord(0)  ), &
            ncell(0,rcoord(0)+1))

       do n=1,finest_level
          min_lev = min(ncell(n,rcoord(n)-1), &
               ncell(n,rcoord(n)  ), &
               ncell(n,rcoord(n)+1))

          if (min_lev .gt. min_all) then
             min_all = min_lev
             which_lev(r) = n
          end if
       end do

       ! if the min hit count at all levels is zero, we expand the search
       ! to find the closest instance of where the hitcount becomes nonzero
       j = 1
       do while (min_all .eq. 0)
          j = j+1
          do n=0,finest_level
             min_lev = max(ncell(n,max(1,rcoord(n)-j)), &
                  ncell(n,min(rcoord(n)+j,nr_irreg-1)))
             if (min_lev .ne. 0) then
                which_lev(r) = n
                min_all = min_lev
                exit
             end if
          end do
       end do

    end do
    !$OMP END PARALLEL DO

    ! squish the list at each level down to exclude points with no contribution
    do n=0,finest_level
       j=0
       do r=0,nr_irreg
          do while(ncell(n,j) .eq. 0)
             j = j+1
             if (j .gt. nr_irreg) then
                exit
             end if
          end do
          if (j .gt. nr_irreg) then
             phisum(n,r:nr_irreg)   = 1.d99
             radii (n,r:nr_irreg+1) = 1.d99
             max_rcoord(n) = r-1
             exit
          end if
          phisum(n,r) = phisum(n,j)
          radii (n,r) = radii (n,j)
          ncell (n,r) = ncell (n,j)
          j = j+1
          if (j .gt. nr_irreg) then
             max_rcoord(n) = r
             exit
          end if
       end do
    end do

    ! compute phibar
    stencil_coord = 0

    !$OMP PARALLEL DO PRIVATE(r,radius,j,limit) FIRSTPRIVATE(stencil_coord)
    do r=0,nr_fine-1

       radius = (dble(r)+0.5d0)*dr(0)

       ! find the closest coordinate
       do j=stencil_coord(which_lev(r)),max_rcoord(which_lev(r))
          if (abs(radius-radii(which_lev(r) ,j  )) .lt. &
               abs(radius-radii(which_lev(r),j+1))) then
             stencil_coord(which_lev(r)) = j
             exit
          end if
       end do

       ! make sure the interpolation points will be in bounds
       if (which_lev(r) .ne. finest_level) then
          stencil_coord(which_lev(r)) = max(stencil_coord(which_lev(r)),1)
       end if
       stencil_coord(which_lev(r)) = min(stencil_coord(which_lev(r)), &
            max_rcoord(which_lev(r))-1)

       if (r > nr_fine - 1 - drdxfac*2.d0**(finest_level-1)) then
          limit = .false.
       else
          limit = .true.
       end if

       call quad_interp(radius, &
            radii(which_lev(r),stencil_coord(which_lev(r))-1), &
            radii(which_lev(r),stencil_coord(which_lev(r))  ), &
            radii(which_lev(r),stencil_coord(which_lev(r))+1), &
            phibar(0,r), &
            phisum(which_lev(r),stencil_coord(which_lev(r))-1), &
            phisum(which_lev(r),stencil_coord(which_lev(r))  ), &
            phisum(which_lev(r),stencil_coord(which_lev(r))+1), limit)

    end do
    !$OMP END PARALLEL DO

  contains

    subroutine cubic_interp(x,x0,x1,x2,x3,y,y0,y1,y2,y3)

      double precision, intent(in   ) :: x,x0,x1,x2,x3,y0,y1,y2,y3
      double precision, intent(  out) :: y

      y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1) &
           + ( ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0) &
           -((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1) ) / (x3-x0) &
           *(x-x0)*(x-x1)*(x-x2)

      if (y .gt. max(y0,y1,y2,y3)) y = max(y0,y1,y2,y3)
      if (y .lt. min(y0,y1,y2,y3)) y = min(y0,y1,y2,y3)

    end subroutine cubic_interp

    subroutine quad_interp(x,x0,x1,x2,y,y0,y1,y2,limit)

      double precision, intent(in   ) :: x,x0,x1,x2,y0,y1,y2
      double precision, intent(  out) :: y
      logical,          intent(in   ) :: limit

      y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)


      if (limit) then
         if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
         if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)
      end if

    end subroutine quad_interp

    subroutine lin_interp(x,x0,x1,y,y0,y1)

      double precision, intent(in   ) :: x,x0,x1,y0,y1
      double precision, intent(  out) :: y

      y = y0 + (y1-y0)/(x1-x0)*(x-x0)

    end subroutine lin_interp

  end subroutine average_sphr

  subroutine sum_phi_3d_sphr(lev,lo,hi,phi,p_lo,p_hi,phisum, &
                             radii, finest_level, dx, ncell, &
                             mask,     m_lo, m_hi, use_mask) &
                             bind (C,name="sum_phi_3d_sphr")

    integer         , intent (in   ) :: lev, lo(3), hi(3), finest_level
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    double precision, intent (in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision, intent (inout) :: phisum(0:finest_level,-1:nr_irreg)
    double precision, intent (in   ) :: radii(0:finest_level,-1:nr_irreg+1)
    double precision, intent (in   ) :: dx(3)
    integer         , intent (inout) :: ncell(0:finest_level,-1:nr_irreg)
    integer         , intent (in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    integer         , intent (in   ) :: use_mask

    ! local
    integer          :: i,j,k,index
    double precision :: x,y,z,radius
    logical          :: cell_valid

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k) + 0.5d0)*dx(3) - center(3)

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j) + 0.5d0)*dx(2) - center(2)

          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i) + 0.5d0)*dx(1) - center(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( use_mask .eq. 1 ) then
                if ( (mask(i,j,k).eq.1) ) cell_valid = .false.
             end if

             if (cell_valid) then

                ! compute distance to center
                radius = sqrt(x**2 + y**2 + z**2)

                ! figure out which radii index this point maps into
                index = ((radius / dx(1))**2 - 0.75d0) / 2.d0

                ! due to roundoff error, need to ensure that we are in the proper radial bin
                if (index .lt. nr_irreg) then
                   if (abs(radius-radii(lev,index)) .gt. abs(radius-radii(lev,index+1))) then
                      index = index+1
                   end if
                end if

                phisum(lev,index) = phisum(lev,index) + phi(i,j,k)
                ncell(lev,index)  = ncell(lev,index) + 1

             end if

          end do
       end do
    end do

  end subroutine sum_phi_3d_sphr

  subroutine compute_radii_sphr(lev,radii,finest_level,dx) bind (C,name="compute_radii_sphr")

    integer         , intent(in   ) :: lev, finest_level
    double precision, intent(inout) :: radii(0:finest_level,-1:nr_irreg+1)
    double precision, intent(in   ) :: dx(3)

    integer :: r

    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_irreg
       radii(lev,r) = sqrt(0.75d0+2.d0*r)*dx(1)
    end do
    !$OMP END PARALLEL DO

    radii(lev,nr_irreg+1) = 1.d99
    radii(lev,-1) = 0.d0

  end subroutine compute_radii_sphr

end module average_module
