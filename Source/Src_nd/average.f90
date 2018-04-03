module average_module

  use base_state_geometry_module, only: max_radial_level, finest_radial_level, nr_fine, &
                                           restrict_base, fill_ghost_base, center
  use amrex_fort_module, only: amrex_spacedim
  use meth_params_module, only: spherical, prob_lo

  implicit none

  private

contains

  subroutine average(lev,lo,hi,phi,p_lo,p_hi,phisum, &
                      radii, nr_irreg, dx, ncell) bind (C,name="average")


    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    double precision, intent (in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision, intent (inout) :: phisum(0:max_radial_level,0:nr_fine-1)
    integer         , intent (in   ) :: nr_irreg
    double precision, intent (in   ) :: radii(0:max_radial_level,0:nr_irreg)
    double precision, intent (in   ) :: dx(3)
    integer         , intent (inout) :: ncell(0:max_radial_level,0:nr_irreg)

    ! local
    integer          :: i,j,k,index
    double precision :: x,y,z,radius

    if (spherical .eq. 0) then
       if (amrex_spacedim .eq. 1) then
          phisum(lev,:) = sum(phi(lo(1):hi(1),0,0));
       else if (amrex_spacedim .eq. 2) then
          do j=lo(2),hi(2)
             phisum(lev,j) = phisum(lev,j) + sum(phi(lo(1):hi(1),j,0))
          end do
       else if (amrex_spacedim .eq. 3) then
          do k=lo(3),hi(3)
             phisum(lev,k) = phisum(lev,k) + sum(phi(lo(1):hi(1),lo(2):hi(2),k))
          end do
       end if
    else
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k) + 0.5d0)*dx(3) - center(3)
       
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j) + 0.5d0)*dx(2) - center(2)
          
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i) + 0.5d0)*dx(1) - center(1)
             
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
                
             end do
          end do
       end do
    end if

  end subroutine average


  subroutine divide_phisum_by_ncell(phisum,ncell) bind (C,name="divide_phisum_by_ncell")

    double precision, intent(inout) :: phisum(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) ::  ncell(0:max_radial_level)

    integer :: n

    do n=0,max_radial_level
       phisum(n,:) = phisum(n,:) / ncell(n)
    end do

    call restrict_base(phisum,1)
    call fill_ghost_base(phisum,1)

  end subroutine divide_phisum_by_ncell

end module average_module
