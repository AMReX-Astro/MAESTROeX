! place problem-specific variables that will go into the fortin namelist
! create a local copy probdata.f90 in your build directory

module probdata_module

  use parallel, only: parallel_IOProcessor

  implicit none

  private

  ! variables for the probdata namelist
  !
  !

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none
  
  integer :: namlen
  integer :: name(namlen)

  if (parallel_IOProcessor()) then
     print*,"Create a local copy of probdata.f90 in your build directory"
  end if

  ! abort program
  call bl_error()

end subroutine probdata_init
