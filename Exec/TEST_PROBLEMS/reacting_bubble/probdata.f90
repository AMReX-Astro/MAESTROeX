! place problem-specific variables that will go into the probdata namelist
! create a local copy probdata.f90 in your build directorymodule probdata_module

module probdata_module

  implicit none

  private

  ! variables for the probdata namelist
  double precision, save, public :: pert_temp_factor = 1.d0
  double precision, save, public :: pert_rad_factor = 1.d0
  logical         , save, public :: do_small_domain = .false.

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none
  
  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /probdata/ pert_temp_factor
  namelist /probdata/ pert_rad_factor
  namelist /probdata/ do_small_domain

  ! default values
  pert_temp_factor = 1.d0
  pert_rad_factor =  1.d0
  do_small_domain = .false.

  ! create the filename
  if (namlen > maxlen) then
     print *, 'probin file name too long'
     stop
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=probdata, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     print *, 'ERROR: problem in the probdata namelist'
     stop
  endif

  close (unit=un)

  !$acc update &
  !$acc device(pert_temp_factor, pert_rad_factor, do_small_domain)

end subroutine probdata_init
