! place problem-specific variables that will go into the probdata namelist
! create a local copy probdata.f90 in your build directory

module probdata_module

  implicit none

  private

  ! variables for the probdata namelist
  logical         , save, public :: apply_vel_field = .false.
  double precision, save, public :: velpert_scale = 2.5d6
  double precision, save, public :: velpert_amplitude = 1.d2
  double precision, save, public :: velpert_height_loc = 1.2d8
  integer         , save, public :: num_vortices = 2

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none
  
  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /probdata/ apply_vel_field
  namelist /probdata/ velpert_scale
  namelist /probdata/ velpert_amplitude
  namelist /probdata/ velpert_height_loc
  namelist /probdata/ num_vortices

  ! default values
  apply_vel_field = .false.
  velpert_scale = 2.5d6
  velpert_amplitude = 1.d2
  velpert_height_loc = 1.2d8
  num_vortices = 2

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
  !$acc device(apply_vel_field, velpert_scale, velpert_amplitude) &
  !$acc device(velpert_height_loc, num_vortices)

end subroutine probdata_init
