! place problem-specific variables that will go into the probdata namelist
! create a local copy probdata.f90 in your build directorymodule probdata_module

module probdata_module

  implicit none

  private

  ! variables for the probdata namelist
  double precision, save, public :: rho_1 = 1.d0
  double precision, save, public :: rho_2 = 2.d0
  double precision, save, public :: vel_amplitude = 1.d0
  double precision, save, public :: vel_width = 1.d0
  double precision, save, public :: p0_base = 1.0d0
  integer, save, public :: nmodes = 1

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none

  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /probdata/ rho_1
  namelist /probdata/ rho_2
  namelist /probdata/ vel_amplitude
  namelist /probdata/ vel_width
  namelist /probdata/ p0_base
  namelist /probdata/ nmodes

  ! default values
  rho_1 = 1.d0
  rho_2 =  2.d0
  vel_amplitude = 1.0d0
  vel_width = 1.0d0
  p0_base = 1.0d0
  nmodes = 1

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
  !$acc device(rho_1, rho_2, vel_amplitude, vel_width, nmodes)

end subroutine probdata_init
