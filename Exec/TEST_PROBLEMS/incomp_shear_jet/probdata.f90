! place problem-specific variables that will go into the probdata namelist
! create a local copy probdata.f90 in your build directorymodule probdata_module

module probdata_module

  implicit none

  private

  ! variables for the probdata namelist
  double precision, save, public :: yt = 1.0d0/30.0d0
  double precision, save, public :: delx = 0.5d0
  double precision, save, public :: rho_base = 1.0d0
  double precision, save, public :: p_base = 1.0d0

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none

  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /probdata/ yt
  namelist /probdata/ delx
  namelist /probdata/ rho_base
  namelist /probdata/ p_base

  ! default values
  yt = 1.0d0/30.0d0
  delx = 0.5d0
  rho_base = 1.0d0
  p_base = 1.0d0

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
