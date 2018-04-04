! place problem-specific variables that will go into the probdata namelist
! create a local copy probdata.f90 in your build directorymodule probdata_module

module probdata_module

  implicit none

  private

  ! variables for the probdata namelist
  double precision, save, public :: velpert_amplitude = 0.d0
  double precision, save, public :: velpert_radius = 0.75d8
  double precision, save, public :: velpert_scale = 0.8d8 
  double precision, save, public :: velpert_steep = 1.d0
  double precision, save, public :: tag_density_1 = 5.d7
  double precision, save, public :: tag_density_2 = 1.d8
  double precision, save, public :: tag_density_3 = 1.d8
  double precision, save, public :: particle_temp_cutoff = 6.e8
  double precision, save, public :: particle_tpert_threshold = 2.e7

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none
  
  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /probdata/ velpert_amplitude 
  namelist /probdata/ velpert_radius 
  namelist /probdata/ velpert_scale     
  namelist /probdata/ velpert_steep 
  namelist /probdata/ tag_density_1 
  namelist /probdata/ tag_density_2 
  namelist /probdata/ tag_density_3 
  namelist /probdata/ particle_temp_cutoff 
  namelist /probdata/ particle_tpert_threshold

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

end subroutine probdata_init
