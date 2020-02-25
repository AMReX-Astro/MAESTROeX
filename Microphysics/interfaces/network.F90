! the network module provides the information about the species we are
! advecting:
!
! This module contains the following routines:
!
!  network_init          -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!
!  network_finalize      -- do any network cleanup

module network

  use amrex_fort_module, only : rt => amrex_real
  use actual_network

  implicit none

  logical :: network_initialized = .false.

contains

  subroutine network_init()

    use amrex_error_module, only : amrex_error
    use amrex_constants_module, only : ONE

    implicit none

    ! First, we call the specific network initialization.
    ! This should set the number of species and number of
    ! aux variables, and the components of the species.
    ! Note that the network MUST define nspec and naux
    ! as parameters, or else the compiler will throw an error.

    call actual_network_init()

    ! Check to make sure, and if not, throw an error.

    if ( nspec .le. 0 ) then
       call amrex_error("Network cannot have a negative number of species.")
    endif

    if ( naux .lt. 0 ) then
       call amrex_error("Network cannot have a negative number of auxiliary variables.")
    endif

    network_initialized = .true.

  end subroutine network_init


  function get_network_name() result(name)

    use actual_network, only: network_name

    character(len=128) :: name

    name = network_name

  end function get_network_name


  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          r = n
          return
       endif
    enddo

  end function network_species_index


  function get_network_species_name(index) result(name)

    character(len=128) :: name
    integer :: index

    if (index < 1 .or. index > nspec) then
       name = ""
    else
       name = spec_names(index)
    endif

  end function get_network_species_name


  function get_network_short_species_name(index) result(name)

    character(len=128) :: name
    integer :: index

    if (index < 1 .or. index > nspec) then
       name = ""
    else
       name = short_spec_names(index)
    endif

  end function get_network_short_species_name


  subroutine network_finalize()
    implicit none

    call actual_network_finalize()

  end subroutine network_finalize

end module network
