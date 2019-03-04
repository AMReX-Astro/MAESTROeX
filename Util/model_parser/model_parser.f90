module model_parser_module

  ! read in an initial model and return arrays with the model data.
  ! take care to match the species available in the model file to
  ! those defined by the network

  ! the model file is assumed to be of the follow form:
  ! # npts = 896                                     
  ! # num of variables = 6
  ! # density 
  ! # temperature 
  ! # pressure
  ! # carbon-12 
  ! # oxygen-16
  ! # magnesium-24
  ! 195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
  ! 585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0
  !
  ! we read in the number of variables and their order and use this to
  ! map them into the model_state array.  We ignore anything other than
  ! density, temperature, pressure and composition.
  !
  ! composition is assumed to be in terms of mass fractions     

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use network

  implicit none

  ! integer keys for indexing the model_state array
  integer, parameter :: nvars_model = 3 + nspec
  integer, parameter :: idens_model = 1
  integer, parameter :: itemp_model = 2
  integer, parameter :: ipres_model = 3
  integer, parameter :: ispec_model = 4

  ! number of points in the model file
  integer, save :: npts_model

  ! arrays for storing the model data
  double precision, allocatable, save :: model_state(:,:)
  double precision, allocatable, save :: model_r(:)

  ! model_initialized will be .true. once the model is read in and the
  ! model data arrays are initialized and filled
  logical, save :: model_initialized = .false.

  integer, parameter :: MAX_VARNAME_LENGTH=80

  public :: read_model_file, interpolate, close_model_file

contains

  subroutine read_model_file(model_file)

    use amrex_constants_module
    use amrex_error_module

    character(len=*), intent(in   ) :: model_file

    ! local variables
    integer :: nvars_model_file
    integer :: ierr

    integer :: i, j, comp

    double precision, allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found_model, found_dens, found_temp, found_pres
    logical :: found_spec(nspec)
    integer :: ipos
    character (len=256) :: header_line


    ! open the model file
    open(99,file=trim(model_file),status='old',iostat=ierr)

    if (ierr .ne. 0) then
       print *,'Couldnt open model_file: ',model_file
       call amrex_error('Aborting now -- please supply model_file')
    end if

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model

    ! now read in the number of variables
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

    allocate (vars_stored(nvars_model_file))
    allocate (varnames_stored(nvars_model_file))

    ! now read in the names of the variables
    do i = 1, nvars_model_file
       read (99, '(a256)') header_line
       ipos = index(header_line, '#') + 1
       varnames_stored(i) = trim(adjustl(header_line(ipos:)))
    enddo

    ! allocate storage for the model data
    allocate (model_state(npts_model, nvars_model))
    allocate (model_r(npts_model))

887 format(78('-'))
889 format(a60)

    if ( parallel_IOProcessor()) then
       write (*,889) ' '
       write (*,887)
       write (*,*)   'reading initial model'
       write (*,*)   npts_model, 'points found in the initial model file'
       write (*,*)   nvars_model_file, ' variables found in the initial model file'
    endif


    ! start reading in the data
    do i = 1, npts_model
       read(99,*) model_r(i), (vars_stored(j), j = 1, nvars_model_file)

       model_state(i,:) = ZERO

       ! make sure that each of the variables that MAESTRO cares about
       ! are found
       found_dens = .false.
       found_temp = .false.
       found_pres = .false.
       found_spec(:) = .false.

       do j = 1,nvars_model_file

          ! keep track of whether the current variable from the model
          ! file is one that MAESTRO cares about
          found_model = .false.


          if (varnames_stored(j) == "density") then
             model_state(i,idens_model) = vars_stored(j)
             found_model = .true.
             found_dens  = .true.

          else if (varnames_stored(j) == "temperature") then
             model_state(i,itemp_model) = vars_stored(j)
             found_model = .true.
             found_temp  = .true.

          else if (varnames_stored(j) == "pressure") then
             model_state(i,ipres_model) = vars_stored(j)
             found_model = .true.
             found_pres  = .true.

          else
             do comp = 1, nspec
                if (varnames_stored(j) == spec_names(comp)) then
                   model_state(i,ispec_model-1+comp) = vars_stored(j)
                   found_model = .true.
                   found_spec(comp) = .true.
                   exit
                endif
             enddo
          endif

          ! is the current variable from the model file one that we
          ! care about?
          if (.NOT. found_model .and. i == 1) then
             if ( parallel_IOProcessor() ) then
                print *, 'WARNING: variable not found: ', &
                     trim(varnames_stored(j))
             end if
          endif

       enddo   ! end loop over nvars_model_file

       ! were all the variables we care about provided?
       if (i == 1) then
          if (.not. found_dens) then
             if ( parallel_IOProcessor() ) then
                print *, 'WARNING: density not provided in inputs file'
             end if
          endif

          if (.not. found_temp) then
             if ( parallel_IOProcessor() ) then
                print *, 'WARNING: temperature not provided in inputs file'
             end if
          endif

          if (.not. found_pres) then
             if ( parallel_IOProcessor() ) then
                print *, 'WARNING: pressure not provided in inputs file'
             end if
          endif

          do comp = 1, nspec
             if (.not. found_spec(comp)) then
                if ( parallel_IOProcessor() ) then
                   print *, 'WARNING: ', trim(spec_names(comp)), &
                        ' not provided in inputs file'
                end if
             endif
          enddo
       endif

    end do   ! end loop over npts_model

    close(99)

    model_initialized = .true.

    deallocate(vars_stored,varnames_stored)

  end subroutine read_model_file


  function get_model_npts(model_file)

    integer :: get_model_npts
  
    ! look in the model file and return the number of points
    character(len=256), intent(in   ) :: model_file

    character (len=256) :: header_line
    integer :: ipos

    open(99,file=model_file)

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) get_model_npts

    close(99)

    return

  end function get_model_npts

  function interpolate(r, ivar, interpolate_top)

    ! use the module's array of model coordinates (model_r), and
    ! variables (model_state), to find the value of model_var at point
    ! r using linear interpolation.
    !
    ! Eventually, this should all be made a class.

    double precision :: interpolate
    double precision, intent(in) :: r
    integer, intent(in) :: ivar
    logical, intent(in), optional :: interpolate_top

    double precision :: slope
    double precision :: minvar, maxvar

    logical :: interp_top

    integer :: i, id

    if (present(interpolate_top)) then
       interp_top = interpolate_top
    else
       interp_top = .false.
    endif

    ! find the location in the coordinate array where we want to interpolate
    do i = 1, npts_model
       if (model_r(i) >= r) exit
    enddo
    if (i > 1 .and. i < npts_model+1) then
       if(abs(r-model_r(i-1)) < abs(r-model_r(i))) then
          i = i-1
       end if
    end if
    if (i == npts_model+1) then
       i = npts_model
    end if

    id = i

    if (id == 1) then

       slope = (model_state(id+1,ivar) - model_state(id,ivar))/(model_r(id+1) - model_r(id))
       interpolate = slope*(r - model_r(id)) + model_state(id,ivar)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_state(id+1,ivar), model_state(id,ivar))
       maxvar = max(model_state(id+1,ivar), model_state(id,ivar))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    else if (id == npts_model) then

       slope = (model_state(id,ivar) - model_state(id-1,ivar))/(model_r(id) - model_r(id-1))
       interpolate = slope*(r - model_r(id)) + model_state(id,ivar)

       ! safety check to make sure interpolate lies within the bounding points
       if (.not. interp_top) then
          minvar = min(model_state(id,ivar),model_state(id-1,ivar))
          maxvar = max(model_state(id,ivar),model_state(id-1,ivar))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
       endif

    else

       if (r >= model_r(id)) then

          ! we should not wind up in here

          slope = (model_state(id+1,ivar) - model_state(id,ivar))/(model_r(id+1) - model_r(id))
          interpolate = slope*(r - model_r(id)) + model_state(id,ivar)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_state(id+1,ivar),model_state(id,ivar))
          maxvar = max(model_state(id+1,ivar),model_state(id,ivar))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       else

          slope = (model_state(id,ivar) - model_state(id-1,ivar))/(model_r(id) - model_r(id-1))
          interpolate = slope*(r - model_r(id)) + model_state(id,ivar)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_state(id,ivar),model_state(id-1,ivar))
          maxvar = max(model_state(id,ivar),model_state(id-1,ivar))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       end if

    end if

    return

  end function interpolate

  subroutine close_model_file
    
    if (model_initialized) then
       deallocate(model_r)
       deallocate(model_state)
       npts_model = -1
       model_initialized = .false.
    endif
  end subroutine close_model_file

end module model_parser_module

