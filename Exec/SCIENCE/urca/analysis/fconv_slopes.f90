! This program reads the contents of a model file, calculates the
! various thermodynamic gradients (actual, adiabatic, & Ledoux) used
! in studying convective instability, and then dumps the output to a
! file.
!
! The format of the model file is the same as in MAESTRO.

program fconv_slopes

  use eos_module, only: eos_init, eos
  use eos_type_module, only: eos_input_rt, eos_t
  use network
  use bl_constants_module
  use bl_error_module
  use bl_types
  
  implicit none

  !f2kcli stuff
  integer :: narg, farg
  character(len=256) :: fname

  ! runtime variables
  character(len=256) :: inputfile, outputfile

  ! local variables
  integer :: u_in, u_out, ipos
  character(len=256) :: line

  integer :: idens, itemp, ipres, ispec, nvars, npoints
  integer :: nvars_file
  real(kind=dp_t), allocatable :: r(:), model_data(:,:)
  real(kind=dp_t), allocatable :: actual(:), adiabatic(:), ledoux(:)

  ! the chi_* variables are (\partial p / \partial *) at constant everything 
  ! else
  real(kind=dp_t) :: chi_rho, chi_t, chi_X(nspec)
  real(kind=dp_t) :: dT, dP, dcomp(nspec), dXdP(nspec)

  type(eos_t) :: eos_state
  
  character(len=80), allocatable :: varnames_stored(:)
  real(kind=dp_t), allocatable :: vars_stored(:)

  logical :: found

  integer :: i, j, n

  real(kind=dp_t), parameter :: small = 1e-10

  ! set defaults
  inputfile = ''
  outputfile = ''

  ! variables needed are density, temperature, pressure and composition
  idens = 1
  itemp = idens + 1
  ipres = itemp + 1
  ispec = ipres + 1
  nvars = ispec + nspec - 1

  ! parse the runtime variables
  narg = command_argument_count()
  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value=fname)
     
     select case(fname)
     case('-i','--input')
        farg = farg + 1
        call get_command_argument(farg, value = inputfile)

     case('-o','--output')
        farg = farg + 1
        call get_command_argument(farg, value = outputfile)

     case default
        print *, 'Unknown command line option.'
        exit

     end select

     farg = farg + 1
  enddo

  ! sanity check
  if (trim(inputfile) == '') then
     call print_usage()
     stop
  endif

  if (trim(outputfile) == '') outputfile = trim(inputfile) // '.out'

  ! initialize
  call network_init()
  call eos_init()
  
  ! open the data file
  u_in = 11
  open(u_in, file=inputfile)

  ! number of points in the model
  read(u_in, 100) line
  ipos = index(line, '=') + 1
  read(line(ipos:),*) npoints

  allocate(r(npoints), model_data(npoints,nvars))
  allocate(actual(npoints), adiabatic(npoints), ledoux(npoints))

  ! number of variables in the model file
  read(u_in, 100) line
  ipos = index(line, '=') + 1
  read(line(ipos:),*) nvars_file
  allocate(varnames_stored(nvars_file), vars_stored(nvars_file))

  ! read the variable names
  do i = 1, nvars_file
     read(u_in,100) line
     ipos = index(line, '#') + 1
     varnames_stored(i) = trim(adjustl(line(ipos:)))
  enddo

  ! read in the data
  do i = 1, npoints
     read(u_in,*) r(i), (vars_stored(j), j = 1, nvars_file)

     do j = 1, nvars_file
        found = .false.

        select case(trim(varnames_stored(j)))
        case("density")
           model_data(i,idens) = vars_stored(j)
           found = .true.
        case("temperature")
           model_data(i,itemp) = vars_stored(j)
           found = .true.
        case("pressure") 
           model_data(i,ipres) = vars_stored(j)
           found = .true. 
        case default
           do n = 1, nspec
              if (trim(varnames_stored(j)) == spec_names(n)) then
                 model_data(i,ispec+n-1) = vars_stored(j)
                 found = .true.
                 exit
              endif
           enddo
        end select

        if (.not. found) call bl_error("ERROR: variable not found: ", &
             trim(varnames_stored(j)))
     enddo
  enddo
  close(unit=u_in)

  ! build the gradients
  do i = 1, npoints

     eos_state % rho = model_data(i,idens)
     eos_state % T = model_data(i,itemp)
     eos_state % xn(:) = model_data(i,ispec:nvars)

     call eos(eos_input_rt, eos_state)

     chi_rho = eos_state % rho * eos_state % dPdr / eos_state % p
     chi_t = eos_state % T * eos_state % dPdT / eos_state % p
     chi_X(:) = eos_state % xn(:) * eos_state % dPdX(:) / (eos_state % p * chi_t)

     ! adiabatic gradient
     adiabatic(i) = (eos_state % gam1 - chi_rho) / (eos_state % gam1 * chi_t)

     ! forward difference
     if (i == 1) then
        dT = model_data(i+1,itemp) - model_data(i  ,itemp)
        dP = model_data(i+1,ipres) - model_data(i  ,ipres)
        dcomp(:) = model_data(i+1,ispec:nvars) - model_data(i  ,ispec:nvars)
     ! backward difference
     else if (i == npoints) then
        dT = model_data(i   ,itemp) - model_data(i-1,itemp)
        dP = model_data(i   ,ipres) - model_data(i-1,ipres)
        dcomp(:) = model_data(i  ,ispec:nvars) - model_data(i-1,ispec:nvars)
     ! centered difference
     else
        dT = model_data(i+1,itemp) - model_data(i-1,itemp)
        dP = model_data(i+1,ipres) - model_data(i-1,ipres)
        dcomp(:) = model_data(i+1,ispec:nvars) - model_data(i-1,ispec:nvars)
     endif

     ! actual gradient; if we are in a constant pressure region, set
     ! gradient to ZERO to prevent div by 0
     if (abs(dp) < small) then
        actual(i) = ZERO
        dXdP = ZERO
     else
        actual(i) = model_data(i,ipres) * dT / (dP * model_data(i,itemp))
        dXdP(:) = dcomp(:)/dP
     endif

     ! ledoux gradient
     ledoux(i) = adiabatic(i)
     do n = 1, nspec
        if (model_data(i,ispec+n-1) > ZERO) then
           ledoux(i) = ledoux(i) - &
                chi_X(n) * model_data(i,ipres)*dXdP(n) &
                / model_data(i,ispec+n-1)
        endif
     enddo

  enddo

  u_out = 22
  open(u_out, file=trim(outputfile))
  write(u_out,*) "# r, actual, adiabatic, ledoux"
  do i = 1, npoints
     write(u_out,*) r(i), actual(i), adiabatic(i), ledoux(i)
  enddo
  
  close(unit=u_out)
     
100 format(a256)

contains
  subroutine print_usage()
    implicit none

    print *, 'This program reads the contents of a model file, calculates the'
    print *, 'various thermodynamic gradients (actual, adiabatic, & Ledoux) used'
    print *, 'in studying convective instability, and then dumps the output to a'
    print *, 'file. The format of the model file is the same as in MAESTRO.'
    print *, ''
    print *, ' calling sequence: '
    print *, ' fconv_slopes -i <inputfile> [-o <outputfile>]'
    print *, ''
    print *, ' arguments:'
    print *, '-i | --input:'
    print *, '     specify the input plotfile used to calculate the'
    print *, '     thermodynamic gradients.  required.'
    print *, '-o | --output:'
    print *, '     specify the output filename for the thermodynamic gradients.'
    print *, '     defaults to "<inputfile>.out"'
    print *, ''

  end subroutine print_usage
  
end program fconv_slopes
