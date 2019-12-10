program fneutrinos

  use network
  use physical_constants, only: erg_per_MeV, N_avo
  use bl_space, only: MAX_SPACEDIM
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use bl_types
  use plotfile_module
  use multifab_module
  use omp_module

  implicit none

  ! argument variables
  character(len=256) :: params_file
  
  ! parameter variables
  character(len=256) :: pltfile
  integer :: nbins
  real(kind=dp_t) :: energy_delta, energy_minimum

  ! f2kcli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local variables
  integer :: uin, dim, nlevs, i, j, ii, jj, kk, ibin
  type(plotfile) :: pf
  integer :: refinement_ratio, rr_rel_fine
  real(kind=dp_t), dimension(MAX_SPACEDIM) :: dx

  integer :: ecap23_comp, beta23_comp, dens_comp
  integer :: xna23_comp, xne23_comp
  integer :: enu_ecap23_comp, enu_beta23_comp

  real(kind=dp_t) :: ecap23, beta23, yna23, yne23, enu_ecap23, enu_beta23
  real(kind=dp_t) :: density, lambda, dvol
  
  integer, dimension(MAX_SPACEDIM) :: lo, hi
  integer, dimension(MAX_SPACEDIM) :: flo, fhi
  real(kind=dp_t), pointer :: p(:,:,:,:)
  logical, allocatable :: imask(:,:,:)

  ! energy bins and weights
  real(kind=dp_t), allocatable :: energy_bins(:), eneut_weights(:), &
                                  ecap23_weights(:), beta23_weights(:)

  ! variables for output
  character(len=256) :: filename
  character(len=50) :: fmt_header, fmt_data

  namelist /params/ pltfile, nbins, energy_minimum, energy_delta

  ! For AMReX, disable nested parallel regions
  if (omp_get_max_threads() > 1) call omp_set_nested(.false.)

  ! defaults
  params_file = ''
  pltfile = ''
  nbins = 1000
  energy_minimum = 1.0_dp_t
  energy_delta = 1.0_dp_t

  ! parse arguments
  narg = command_argument_count()

  farg = 1
  do while (farg<=narg)
     call get_command_argument(farg,value=fname)

     select case(fname)

     case ('-p', '--parameters')
        farg = farg + 1
        call get_command_argument(farg,value=params_file)
     case default
        exit
     end select
     farg = farg + 1
  end do

  ! sanity check
  if (params_file == '') then
     call print_usage()
     stop
  end if

  ! Read parameters file
  uin = unit_new()
  open(unit=uin, file=trim(params_file))
  read(uin, params)
  close(unit=uin)

  print *, 'working on pltfile: ', trim(pltfile)

  call network_init()

  ! build the input plotfile
  uin = unit_new()
  call build(pf,pltfile,uin)

  nlevs = plotfile_nlevels(pf)
  dim = plotfile_dim(pf)

  ! get dx for the coarse level.
  dx = plotfile_get_dx(pf, 1)

  ! get the index bounds for the finest level.
  ! Note, lo and hi are ZERO-based indicies
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  ! imask is a flag for each zone at the finest possible level.
  ! imask = .true. means that we have not output data at this
  ! physical location from any level of data.  Once we output data
  ! at this location, imask is set to .false.  As we loop over
  ! levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)))
  imask(:,:,:) = .true.

  ! Sanity check the variables we need are present
  dens_comp = plotfile_var_index(pf,"density")
  ecap23_comp = plotfile_var_index(pf,"ecap23")
  beta23_comp = plotfile_var_index(pf,"beta23")
  xna23_comp = plotfile_var_index(pf,"X(na23)")
  xne23_comp = plotfile_var_index(pf,"X(ne23)")
  enu_ecap23_comp = plotfile_var_index(pf,"enu_ecap23")
  enu_beta23_comp = plotfile_var_index(pf,"enu_beta23")

  if (dens_comp < 0 .or. ecap23_comp < 0 .or. beta23_comp < 0 .or. &
      xna23_comp < 0 .or. xne23_comp < 0 .or. &
      enu_ecap23_comp < 0 .or. enu_beta23_comp < 0) then
     print *, ecap23_comp, beta23_comp
     print *, xna23_comp, xne23_comp
     print *, enu_ecap23_comp, enu_beta23_comp
     call bl_error("Variables not found")
  endif

  ! Allocate memory for the bins and weights
  allocate(energy_bins(nbins))
  allocate(eneut_weights(nbins))
  allocate(ecap23_weights(nbins))
  allocate(beta23_weights(nbins))

  eneut_weights  = ZERO
  ecap23_weights = ZERO
  beta23_weights = ZERO

  ! Set the values of energy_bins to the lower bound
  ! of the energy for each respective bin.
  do ibin = 1, nbins
     energy_bins(ibin) = energy_minimum + real(ibin-1, kind=dp_t)*energy_delta
  end do

  ! loop over the plotfile data starting at the finest
  rr_rel_fine = 1

  do i = nlevs, 1, -1

     ! refinement_ratio is the factor between the COARSEST level
     ! grid spacing and the current level.
     refinement_ratio = real(product(pf%refrat(1:i-1,1)), kind=dp_t)
     dvol = product(dx(1:MAX_SPACEDIM))/(refinement_ratio**MAX_SPACEDIM)

     ! loop over each box at this level
     do j = 1, nboxes(pf,i)
        ! read in the data 1 patch at a time
        call fab_bind(pf,i,j)

        lo(1:dim) = lwb(get_box(pf,i,j))
        hi(1:dim) = upb(get_box(pf,i,j))

        p => dataptr(pf,i,j)

        !$OMP PARALLEL DO PRIVATE(kk, jj, ii, ibin) &
        !$OMP PRIVATE(ecap23, beta23, yna23, yne23, enu_ecap23, enu_beta23) &
        !$OMP PRIVATE(density, lambda) &
        !$OMP REDUCTION(+:ecap23_weights) &
        !$OMP REDUCTION(+:beta23_weights) &
        !$OMP SCHEDULE(DYNAMIC,1)
        do kk = lo(3), hi(3)
           do jj = lo(2), hi(2)
              do ii = lo(1), hi(1)

                 ! Convert the cell-centered indices at the current
                 ! level into the corresponding RANGE on the finest
                 ! level, and test if we've used data in any of those
                 ! locations.  If we haven't then we use this level's
                 ! data.
                 if ( any(imask(ii*rr_rel_fine:(ii+1)*rr_rel_fine-1, &
                                jj*rr_rel_fine:(jj+1)*rr_rel_fine-1, &
                                kk*rr_rel_fine:(kk+1)*rr_rel_fine-1) ) ) then

                    density = p(ii,jj,kk,dens_comp)

                    ecap23 = p(ii,jj,kk,ecap23_comp)
                    beta23 = p(ii,jj,kk,beta23_comp)

                    yna23  = p(ii,jj,kk,xna23_comp)/23.0_dp_t
                    yne23  = p(ii,jj,kk,xne23_comp)/23.0_dp_t

                    ! Convert neutrino energies from erg to MeV
                    enu_ecap23 = p(ii,jj,kk,enu_ecap23_comp)/erg_per_MeV
                    enu_beta23 = p(ii,jj,kk,enu_beta23_comp)/erg_per_MeV

                    ! Update the ecap23 weights
                    ibin = floor((enu_ecap23 - energy_minimum)/energy_delta) + 1
                    write(*,*) 'enu_ecap23 = ', enu_ecap23
                    write(*,*) 'ibin = ', ibin
                    if (ibin .le. nbins) then
                       lambda = ecap23 * yna23 * N_avo * density * dvol
                       ecap23_weights(ibin) = ecap23_weights(ibin) + max(lambda, ZERO)
                    end if

                    ! Update the beta23 weights
                    ibin = floor((enu_beta23 - energy_minimum)/energy_delta) + 1
                    write(*,*) 'enu_beta23 = ', enu_beta23
                    write(*,*) 'ibin = ', ibin
                    if (ibin .le. nbins) then
                       lambda = beta23 * yne23 * N_avo * density * dvol
                       beta23_weights(ibin) = beta23_weights(ibin) + max(lambda, ZERO)
                    end if

                    ! mark this range of the domain as used in the mask
                    imask(ii*rr_rel_fine:(ii+1)*rr_rel_fine-1, &
                          jj*rr_rel_fine:(jj+1)*rr_rel_fine-1, &
                          kk*rr_rel_fine:(kk+1)*rr_rel_fine-1) = .false.

                 end if

              end do
           end do
        end do
        !$OMP END PARALLEL DO

        call fab_unbind(pf,i,j)
        
     end do

     ! adjust rr_rel_fine (refinement factor between the current level grid
     ! spacing and the FINEST level) for the next lowest level
     if ( i /= 1 ) rr_rel_fine = rr_rel_fine*pf%refrat(i-1,1)

  end do

  call destroy(pf)

  ! Calculate total weight
  eneut_weights = ecap23_weights + beta23_weights

  ! Save weights and bins
  uin = unit_new()
  filename = trim(pltfile) // "/neutrino_spectrum.dat"
  fmt_header = trim("(4(A25))")
  fmt_data = trim("(4(ES25.14))")
  open(unit=uin, file=filename)
  write(uin, *) "Energy in (MeV), Lambda in (#/s), Energy Spacing (MeV):"
  write(uin, fmt=trim("(ES25.14)")) energy_delta
  write(uin, fmt=fmt_header) "Energy", "Lambda_ecap23", "Lambda_beta23", "Lambda_total"
  do ibin = 1, nbins
     write(uin, fmt=fmt_data) energy_bins(ibin), ecap23_weights(ibin), &
                              beta23_weights(ibin), eneut_weights(ibin)
  end do
  close(unit=uin)

contains

  subroutine print_usage()
    implicit none
    
    print *,""
    print *, "This program takes a 3D plotfile generated by furcashell and "
    print *, " generates rate counts of neutrino energies from "
    print *, " electron capture, beta decay, and the total. "
    print *, ""
    print *, " The rate count is calculated as the total number of "
    print *, " neutrinos per second in a given energy bin. "
    print *, ""
    print *, " This can be used to approximate the neutrino energy spectrum. "
    print *, ""
    print *, " The output will be written in the plotfile directory. "
    print *, ""
    print *, "This is set up for the URCA-simple network in StarKiller Microphysics."
    print *, ""
    print *, "usage: "
    print *, " *fneutrinos* -p|--parameters <name of parameters file> "
    print *, ""
    print *, "The parameters file should contain a 'params' namelist "
    print *, "  defining the following variables. "
    print *, ""
    print *, "    pltfile: (Character string)"
    print *, "        Specify which plotfile to work on. (Required)"
    print *, ""
    print *, "    nbins: (Integer)"
    print *, "        Specify the number of energy bins to use. Defaults to 1000. "
    print *, ""
    print *, "    energy_minimum: (Real)"
    print *, "        Specify the minimum energy for binning in MeV. Defaults to 1."
    print *, ""
    print *, "    energy_delta: (Real)"
    print *, "        Specify the energy bin spacing in MeV. Defaults to 1. "
    print *, ""

  end subroutine print_usage

end program fneutrinos
