!
! This program takes a 3-d cartesian plotfile and calculates
! the actual, adiabatic, and ledoux thermodynamic gradients.
!
! The actual thermodynamic gradient is computed along the
! radial direction, constructed from gradients along
! the cartesian coordinate directions.
!
! See MAESTRO/docs/thermo_notes for details.
!

program fconv_radial

  use bl_space, only: MAX_SPACEDIM
  use bl_error_module
  use bl_types
  use bl_constants_module, only: ZERO, HALF
  use bl_IO_module
  use plotfile_module
  use multifab_module
  use network
  use eos_module
  use eos_type_module
  use omp_module

  implicit none

  ! argument variables
  character(len=256) :: pltfile, outputfile
  real(kind=dp_t) :: low_cutoff

  ! f2kcli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local variables
  type(plotfile) :: pf
  type(layout) :: la
  type(boxarray) :: ba
  type(list_box) :: bl
  type(box) :: pd
  integer :: rr
  integer, allocatable :: ref_ratio(:)
  real(kind=dp_t),dimension(MAX_SPACEDIM) :: prob_lo, prob_hi
  real(kind=dp_t) :: time
  type(multifab), allocatable :: conv_grad(:)

  integer :: uin
  integer :: dim, nlevs
  integer :: i, j, n, ii, jj, kk
  real(kind=dp_t) :: xctr, yctr, zctr ! center coordinates
  real(kind=dp_t) :: xpos, ypos, zpos, rpos ! relative to center
  real(kind=dp_t) :: xx, yy, zz ! relative to domain
  real(kind=dp_t) :: dxx, dyy, dzz, drr
  real(kind=dp_t), dimension(MAX_SPACEDIM) :: dx
  integer, dimension(MAX_SPACEDIM) :: flo, fhi, lo, hi
  real(kind=dp_t), pointer :: p(:,:,:,:), ap(:,:,:,:)
  real(kind=dp_t), allocatable :: pres(:,:,:)
  real(kind=dp_t) :: chi_rho, chi_t, nabla
  real(kind=dp_t) :: chi_X(nspec), dXdP(nspec)
  real(kind=dp_t) :: dtdxx, dtdyy, dtdzz, dtdrr
  real(kind=dp_t) :: dpdxx, dpdyy, dpdzz, dpdrr
  real(kind=dp_t) :: dxdxx(nspec), dxdyy(nspec), dxdzz(nspec), dxdrr(nspec)

  type(eos_t) :: eos_state

  integer :: dens_comp, spec_comp, temp_comp
  integer :: actl_comp, adic_comp, ledx_comp

  character(len=20) :: component_names(3)

  real(kind=dp_t), parameter :: small = 1.e-14

  ! For AMReX, disable nested parallel regions
  if (omp_get_max_threads() > 1) call omp_set_nested(.false.)

  ! Set unit for input plotfile
  uin = unit_new()

  ! defaults
  pltfile = ''
  outputfile = 'conv_radial'
  low_cutoff = 1.e4

  component_names(1) = "actual gradient"
  component_names(2) = "adiabatic gradient"
  component_names(3) = "ledoux gradient"

  ! parse the arguements
  narg = command_argument_count()

  farg = 1
  do while (farg <= narg) 
     call get_command_argument(farg, value = fname)

     select case(fname)
        
     case ('-i', '--input')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-o', '--output')
        farg = farg + 1
        call get_command_argument(farg, value = outputfile)

     case ('-l', '--low_cutoff')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) low_cutoff

     case default
        exit

     end select
     farg = farg + 1
  enddo

  ! sanity check
  if (pltfile == '') then
     call print_usage()
     stop
  endif

  print *, 'working on pltfile:', pltfile
  print *, 'dumping data to file: ', outputfile
  print *, 'using low density cutoff:', low_cutoff

  call network_init()
  call eos_init()

  ! build the input plotfile
  call build(pf, pltfile, uin)

  nlevs = pf%flevel
  dim = pf%dim

  dens_comp = plotfile_var_index(pf, "density")
  ! For the URCA-simple network
  spec_comp = plotfile_var_index(pf, "X(n)")
  ! For the Chamulak network
  ! spec_comp = plotfile_var_index(pf, "X(C12)")  
  temp_comp = plotfile_var_index(pf, "tfromp")

  actl_comp = 1
  adic_comp = 2
  ledx_comp = 3

  if (dens_comp < 0 .or. spec_comp < 0 .or. temp_comp < 0) &
       call bl_error("Variables not found")

  allocate(ref_ratio(nlevs-1),conv_grad(nlevs))

  ! get the index bounds for the finest level
  flo(1:dim) = lwb(plotfile_get_pd_box(pf, nlevs))
  fhi(1:dim) = upb(plotfile_get_pd_box(pf, nlevs))

  ! get dx for the coarse level.  
  dx = plotfile_get_dx(pf, 1)

  ! the default for the center of the star will be the geometric center
  ! of the domain
  xctr = HALF*(pf%phi(1) + pf%plo(1))
  yctr = HALF*(pf%phi(2) + pf%plo(2))
  zctr = HALF*(pf%phi(3) + pf%plo(3))

  allocate(pres(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3)))

  ! define the ref ratios
  do i = 1, nlevs-1
     ref_ratio(i) = pf%refrat(i,1)
  enddo

  do i = 1, pf%flevel
     do j = 1, nboxes(pf,i)
        call push_back(bl,get_box(pf,i,j))
     enddo

     call build(ba,bl)

     call layout_build_ba(la, ba, boxarray_bbox(ba))
     
     call destroy(bl)
     call destroy(ba)

     ! build the mutifab with 0 ghost cells and 3 components
     ! component 1: actual convective gradient
     ! component 2: adiabatic gradient condition
     ! component 3: ledoux gradient condition
     call multifab_build(conv_grad(i),la,3,0)

  enddo
     
  ! loop over the plotfile data starting at the first
  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))
     
     ! loop over each box at this level
     do j = 1, nboxes(pf, i)

        ! read in the data 1 patch at a time
        call fab_bind(pf, i, j)

        ! get the integer bounds of the current box, in terms of this 
        ! level's index space
        lo(1:dim) = lwb(get_box(pf,i,j))
        hi(1:dim) = upb(get_box(pf,i,j))

        ! pointers to the plotfile data and the conv_grad data
        p => dataptr(pf, i, j)
        ap => dataptr(conv_grad(i), j)

        ! do a first 3D loop through the data to initialize pressure array
        ! for finite differencing in the next 3D loop
        !$OMP PARALLEL DO PRIVATE(kk, jj, ii, eos_state) &
        !$OMP SCHEDULE(DYNAMIC,1)
        do kk = lo(3), hi(3)
           do jj = lo(2), hi(2)
              do ii = lo(1), hi(1)
        
                 eos_state % rho = p(ii,jj,kk,dens_comp)
                 eos_state % T = p(ii,jj,kk,temp_comp)
                 eos_state % xn(:) = p(ii,jj,kk,spec_comp:spec_comp+nspec-1)

                 call eos(eos_input_rt, eos_state)

                 pres(ii,jj,kk) = eos_state % p

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(kk, dzz, zz, jj, dyy, yy, ii, dxx, xx) &
        !$OMP PRIVATE(eos_state, chi_rho, chi_t, chi_X, dXdP, n, nabla) &
        !$OMP PRIVATE(dtdxx, dpdxx, dxdxx) &
        !$OMP PRIVATE(dtdyy, dpdyy, dxdyy) &
        !$OMP PRIVATE(dtdzz, dpdzz, dxdzz) &
        !$OMP PRIVATE(dtdrr, dpdrr, dxdrr) &
        !$OMP PRIVATE(xpos, ypos, zpos, rpos, drr) &
        !$OMP SCHEDULE(DYNAMIC,1)
        do kk = lo(3), hi(3)
           dzz = dx(3)/rr
           zz = (kk + HALF)*dzz
           
           do jj = lo(2), hi(2)
              dyy = dx(2)/rr
              yy = (jj + HALF)*dyy
              
              do ii = lo(1), hi(1)
                 dxx = dx(1)/rr
                 xx = (ii + HALF)*dxx
                 
                 eos_state % rho = p(ii,jj,kk,dens_comp)
                 eos_state % T = p(ii,jj,kk,temp_comp)
                 eos_state % xn(:) = p(ii,jj,kk,spec_comp:spec_comp+nspec-1)

                 call eos(eos_input_rt, eos_state)

                 ! The previous 3D loop has already stored the pressure

                 chi_rho = eos_state % rho * eos_state % dPdr / eos_state % p
                 chi_t   = eos_state % T * eos_state % dPdT / eos_state % p

                 ! calculate adiabatic gradient
                 ap(ii,jj,kk,adic_comp) = (eos_state % gam1 - chi_rho) / (eos_state % gam1 * chi_t)

                 ! calculate actual gradient
                 ! Here the derivative dt/dp will be calculated radially
                 ! by determining dt/dr and dp/dr and dividing (dt/dr)/(dp/dr)
                 ! The derivative dX/dp will be similarly computed.

                 ! calculate gradients along x
                 if (ii==lo(1)) then
                    ! forward difference
                    dtdxx = p(ii+1,jj,kk,temp_comp) - p(ii,jj,kk,temp_comp)
                    dpdxx = pres(ii+1,jj,kk) - pres(ii,jj,kk)
                    dxdxx(:) = p(ii+1,jj,kk,spec_comp:spec_comp+nspec-1) - p(ii,jj,kk,spec_comp:spec_comp+nspec-1)
                    dpdxx = dpdxx/dxx                    
                    dtdxx = dtdxx/dxx                    
                    dxdxx(:) = dxdxx(:)/dxx
                 else if (ii == hi(1)) then
                    ! backward difference
                    dtdxx = p(ii,jj,kk,temp_comp) - p(ii-1,jj,kk,temp_comp)
                    dpdxx = pres(ii,jj,kk) - pres(ii-1,jj,kk)
                    dxdxx(:) = p(ii,jj,kk,spec_comp:spec_comp+nspec-1) - p(ii-1,jj,kk,spec_comp:spec_comp+nspec-1)
                    dpdxx = dpdxx/dxx                    
                    dtdxx = dtdxx/dxx                    
                    dxdxx(:) = dxdxx(:)/dxx
                 else
                    ! centered difference
                    dtdxx = p(ii+1,jj,kk,temp_comp) - p(ii-1,jj,kk,temp_comp)
                    dpdxx = pres(ii+1,jj,kk) - pres(ii-1,jj,kk)
                    dxdxx(:) = p(ii+1,jj,kk,spec_comp:spec_comp+nspec-1) - p(ii-1,jj,kk,spec_comp:spec_comp+nspec-1)
                    dpdxx = HALF*dpdxx/dxx
                    dtdxx = HALF*dtdxx/dxx
                    dxdxx(:) = HALF*dxdxx(:)/dxx
                 endif
                 
                 ! calculate gradients along y
                 if (jj==lo(2)) then
                    ! forward difference
                    dtdyy = p(ii,jj+1,kk,temp_comp) - p(ii,jj,kk,temp_comp)
                    dpdyy = pres(ii,jj+1,kk) - pres(ii,jj,kk)
                    dxdyy(:) = p(ii,jj+1,kk,spec_comp:spec_comp+nspec-1) - p(ii,jj,kk,spec_comp:spec_comp+nspec-1)
                    dpdyy = dpdyy/dyy                    
                    dtdyy = dtdyy/dyy                    
                    dxdyy(:) = dxdyy(:)/dyy
                 else if (jj == hi(2)) then
                    ! backward difference
                    dtdyy = p(ii,jj,kk,temp_comp) - p(ii,jj-1,kk,temp_comp)
                    dpdyy = pres(ii,jj,kk) - pres(ii,jj-1,kk)
                    dxdyy(:) = p(ii,jj,kk,spec_comp:spec_comp+nspec-1) - p(ii,jj-1,kk,spec_comp:spec_comp+nspec-1)
                    dpdyy = dpdyy/dyy                    
                    dtdyy = dtdyy/dyy                    
                    dxdyy(:) = dxdyy(:)/dyy
                 else
                    ! centered difference
                    dtdyy = p(ii,jj+1,kk,temp_comp) - p(ii,jj-1,kk,temp_comp)
                    dpdyy = pres(ii,jj+1,kk) - pres(ii,jj-1,kk)
                    dxdyy(:) = p(ii,jj+1,kk,spec_comp:spec_comp+nspec-1) - p(ii,jj-1,kk,spec_comp:spec_comp+nspec-1)
                    dpdyy = HALF*dpdyy/dyy
                    dtdyy = HALF*dtdyy/dyy
                    dxdyy(:) = HALF*dxdyy(:)/dyy
                 endif

                 ! calculate gradients along z
                 if (kk==lo(3)) then
                    ! forward difference
                    dtdzz = p(ii,jj,kk+1,temp_comp) - p(ii,jj,kk,temp_comp)
                    dpdzz = pres(ii,jj,kk+1) - pres(ii,jj,kk)
                    dxdzz(:) = p(ii,jj,kk+1,spec_comp:spec_comp+nspec-1) - p(ii,jj,kk,spec_comp:spec_comp+nspec-1)
                    dpdzz = dpdzz/dzz                    
                    dtdzz = dtdzz/dzz                    
                    dxdzz(:) = dxdzz(:)/dzz
                 else if (kk == hi(3)) then
                    ! backward difference
                    dtdzz = p(ii,jj,kk,temp_comp) - p(ii,jj,kk-1,temp_comp)
                    dpdzz = pres(ii,jj,kk) - pres(ii,jj,kk-1)
                    dxdzz(:) = p(ii,jj,kk,spec_comp:spec_comp+nspec-1) - p(ii,jj,kk-1,spec_comp:spec_comp+nspec-1)
                    dpdzz = dpdzz/dzz                    
                    dtdzz = dtdzz/dzz                    
                    dxdzz(:) = dxdzz(:)/dzz
                 else
                    ! centered difference
                    dtdzz = p(ii,jj,kk+1,temp_comp) - p(ii,jj,kk-1,temp_comp)
                    dpdzz = pres(ii,jj,kk+1) - pres(ii,jj,kk-1)
                    dxdzz(:) = p(ii,jj,kk+1,spec_comp:spec_comp+nspec-1) - p(ii,jj,kk-1,spec_comp:spec_comp+nspec-1)
                    dpdzz = HALF*dpdzz/dzz
                    dtdzz = HALF*dtdzz/dzz
                    dxdzz(:) = HALF*dxdzz(:)/dzz
                 endif

                 ! calculate position relative to center of star
                 xpos = xx - xctr
                 ypos = yy - yctr
                 zpos = zz - zctr
                 rpos = sqrt(xpos**2 + ypos**2 + zpos**2)
                 drr = sqrt(dxx**2 + dyy**2 + dzz**2)

                 ! now calculate radial gradients
                 dtdrr = dtdxx * xpos/rpos + dtdyy * ypos/rpos + dtdzz * zpos/rpos
                 dpdrr = dpdxx * xpos/rpos + dpdyy * ypos/rpos + dpdzz * zpos/rpos
                 dxdrr = dxdxx(:) * xpos/rpos + dxdyy(:) * ypos/rpos + dxdzz(:) * zpos/rpos                 

                 ! actual gradient
                 if (p(ii,jj,kk,dens_comp) <= low_cutoff .or. abs(dpdrr*drr) <= small) then
                    nabla = ZERO
                    dXdP  = ZERO
                 else
                    nabla = pres(ii,jj,kk) * dtdrr / (dpdrr * p(ii,jj,kk,temp_comp))
                    dXdP(:) = dxdrr(:)/dpdrr
                 endif

                 ap(ii,jj,kk,actl_comp) = nabla

                 ! calculate ledoux gradient
                 
                 ! initialize ledoux nabla to adiabatic nabla
                 nabla = ap(ii,jj,kk,adic_comp)

                 ! add species contributions
                 chi_X(:) = eos_state % xn(:) * eos_state % dPdX(:) / (eos_state % p * chi_t)

                 do n = 1, nspec
                    if (p(ii,jj,kk,spec_comp+n-1) > ZERO) then
                       nabla = nabla - &
                            chi_X(n) * pres(ii,jj,kk) * dXdP(n) &
                            / p(ii,jj,kk,spec_comp+n-1)
                    endif
                 enddo

                 ap(ii,jj,kk,ledx_comp) = nabla
                 
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

        call fab_unbind(pf, i,j)

     enddo
  enddo
  
  pd = plotfile_get_pd_box(pf,1)
  prob_lo(1:dim) = pf%plo
  prob_hi(1:dim) = pf%phi
  time = plotfile_time(pf)

  call fabio_ml_multifab_write_d(conv_grad, ref_ratio, trim(outputfile), &
                                 component_names, pd, &
                                 prob_lo(1:dim), prob_hi(1:dim), time, &
                                 dx(1:dim))

  call destroy(pf)
  deallocate(pres)
  deallocate(ref_ratio)
  deallocate(conv_grad)

contains
  subroutine print_usage()
    implicit none

    print *, 'This program takes a 3-d cartesian plotfile and calculates'
    print *, 'the actual, adiabatic, and ledoux thermodynamic gradients.'
    print *, ''
    print *, 'The actual thermodynamic gradient is computed along the'
    print *, 'radial direction, constructed from gradients along'
    print *, 'the cartesian coordinate directions.'
    print *, ''
    print *, ' calling sequence: '
    print *, ' fconv_radial -i <inputfile> [-o <outputfile>, -l <cutoff>]'
    print *, ''
    print *, ' arguments:'
    print *, '-i | --input:'
    print *, '     specify the input plotfile used to calculate the adiabatic'
    print *, '     excess.  required.'
    print *, '-o | --output:'
    print *, '     specify the output filename for the adiabatic excess.'
    print *, '     defaults to "conv_radial"'
    print *, '-l | --low_cutoff:'
    print *, '     specifies the low density cutoff.  for densities below the'
    print *, '     low density cutoff, the actual thermal gradient is set to'
    print *, '     ZERO.  defaults to 1.e4.'
    print *, ''

  end subroutine print_usage

end program fconv_radial
