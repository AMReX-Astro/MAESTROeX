program frates

  use network
  use actual_rhs_module, only: rate_eval_t, evaluate_rates
  use burner_module, only: burner_init
  use burn_type_module, only: burn_t, burn_to_eos, eos_to_burn
  use eos_module, only: eos, eos_init
  use eos_type_module, only: eos_t, eos_input_rt
  use rpar_indices
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
  character(len=256) :: pltfile, outputfile
  logical :: use_tfromp

  ! f2kcli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local variables
  integer :: uin, dim, nlevs, i, j, ii, jj, kk
  type(plotfile) :: pf
  type(layout) :: la
  type(boxarray) :: ba
  type(list_box) :: bl
  type(box) :: bx,pd
  type(multifab), allocatable :: rates(:)
  integer :: dens_comp, temp_comp, spec_comp
  integer, dimension(MAX_SPACEDIM) :: lo, hi
  integer, allocatable :: rr(:,:)
  real(kind=dp_t), pointer :: p(:,:,:,:), r(:,:,:,:)

  type (burn_t) :: burn_state
  type (eos_t)  :: eos_state
  type (rate_eval_t) :: rate_state

  integer, parameter :: size_rate_eval = 5*nrates + 2*nrat_tabular

  character(len=20) :: plot_names(size_rate_eval)

  ! For AMReX, disable nested parallel regions
  if (omp_get_max_threads() > 1) call omp_set_nested(.false.)

  uin = unit_new()

  ! defaults
  pltfile =''
  outputfile = 'rates'
  use_tfromp = .false.
  
  ! parse arguments
  narg = command_argument_count()

  farg = 1
  do while (farg<=narg)
     call get_command_argument(farg,value=fname)

     select case(fname)

     case ('-i', '--input')
        farg = farg + 1
        call get_command_argument(farg,value=pltfile)
     case ('-o', '--output')
        farg = farg + 1
        call get_command_argument(farg,value=outputfile)
     case ('--tfromp')
        use_tfromp = .true.
     case default
        exit
     end select
     farg = farg + 1
  end do

  ! sanity check
  if (pltfile == '') then
     call print_usage()
     stop
  end if

  print *, 'working on pltfile: ', trim(pltfile)
  print *, 'saving to pltfile: ', trim(outputfile)
  if (use_tfromp) print *, 'using tfromp instead of tfromh'

  call network_init()
  call eos_init()
  call burner_init()  
  
  ! build the input plotfile
  call build(pf,pltfile,uin)

  nlevs = plotfile_nlevels(pf)
  dim = plotfile_dim(pf)

  allocate(rr(nlevs,dim),rates(nlevs))
  rr = plotfile_refrat(pf)

  do i = 1, nrates
     plot_names(i) = "screened_rate_" // int_to_str(i)
     plot_names(i+nrates) = "rate_" // int_to_str(i)
     plot_names(i+2*nrates) = "drate_dt_" // int_to_str(i)
     plot_names(i+3*nrates) = "screen_" // int_to_str(i)
     plot_names(i+4*nrates) = "dscreen_dt_" // int_to_str(i)
  end do

  do i = 1, nrat_tabular
     plot_names(i+5*nrates) = "add_energy_" // int_to_str(i)
     plot_names(i+5*nrates+nrat_tabular) = "add_energy_rate_" // int_to_str(i)
  end do

  dens_comp = plotfile_var_index(pf,"density")
  if (use_tfromp) then
     temp_comp = plotfile_var_index(pf,"tfromp")
  else
     temp_comp = plotfile_var_index(pf,"tfromh")
  end if
  spec_comp = plotfile_var_index(pf,"X(" // trim(short_spec_names(1)) // ")")

  if (dens_comp < 0 .or. spec_comp < 0 .or. temp_comp < 0) then
     print *, dens_comp, temp_comp, spec_comp
     call bl_error("Variables not found")
  endif

  do i = 1, nlevs
     do j = 1, nboxes(pf,i)
        call push_back(bl,get_box(pf,i,j))
     end do

     call build(ba,bl)
     call build(la,ba,plotfile_get_pd_box(pf,i))
     call destroy(bl)
     call destroy(ba)

     ! build the multifab with 0 ghost cells and size_rate_eval components
     call multifab_build(rates(i),la,size_rate_eval,0)
  end do

  ! loop over the plotfile data starting at the finest
  do i = nlevs, 1, -1
     ! loop over each box at this level
     do j = 1, nboxes(pf,i)
        ! read in the data 1 patch at a time
        call fab_bind(pf,i,j)

        lo(1:dim) = lwb(get_box(pf,i,j))
        hi(1:dim) = upb(get_box(pf,i,j))

        p => dataptr(pf,i,j)
        r => dataptr(rates(i),j)

        !$OMP PARALLEL DO PRIVATE(kk, jj, ii, eos_state, burn_state, rate_state) &
        !$OMP SCHEDULE(DYNAMIC,1)
        do kk = lo(3), hi(3)
           do jj = lo(2), hi(2)
              do ii = lo(1), hi(1)

                 eos_state % rho = p(ii,jj,kk,dens_comp)
                 eos_state % T   = p(ii,jj,kk,temp_comp)
                 eos_state % xn  = p(ii,jj,kk,spec_comp:spec_comp+nspec-1)

                 call eos(eos_input_rt, eos_state)
                 call eos_to_burn(eos_state, burn_state)

                 call evaluate_rates(burn_state, rate_state)

                 call flatten_rate_state(r(ii,jj,kk,:), rate_state)

              end do
           end do
        end do
        !$OMP END PARALLEL DO

        call fab_unbind(pf,i,j)
        
     end do
  end do

  call fabio_ml_multifab_write_d(rates, rr(:,1), &
                                 trim(outputfile), &
                                 plot_names, plotfile_get_pd_box(pf,1), &
                                 pf%plo, pf%phi, plotfile_time(pf), &
                                 plotfile_get_dx(pf,1))
  call destroy(pf)


contains

  subroutine print_usage()
    implicit none
    
    print *,""
    print *, "This program takes a 3D plotfile, calls evaluate_rates in actual_rhs_module, and "
    print *, "dumps a new plotfile containing the reaction rates and screening factors in each zone."
    print *, ""
    print *, "This is set up for the URCA-simple network in StarKiller Microphysics."
    print *, ""
    print *, "usage: "
    print *, " *frates* -i|--input <pltfile in> [-o|--output <pltfile out>]", &
         " [--tfromp]"    
    print *, ""
    print *, "    -i|--input: <pltfile in>"
    print *, "        Specify which plotfile to work on. (Required)"
    print *, "    -o|--output:"
    print *, "        Name of the out new plotfile to create. (Default: 'rates')"
    print *, "    --tfromp:"
    print *, "        Toggles the use of 'temperature' to be tfromp instead", &
         " of tfromh."
    print *, "        (Default: use tfromh)"
    print *, ""

  end subroutine print_usage

  function int_to_str(i) result(s)
    implicit none
    
    character(:), allocatable :: s
    integer, intent(in) :: i
    character(range(i)+2) :: scratch

    write(scratch,'(i0)') i
    s = trim(scratch)
  end function int_to_str

  subroutine flatten_rate_state(rflat, rate_state)
    implicit none

    real (kind=dp_t), intent(inout) :: rflat(:)
    type (rate_eval_t), intent(in) :: rate_state

    integer :: i

    do i = 1, nrates
       rflat(i) = rate_state % screened_rates(i)
       rflat(i +   nrates) = rate_state % unscreened_rates(1, i)
       rflat(i + 2*nrates) = rate_state % unscreened_rates(2, i)
       rflat(i + 3*nrates) = rate_state % unscreened_rates(3, i)
       rflat(i + 4*nrates) = rate_state % unscreened_rates(4, i)
    end do

    do i = 1, nrat_tabular
       rflat(i + 5*nrates) = rate_state % add_energy(i)
       rflat(i + 5*nrates + nrat_tabular) = rate_state % add_energy_rate(i)
    end do
  end subroutine flatten_rate_state
  
end program frates
