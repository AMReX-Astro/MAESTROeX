program furcashell

  use network
  use table_rates, only: j_na23_ne23, j_ne23_na23, get_entries, &
                         table_meta, jtab_rate, jtab_nuloss
  use actual_rhs_module, only: rate_eval_t, evaluate_rates, rhs_nuc, actual_rhs
  use physical_constants, only: N_AVO
  use sneut_module, only: sneut5
  use burner_module, only: burner_init
  use burn_type_module, only: burn_t, burn_to_eos, eos_to_burn, net_ienuc
  use eos_module, only: eos, eos_init
  use eos_type_module, only: eos_t, eos_input_rt

  use extern_probin_module, only: table_interpolation_scheme
  
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
  integer :: uin, dim, nlevs, i, j, ii, jj, kk, itab
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
  double precision :: y_nuc(nspec), ydot_nuc(nspec), enuc, enucion, enucdiff
  double precision :: table_entries(7) ! Hard-coded for Urca tabulated rates

  integer, parameter :: size_rate_eval = 14

  ! For sneut5
  double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz

  character(len=20) :: plot_names(size_rate_eval)

  ! For AMReX, disable nested parallel regions
  if (omp_get_max_threads() > 1) call omp_set_nested(.false.)

  uin = unit_new()

  ! defaults
  pltfile =''
  outputfile = 'urca'
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
     case ('--bicubic')
        table_interpolation_scheme = 3
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

  call network_init()
  call burner_init()
  call eos_init()
  
  ! build the input plotfile
  call build(pf,pltfile,uin)

  nlevs = plotfile_nlevels(pf)
  dim = plotfile_dim(pf)

  allocate(rr(nlevs,dim),rates(nlevs))
  rr = plotfile_refrat(pf)

  plot_names(1) = "ecap23"
  plot_names(2) = "beta23"
  plot_names(3) = "epart_ecap23"
  plot_names(4) = "epart_beta23"
  plot_names(5) = "dqweak_ecap23"
  plot_names(6) = "dqweak_beta23"
  plot_names(7) = "X(na23)"
  plot_names(8) = "X(ne23)"
  plot_names(9) = "density"
  plot_names(10) = "ionenuc"
  plot_names(11) = "sneut"
  plot_names(12) = "enucrdiff"
  plot_names(13) = "enu_ecap23"
  plot_names(14) = "enu_beta23"

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
        !$OMP PRIVATE(enucion, enuc, enucdiff, y_nuc, ydot_nuc) &
        !$OMP PRIVATE(sneut, dsneutdt, dsneutdd, snuda, snudz) &
        !$OMP PRIVATE(itab, table_entries) &
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

                 y_nuc(:) = burn_state % xn(:) * aion_inv(:)
                 call rhs_nuc(ydot_nuc, y_nuc, rate_state % screened_rates, burn_state % rho)
                 burn_state % ydot(1:nspec) = ydot_nuc
                 call ener_gener_rate(ydot_nuc, enucion)
                 enuc = enucion
                 
                 ! Electron capture rate (A=23)
                 r(ii,jj,kk,1) = rate_state % screened_rates(k_na23__ne23)

                 ! Beta decay rate (A=23)
                 r(ii,jj,kk,2) = rate_state % screened_rates(k_ne23__na23)

                 ! Particle energy from electron capture (A=23) (erg/g/s)
                 r(ii,jj,kk,3) = rate_state % epart(j_na23_ne23) * y_nuc(jna23) * N_AVO
                 enuc = enuc + r(ii,jj,kk,3)

                 ! Particle energy from beta decay (A=23) (erg/g/s)
                 r(ii,jj,kk,4) = rate_state % epart(j_ne23_na23) * y_nuc(jne23) * N_AVO
                 enuc = enuc + r(ii,jj,kk,4)
                 
                 ! dQ energy correction from electron capture (A=23) (erg/g/s)
                 r(ii,jj,kk,5) = rate_state % dqweak(j_na23_ne23) * ydot_nuc(jna23) * N_AVO
                 enuc = enuc + r(ii,jj,kk,5)
                 
                 ! dQ energy correction from beta decay (A=23) (erg/g/s)
                 r(ii,jj,kk,6) = rate_state % dqweak(j_ne23_na23) * ydot_nuc(jne23) * N_AVO
                 enuc = enuc + r(ii,jj,kk,6)
                 
                 ! Mass fraction of Na-23
                 r(ii,jj,kk,7) = burn_state % xn(jna23)

                 ! Mass fraction of Ne-23
                 r(ii,jj,kk,8) = burn_state % xn(jne23)

                 ! Density
                 r(ii,jj,kk,9) = burn_state % rho

                 ! Raw Ion binding energy enucdot (erg/g/s) (without dQ energy corrections)
                 r(ii,jj,kk,10) = enucion

                 ! non-Urca neutrino losses: sneut
                 call sneut5(burn_state % T, burn_state % rho, &
                             burn_state % abar, burn_state % zbar, &
                             sneut, dsneutdt, dsneutdd, snuda, snudz)
                 r(ii,jj,kk,11) = sneut
                 enuc = enuc - sneut

                 call actual_rhs(burn_state)

                 enucdiff = burn_state % ydot(net_ienuc) - enuc
                 r(ii,jj,kk,12) = abs(enucdiff/burn_state % ydot(net_ienuc))

                 ! Average neutrino energy for electron capture (A=23)
                 call get_entries(table_meta(j_na23_ne23), &
                                  burn_state % rho * burn_state % y_e, &
                                  burn_state % T, table_entries)
                 if (table_entries(jtab_rate) == 0.0d0) then
                    r(ii,jj,kk,13) = 0.0d0
                 else
                    r(ii,jj,kk,13) = table_entries(jtab_nuloss)/table_entries(jtab_rate)
                 end if
                 if (isnan(r(ii,jj,kk,13))) then
                    print *, "got NaN in enu_ecap23 with energy loss ", table_entries(jtab_nuloss), &
                             " and rate ", table_entries(jtab_rate)
                 else if (r(ii,jj,kk,13) < 0.0d0) then
                    print *, "got negative energy in enu_ecap23 with energy loss ", table_entries(jtab_nuloss), &
                             " and rate ", table_entries(jtab_rate)
                 end if
                 
                 ! Average neutrino energy for beta decay (A=23)
                 call get_entries(table_meta(j_ne23_na23), &
                                  burn_state % rho * burn_state % y_e, &
                                  burn_state % T, table_entries)
                 if (table_entries(jtab_rate) == 0.0d0) then
                    r(ii,jj,kk,14) = 0.0d0
                 else
                    r(ii,jj,kk,14) = table_entries(jtab_nuloss)/table_entries(jtab_rate)
                 end if
                 if (isnan(r(ii,jj,kk,14))) then
                    print *, "got NaN in enu_beta23 with energy ", table_entries(jtab_nuloss), &
                             " and rate ", table_entries(jtab_rate)
                 else if (r(ii,jj,kk,14) < 0.0d0) then
                    print *, "got negative energy in enu_beta23 with energy loss ", table_entries(jtab_nuloss), &
                             " and rate ", table_entries(jtab_rate)
                 end if

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
    print *, "This program takes a 3D plotfile and extracts the electron capture, "
    print *, " beta decay rate, mass fractions, ... "
    print *, " particle energy, and other neutrino losses -- "
    print *, " then dumps a new plotfile containing these quantities."
    print *, ""
    print *, "This is set up for the URCA-simple network in StarKiller Microphysics."
    print *, ""
    print *, "usage: "
    print *, " *furcashell* -i|--input <pltfile in> [-o|--output <pltfile out>]"
    print *, ""
    print *, "    -i|--input: <pltfile in>"
    print *, "        Specify which plotfile to work on. (Required)"
    print *, "    -o|--output:"
    print *, "        Name of the out new plotfile to create. (Default: 'urcashell')"
    print *, "    --tfromp:"
    print *, "        Toggles the use of 'temperature' to be tfromp instead", &
         " of tfromh."
    print *, "        (Default: use tfromh)"
    print *, "    --bicubic:"
    print *, "        Toggles the use of bicubic interpolation in weak rate tables. "
    print *, "        (Default: bilinear interpolation)"
    print *, ""

  end subroutine print_usage

end program furcashell
