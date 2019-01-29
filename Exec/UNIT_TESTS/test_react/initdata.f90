module initdata_module

  use parallel, only: parallel_IOProcessor
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, &
       spec_comp, pi_comp, prob_lo, prob_hi
  use probin_module, only:  dens_max, temp_max, dens_min, temp_min
  use eos_module
  use eos_type_module
  use amrex_constants_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  ! use fill_3d_data_module, only: put_1d_array_on_cart

  implicit none

  private
  integer, save :: rho_c, h_c, spec_c, t_c, omegadot_c, hnuc_c, &
       lhnuc_c, hext_c, dxn_con_c, h_con_c, ncomps
  integer, save :: nlevs

  character(len=20), save, allocatable :: varnames(:)

  logical, save :: react_is_init = .false.

contains

  subroutine initdata(lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       domlo, domhi, &
       s0_init, p0_init, dx) bind(C, name="initdata_thermal")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: domlo(3), domhi(3)
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), 1:nc_s)
    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k,xn_hi
    double precision :: temp_zone, dens_zone
    double precision :: dlogrho, dlogT
    double precision, pointer :: xn_zone(:, :)
    type (eos_t) :: eos_state

    scal(scal_lo(1):scal_hi(1), scal_lo(2):scal_hi(2), scal_lo(3):scal_hi(3), 1:nc_s) = ZERO

    xn_hi = domhi(3)-domlo(3)

    call bl_allocate(xn_zone, 1, nspec, 0, xn_hi)

    dlogrho = (log10(dens_max) - log10(dens_min)) / (domhi(1) - domlo(1))
    dlogT = (log10(temp_max) - log10(temp_min)) / (domhi(2) - domlo(2))
    call get_xn(xn_zone, 0, xn_hi)

    ! initialize the scalars
    do k=lo(3),hi(3)
       do j = lo(2), hi(2)
          ! Set the temperature
          temp_zone = 10.d0**(log10(temp_min) + dble(j)*dlogT)

          do i = lo(1), hi(1)
             ! Set the density
             dens_zone = 10.d0**(log10(dens_min) + dble(i)*dlogrho)

             ! Call the EoS w/ rho, temp, & X as inputs
             eos_state%T = temp_zone
             eos_state%rho = dens_zone
             eos_state%xn(1:nspec) = xn_zone(1:nspec,k)

             call eos(eos_input_rt, eos_state)

             ! Initialize this element of the state
             scal(i,j,k,rho_comp)  = dens_zone
             scal(i,j,k,rhoh_comp) = dens_zone * eos_state%h
             scal(i,j,k,temp_comp) = temp_zone
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  dens_zone * xn_zone(1:nspec,k)

          enddo
       enddo
    enddo

    call bl_deallocate(xn_zone)

  end subroutine initdata

  subroutine get_xn(xn_zone, x_lo, x_hi)
    !=== Data ===
    !Modules
    use network,       only: nspec, spec_names
    use probin_module, only: xin_file

    !Args
    integer,         intent(in   ) :: x_lo, x_hi
    double precision, intent(  out) :: xn_zone(1:nspec,x_lo:x_hi)

    !Local data
    integer         :: un
    integer         :: i
    double precision :: summ, usr_in
    character (len=1024) :: line

    !=== Execution ===
    if (xin_file .eq. 'uniform') then
       summ = ZERO
       do i=1, nspec - 1
          print *, trim(adjustl(spec_names(i))) // ' mass fraction: '
          read(*,*) usr_in
          xn_zone(i,:) = usr_in
          summ = summ + usr_in
       enddo

       if(summ > ONE) then
          print *, 'ERROR: Mass fraction sum exceeds 1.0!'
          stop
       else
          xn_zone(nspec,:) = ONE - summ
       endif
    else
       ! read in an inputs file containing the mass fractions.
       ! each species is on its own line.
       ! Allow for comment lines with '#' in the first column
       ! un = unit_new()
       un = 1
       open(unit=un, file=xin_file, status='old', action='read')

       summ = ZERO

       !TODO: Add xn <= 1.0 error checking
       !TODO: Add proper cell count error checking
       i = 1
       print *, 'nspec: ', nspec
       do while (i <= nspec)
          ! read the line into a character buffer
          read (un,'(a)') line

          ! skip comments
          if (index(line, '#') == 1) cycle

          ! parse the line
          read(line,*) xn_zone(i,:)

          i = i + 1
       enddo

       do i = 1, size(xn_zone,1)
          print *, i, sum(xn_zone(:,i))
       enddo

       close(un)
    endif
  end subroutine get_xn

end module initdata_module
