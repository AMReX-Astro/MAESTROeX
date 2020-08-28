module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use model_parser_module
  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level, center
  use meth_params_module, only: nscal, model_file, rho_comp, prob_lo, spherical, &
       heating_cutoff_density_lo, heating_cutoff_density_hi
  use probin_module, only: use_analytic_heating
  
  implicit none

  private

contains

  subroutine make_heating(lo, hi, &
                          rho_Hext, r_lo, r_hi, &
                          scal,     s_lo, s_hi, nc_s, &
                          dx, time ) bind (C,name="make_heating")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3), nc_s
    double precision, intent (inout) :: rho_Hext(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)     )
    double precision, intent (in   ) ::     scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in   ) :: dx(3), time

    integer :: n, r
    double precision :: rloc, starting_rad, xloc(3), rho
    integer :: i, j, k

    if (.not. model_initialized) then 
        call read_model_file(model_file)
    endif

    ! if (use_analytic_heating) then

    if (.not. spherical) then
        starting_rad = prob_lo(AMREX_SPACEDIM)
    else
        starting_rad = ZERO
    endif

    do k=lo(3),hi(3)
        do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                xloc(1) = prob_lo(1) + (dble(i)+0.5d0)*dx(1) - center(1)
                xloc(2) = prob_lo(2) + (dble(j)+0.5d0)*dx(2) - center(2)
                xloc(3) = prob_lo(3) + (dble(k)+0.5d0)*dx(3) - center(3)

                if (AMREX_SPACEDIM .eq. 2) then
                    rloc = xloc(2)
                else if (AMREX_SPACEDIM .eq. 3) then
                    if (.not. spherical) then
                        rloc = xloc(3)
                    else
                        ! compute distance to the center of the star
                        rloc = ZERO
                        do n=1,3
                            rloc = rloc + xloc(n)**2
                        enddo
                        rloc = sqrt(rloc)
                    end if
                end if


                rho = scal(i,j,k,rho_comp)
                if (rho > heating_cutoff_density_lo .and. rho < heating_cutoff_density_hi) then
                    rho_Hext(i,j,k)  = interpolate(rloc, ienuc_model) * scal(i,j,k,rho_comp)
                else
                    rho_Hext(i,j,k) = 0.d0
                end if
                
            end do
        end do
    end do

    ! endif

  end subroutine make_heating

end module make_heating_module

