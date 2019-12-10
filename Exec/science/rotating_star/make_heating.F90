module make_heating_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use model_parser_module
  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use meth_params_module, only: nscal, model_file, rho_comp, prob_lo, spherical
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
    double precision :: rloc,model_dr,rmax,starting_rad,mod_dr
    integer :: i, j, k
    double precision :: enuc(0:max_radial_level,0:nr_fine-1)

    if (use_analytic_heating) then

        if (spherical .eq. 0) then
            starting_rad = prob_lo(AMREX_SPACEDIM)
        else
            starting_rad = ZERO
        endif

        do n=0,max_radial_level
            do r=0,nr(n)-1
                rloc = starting_rad + (dble(r) + HALF)*dr(n)
                enuc(n,r) = interpolate(rloc, ienuc_model)
            end do
        enddo

        do k=lo(3),hi(3)
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)

                        if (AMREX_SPACEDIM .eq. 2) then
                            r = j
                        else if (AMREX_SPACEDIM .eq. 3) then
                            r = k
                        end if

                        rho_Hext(i,j,k)  = enuc(0, r) * scal(i,j,k,rho_comp)
                end do
            end do
        end do  

    endif
    
  end subroutine make_heating

end module make_heating_module

