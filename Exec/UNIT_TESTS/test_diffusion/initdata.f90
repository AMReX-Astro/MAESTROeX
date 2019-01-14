module initdata_module

  use parallel, only: parallel_IOProcessor
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, &
       spec_comp, pi_comp, prob_lo
  use probin_module, only: ambient_h, ambient_dens, &
       ambient_he4, ambient_c12, ambient_fe56, peak_h, t0, diffusion_coefficient
  use eos_module
  use eos_type_module
  use extern_probin_module, only: const_conductivity
  use amrex_constants_module
  use fill_3d_data_module, only: put_1d_array_on_cart

  implicit none

  private

contains

  subroutine initdata(lev, time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, dx) bind(C, name="initdata")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), 1:nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), 1:nc_v)
    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k,r,iter
    integer, parameter :: max_iter = 50
    double precision, parameter :: tol = 1.d-12
    double precision :: x,y,dist2,dens_zone,temp_zone,del_dens, del_temp
    double precision :: pres_zone, del_pres
    double precision :: h_zone, dhdt
    logical :: converged
    double precision :: dens_pert, rhoh_pert, temp_pert
    double precision :: rhoX_pert(nspec)
    type (eos_t) :: eos_state
    double precision :: s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nscal)

    ! set velocity to zero
    vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3),1:nc_v) = ZERO

    scal(scal_lo(1):scal_hi(1), scal_lo(2):scal_hi(2), scal_lo(3):scal_hi(3), 1:nc_s) = ZERO

    do i = 1, nscal
        call put_1d_array_on_cart(lev, lo, hi, s0_cart(:,:,:,i), lo, hi, 1, s0_init(:,:,i), 0, 0)
    enddo

    ! initialize the scalars
    do k=lo(3),hi(3)
       do j = lo(2), hi(2)

          y = prob_lo(2) + (dble(j)+HALF) * dx(2)

          do i = lo(1), hi(1)

             x = prob_lo(1) + (dble(i)+HALF) * dx(1)

             ! apply the guassian enthalpy pulse at constant density
             dist2 = (center(1) - x)**2 + (center(2) - y)**2

             h_zone = (peak_h - ambient_h) * &
                  exp(-dist2/(4.0d0*diffusion_coefficient*t0)) + ambient_h

             temp_zone = s0_cart(i,j,k,temp_comp)

             eos_state%xn(1:nspec) = s0_cart(i,j,k,spec_comp:spec_comp+nspec-1) / &
                  s0_cart(i,j,k,rho_comp)

             eos_state%rho = s0_cart(i,j,k,rho_comp)

             converged = .false.

             eos_state%T = temp_zone

             do iter = 1, max_iter
                eos_state%T = temp_zone

                call eos(eos_input_rt, eos_state)

                dhdt = eos_state%cv + eos_state%dpdt/eos_state%rho

                del_temp = -(eos_state%h - h_zone) / dhdt

                temp_zone = temp_zone + del_temp

                if (abs(del_temp) < tol*temp_zone) then
                   converged = .true.
                   exit
                endif
             enddo

             if (.not. converged) &
                  call bl_error("iters did not converge in initscalars")

             ! call eos one last time
             eos_state%T = temp_zone

             call eos(eos_input_rt, eos_state)

             scal(i,j,k,rho_comp)  = eos_state%rho
             scal(i,j,k,rhoh_comp) = eos_state%rho * eos_state%h
             scal(i,j,k,temp_comp) = eos_state%T
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  eos_state%xn(1:nspec) * eos_state%rho

          enddo
       enddo
    enddo

  end subroutine initdata

end module initdata_module
