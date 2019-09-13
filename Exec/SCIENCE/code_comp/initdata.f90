
module initdata_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use amrex_constants_module
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim, amrex_random
  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       prob_lo, prob_hi
  use eos_module
  use eos_type_module
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr
  use probin_module, only: rho_0

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
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer :: i, j, k, r
    double precision :: x, y, z, fheat, rhopert

    type(eos_t) :: eos_state

    ! set velocity to zero
    vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_v) = 0.d0

    do k=lo(3),hi(3)
       z = prob_lo(3) + dx(3)*(dble(k) + HALF)
       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j) + HALF)

          if (amrex_spacedim .eq. 2) then
            if (y < 1.125d0 * 4.d8) then
                fheat = sin(8.d0 * M_PI * (y/ 4.d8 - ONE))
            else
                fheat = ZERO
            endif
          else
            if (z < 1.125d0 * 4.d8) then
                fheat = sin(8.d0 * M_PI * (z/ 4.d8 - ONE))
            else
                fheat = ZERO
            endif
          endif

          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1)*(dble(i) + HALF) 

             rhopert = 5.d-5 * rho_0 * fheat * (sin(3.d0 * M_PI * x / 4.d8) + &
                                                cos(M_PI * x / 4.d8)) * &
                     (sin(3 * M_PI * y/4.d8) - cos(M_PI * y/4.d8))

             if (amrex_spacedim .eq. 2) then
                r = j
             else if (amrex_spacedim .eq. 3) then
                r = k
             end if

             eos_state % rho = s0_init(lev,r,rho_comp) + rhopert 
             eos_state % p = p0_init(lev, r)
             eos_state % xn(:) = s0_init(lev,r,spec_comp:spec_comp+nspec-1)

             call eos(eos_input_rp, eos_state)

             ! set scalars using s0
             scal(i,j,k,rho_comp)  = eos_state % rho
             scal(i,j,k,rhoh_comp) = eos_state % rho * eos_state % h
             scal(i,j,k,temp_comp) = eos_state % t
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = eos_state % xn(:)

             ! initialize pi to zero for now
             scal(i,j,k,pi_comp) = 0.d0

          end do
       end do
    end do

  end subroutine initdata

  subroutine initdata_sphr(time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, &
       dx, r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind(C, name="initdata_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc  (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: cc_to_r   (ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

  end subroutine initdata_sphr

end module initdata_module
