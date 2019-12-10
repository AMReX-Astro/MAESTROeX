
module initdata_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
                                prob_lo, prob_hi, grav_const
  use eos_module
  use eos_type_module
  use amrex_constants_module
  use probin_module, only : pert_amp, scale_height, pres_base, k_hoz, k_vert

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

    ! local variables

    integer          :: i,j,k,r
    type (eos_t) :: eos_state
    double precision :: x, y, rho0, rho_local, rho_base 

!    ! abort program
!    call amrex_error()

    ! set velocity to zero
    vel(vel_lo(1):vel_hi(1),vel_lo(2):vel_hi(2),vel_lo(3):vel_hi(3),1:nc_v) = 0.d0

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (amrex_spacedim .eq. 1) then
                r = i
             else if (amrex_spacedim .eq. 2) then
                r = j
             else if (amrex_spacedim .eq. 3) then
                r = k
             end if

             ! set scalars using s0
             scal(i,j,k,rho_comp)  = s0_init(lev,r,rho_comp)
             scal(i,j,k,rhoh_comp) = s0_init(lev,r,rhoh_comp)
             scal(i,j,k,temp_comp) = s0_init(lev,r,temp_comp)
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  s0_init(lev,r,spec_comp:spec_comp+nspec-1)

             ! initialize pi to zero for now
             scal(i,j,k,pi_comp) = 0.d0

          end do
       end do
    end do

    ! introduce density fluctuations - only implemented for 2D

    if ((amrex_spacedim == 1) .or. amrex_spacedim == 3) then 
      call amrex_error("density perturbation only implemented for 2d")
    endif

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      r = j ! this loop is currently 2d only
      x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
      y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)


      ! This seems to work ok with sealed box? 

!      rho0 = s0_init(lev,r,rho_comp)
!      rho_local = rho0 * (1.0 + &
!        pert_amp * exp(-y/(2.0_num * scale_height)) * &
!        cos(x * k_hoz * M_PI / (prob_hi(1) - prob_lo(1)) ) * &
!        sin(y * k_vert * M_PI / (prob_hi(2) - prob_lo(2)) )  &
!      )

      ! What about this ?
 
      rho_base = pres_base / scale_height / abs(grav_const)
      rho0 = s0_init(lev,r,rho_comp)
      rho_local = rho0 &
        + rho_base * pert_amp * exp(-y/(2d0 * scale_height)) * &
        cos(x * k_hoz * m_pi / (prob_hi(1) - prob_lo(1)) ) * &
        sin(y * k_vert * m_pi / (prob_hi(2) - prob_lo(2)) ) 

      ! if k_vert is 0, dont multiply by sin(0) = 0
      ! (make this all more elegant later)
      
      if (k_vert == 0d0) then
        rho_local = rho0 &
          + rho_base * pert_amp * exp(-y/(2d0 * scale_height)) * &
          cos(x * k_hoz * m_pi / (prob_hi(1) - prob_lo(1)) )
      endif

      !
      eos_state%rho   = rho_local 
      eos_state%p     = p0_init(lev,r)
      eos_state%T     = s0_init(lev,r,temp_comp) 
      eos_state%xn(:) = 1.0d0 ! single fluid
      
      ! (rho,p) --> T, h
      call eos(eos_input_rp, eos_state)
      
      scal(i,j,k,rho_comp)  = eos_state%rho
      scal(i,j,k,rhoh_comp) = eos_state%rho * eos_state%h
      scal(i,j,k,temp_comp) = eos_state%T
      scal(i,j,k,spec_comp:spec_comp+nspec-1) = eos_state%rho * eos_state%xn ! single fluid

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

    call amrex_error("error: called initdata_sphr, but is a stub")

  end subroutine initdata_sphr

end module initdata_module
