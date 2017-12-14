
module initdata_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
                                perturb_model, pert_temp_factor, pert_rad_factor,&
                                do_small_domain, prob_lo
  
  implicit none

  private

contains

  subroutine initdata(lev, time, lo, hi, &
                      scal, scal_lo, scal_hi, nc_s, &
                      vel, vel_lo, vel_hi, nc_v, &
                      s0_init, p0_init, &
                      dx) bind(C, name="initdata")
    
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
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k,r
    double precision :: x,y,z

    double precision :: dens_pert, rhoh_pert, temp_pert
    double precision :: rhoX_pert(nspec)

    ! set velocity to zero
    vel = 0.d0

    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)

          if (amrex_spacedim .eq. 2) then
             r = j
          else
             r = k
          end if

          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

             scal(i,j,k,rho_comp)  = s0_init(lev,r,rho_comp)
             scal(i,j,k,rhoh_comp) = s0_init(lev,r,rhoh_comp)
             scal(i,j,k,temp_comp) = s0_init(lev,r,temp_comp)
             scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  s0_init(lev,r,spec_comp:spec_comp+nspec-1)

             ! initialize pi_comp to zero for now
             scal(i,j,k,pi_comp) = 0.d0

          end do
       end do
    end do

    ! add an optional perturbation
    if (perturb_model) then

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
          do j = lo(2), hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
             do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                if (amrex_spacedim .eq. 2) then
                   call perturb_2d(x, y, p0_init(lev,j), s0_init(lev,j,:), &
                        dens_pert, rhoh_pert, rhoX_pert, temp_pert)
                else
                   call perturb_3d(x, y, z, p0_init(lev,k), s0_init(lev,k,:), &
                                   dens_pert, rhoh_pert, rhoX_pert, temp_pert)
                end if

                scal(i,j,k,rho_comp) = dens_pert
                scal(i,j,k,rhoh_comp) = rhoh_pert
                scal(i,j,k,temp_comp) = temp_pert
                scal(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)

             end do
          end do
       end do

    end if

  end subroutine initdata

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: x, y
    double precision, intent(in ) :: p0_init, s0_init(1:nscal)
    double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
    double precision, intent(out) :: rhoX_pert(1:nspec)

    double precision :: temp,t0
    double precision :: x1, y1, r1, x2, y2, r2, x3, y3, r3, x4, y4, r4

    type (eos_t) :: eos_state

    t0 = s0_init(temp_comp)

    x1 = 5.0d7
    y1 = 6.5d7
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / (2.5d6*pert_rad_factor)

    x2 = 1.2d8
    y2 = 8.5d7
    r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / (2.5d6*pert_rad_factor)

    x3 = 2.0d8
    y3 = 7.5d7
    r3 = sqrt( (x-x3)**2 +(y-y3)**2 ) / (2.5d6*pert_rad_factor)

    ! this is a tiny bubble for inputs_2d_smalldomain
    x4 = 0.5d0
    y4 = 85218750.25d0
    r4 = sqrt( (x-x4)**2 +(y-y4)**2 ) / (2.5d-2*pert_rad_factor)

    if (do_small_domain) then
       temp = t0 * (1.d0 + pert_temp_factor* &
            (0.150d0 * (1.d0 + tanh(2.d0-r1)) + &
             0.300d0 * (1.d0 + tanh(2.d0-r2)) + &
             0.225d0 * (1.d0 + tanh(2.d0-r3)) + &
             0.300d0 * (1.d0 + tanh(2.d0-r4))))
    else
       temp = t0 * (1.d0 + pert_temp_factor* &
            (0.150d0 * (1.d0 + tanh(2.d0-r1)) + &
             0.300d0 * (1.d0 + tanh(2.d0-r2)) + &
             0.225d0 * (1.d0 + tanh(2.d0-r3))))
    end if

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    double precision, intent(in ) :: x, y, z
    double precision, intent(in ) :: p0_init, s0_init(1:nscal)
    double precision, intent(out) :: dens_pert, rhoh_pert, temp_pert
    double precision, intent(out) :: rhoX_pert(1:nspec)

    double precision :: temp, t0
    double precision :: x0, y0, z0, r0

    type (eos_t) :: eos_state

    t0 = s0_init(temp_comp)

    x0 = 3.6d7
    y0 = 3.6d7
    z0 = 8.5d7

    r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 2.5d6

    temp = t0 * (1.d0 + 2.d0 * (0.15d0 * (1.d0 + tanh((2.d0-r0)))))

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

  end subroutine perturb_3d

end module initdata_module
