
subroutine initdata(level, max_levs, time, lo, hi, &
                    scal, scal_lo, scal_hi, &
                    vel, vel_lo, vel_hi, &
                    s0_init, p0_init, &
                    dx) bind(C, name="initdata")

  use amrex_fort_module, only : amrex_spacedim
  use base_state_geometry_module, only: nr_fine
  use meth_params_module, only: nscal, prob_lo, rho_comp, rhoh_comp, temp_comp, spec_comp
  use network, only: nspec

  implicit none
  integer, intent(in) :: level, max_levs, lo(3), hi(3)
  integer, intent(in) :: scal_lo(3), scal_hi(3)
  integer, intent(in) :: vel_lo(3), vel_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
       &                                  scal_lo(2):scal_hi(2), &
       &                                  scal_lo(3):scal_hi(3), 1:nscal)
  double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
       &                                 vel_lo(2):vel_hi(2), &
       &                                 vel_lo(3):vel_hi(3), 1:amrex_spacedim)
  double precision, intent(inout) :: s0_init(1:max_levs,0:nr_fine-1,1:nscal)
  double precision, intent(inout) :: p0_init(1:max_levs,0:nr_fine-1)
  double precision, intent(in) :: dx(3)

  integer          :: i,j,k,r
  double precision :: x,y,z

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
           
           scal(i,j,k,rho_comp)  = s0_init(level+1,r,rho_comp)
           scal(i,j,k,rhoh_comp) = s0_init(level+1,r,rhoh_comp)
           scal(i,j,k,temp_comp) = s0_init(level+1,r,temp_comp)
           scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                s0_init(level+1,r,spec_comp:spec_comp+nspec-1)

        end do
     end do
  end do

end subroutine initdata
