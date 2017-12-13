module rhoh_vs_t_module

  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
                                use_eos_e_instead_of_h, use_pprime_in_tfromp
  use base_state_geometry_module, only:  max_radial_level, nr_fine

  implicit none

  private

  public :: makeTfromRhoH, makeTfromRhoP

contains

  subroutine makeTfromRhoH(lev,lo,hi,state,s_lo,s_hi,nc_s,p0) bind(C,name="makeTfromRhoH")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3), nc_s
    double precision, intent (inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    if (use_eos_e_instead_of_h) then

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 1)
       r = i
#elif (AMREX_SPACEDIM == 2)
       r = j
#elif (AMREX_SPACEDIM == 3)
       r = k
#endif
             
          ! (rho, (h->e)) --> T, p
            
          eos_state%rho   = state(i,j,k,rho_comp)
          eos_state%T     = state(i,j,k,temp_comp)
          eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

          ! e = h - p/rho
          eos_state%e = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp) - &
               p0(lev,k) / state(i,j,k,rho_comp)
          
          pt_index(:) = (/i, j, k/)
          
          call eos(eos_input_re, eos_state, pt_index)
          
          state(i,j,k,temp_comp) = eos_state%T
          
       enddo
       enddo
       enddo

    else

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 1)
       r = i
#elif (AMREX_SPACEDIM == 2)
       r = j
#elif (AMREX_SPACEDIM == 3)
       r = k
#endif
             
          ! (rho, h) --> T, p
             
          eos_state%rho   = state(i,j,k,rho_comp)
          eos_state%T     = state(i,j,k,temp_comp)
          eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

          eos_state%h = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
          
          pt_index(:) = (/i, j, k/)
                
          call eos(eos_input_rh, eos_state, pt_index)
                
          state(i,j,k,temp_comp) = eos_state%T
             
       enddo
       enddo
       enddo

    endif

  end subroutine makeTfromRhoH

  subroutine makeTfromRhoP(lev,lo,hi,state,s_lo,s_hi,nc_s,p0) bind(C,name="makeTfromRhoP")

    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3), nc_s
    double precision, intent (inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
    do j = lo(2), hi(2)
    do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 1)
       r = i
#elif (AMREX_SPACEDIM == 2)
       r = j
#elif (AMREX_SPACEDIM == 3)
       r = k
#endif
             
       ! (rho, p) --> T
             
       eos_state%rho   = state(i,j,k,rho_comp)
       eos_state%T     = state(i,j,k,temp_comp)
       if (use_pprime_in_tfromp) then
          eos_state%p     = p0(lev,r) + state(i,j,k,pi_comp)
       else
          eos_state%p     = p0(lev,r)
       endif

       eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

       pt_index(:) = (/i, j, k/)
             
       call eos(eos_input_rp, eos_state, pt_index)
             
       state(i,j,k,temp_comp) = eos_state%T
       
    enddo
    enddo
    enddo

  end subroutine makeTfromRhoP

end module rhoh_vs_t_module
