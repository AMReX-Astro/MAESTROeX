module plot_variables_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: spherical, rho_comp, rhoh_comp, temp_comp, spec_comp, base_cutoff_density
  ! , pi_comp, &
  !      use_eos_e_instead_of_h, use_pprime_in_tfromp
  use base_state_geometry_module, only:  max_radial_level, nr_fine
  ! use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

contains

  subroutine make_ad_excess(lo,hi,state,s_lo,s_hi,nc_s,&
       ad_excess,a_lo,a_hi) bind(C, name="make_ad_excess")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3), nc_s
    integer, intent(in) :: a_lo(3), a_hi(3)
    double precision, intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(inout) :: ad_excess(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    integer :: i, j, k
    integer :: pt_index(3)
    double precision :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: chi_rho, chi_t, dt, dp, nabla

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)

             pres(i,j,k) = eos_state%p

             chi_rho = eos_state%rho * eos_state%dpdr / eos_state%p
             chi_t = eos_state%T * eos_state%dpdt / eos_state%p
             nabla_ad(i,j,k) = (eos_state%gam1 - chi_rho) / (chi_t * eos_state%gam1)

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = 0.0d0
             else
                ! forward difference
                if (k == lo(3)) then
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                ! prevent Inf
                if (dp == 0.0d0) then
                   nabla = -huge(0.0d0)
                else
                   nabla = pres(i,j,k) * dt / (dp * state(i,j,k,temp_comp))
                endif
             endif

             ad_excess(i,j,k) = nabla - nabla_ad(i,j,k)

          enddo
       enddo
    enddo


  end subroutine make_ad_excess


  subroutine make_ad_excess_sphr(lo,hi,state,s_lo,s_hi,nc_s,normal,n_lo,n_hi,&
       ad_excess,a_lo,a_hi) bind(C, name="make_ad_excess_sphr")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3), nc_s
    integer, intent(in) :: n_lo(3), n_hi(3)
    integer, intent(in) :: a_lo(3), a_hi(3)
    double precision, intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in) :: normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),3)
    double precision, intent(inout) :: ad_excess(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    integer :: i, j, k, c
    integer :: pt_index(3)
    double precision :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: chi_rho, chi_t, dp(4), dt(4), nabla

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)

             pres(i,j,k) = eos_state%p

             chi_rho = eos_state%rho * eos_state%dpdr / eos_state%p
             chi_t = eos_state%T * eos_state%dpdt / eos_state%p
             nabla_ad(i,j,k) = (eos_state%gam1 - chi_rho) / (chi_t * eos_state%gam1)

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = 0.0d0
             else
                ! compute gradient

                ! forward difference
                if (k == lo(3)) then
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt(3) = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                if (j == lo(2)) then
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j,k)
                   ! backward difference
                else if (j == hi(2)) then
                   dt(2) = state(i,j,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j,k) - pres(i,j-1,k)
                   ! centered difference
                else
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j-1,k)
                endif

                if (i == lo(1)) then
                   dt(1) = state(i+1,j,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i,j,k)
                   ! backward difference
                else if (i == hi(1)) then
                   dt(1) = state(i,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i,j,k) - pres(i-1,j,k)
                   ! centered difference
                else
                   dt(1) = state(i+1,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i-1,j,k)
                endif

                ! dot into normal to get d/dr
                dp(4) = 0.d0
                dt(4) = 0.d0
                do c = 1,3
                   dp(4) = dp(4) + dp(c)*normal(i,j,k,c)
                   dt(4) = dt(4) + dt(c)*normal(i,j,k,c)
                enddo

                ! prevent Inf
                if (dp(4) == 0.0d0) then
                   nabla = -huge(0.0d0)
                else
                   nabla = pres(i,j,k)*dt(4) / (dp(4)*state(i,j,k,temp_comp))
                endif
             endif

             ad_excess(i,j,k) = nabla - nabla_ad(i,j,k)

          enddo
       enddo
    enddo


  end subroutine make_ad_excess_sphr


  subroutine make_magvel(lev,lo,hi,vel,v_lo,v_hi,w0,magvel,m_lo,m_hi) bind(C, name="make_magvel")

    integer, intent(in) :: lev, lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: m_lo(3), m_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent (in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(inout) :: magvel(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 1)
             magvel(i,j,k) = sqrt( (vel(i,j,k,1) + 0.5d0*(w0(lev,i) + w0(lev,i+1)) )**2 )
#elif (AMREX_SPACEDIM == 2)
             magvel(i,j,k) = sqrt(  vel(i,j,k,1)**2 + &
                  ( vel(i,j,k,2) + 0.5d0*(w0(lev,j) + w0(lev,j+1)) )**2 )
#elif (AMREX_SPACEDIM == 3)
             magvel(i,j,k) = sqrt(  vel(i,j,k,1)**2 + &
                  vel(i,j,k,2)**2 + &
                  ( vel(i,j,k,3) + 0.5d0*(w0(lev,k) + w0(lev,k+1)) )**2 )
#endif
          enddo
       enddo
    enddo


  end subroutine make_magvel

  subroutine make_magvel_sphr(lo,hi,vel,v_lo,v_hi,w0cart,w_lo,w_hi, &
       magvel,m_lo,m_hi) bind(C, name="make_magvel_sphr")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: m_lo(3), m_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent (in) :: w0cart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),1)
    double precision, intent(inout) :: magvel(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 1)
             magvel(i,j,k) = sqrt((vel(i,j,k,1) + w0cart(i,j,k,1) )**2)
#elif (AMREX_SPACEDIM == 2)
             magvel(i,j,k) = sqrt(vel(i,j,k,1)**2 + (vel(i,j,k,2) + w0cart(i,j,k,1))**2 )
#elif (AMREX_SPACEDIM == 3)
             magvel(i,j,k) = sqrt(vel(i,j,k,1)**2 + vel(i,j,k,2)**2 + &
                  (vel(i,j,k,3) + w0cart(i,j,k,1))**2 )
#endif
          enddo
       enddo
    enddo


  end subroutine make_magvel_sphr

  subroutine make_velrc(lo,hi,vel,v_lo,v_hi,w0rcart,w_lo,w_hi,normal,n_lo,n_hi,&
       rad_vel,r_lo,r_hi,circ_vel,c_lo,c_hi) bind(C, name="make_velrc")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: n_lo(3), n_hi(3)
    integer, intent(in) :: r_lo(3), r_hi(3)
    integer, intent(in) :: c_lo(3), c_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent (in) :: w0rcart(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),1)
    double precision, intent (in) :: normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),3)
    double precision, intent(inout) :: rad_vel(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    double precision, intent(inout) :: circ_vel(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

    integer :: i,j,k,n

    circ_vel(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.0d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rad_vel(i,j,k) = vel(i,j,k,1) * normal(i,j,k,1) + &
                  vel(i,j,k,2) * normal(i,j,k,2) + &
                  vel(i,j,k,3) * normal(i,j,k,3)

             do n = 1,3
                circ_vel(i,j,k) = circ_vel(i,j,k) + (vel(i,j,k,n) - rad_vel(i,j,k) * normal(i,j,k,n))**2
             enddo

             circ_vel(i,j,k) = sqrt(circ_vel(i,j,k))

             ! add base state vel to get full radial velocity
             rad_vel(i,j,k) = rad_vel(i,j,k) + w0rcart(i,j,k,1)

          enddo
       enddo
    enddo

  end subroutine make_velrc

end module plot_variables_module
