module plot_variables_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use eos_type_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only: spherical
  !rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
  !      use_eos_e_instead_of_h, use_pprime_in_tfromp
  use base_state_geometry_module, only:  max_radial_level, nr_fine
  ! use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

contains

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

end module plot_variables_module
