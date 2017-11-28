
module fill_3d_data_module

  use base_state_geometry_module, only: nr_fine, max_radial_level
  
  implicit none

  private

contains

  subroutine put_1d_array_on_cart(lev, lo, hi, &
                                  s0_cart, s0_cart_lo, s0_cart_hi, nc_s, &
                                  s0, is_input_edge_centered, is_output_a_vector) &
                                  bind(C, name="put_1d_array_on_cart")
    
    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: s0_cart_lo(3), s0_cart_hi(3), nc_s
    double precision, intent(inout) :: s0_cart(s0_cart_lo(1):s0_cart_hi(1), &
                                               s0_cart_lo(2):s0_cart_hi(2), &
                                               s0_cart_lo(3):s0_cart_hi(3), 1:nc_s)
    double precision, intent(inout) :: s0(0:max_radial_level,0:nr_fine-1+is_input_edge_centered)
    integer         , intent(in   ) :: is_input_edge_centered, is_output_a_vector

    ! local
    integer i,j,k,r
    integer outcomp

    ! zero s0_cart, then fill in the non-zero values
    s0_cart = 0.d0

    if (is_input_edge_centered .eq. 1) then

       if (is_output_a_vector .eq. 1) then
          outcomp = AMREX_SPACEDIM
       else
          outcomp = 1
       end if

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 1)
          r = i
#elif (AMREX_SPACEDIM == 2)
          r = j
#elif (AMREX_SPACEDIM == 3)
          r = k
#endif
          s0_cart(i,j,k,outcomp) = 0.5d0*( s0(lev,r) + s0(lev,r+1) )

       end do
       end do
       end do

    else

       if (is_output_a_vector .eq. 1) then
          outcomp = AMREX_SPACEDIM
       else
          outcomp = 1
       end if

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 1)
          r = i
#elif (AMREX_SPACEDIM == 2)
          r = j
#elif (AMREX_SPACEDIM == 3)
          r = k
#endif
          s0_cart(i,j,k,outcomp) = s0(lev,r)

       end do
       end do
       end do

    end if

  end subroutine put_1d_array_on_cart

end module fill_3d_data_module
