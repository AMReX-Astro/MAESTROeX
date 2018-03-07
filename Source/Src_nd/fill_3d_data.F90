
module fill_3d_data_module

  use base_state_geometry_module, only: nr_fine, max_radial_level, center
  use bl_constants_module
  use meth_params_module, only: prob_lo, spherical
  
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

    if (is_output_a_vector .eq. 1) then
       outcomp = AMREX_SPACEDIM
    else
       outcomp = 1
    end if

    if (is_input_edge_centered .eq. 1) then

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

  subroutine addw0(lev, lo, hi, &
                   uedge, u_lo, u_hi, &
#if (AMREX_SPACEDIM >= 2)
                   vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                   wedge, w_lo, w_hi, &
#endif
#endif
                   w0,mult) bind(C, name="addw0")
    
    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(inout) :: uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
#if (AMREX_SPACEDIM >= 2)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(inout) :: vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(inout) :: wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
#endif
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: mult

    ! local
    integer i,j,k

#if (AMREX_SPACEDIM == 1)
    j = lo(2)
    k = lo(3)
    do i = lo(1),hi(1)+1
       uedge(i,j,k) = uedge(i,j,k) + mult * w0(lev,i)
    end do

#elif (AMREX_SPACEDIM == 2)
    k = lo(3)
    do j = lo(2),hi(2)+1
    do i = lo(1)-1,hi(1)+1
       vedge(i,j,k) = vedge(i,j,k) + mult * w0(lev,j)
    end do
    end do
#elif (AMREX_SPACEDIM == 3)
    do k = lo(3),hi(3)+1
    do j = lo(2)-1,hi(2)+1
    do i = lo(1)-1,hi(1)+1
       wedge(i,j,k) = wedge(i,j,k) + mult * w0(lev,k)
    end do
    end do
    end do

#endif

  end subroutine addw0

  subroutine make_normal(normal,n_lo,n_hi,dx) bind(C, name="make_normal")
    
    integer         , intent(in   ) :: n_lo(3), n_hi(3)
    double precision, intent(  out) :: normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),3)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k
    double precision :: x,y,z,radius

    ! normal is the unit vector in the radial direction (e_r) in spherical
    ! coordinates.
    !
    ! in terms of Cartesian coordinates, with unit vectors e_x, e_y, e_z,
    !    e_r = sin(theta)cos(phi) e_x + sin(theta)sin(phi) e_y + cos(theta) e_z
    ! or
    !    e_r = (x/R) e_x + (y/R) e_y + (z/R) e_z

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius)
       do k = n_lo(3),n_hi(3)
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = n_lo(2),n_hi(2)
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = n_lo(1),n_hi(1)
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)

                normal(i,j,k,1) = x * (ONE / radius)
                normal(i,j,k,2) = y * (ONE / radius)
                normal(i,j,k,3) = z * (ONE / radius)

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else 
       call bl_error('SHOULDNT CALL MAKE_3D_NORMAL WITH SPHERICAL = 0')
    end if

  end subroutine make_normal

end module fill_3d_data_module
