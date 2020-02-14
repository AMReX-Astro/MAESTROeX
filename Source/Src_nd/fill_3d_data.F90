
module fill_3d_data_module

use amrex_error_module
use amrex_mempool_module, only : bl_allocate, bl_deallocate
use base_state_geometry_module, only: nr_fine, max_radial_level, center, dr
use amrex_constants_module
use meth_params_module, only: prob_lo, spherical, s0_interp_type, w0_interp_type, &
        w0mac_interp_type, s0mac_interp_type, &
        use_exact_base_state

implicit none

contains

subroutine put_1d_array_on_cart_sphr(lo, hi, &
        s0_cart, s0_cart_lo, s0_cart_hi, nc_s, &
        s0, dx, &
        is_input_edge_centered, &
        is_output_a_vector, &
        r_cc_loc, r_edge_loc, &
        cc_to_r, ccr_lo, ccr_hi) &
        bind(C, name="put_1d_array_on_cart_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s0_cart_lo(3), s0_cart_hi(3)
    integer  , value, intent(in   ) :: nc_s
    double precision, intent(inout) :: s0_cart(s0_cart_lo(1):s0_cart_hi(1), &
        s0_cart_lo(2):s0_cart_hi(2), &
        s0_cart_lo(3):s0_cart_hi(3), nc_s)
    double precision, intent(in   ) :: s0(0:max_radial_level,0:nr_fine-1+is_input_edge_centered)
    double precision, intent(in   ) :: dx(3)
    integer  , value, intent(in   ) :: is_input_edge_centered, is_output_a_vector
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
        ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    ! Local variables
    integer          :: i,j,k,index
    double precision :: x,y,z
    double precision :: radius,rfac,s0_cart_val

    !$gpu

    if (use_exact_base_state) then

        if (is_input_edge_centered .eq. 1) then

        ! we implemented three different ideas for computing s0_cart,
        ! where s0 is edge-centered.
        ! 1.  Piecewise constant
        ! 2.  Piecewise linear
        ! 3.  Quadratic
        if (w0_interp_type .eq. 1) then
            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index = cc_to_r(i,j,k)

                    if (index .lt. nr_fine) then 
                        rfac = (radius - r_edge_loc(0,index+1)) / (r_cc_loc(0,index+1) - r_cc_loc(0,index))
                    else 
                        rfac = (radius - r_edge_loc(0,index+1)) / (r_cc_loc(0,index) - r_cc_loc(0,index-1))
                    endif

                    if (rfac .gt. 0.5d0) then
                        s0_cart_val = s0(0,index+1)
                    else
                        s0_cart_val = s0(0,index)
                    end if

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do
        else if (w0_interp_type .eq. 2) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) +(dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) +(dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = cc_to_r(i,j,k)

                    if (index .lt. nr_fine-1) then 
                        rfac = (radius - r_edge_loc(0,index+1)) / (r_cc_loc(0,index+1) - r_cc_loc(0,index))
                    else 
                        rfac = (radius - r_edge_loc(0,index+1)) / (r_cc_loc(0,index) - r_cc_loc(0,index-1))
                    endif

                    if (index .lt. nr_fine) then
                        s0_cart_val = rfac * s0(0,index+1) + (ONE-rfac) * s0(0,index)
                    else
                        s0_cart_val = s0(0,nr_fine)
                    end if

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do

        else if (w0_interp_type .eq. 3) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = cc_to_r(i,j,k) + 1

                    ! index refers to the lo point in the quadratic stencil
                    if (index .le. 0) then
                        index = 0
                    else if (index .ge. nr_fine-1) then
                        index = nr_fine-2
                    else if (radius-r_edge_loc(0,index) .lt. r_edge_loc(0,index+1)) then
                        index = index-1
                    end if

                    call quad_interp(radius, &
                            r_edge_loc(0,index),r_edge_loc(0,index+1), &
                            r_edge_loc(0,index+2), &
                            s0_cart_val, &
                            s0(0,index),s0(0,index+1),s0(0,index+2))

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do
#ifndef AMREX_USE_CUDA
        else
            call amrex_error('Error: w0_interp_type not defined')
#endif
        end if

        else

        ! we directly inject the spherical values into each cell center
        ! because s0 is also bin-centered.

        do k = lo(3),hi(3)
            z = prob_lo(3) +(dble(k)+HALF)*dx(3) - center(3)
            do j = lo(2),hi(2)
                y = prob_lo(2) +(dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = cc_to_r(i,j,k)

                    s0_cart_val = s0(0,index)

                    if (is_output_a_vector .eq. 1) then
                    s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                    s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                    s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                    s0_cart(i,j,k,1) = s0_cart_val
                    end if

                end do
            end do
        end do

        end if  ! is_input_edge_centered

    else

        if (is_input_edge_centered .eq. 1) then

        ! we currently have three different ideas for computing s0_cart,
        ! where s0 is edge-centered.
        ! 1.  Piecewise constant
        ! 2.  Piecewise linear
        ! 3.  Quadratic

        if (w0_interp_type .eq. 1) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = int(radius / dr(0))

                    rfac = (radius - dble(index)*dr(0)) / dr(0)

                    if (rfac .gt. 0.5d0) then
                        s0_cart_val = s0(0,index+1)
                    else
                        s0_cart_val = s0(0,index)
                    end if

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do

        else if (w0_interp_type .eq. 2) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) +(dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) +(dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = int(radius / dr(0))

                    rfac = (radius - dble(index)*dr(0)) / dr(0)

                    if (index .lt. nr_fine) then
                        s0_cart_val = rfac * s0(0,index+1) + (ONE-rfac) * s0(0,index)
                    else
                        s0_cart_val = s0(0,nr_fine)
                    end if

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do

        else if (w0_interp_type .eq. 3) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = int(radius / dr(0))

                    ! index refers to the lo point in the quadratic stencil
                    if (index .le. 0) then
                        index = 0
                    else if (index .ge. nr_fine-1) then
                        index = nr_fine-2
                    else if (radius-r_edge_loc(0,index) .lt. r_edge_loc(0,index+1)) then
                        index = index-1
                    end if

                    call quad_interp(radius, &
                            r_edge_loc(0,index),r_edge_loc(0,index+1), &
                            r_edge_loc(0,index+2), &
                            s0_cart_val, &
                            s0(0,index),s0(0,index+1),s0(0,index+2))

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do
#ifndef AMREX_USE_CUDA
        else
            call amrex_error('Error: w0_interp_type not defined')
#endif
        end if

        else

        ! we currently have three different ideas for computing s0_cart,
        ! where s0 is bin-centered.
        ! 1.  Piecewise constant
        ! 2.  Piecewise linear
        ! 3.  Quadratic

        if (s0_interp_type .eq. 1) then

            do k = lo(3),hi(3)
                z = prob_lo(3) +(dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) +(dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = int(radius / dr(0))

                    s0_cart_val = s0(0,index)

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do

        else if (s0_interp_type .eq. 2) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = int(radius / dr(0))

                    if (radius .ge. r_cc_loc(0,index)) then
                        if (index .ge. nr_fine-1) then
                            s0_cart_val = s0(0,nr_fine-1)
                        else
                            s0_cart_val = s0(0,index+1)*(radius-r_cc_loc(0,index))/dr(0) &
                                + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dr(0)
                        endif
                    else
                        if (index .eq. 0) then
                            s0_cart_val = s0(0,index)
                        else if (index .gt. nr_fine-1) then
                            s0_cart_val = s0(0,nr_fine-1)
                        else
                            s0_cart_val = s0(0,index)*(radius-r_cc_loc(0,index-1))/dr(0) &
                                + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dr(0)
                        end if
                    end if

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do

        else if (s0_interp_type .eq. 3) then

            do k = lo(3),hi(3)
                z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
                do j = lo(2),hi(2)
                    y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                    do i = lo(1),hi(1)
                    x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                    radius = sqrt(x**2 + y**2 + z**2)
                    index  = int(radius / dr(0))

                    ! index refers to the center point in the quadratic stencil.
                    ! we need to modify this if we're too close to the edge
                    if (index .eq. 0) then
                        index = 1
                    else if (index .ge. nr_fine-1) then
                        index = nr_fine-2
                    end if

                    call quad_interp(radius, &
                            r_cc_loc(0,index-1),r_cc_loc(0,index), &
                            r_cc_loc(0,index+1), &
                            s0_cart_val, &
                            s0(0,index-1),s0(0,index),s0(0,index+1))

                    if (is_output_a_vector .eq. 1) then
                        s0_cart(i,j,k,1) = s0_cart_val * x * (ONE / radius)
                        s0_cart(i,j,k,2) = s0_cart_val * y * (ONE / radius)
                        s0_cart(i,j,k,3) = s0_cart_val * z * (ONE / radius)
                    else
                        s0_cart(i,j,k,1) = s0_cart_val
                    end if

                    end do
                end do
            end do
        else
#ifndef AMREX_USE_CUDA
            call amrex_error('Error: s0_interp_type not defined')
#endif
        end if

        end if  ! is_input_edge_centered

    end if  ! use_exact_base_state

end subroutine put_1d_array_on_cart_sphr

    subroutine quad_interp(x,x0,x1,x2,y,y0,y1,y2)

      double precision, value, intent(in   ) :: x,x0,x1,x2,y0,y1,y2
      double precision, intent(  out) :: y

      !$gpu

      y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)

      if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
      if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)

    end subroutine quad_interp

end module fill_3d_data_module
