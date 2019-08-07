
module fill_3d_data_module

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use base_state_geometry_module, only: nr_fine, max_radial_level, center, dr
  use amrex_constants_module
  use meth_params_module, only: prob_lo, spherical, s0_interp_type, w0_interp_type, &
       w0mac_interp_type, s0mac_interp_type, &
       use_exact_base_state

  implicit none

  ! private :: addw0, addw0_sphr, make_w0mac_sphr, make_s0mac_sphr, &
  !  make_s0mac_sphr_irreg, make_normal, put_data_on_faces
  !
  ! public

contains

  subroutine put_1d_array_on_cart(lo, hi, lev, &
       s0_cart, s0_cart_lo, s0_cart_hi, nc_s, &
       s0, is_input_edge_centered, is_output_a_vector) &
       bind(C, name="put_1d_array_on_cart")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: s0_cart_lo(3), s0_cart_hi(3)
    integer  , value, intent(in   ) :: nc_s
    double precision, intent(inout) :: s0_cart(s0_cart_lo(1):s0_cart_hi(1), &
         s0_cart_lo(2):s0_cart_hi(2), &
         s0_cart_lo(3):s0_cart_hi(3), 1:nc_s)
    double precision, intent(inout) :: s0(0:max_radial_level,0:nr_fine-1+is_input_edge_centered)
    integer  , value, intent(in   ) :: is_input_edge_centered, is_output_a_vector

    ! local
    integer i,j,k,r
    integer outcomp

    !$gpu

    ! zero s0_cart, then fill in the non-zero values
    s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:nc_s) = 0.d0

    if (is_output_a_vector .eq. 1) then
       outcomp = AMREX_SPACEDIM
    else
       outcomp = 1
    end if

    if (is_input_edge_centered .eq. 1) then

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 2)
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

#if (AMREX_SPACEDIM == 2)
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

                      rfac = (radius - r_edge_loc(0,index+1)) / (r_cc_loc(0,index+1) - r_cc_loc(0,index))

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

                      rfac = (radius - r_edge_loc(0,index+1)) / (r_cc_loc(0,index+1) - r_cc_loc(0,index))

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

    double precision, intent(in   ) :: x,x0,x1,x2,y0,y1,y2
    double precision, intent(  out) :: y

    !$gpu

    y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
         + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)

    if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
    if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)

  end subroutine quad_interp

  subroutine addw0(lev, lo, hi, &
       uedge, u_lo, u_hi, &
       vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wedge, w_lo, w_hi, &
#endif
       w0,mult) bind(C, name="addw0")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(inout) :: uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(inout) :: vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(inout) :: wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: mult

    ! local
    integer i,j,k

#if (AMREX_SPACEDIM == 2)
    k = lo(3)
    do j = lo(2),hi(2)+1
       do i = lo(1)-1,hi(1)+1
          vedge(i,j,k) = vedge(i,j,k) + mult * w0(lev,j)
       end do
    end do
#elif (AMREX_SPACEDIM == 3)
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)+1
       do j = lo(2)-1,hi(2)+1
          do i = lo(1)-1,hi(1)+1
             wedge(i,j,k) = wedge(i,j,k) + mult * w0(lev,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

#endif

  end subroutine addw0

  subroutine addw0_sphr(lo, hi, &
       umac, u_lo, u_hi, &
       vmac, v_lo, v_hi, &
       wmac, w_lo, w_hi, &
       w0macx, x_lo, x_hi, &
       w0macy, y_lo, y_hi, &
       w0macz, z_lo, z_hi, &
       mult) bind(C, name="addw0_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(inout) ::   umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(inout) ::   vmac(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(inout) ::   wmac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) :: w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) :: w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) :: w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    double precision, intent(in   ) :: mult

    ! local variable
    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             umac(i,j,k) = umac(i,j,k) + mult * w0macx(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             vmac(i,j,k) = vmac(i,j,k) + mult * w0macy(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             wmac(i,j,k) = wmac(i,j,k) + mult * w0macz(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine addw0_sphr

  subroutine make_w0mac_sphr(lo, hi, w0, &
       w0macx, x_lo, x_hi, &
       w0macy, y_lo, y_hi, &
       w0macz, z_lo, z_hi, &
       w0_cart, w0_lo, w0_hi, nc_w0, &
       dx, &
       r_edge_loc) bind(C, name="make_w0mac_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(inout) ::  w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(inout) ::  w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(inout) ::  w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    integer         , intent(in   ) :: w0_lo(3), w0_hi(3), nc_w0
    double precision, intent(inout) :: w0_cart(w0_lo(1):w0_hi(1),w0_lo(2):w0_hi(2), &
         w0_lo(3):w0_hi(3),nc_w0)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)

    ! Local variables
    integer          :: i,j,k,index
    double precision :: x,y,z
    double precision :: radius,w0_cart_val,rfac
    double precision, pointer :: w0_nodal(:,:,:,:)

    ! we currently have three different ideas for computing w0mac
    ! 1.  Interpolate w0 to cell centers, then average to edges
    ! 2.  Interpolate w0 to edges directly using linear interpolation
    ! 3.  Interpolate w0 to edges directly using quadratic interpolation
    ! 4.  Interpolate w0 to nodes, then average to edges

    if (w0mac_interp_type .eq. 1) then

       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                w0macx(i,j,k) = HALF* (w0_cart(i-1,j,k,1) + w0_cart(i,j,k,1))
             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                w0macy(i,j,k) = HALF* (w0_cart(i,j-1,k,2) + w0_cart(i,j,k,2))
             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                w0macz(i,j,k) = HALF* (w0_cart(i,j,k-1,3) + w0_cart(i,j,k,3))
             end do
          end do
       end do

    else if (w0mac_interp_type .eq. 2) then

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                rfac = (radius - dble(index)*dr(0)) / dr(0)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(0,index+1) + (ONE-rfac) * w0(0,index)
                else
                   w0_cart_val = w0(0,nr_fine)
                end if

                w0macx(i,j,k) = w0_cart_val * x / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                rfac = (radius - dble(index)*dr(0)) / dr(0)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(0,index+1) + (ONE-rfac) * w0(0,index)
                else
                   w0_cart_val = w0(0,nr_fine)
                end if

                w0macy(i,j,k) = w0_cart_val * y / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                rfac = (radius - dble(index)*dr(0)) / dr(0)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(0,index+1) + (ONE-rfac) * w0(0,index)
                else
                   w0_cart_val = w0(0,nr_fine)
                end if

                w0macz(i,j,k) = w0_cart_val * z / radius

             end do
          end do
       end do

    else if (w0mac_interp_type .eq. 3) then

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
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
                     w0_cart_val, &
                     w0(0,index),w0(0,index+1),w0(0,index+2))

                w0macx(i,j,k) = w0_cart_val * x / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
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
                     w0_cart_val, &
                     w0(0,index),w0(0,index+1),w0(0,index+2))

                w0macy(i,j,k) = w0_cart_val * y / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
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
                     w0_cart_val, &
                     w0(0,index),w0(0,index+1),w0(0,index+2))

                w0macz(i,j,k) = w0_cart_val * z / radius

             end do
          end do
       end do

    else if (w0mac_interp_type .eq. 4) then

       call bl_allocate(w0_nodal,lo(1)-1,hi(1)+2,lo(2)-1,hi(2)+2,lo(3)-1,hi(3)+2,1,3)

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k))*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j))*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i))*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                rfac = (radius - dble(index)*dr(0)) / dr(0)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(0,index+1) + (ONE-rfac) * w0(0,index)
                else
                   w0_cart_val = w0(0,nr_fine)
                end if

                w0_nodal(i,j,k,1) = w0_cart_val * x * (ONE / radius)
                w0_nodal(i,j,k,2) = w0_cart_val * y * (ONE / radius)
                w0_nodal(i,j,k,3) = w0_cart_val * z * (ONE / radius)

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                w0macx(i,j,k) = FOURTH*( w0_nodal(i,j,k  ,1) + w0_nodal(i,j+1,k  ,1) &
                     +w0_nodal(i,j,k+1,1) + w0_nodal(i,j+1,k+1,1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                w0macy(i,j,k) = FOURTH*( w0_nodal(i,j,k  ,2) + w0_nodal(i+1,j,k  ,2) &
                     +w0_nodal(i,j,k+1,2) + w0_nodal(i+1,j,k+1,2))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                w0macz(i,j,k) = FOURTH*( w0_nodal(i,j  ,k,3) + w0_nodal(i+1,j  ,k,3) &
                     +w0_nodal(i,j+1,k,3) + w0_nodal(i+1,j+1,k,3))
             end do
          end do
       end do

       call bl_deallocate(w0_nodal)

    else
       call amrex_error('Error: w0mac_interp_type not defined')
    end if

  end subroutine make_w0mac_sphr

  subroutine make_s0mac_sphr(lo, hi, s0, &
       s0macx, x_lo, x_hi, &
       s0macy, y_lo, y_hi, &
       s0macz, z_lo, z_hi, &
       s0_cart, s0_lo, s0_hi, &
       dx, &
       r_cc_loc) bind(C, name="make_s0mac_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: s0(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(inout) ::  s0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(inout) ::  s0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(inout) ::  s0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    integer         , intent(in   ) :: s0_lo(3), s0_hi(3)
    double precision, intent(inout) :: s0_cart(s0_lo(1):s0_hi(1),s0_lo(2):s0_hi(2), &
         s0_lo(3):s0_hi(3))
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer          :: i,j,k,index
    double precision :: x,y,z
    double precision :: radius

    ! we currently have three different ideas for computing s0mac
    ! 1.  Interpolate s0 to cell centers, then average to edges
    ! 2.  Interpolate s0 to edges directly using linear interpolation
    ! 3.  Interpolate s0 to edges directly using quadratic interpolation
    ! 4.  Interpolate s0 to nodes, then average to edges

    if (s0mac_interp_type .eq. 1) then

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                s0macx(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i-1,j,k))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                s0macy(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i,j-1,k))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                s0macz(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i,j,k-1))
             end do
          end do
       end do

    else if (s0mac_interp_type .eq. 2) then

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                if (radius .ge. r_cc_loc(0,index)) then
                   if (index .ge. nr_fine-1) then
                      s0macx(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macx(i,j,k) = s0(0,index+1)*(radius-r_cc_loc(0,index))/dr(0) &
                           + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dr(0)
                   endif
                else
                   if (index .eq. 0) then
                      s0macx(i,j,k) = s0(0,index)
                   else if (index .gt. nr_fine-1) then
                      s0macx(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macx(i,j,k) = s0(0,index)*(radius-r_cc_loc(0,index-1))/dr(0) &
                           + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dr(0)
                   end if
                end if

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                if (radius .ge. r_cc_loc(0,index)) then
                   if (index .ge. nr_fine-1) then
                      s0macy(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macy(i,j,k) = s0(0,index+1)*(radius-r_cc_loc(0,index))/dr(0) &
                           + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dr(0)
                   endif
                else
                   if (index .eq. 0) then
                      s0macy(i,j,k) = s0(0,index)
                   else if (index .gt. nr_fine-1) then
                      s0macy(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macy(i,j,k) = s0(0,index)*(radius-r_cc_loc(0,index-1))/dr(0) &
                           + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dr(0)
                   end if
                end if

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(0))

                if (radius .ge. r_cc_loc(0,index)) then
                   if (index .ge. nr_fine-1) then
                      s0macz(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macz(i,j,k) = s0(0,index+1)*(radius-r_cc_loc(0,index))/dr(0) &
                           + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dr(0)
                   endif
                else
                   if (index .eq. 0) then
                      s0macz(i,j,k) = s0(0,index)
                   else if (index .gt. nr_fine-1) then
                      s0macz(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macz(i,j,k) = s0(0,index)*(radius-r_cc_loc(0,index-1))/dr(0) &
                           + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dr(0)
                   end if
                end if

             end do
          end do
       end do

    else if (s0mac_interp_type .eq. 3) then

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
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
                     s0macx(i,j,k), &
                     s0(0,index-1),s0(0,index),s0(0,index+1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
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
                     s0macy(i,j,k), &
                     s0(0,index-1),s0(0,index),s0(0,index+1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
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
                     s0macz(i,j,k), &
                     s0(0,index-1),s0(0,index),s0(0,index+1))
             end do
          end do
       end do

    else

       call amrex_error('Error: s0mac_interp_type not defined')

    end if

  end subroutine make_s0mac_sphr

  subroutine make_s0mac_sphr_irreg(lo, hi, s0, &
         s0macx, x_lo, x_hi, &
         s0macy, y_lo, y_hi, &
         s0macz, z_lo, z_hi, &
         s0_cart, s0_lo, s0_hi, &
         dx, r_cc_loc) bind(C, name="make_s0mac_sphr_irreg")

    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: s0(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(inout) ::  s0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(inout) ::  s0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(inout) ::  s0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    integer         , intent(in   ) :: s0_lo(3), s0_hi(3)
    double precision, intent(inout) :: s0_cart(s0_lo(1):s0_hi(1),s0_lo(2):s0_hi(2), &
         s0_lo(3):s0_hi(3))
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)

    ! Local variables
    integer          :: i,j,k,index
    double precision :: x,y,z
    double precision :: radius, dri

    ! we currently have three different ideas for computing s0mac
    ! 1.  Interpolate s0 to cell centers, then average to edges
    ! 2.  Interpolate s0 to edges directly using linear interpolation
    ! 3.  Interpolate s0 to edges directly using quadratic interpolation
    ! 4.  Interpolate s0 to nodes, then average to edges

    if (s0mac_interp_type .eq. 1) then

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                s0macx(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i-1,j,k))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                s0macy(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i,j-1,k))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                s0macz(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i,j,k-1))
             end do
          end do
       end do

    else if (s0mac_interp_type .eq. 2) then

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = nint( (radius/dx(1))**2 - 0.375d0 )  ! closest radial index to edge-centered point

                if (radius .ge. r_cc_loc(0,index)) then
                   dri = r_cc_loc(0,index+1) - r_cc_loc(0,index)
                   if (index .ge. nr_fine-1) then
                      s0macx(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macx(i,j,k) = s0(0,index+1)*(radius-r_cc_loc(0,index))/dri &
                           + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dri
                   endif
                else
                   dri = r_cc_loc(0,index) - r_cc_loc(0,index-1)
                   if (index .eq. 0) then
                      s0macx(i,j,k) = s0(0,index)
                   else if (index .gt. nr_fine-1) then
                      s0macx(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macx(i,j,k) = s0(0,index)*(radius-r_cc_loc(0,index-1))/dri &
                           + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dri
                   end if
                end if

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = nint( (radius/dx(2))**2 - 0.375d0 )  ! closest radial index to edge-centered point

                if (radius .ge. r_cc_loc(0,index)) then
                   dri = r_cc_loc(0,index+1) - r_cc_loc(0,index)
                   if (index .ge. nr_fine-1) then
                      s0macy(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macy(i,j,k) = s0(0,index+1)*(radius-r_cc_loc(0,index))/dri &
                           + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dri
                   endif
                else
                   dri = r_cc_loc(0,index) - r_cc_loc(0,index-1)
                   if (index .eq. 0) then
                      s0macy(i,j,k) = s0(0,index)
                   else if (index .gt. nr_fine-1) then
                      s0macy(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macy(i,j,k) = s0(0,index)*(radius-r_cc_loc(0,index-1))/dri &
                           + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dri
                   end if
                end if

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = nint( (radius/dx(3))**2 - 0.375d0 )  ! closest radial index to edge-centered point

                if (radius .ge. r_cc_loc(0,index)) then
                   dri = r_cc_loc(0,index+1) - r_cc_loc(0,index)
                   if (index .ge. nr_fine-1) then
                      s0macz(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macz(i,j,k) = s0(0,index+1)*(radius-r_cc_loc(0,index))/dri &
                           + s0(0,index)*(r_cc_loc(0,index+1)-radius)/dri
                   endif
                else
                   dri = r_cc_loc(0,index) - r_cc_loc(0,index-1)
                   if (index .eq. 0) then
                      s0macz(i,j,k) = s0(0,index)
                   else if (index .gt. nr_fine-1) then
                      s0macz(i,j,k) = s0(0,nr_fine-1)
                   else
                      s0macz(i,j,k) = s0(0,index)*(radius-r_cc_loc(0,index-1))/dri &
                           + s0(0,index-1)*(r_cc_loc(0,index)-radius)/dri
                   end if
                end if

             end do
          end do
       end do

    else if (s0mac_interp_type .eq. 3) then

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = nint( (radius/dx(1))**2 - 0.375d0 )  ! closest radial index to edge-centered point

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
                     s0macx(i,j,k), &
                     s0(0,index-1),s0(0,index),s0(0,index+1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = nint( (radius/dx(2))**2 - 0.375d0 )  ! closest radial index to edge-centered point

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
                     s0macy(i,j,k), &
                     s0(0,index-1),s0(0,index),s0(0,index+1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = nint( (radius/dx(3))**2 - 0.375d0 )  ! closest radial index to edge-centered point

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
                     s0macz(i,j,k), &
                     s0(0,index-1),s0(0,index),s0(0,index+1))
             end do
          end do
       end do

    else

       call amrex_error('Error: s0mac_interp_type not defined')

    end if

  end subroutine make_s0mac_sphr_irreg

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

    else
       call amrex_error('SHOULDNT CALL MAKE_3D_NORMAL WITH SPHERICAL = 0')
    end if

  end subroutine make_normal

  subroutine put_data_on_faces(lo, hi, &
       scc, cc_lo, cc_hi, &
       facex, x_lo, x_hi, &
       facey, y_lo, y_hi, &
#if (AMREX_SPACEDIM == 3)
       facez, z_lo, z_hi, &
#endif
       harmonic_avg) bind(C, name="put_data_on_faces")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: cc_lo(3), cc_hi(3)
    double precision, intent(in   ) :: scc(cc_lo(1):cc_hi(1),cc_lo(2):cc_hi(2),cc_lo(3):cc_hi(3))
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(inout) :: facex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(inout) :: facey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(inout) :: facez(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
#endif
    integer         , intent(in   ) :: harmonic_avg

    ! local
    integer :: i,j,k
    double precision :: denom

    if (harmonic_avg .eq. 1) then

       k = lo(3)
       j = lo(2)
#if (AMREX_SPACEDIM == 3)
       do k = lo(3),hi(3)
#endif
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                denom = (scc(i,j,k) + scc(i-1,j,k))
                if (denom .ne. 0.d0) then
                   facex(i,j,k) = TWO*(scc(i,j,k) * scc(i-1,j,k)) / denom
                else
                   facex(i,j,k) = HALF*denom
                end if
             end do
          end do
#if (AMREX_SPACEDIM == 3)
       end do
#endif

#if (AMREX_SPACEDIM == 3)
       do k = lo(3),hi(3)
#endif
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                denom = (scc(i,j,k) + scc(i,j-1,k))
                if (denom .ne. 0.d0) then
                   facey(i,j,k) = TWO*(scc(i,j,k) * scc(i,j-1,k)) / denom
                else
                   facey(i,j,k) = HALF*denom
                end if
             end do
          end do
#if (AMREX_SPACEDIM == 3)
       end do
#endif

#if (AMREX_SPACEDIM == 3)
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                denom = (scc(i,j,k) + scc(i,j,k-1))
                if (denom .ne. 0.d0) then
                   facez(i,j,k) = TWO*(scc(i,j,k) * scc(i,j,k-1)) / denom
                else
                   facez(i,j,k) = HALF*denom
                end if
             end do
          end do
       end do
#endif

    else

#if (AMREX_SPACEDIM == 3)
       do k = lo(3),hi(3)
#endif
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                facex(i,j,k) = HALF*(scc(i,j,k)+scc(i-1,j,k))
             end do
          end do
#if (AMREX_SPACEDIM == 3)
       end do
#endif

#if (AMREX_SPACEDIM == 3)
       do k = lo(3),hi(3)
#endif
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                facey(i,j,k) = HALF*(scc(i,j,k)+scc(i,j-1,k))
             end do
          end do
#if (AMREX_SPACEDIM == 3)
       end do
#endif

#if (AMREX_SPACEDIM == 3)
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                facez(i,j,k) = HALF*(scc(i,j,k)+scc(i,j,k-1))
             end do
          end do
       end do
#endif

    end if

  end subroutine put_data_on_faces

end module fill_3d_data_module
