
module mac_solver_module

  use base_state_geometry_module, only: nr_fine, max_radial_level

  implicit none

  private

contains

  subroutine mac_solver_rhs(lo, hi, lev, &
       newrhs, nrhs_lo, nrhs_hi, &
       oldrhs, orhs_lo, orhs_hi, &
       uedge, u_lo, u_hi, &
       vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wedge, w_lo, w_hi, &
#endif
       dx) bind(C, name="mac_solver_rhs")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: nrhs_lo(3), nrhs_hi(3)
    double precision, intent(inout) :: newrhs(nrhs_lo(1):nrhs_hi(1),nrhs_lo(2):nrhs_hi(2),nrhs_lo(3):nrhs_hi(3))
    integer         , intent(in   ) :: orhs_lo(3), orhs_hi(3)
    double precision, intent(in   ) :: oldrhs(orhs_lo(1):orhs_hi(1),orhs_lo(2):orhs_hi(2),orhs_lo(3):orhs_hi(3))
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent(in   ) :: dx(3)

    ! local
    integer i,j,k

    !$gpu

    ! Compute newrhs = oldrhs - div(Uedge)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             newrhs(i,j,k) = oldrhs(i,j,k) &
                  - ( (uedge(i+1,j,k)-uedge(i,j,k))/dx(1) &
                  + (vedge(i,j+1,k)-vedge(i,j,k))/dx(2) &
#if (AMREX_SPACEDIM == 3)
                  + (wedge(i,j,k+1)-wedge(i,j,k))/dx(3) &
#endif
                  )

          end do
       end do
    end do

  end subroutine mac_solver_rhs

  subroutine mult_beta0(lev, lo, hi, &
       uedge, u_lo, u_hi, &
       vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
       wedge, w_lo, w_hi, &
#endif
       beta0_edge, &
       beta0, &
       mult_or_div) bind(C, name="mult_beta0")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(inout) :: uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(inout) :: vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(inout) :: wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    double precision, intent(in   ) :: beta0_edge(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    integer, intent(in   ) :: mult_or_div

    ! local
    integer i,j,k

    ! Neglecting ghost cells
    if (mult_or_div .eq. 1) then

       ! Multiply
#if (AMREX_SPACEDIM == 2)
       k = lo(3)
       ! Use cell-centered beta0 to update u-velocity
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             uedge(i,j,k) = uedge(i,j,k) * beta0(lev,j)
          end do
       end do
       ! Use edge beta0 to update v-velocity
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             vedge(i,j,k) = vedge(i,j,k) * beta0_edge(lev,j)
          end do
       end do

#elif (AMREX_SPACEDIM == 3)
       ! Use cell-centered beta0 to update u-velocity
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                uedge(i,j,k) = uedge(i,j,k) * beta0(lev,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       ! Use edge beta0 to update v-velocity
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                vedge(i,j,k) = vedge(i,j,k) * beta0(lev,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       ! Use cell-centered beta0 to update w-velocity
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                wedge(i,j,k) = wedge(i,j,k) * beta0_edge(lev,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
#endif

    else if (mult_or_div .eq. 0) then

       ! Divide
#if (AMREX_SPACEDIM == 2)
       k = lo(3)
       ! Use cell-centered beta0 to update u-velocity
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             uedge(i,j,k) = uedge(i,j,k) / beta0(lev,j)
          end do
       end do
       ! Use edge beta0 to update v-velocity
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             vedge(i,j,k) = vedge(i,j,k) / beta0_edge(lev,j)
          end do
       end do

#elif (AMREX_SPACEDIM == 3)
       ! Use cell-centered beta0 to update u-velocity
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                uedge(i,j,k) = uedge(i,j,k) / beta0(lev,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       ! Use edge beta0 to update v-velocity
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                vedge(i,j,k) = vedge(i,j,k) / beta0(lev,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       ! Use cell-centered beta0 to update w-velocity
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                wedge(i,j,k) = wedge(i,j,k) / beta0_edge(lev,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
#endif

    end if

  end subroutine mult_beta0

  subroutine mac_bcoef_face(lev, lo, hi, &
       xface, x_lo, x_hi, &
       yface, y_lo, y_hi, &
#if (AMREX_SPACEDIM == 3)
       zface, z_lo, z_hi, &
#endif
       rhocc, r_lo, r_hi) bind(C, name="mac_bcoef_face")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(inout) :: xface(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(inout) :: yface(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(inout) :: zface(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
#endif
    integer         , intent(in   ) :: r_lo(3), r_hi(3)
    double precision, intent(in   ) :: rhocc(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))

    ! local
    integer i,j,k

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             xface(i,j,k) = 2.0/(rhocc(i,j,k) + rhocc(i-1,j,k))
          end do
       end do
    end do

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             yface(i,j,k) = 2.0/(rhocc(i,j,k) + rhocc(i,j-1,k))
          end do
       end do
    end do


#if (AMREX_SPACEDIM == 3)
    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             zface(i,j,k) = 2.0/(rhocc(i,j,k) + rhocc(i,j,k-1))
          end do
       end do
    end do
#endif

  end subroutine mac_bcoef_face

end module mac_solver_module
