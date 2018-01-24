
module mac_solver_module

  use base_state_geometry_module, only: nr_fine, max_radial_level
  
  implicit none

  private

contains

  subroutine mac_solver_rhs(lev, lo, hi, &
                             newrhs, nrhs_lo, nrhs_hi, &
                             oldrhs, orhs_lo, orhs_hi, &
                             uedge, u_lo, u_hi, &
#if (AMREX_SPACEDIM >= 2)
                             vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                             wedge, w_lo, w_hi, &
#endif
#endif
                             dx) bind(C, name="mac_solver_rhs")
    
    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: nrhs_lo(3), nrhs_hi(3) 
    double precision, intent(inout) :: newrhs(nrhs_lo(1):nrhs_hi(1),nrhs_lo(2):nrhs_hi(2),nrhs_lo(3):nrhs_hi(3))
    integer         , intent(in   ) :: orhs_lo(3), orhs_hi(3)
    double precision, intent(in   ) :: oldrhs(orhs_lo(1):orhs_hi(1),orhs_lo(2):orhs_hi(2),orhs_lo(3):orhs_hi(3))
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: uedge(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
#if (AMREX_SPACEDIM >= 2)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vedge(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wedge(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
#endif
    double precision, intent(in   ) :: dx(3)

    ! local
    integer i,j,k

    ! Compute newrhs = oldrhs - div(Uedge)
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
       newrhs(i,j,k) = oldrhs(i,j,k) & 
                                     - ((uedge(i+1,j,k)-uedge(i,j,k))/dx(1) & 
#if (AMREX_SPACEDIM >= 2)
                                        + (vedge(i,j+1,k)-vedge(i,j,k))/dx(2) & 
#if (AMREX_SPACEDIM == 3)
                                        + (wedge(i,j,k+1)-wedge(i,j,k))/dx(3) &
#endif
#endif
                                        )
    end do
    end do
    end do

  end subroutine mac_solver_rhs

  subroutine mult_beta0(lev, lo, hi, &
                   uedge, u_lo, u_hi, &
#if (AMREX_SPACEDIM >= 2)
                   vedge, v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                   wedge, w_lo, w_hi, &
#endif
                   beta0_edge, &
#endif
                   beta0, & 
                   mult_or_div) bind(C, name="mult_beta0")
    
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
    double precision, intent(in   ) :: beta0_edge(0:max_radial_level,0:nr_fine)
#endif
    double precision, intent(in   ) :: beta0(0:max_radial_level,0:nr_fine-1)
    integer, intent(in   ) :: mult_or_div

    ! local
    integer i,j,k

    ! Neglecting ghost cells 
    if (mult_or_div .eq. 1) then

    ! Multiply
#if (AMREX_SPACEDIM == 1)
       j = lo(2)
       k = lo(3)
       do i = lo(1),hi(1)+1
          uedge(i,j,k) = uedge(i,j,k) * beta0(lev,j)
       end do
    
#elif (AMREX_SPACEDIM == 2)
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
       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          uedge(i,j,k) = uedge(i,j,k) * beta0(lev,j)
       end do
       end do
       end do
       ! Use edge beta0 to update v-velocity
       do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          vedge(i,j,k) = vedge(i,j,k) * beta0_edge(lev,j)
       end do
       end do
       end do
       ! Use cell-centered beta0 to update w-velocity
       do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          wedge(i,j,k) = wedge(i,j,k) * beta0(lev,j)
       end do
       end do
       end do
#endif

    elseif (mult_or_div .eq. 0) then 

   ! Divide
#if (AMREX_SPACEDIM == 1)
       j = lo(2)
       k = lo(3)
       do i = lo(1),hi(1)+1
          uedge(i,j,k) = uedge(i,j,k) / beta0(lev,j)
       end do
    
#elif (AMREX_SPACEDIM == 2)
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
       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          uedge(i,j,k) = uedge(i,j,k) / beta0(lev,j)
       end do
       end do
       end do
       ! Use edge beta0 to update v-velocity
       do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          vedge(i,j,k) = vedge(i,j,k) / beta0_edge(lev,j)
       end do
       end do
       end do
       ! Use cell-centered beta0 to update w-velocity
       do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          wedge(i,j,k) = wedge(i,j,k) / beta0(lev,j)
       end do
       end do
       end do
#endif

    end if

  end subroutine mult_beta0

end module mac_solver_module
