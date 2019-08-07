#include "AMReX_BC_TYPES.H"

module plot_variables_module

  use eos_type_module
  use eos_module
  use network, only: nspec, aion
  use meth_params_module, only: spherical, rho_comp, rhoh_comp, temp_comp, spec_comp, &
       pi_comp, use_pprime_in_tfromp, base_cutoff_density
  use base_state_geometry_module, only:  max_radial_level, nr_fine

  implicit none

contains

  subroutine make_ad_excess(lo,hi,state,s_lo,s_hi,nc_s,&
       ad_excess,a_lo,a_hi) bind(C, name="make_ad_excess")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, value, intent(in) :: nc_s
    integer, intent(in) :: a_lo(3), a_hi(3)
    double precision, intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent(inout) :: ad_excess(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    integer :: i, j, k
    integer :: pt_index(3)
    double precision :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: chi_rho, chi_t, dt, dp, nabla

    type (eos_t) :: eos_state

    !$gpu

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
#if (AMREX_SPACEDIM == 2)
                ! forward difference
                if (j == lo(2)) then
                   dt = state(i,j+1,k,temp_comp) - state(i,j,k,temp_comp)
                   dp = pres(i,j+1,k) - pres(i,j,k)
                   ! backward difference
                else if (j == hi(2)) then
                   dt = state(i,j,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp = pres(i,j,k) - pres(i,j-1,k)
                   ! centered difference
                else
                   dt = state(i,j+1,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp = pres(i,j+1,k) - pres(i,j-1,k)
                endif
#else
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
#endif

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
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, value, intent(in) :: nc_s
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

    !$gpu

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


  subroutine make_vorticity(lo,hi,vel,v_lo,v_hi,dx,vort,d_lo,d_hi,bc) bind(C, name="make_vorticity")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: vort(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    integer, intent(in) :: bc(AMREX_SPACEDIM,2)

    integer :: i, j, k
    logical :: fix_lo_x,fix_hi_x,fix_lo_y,fix_hi_y,fix_lo_z,fix_hi_z
    double precision :: wy,vz,uz,wx,vx,uy

#if (AMREX_SPACEDIM == 2)
    call make_vorticity_2d(lo,hi,vel,v_lo,v_hi,dx,vort,d_lo,d_hi,bc)
#else
    call make_vorticity_3d(lo,hi,vel,v_lo,v_hi,dx,vort,d_lo,d_hi,bc)
#endif

  end subroutine make_vorticity

  subroutine make_vorticity_2d(lo,hi,vel,v_lo,v_hi,dx,vort,d_lo,d_hi,bc)

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: vort(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    integer, intent(in) :: bc(2,2)

    integer :: i, j, k
    double precision :: vx, uy

    k = lo(3)

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          vx = (vel(i+1,j,k,2) - vel(i-1,j,k,2)) / (2.d0*dx(1))
          uy = (vel(i,j+1,k,1) - vel(i,j-1,k,1)) / (2.d0*dx(2))
          vort(i,j,k) = vx - uy
       enddo
    enddo

    if (bc(1,1) .eq. Inflow .or. bc(1,1) .eq. SlipWall .or. bc(1,1) .eq. NoSlipWall) then
       i = lo(1)
       do j = lo(2), hi(2)
          vx = (vel(i+1,j,k,2) + 3.d0*vel(i,j,k,2) - 4.d0*vel(i-1,j,k,2)) / dx(1)
          uy = (vel(i,j+1,k,1) - vel(i,j-1,k,1)) / (2.d0*dx(2))
          vort(i,j,k) = vx - uy
       end do
    end if

    if (bc(1,2) .eq. Inflow .or. bc(1,2) .eq. SlipWall .or. bc(1,2) .eq. NoSlipWall) then
       i = hi(1)
       do j = lo(2), hi(2)
          vx = -(vel(i-1,j,k,2) + 3.d0*vel(i,j,k,2) - 4.d0*vel(i+1,j,k,2)) / dx(1)
          uy = (vel(i,j+1,k,1) - vel(i,j-1,k,1)) / (2.d0*dx(2))
          vort(i,j,k) = vx - uy
       end do
    end if

    if (bc(2,1) .eq. Inflow .or. bc(2,1) .eq. SlipWall .or. bc(2,1) .eq. NoSlipWall) then
       j = lo(2)
       do i = lo(1), hi(1)
          vx = (vel(i+1,j,k,2) - vel(i-1,j,k,2)) / (2.d0*dx(1))
          uy = (vel(i,j+1,k,1) + 3.d0*vel(i,j,k,1) - 4.d0*vel(i,j-1,k,1)) / dx(2)
          vort(i,j,k) = vx - uy
       end do
    end if

    if (bc(2,2) .eq. Inflow .or. bc(2,2) .eq. SlipWall .or. bc(2,2) .eq. NoSlipWall) then
       j = hi(2)
       do i = lo(1), hi(1)
          vx =  (vel(i+1,j,k,2) - vel(i-1,j,k,2)) / (2.d0*dx(1))
          uy = -(vel(i,j-1,k,1) + 3.d0*vel(i,j,k,1) - 4.d0*vel(i,j+1,k,1)) / dx(2)
          vort(i,j,k) = vx - uy
       end do
    end if


  end subroutine make_vorticity_2d


  subroutine make_vorticity_3d(lo,hi,vel,v_lo,v_hi,dx,vort,d_lo,d_hi,bc)

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: vort(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    integer, intent(in) :: bc(3,2)

    integer :: i, j, k
    logical :: fix_lo_x,fix_hi_x,fix_lo_y,fix_hi_y,fix_lo_z,fix_hi_z
    double precision :: wy,vz,uz,wx,vx,uy

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uy = uycen(i,j,k)
             uz = uzcen(i,j,k)
             vx = vxcen(i,j,k)
             vz = vzcen(i,j,k)
             wx = wxcen(i,j,k)
             wy = wycen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          enddo
       enddo
    enddo

    fix_lo_x = ( bc(1,1) .eq. Inflow .or. bc(1,1) .eq. NoSlipWall )
    fix_hi_x = ( bc(1,2) .eq. Inflow .or. bc(1,2) .eq. NoSlipWall )

    fix_lo_y = ( bc(2,1) .eq. Inflow .or. bc(2,1) .eq. NoSlipWall )
    fix_hi_y = ( bc(2,2) .eq. Inflow .or. bc(2,2) .eq. NoSlipWall )

    fix_lo_z = ( bc(3,1) .eq. Inflow .or. bc(3,1) .eq. NoSlipWall )
    fix_hi_z = ( bc(3,2) .eq. Inflow .or. bc(3,2) .eq. NoSlipWall )

    !
    !     First do all the faces
    !
    if (fix_lo_x) then
       i = lo(1)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxlo(i,j,k)
             wx = wxlo(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
    end if

    if (fix_hi_x) then
       i = hi(1)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxhi(i,j,k)
             wx = wxhi(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
    end if

    if (fix_lo_y) then
       j = lo(2)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uylo(i,j,k)
             wy = wylo(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
    end if

    if (fix_hi_y) then
       j = hi(2)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uyhi(i,j,k)
             wy = wyhi(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
    end if

    if (fix_lo_z) then
       k = lo(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzlo(i,j,k)
             vz = vzlo(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
    end if

    if (fix_hi_z) then
       k = hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzhi(i,j,k)
             vz = vzhi(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
    end if
    !
    !     Next do all the edges
    !
    if (fix_lo_x .and. fix_lo_y) then
       i = lo(1)
       j = lo(2)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_lo_y) then
       i = hi(1)
       j = lo(2)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_hi_y) then
       i = lo(1)
       j = hi(2)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_hi_y) then
       i = hi(1)
       j = hi(2)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_lo_z) then
       i = lo(1)
       k = lo(3)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_lo_z) then
       i = hi(1)
       k = lo(3)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_x .and. fix_hi_z) then
       i = lo(1)
       k = hi(3)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_x .and. fix_hi_z) then
       i = hi(1)
       k = hi(3)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_y .and. fix_lo_z) then
       j = lo(2)
       k = lo(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_y .and. fix_lo_z) then
       j = hi(2)
       k = lo(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_lo_y .and. fix_hi_z) then
       j = lo(2)
       k = hi(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if

    if (fix_hi_y .and. fix_hi_z) then
       j = hi(2)
       k = hi(3)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
    end if
    !
    !     Finally do all the corners
    !
    if (fix_lo_x .and. fix_lo_y .and. fix_lo_z) then
       i = lo(1)
       j = lo(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_lo_z) then
       i = hi(1)
       j = lo(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_lo_z) then
       i = lo(1)
       j = hi(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_lo_z) then
       i = hi(1)
       j = hi(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_lo_y .and. fix_hi_z) then
       i = lo(1)
       j = lo(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_hi_z) then
       i = hi(1)
       j = lo(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_hi_z) then
       i = lo(1)
       j = hi(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_hi_z) then
       i = hi(1)
       j = hi(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

  contains

    function uycen(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = 0.5d0*(vel(i,j+1,k,1)-vel(i,j-1,k,1))/dx(2)
    end function uycen

    function uylo(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = (vel(i,j+1,k,1)+3.0d0*vel(i,j,k,1)-4.0d0*vel(i,j-1,k,1))/(3.0d0*dx(2))
    end function uylo

    function uyhi(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = -(vel(i,j-1,k,1)+3.0d0*vel(i,j,k,1)-4.0d0*vel(i,j+1,k,1))/(3.0d0*dx(2))
    end function uyhi

    function uzcen(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = 0.5d0*(vel(i,j,k+1,1)-vel(i,j,k-1,1))/dx(3)
    end function uzcen

    function uzlo(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = (vel(i,j,k+1,1)+3.0d0*vel(i,j,k,1)-4.0d0*vel(i,j,k-1,1))/(3.0d0*dx(3))
    end function uzlo

    function uzhi(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r =-(vel(i,j,k-1,1)+3.0d0*vel(i,j,k,1)-4.0d0*vel(i,j,k+1,1))/(3.0d0*dx(3))
    end function uzhi

    function vxcen(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = 0.5d0*(vel(i+1,j,k,2)-vel(i-1,j,k,2))/dx(1)
    end function vxcen

    function vxlo(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = (vel(i+1,j,k,2)+3.0d0*vel(i,j,k,2)-4.0d0*vel(i-1,j,k,2))/(3.0d0*dx(1))
    end function vxlo

    function vxhi(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r =-(vel(i-1,j,k,2)+3.0d0*vel(i,j,k,2)-4.0d0*vel(i+1,j,k,2))/(3.0d0*dx(1))
    end function vxhi

    function vzcen(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = 0.5d0*(vel(i,j,k+1,2)-vel(i,j,k-1,2))/dx(3)
    end function vzcen

    function vzlo(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = (vel(i,j,k+1,2)+3.0d0*vel(i,j,k,2)-4.0d0*vel(i,j,k-1,2))/(3.0d0*dx(3))
    end function vzlo

    function vzhi(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r =-(vel(i,j,k-1,2)+3.0d0*vel(i,j,k,2)-4.0d0*vel(i,j,k+1,2))/(3.0d0*dx(3))
    end function vzhi

    function wxcen(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = 0.5d0*(vel(i+1,j,k,3)-vel(i-1,j,k,3))/dx(1)
    end function wxcen

    function wxlo(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = (vel(i+1,j,k,3)+3.0d0*vel(i,j,k,3)-4.0d0*vel(i-1,j,k,3))/(3.0d0*dx(1))
    end function wxlo

    function wxhi(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r =-(vel(i-1,j,k,3)+3.0d0*vel(i,j,k,3)-4.0d0*vel(i+1,j,k,3))/(3.0d0*dx(1))
    end function wxhi

    function wycen(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = 0.5d0*(vel(i,j+1,k,3)-vel(i,j-1,k,3))/dx(2)
    end function wycen

    function wylo(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r = (vel(i,j+1,k,3)+3.0d0*vel(i,j,k,3)-4.0d0*vel(i,j-1,k,3))/(3.0d0*dx(2))
    end function wylo

    function wyhi(i,j,k) result(r)
      integer :: i,j,k
      double precision :: r
      r =-(vel(i,j-1,k,3)+3.0d0*vel(i,j,k,3)-4.0d0*vel(i,j+1,k,3))/(3.0d0*dx(2))
    end function wyhi

    function vorfun(uy,uz,vx,vz,wx,wy) result(r)
      double precision :: uy,uz,vx,vz,wx,wy
      double precision :: r
      r = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
    end function vorfun

  end subroutine make_vorticity_3d


  subroutine make_magvel(lo,hi,lev,vel,v_lo,v_hi,w0,magvel,m_lo,m_hi) bind(C, name="make_magvel")

    integer, intent(in) :: lo(3), hi(3)
    integer, value, intent(in) :: lev
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: m_lo(3), m_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent (in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(inout) :: magvel(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 2)
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

  subroutine make_magvel_sphr(lo,hi,vel,v_lo,v_hi, &
                              w0macx, x_lo, x_hi, &
                              w0macy, y_lo, y_hi, &
                              w0macz, z_lo, z_hi, &
                              magvel,m_lo,m_hi) bind(C, name="make_magvel_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(inout) :: w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(inout) :: w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(inout) :: w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    integer         , intent(in   ) :: m_lo(3), m_hi(3)
    double precision, intent(inout) :: magvel(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             magvel(i,j,k) = sqrt( (vel(i,j,k,1)+0.5d0*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                                   (vel(i,j,k,2)+0.5d0*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                                   (vel(i,j,k,3)+0.5d0*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)
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

    !$gpu

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

  subroutine make_deltagamma(lo,hi,lev,state,s_lo,s_hi,nc_s,p0,gamma1bar,&
       deltagamma,d_lo,d_hi) bind(C,name="make_deltagamma")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer  , value, intent (in   ) :: nc_s
    integer         , intent (in   ) :: d_lo(3), d_hi(3)
    double precision, intent (in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent (inout) :: deltagamma(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
             r = j
#elif (AMREX_SPACEDIM == 3)
             r = k
#endif

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

             deltagamma(i,j,k) = eos_state%gam1 - gamma1bar(lev,r)
          enddo
       enddo
    enddo

  end subroutine make_deltagamma

  subroutine make_deltagamma_sphr(lo,hi,state,s_lo,s_hi,nc_s,&
       p0_cart,p_lo,p_hi,gamma1bar_cart,g_lo,g_hi,&
       deltagamma,d_lo,d_hi) bind(C,name="make_deltagamma_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer  , value, intent (in   ) :: nc_s
    integer         , intent (in   ) :: p_lo(3), p_hi(3)
    integer         , intent (in   ) :: g_lo(3), g_hi(3)
    integer         , intent (in   ) :: d_lo(3), d_hi(3)
    double precision, intent (in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (in   ) :: p0_cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),1)
    double precision, intent (in   ) :: gamma1bar_cart(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),1)
    double precision, intent (inout) :: deltagamma(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i, j, k
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             if (use_pprime_in_tfromp) then
                eos_state%p     = p0_cart(i,j,k,1) + state(i,j,k,pi_comp)
             else
                eos_state%p     = p0_cart(i,j,k,1)
             endif

             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rp, eos_state, pt_index)

             deltagamma(i,j,k) = eos_state%gam1 - gamma1bar_cart(i,j,k,1)
          enddo
       enddo
    enddo

  end subroutine make_deltagamma_sphr

  subroutine make_entropy(lo,hi,lev,state,s_lo,s_hi,nc_s,&
       entropy,d_lo,d_hi) bind(C,name="make_entropy")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer  , value, intent (in   ) :: lev, nc_s
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: d_lo(3), d_hi(3)
    double precision, intent (in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (inout) :: entropy(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i, j, k, r
    integer :: pt_index(3)
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#if (AMREX_SPACEDIM == 2)
             r = j
#elif (AMREX_SPACEDIM == 3)
             r = k
#endif
             eos_state%rho   = state(i,j,k,rho_comp)
             eos_state%T     = state(i,j,k,temp_comp)
             eos_state%xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

             pt_index(:) = (/i, j, k/)

             call eos(eos_input_rt, eos_state, pt_index)

             entropy(i,j,k) = eos_state%s
          enddo
       enddo
    enddo

  end subroutine make_entropy


  subroutine make_divw0(lo,hi,lev,w0,dx,divw0,d_lo,d_hi) bind(C,name="make_divw0")

    integer, intent (in) :: lo(3), hi(3)
    integer, value, intent (in) :: lev
    integer, intent (in) :: d_lo(3), d_hi(3)
    double precision, intent (in) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent (in) :: dx(3)
    double precision, intent (inout) :: divw0(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
#if (AMREX_SPACEDIM == 2)
             divw0(i,j,k) = (w0(lev,j+1) - w0(lev,j)) / dx(2)
#else
             divw0(i,j,k) = (w0(lev,k+1) - w0(lev,k)) / dx(3)
#endif
          enddo
       enddo
    enddo

  end subroutine make_divw0


  subroutine make_divw0_sphr(lo,hi,w0macx,x_lo,x_hi,w0macy,y_lo,y_hi,w0macz,z_lo,z_hi,&
       dx,divw0,d_lo,d_hi) bind(C,name="make_divw0_sphr")

    integer, intent (in) :: lo(3), hi(3)
    integer, intent (in) :: x_lo(3), x_hi(3)
    integer, intent (in) :: y_lo(3), y_hi(3)
    integer, intent (in) :: z_lo(3), z_hi(3)
    integer, intent (in) :: d_lo(3), d_hi(3)
    double precision, intent (in) :: w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    double precision, intent (in) :: w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    double precision, intent (in) :: w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    double precision, intent (in) :: dx(3)
    double precision, intent (inout) :: divw0(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             divw0(i,j,k) = (w0macx(i+1,j,k) - w0macx(i,j,k)) / dx(1) + &
                  (w0macy(i,j+1,k) - w0macy(i,j,k)) / dx(2) + &
                  (w0macz(i,j,k+1) - w0macz(i,j,k)) / dx(3)
          enddo
       enddo
    enddo

  end subroutine make_divw0_sphr


  subroutine make_pidivu(lo,hi,vel,v_lo,v_hi,dx,pi_cc,p_lo,p_hi,nc,&
       pidivu,d_lo,d_hi) bind(C,name="make_pidivu")

    integer, intent (in) :: lo(3), hi(3)
    integer, intent (in) :: v_lo(3), v_hi(3)
    integer, intent (in) :: p_lo(3), p_hi(3)
    integer, value, intent (in) :: nc
    integer, intent (in) :: d_lo(3), d_hi(3)
    double precision, intent(in) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    double precision, intent (in) :: dx(3)
    double precision, intent(in) :: pi_cc(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),nc)
    double precision, intent (inout) :: pidivu(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    ! Local variables
    integer :: i, j, k

    !$gpu

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             pidivu(i,j,k) = pi_cc(i,j,k,pi_comp)*0.5d0*(vel(i+1,j,k,1) - vel(i-1,j,k,1))/dx(1)
             pidivu(i,j,k) = pidivu(i,j,k) + &
                  pi_cc(i,j,k,pi_comp)*0.5d0*(vel(i,j+1,k,2) - vel(i,j-1,k,2))/dx(2)
#if (AMREX_SPACEDIM == 3)
             pidivu(i,j,k) = pidivu(i,j,k) + &
                  pi_cc(i,j,k,pi_comp)*0.5d0*(vel(i,j,k+1,3) - vel(i,j-1,k,3))/dx(3)
#endif
          enddo
       enddo
    enddo



  end subroutine make_pidivu


  subroutine make_abar(lo,hi,state,s_lo,s_hi,nc_s,abar,a_lo,a_hi) &
       bind(C,name="make_abar")
    !
    ! This routine derives the mass fractions of the species.
    !
    integer, intent (in) :: lo(3), hi(3)
    integer, intent (in) :: s_lo(3), s_hi(3)
    integer, value, intent (in) :: nc_s
    integer, intent (in) :: a_lo(3), a_hi(3)
    double precision, intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    double precision, intent (inout) :: abar(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    ! Local variables
    integer :: i, j, k
    double precision :: xn(nspec)

    !$gpu

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             xn(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/state(i,j,k,rho_comp)
             abar(i,j,k) = 1.0d0/(sum(xn(:)/aion(:)))
          enddo
       enddo
    enddo

  end subroutine make_abar

end module plot_variables_module
