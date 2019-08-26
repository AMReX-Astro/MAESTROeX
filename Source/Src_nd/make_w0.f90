module make_w0_module
  ! compute w0 -- the base state velocity.  This is based on the average
  ! heating in a layer (Sbar) and the mixing (the eta quantities).  The
  ! computation of w0 for plane-parallel atmospheres was first described
  ! in paper II, with modifications due to mixing in paper III.  For
  ! spherical geometry, it was first described in paper III.

  use amrex_constants_module
  use make_grav_module
  use amrex_error_module
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use fundamental_constants_module, only: Gconst
  use base_state_geometry_module, only: max_radial_level, finest_radial_level, nr_fine, &
                                        dr, r_start_coord, r_end_coord, restrict_base, nr, &
                                        fill_ghost_base, base_cutoff_density_coord, numdisjointchunks
  use meth_params_module, only: spherical, maestro_verbose, do_planar_invsq_grav, do_2d_planar_octant, &
                               dpdt_factor, base_cutoff_density, grav_const, &
                               use_exact_base_state, average_base_state

  implicit none

  private

  public :: make_w0

contains

  subroutine make_w0(w0,w0_old,w0_force,Sbar_in, &
                     rho0_old,rho0_new,p0_old,p0_new, &
                     gamma1bar_old,gamma1bar_new,p0_minus_peosbar, &
                     etarho_ec,etarho_cc,delta_chi_w0, &
                     r_cc_loc,r_edge_loc, &
                     dt,dtold,is_predictor) bind(C, name="make_w0")
    ! Binds to C function ``make_w0``

    double precision, intent(inout) ::               w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::           w0_old(0:max_radial_level,0:nr_fine  )
    double precision, intent(inout) ::         w0_force(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::          Sbar_in(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           p0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           p0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0_minus_peosbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::        etarho_ec(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::        etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::     delta_chi_w0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::         r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::       r_edge_loc(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: dt,dtold
    integer         , intent(in   ) :: is_predictor

    integer         :: r,n
    double precision :: max_w0

    call bl_proffortfuncstart("Maestro::make_w0")

    w0_force = ZERO

    if (spherical .eq. 0) then

       if (do_planar_invsq_grav .OR. do_2d_planar_octant .eq. 1) then

          call make_w0_planar_var_g(w0,w0_old,Sbar_in, &
                                    rho0_old,rho0_new,p0_old,p0_new, &
                                    gamma1bar_old,gamma1bar_new, &
                                    p0_minus_peosbar, &
                                    etarho_cc,w0_force, &
                                    dt,dtold,r_cc_loc,r_edge_loc)
       else
          call make_w0_planar(w0,w0_old,Sbar_in, &
                              p0_old,p0_new,gamma1bar_old,gamma1bar_new, &
                              p0_minus_peosbar,etarho_cc,w0_force, &
                              dt,dtold,delta_chi_w0,is_predictor)
       endif


    else

       if (use_exact_base_state) then
          call make_w0_sphr_irreg(w0(0,:),w0_old(0,:),Sbar_in(0,:), &
                                 rho0_old(0,:),rho0_new(0,:), &
                                 p0_old(0,:),p0_new(0,:), &
                                 gamma1bar_old(0,:),gamma1bar_new(0,:), &
                                 p0_minus_peosbar(0,:), &
                                 etarho_ec(0,:),etarho_cc(0,:),w0_force(0,:), &
                                 r_cc_loc,r_edge_loc,dt,dtold)
       else
!          if (average_base_state) then
!             call make_w0_spherical_simple(w0(0,:),Sbar_in(0,:), &
!                                           rho0_old(0,:),rho0_new(0,:), &
!                                           p0_old(0,:),p0_new(0,:), &
!                                           gamma1bar_old(0,:),gamma1bar_new(0,:), &
!                                           etarho_cc(0,:), &
!                                           r_cc_loc,r_edge_loc,dt)
!          else
             call make_w0_spherical(w0(0,:),w0_old(0,:),Sbar_in(0,:), &
                                   rho0_old(0,:),rho0_new(0,:), &
                                   p0_old(0,:),p0_new(0,:), &
                                   gamma1bar_old(0,:),gamma1bar_new(0,:), &
                                   p0_minus_peosbar(0,:), &
                                   etarho_ec(0,:),etarho_cc(0,:),w0_force(0,:), &
                                   r_cc_loc,r_edge_loc,dt,dtold)
!          end if
       endif

    end if

    if (maestro_verbose .ge. 2) then
       do n=0,finest_radial_level
          max_w0 = zero
          do r=r_start_coord(n,1),r_end_coord(n,1)+1
             max_w0 = max(max_w0, abs(w0(n,r)))
          end do
          if (parallel_IOProcessor()) then
             write(6,*) '... max CFL of w0: ',max_w0 * dt / dr(n)
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,*) ''
       end if
    end if

    call bl_proffortfuncstop("Maestro::make_w0")

  end subroutine make_w0


  subroutine make_w0_planar(w0,w0_old,Sbar_in,p0_old,p0_new, &
                            gamma1bar_old,gamma1bar_new,p0_minus_peosbar, &
                            etarho_cc,w0_force,dt,dtold,delta_chi_w0,is_predictor)

    double precision, intent(  out) ::               w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::           w0_old(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::          Sbar_in(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           p0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           p0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0_minus_peosbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::        etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::         w0_force(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::     delta_chi_w0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dt,dtold
    integer         , intent(in   ) :: is_predictor

    ! Local variables
    integer         :: r, n, i, j, refrat
    double precision :: w0_old_cen(0:max_radial_level,0:nr_fine-1)
    double precision :: w0_new_cen(0:max_radial_level,0:nr_fine-1)
    double precision :: psi_planar(0:nr_fine-1)
    double precision :: w0_avg, div_avg, dt_avg, gamma1bar_p0_avg
    double precision :: offset

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Multilevel Outline
    !
    ! Compute w0 at level 1 only
    ! Initialize new w0 at bottom of coarse base array to zero.
    ! do n=1,finest_radial_level
    !   Compute w0 on edges at level n
    !   Obtain the starting value of w0 from the coarser grid
    !   if n>1, compare the difference between w0 at top of level n to the
    !           corresponding point on level n-1
    !   do i=n-1,1,-1
    !     Restrict w0 from level n to level i
    !     Offset the w0 on level i above the top of level n
    !   end do
    ! end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w0 = ZERO

    ! Compute w0 on edges at level n
    do n=0,max_radial_level

       psi_planar = ZERO

       do j=1,numdisjointchunks(n)

          if (n .eq. 0) then
             ! Initialize new w0 at bottom of coarse base array to zero.
             w0(0,0) = ZERO
          else
             ! Obtain the starting value of w0 from the coarser grid
             w0(n,r_start_coord(n,j)) = w0(n-1,r_start_coord(n,j)/2)
          end if

          ! compute psi for level n
          do r = r_start_coord(n,j), r_end_coord(n,j)
             if (r .lt. base_cutoff_density_coord(n)) then
                psi_planar(r) = etarho_cc(n,r) * abs(grav_const)
             end if
          end do

          do r=r_start_coord(n,j)+1,r_end_coord(n,j)+1

             gamma1bar_p0_avg = (gamma1bar_old(n,r-1)+gamma1bar_new(n,r-1)) * &
                  (p0_old(n,r-1)+p0_new(n,r-1))/4.0d0

             if (r .lt. base_cutoff_density_coord(n)) then

                if (is_predictor .eq. 1) then
                   delta_chi_w0(n,r-1) = dpdt_factor * &
                        p0_minus_peosbar(n,r-1) / (gamma1bar_old(n,r-1)*p0_old(n,r-1)*dt)
                else
                   delta_chi_w0(n,r-1) = delta_chi_w0(n,r-1) + dpdt_factor * &
                        p0_minus_peosbar(n,r-1) / (gamma1bar_new(n,r-1)*p0_new(n,r-1)*dt)
                end if

             else
                delta_chi_w0(n,r-1) = 0.d0
             end if

             w0(n,r) = w0(n,r-1) + Sbar_in(n,r-1) * dr(n) &
                  - psi_planar(r-1) / gamma1bar_p0_avg * dr(n) &
                  - delta_chi_w0(n,r-1) * dr(n)

          end do

          if (n .gt. 0) then

             ! Compare the difference between w0 at top of level n to
             ! the corresponding point on level n-1
             offset = w0(n,r_end_coord(n,j)+1) - w0(n-1,(r_end_coord(n,j)+1)/2)

             do i=n-1,0,-1

                refrat = 2**(n-i)

                ! Restrict w0 from level n to level i
                do r=r_start_coord(n,j),r_end_coord(n,j)+1
                   if (mod(r,refrat) .eq. 0) then
                      w0(i,r/refrat) = w0(n,r)
                   end if
                end do

                ! Offset the w0 on level i above the top of level n
                do r=(r_end_coord(n,j)+1)/refrat+1,nr(i)
                   w0(i,r) = w0(i,r) + offset
                end do

             end do

          end if

       end do

    end do

       ! zero w0 where there is no corresponding full state array
       do n=1,max_radial_level
          do j=1,numdisjointchunks(n)
             if (j .eq. numdisjointchunks(n)) then
                do r=r_end_coord(n,j)+2,nr(n)
                   w0(n,r) = ZERO
                end do
             else
                do r=r_end_coord(n,j)+2,r_start_coord(n,j+1)-1
                   w0(n,r) = ZERO
                end do
             end if
          end do
       end do

    call restrict_base(w0,0)
    call fill_ghost_base(w0,0)

    do n=0,max_radial_level
       do j=1,numdisjointchunks(n)

          ! Compute the forcing term in the base state velocity
          ! equation, - 1/rho0 grad pi0
          dt_avg = HALF * (dt + dtold)
          do r=r_start_coord(n,j),r_end_coord(n,j)
             w0_old_cen(n,r) = HALF * (w0_old(n,r) + w0_old(n,r+1))
             w0_new_cen(n,r) = HALF * (w0(n,r) + w0(n,r+1))
             w0_avg = HALF * (dt * w0_old_cen(n,r) + dtold *  w0_new_cen(n,r)) / dt_avg
             div_avg = HALF * (dt * (w0_old(n,r+1)-w0_old(n,r)) + &
                  dtold * (w0(n,r+1)-w0(n,r))) / dt_avg
             w0_force(n,r) = (w0_new_cen(n,r)-w0_old_cen(n,r))/dt_avg + w0_avg*div_avg/dr(n)
          end do

       end do
    end do

    call restrict_base(w0_force,1)
    call fill_ghost_base(w0_force,1)

  end subroutine make_w0_planar


  subroutine make_w0_planar_var_g(w0,w0_old,Sbar_in, &
                                  rho0_old,rho0_new,p0_old,p0_new, &
                                  gamma1bar_old,gamma1bar_new, &
                                  p0_minus_peosbar, &
                                  etarho_cc,w0_force, &
                                  dt,dtold,r_cc_loc,r_edge_loc)

    double precision, intent(  out) ::               w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::           w0_old(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::          Sbar_in(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           p0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::           p0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0_minus_peosbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::        etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(  out) ::         w0_force(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::         r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::       r_edge_loc(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: dt,dtold

    ! local variables
    double precision, allocatable :: w0_fine(:)
    double precision, allocatable :: w0bar_fine(:)
    double precision, allocatable :: deltaw0_fine(:)
    double precision, allocatable :: p0_old_fine(:)
    double precision, allocatable :: p0_new_fine(:)
    double precision, allocatable :: p0_nph_fine(:)
    double precision, allocatable :: rho0_old_fine(:)
    double precision, allocatable :: rho0_new_fine(:)
    double precision, allocatable :: rho0_nph_fine(:)
    double precision, allocatable :: gamma1bar_old_fine(:)
    double precision, allocatable :: gamma1bar_new_fine(:)
    double precision, allocatable :: gamma1bar_nph_fine(:)
    double precision, allocatable :: p0_minus_peosbar_fine(:)
    double precision, allocatable :: etarho_cc_fine(:)
    double precision, allocatable :: Sbar_in_fine(:)
    double precision, allocatable :: grav_edge_fine(:)

    integer :: n, r, j

    double precision, allocatable :: A(:), B(:), C(:), u(:), F(:)

    double precision :: gamma1bar_p0_avg, volume_discrepancy, dpdr
    double precision :: dt_avg, w0_avg, div_avg
    double precision :: w0_old_cen(0:finest_radial_level,0:nr(finest_radial_level)-1)
    double precision :: w0_new_cen(0:finest_radial_level,0:nr(finest_radial_level)-1)


    ! The planar 1/r**2 gravity constraint equation is solved
    ! by calling the tridiagonal solver, just like spherical.
    ! This is accomplished by putting all the requisite data
    ! on the finest basestate grid, solving for w0, and then
    ! restricting w0 back down to the coarse grid.


    ! 1) allocate the finely-gridded temporary basestate arrays
    allocate(              w0_fine(0:nr(finest_radial_level)))
    allocate(           w0bar_fine(0:nr(finest_radial_level)))
    allocate(         deltaw0_fine(0:nr(finest_radial_level)))
    allocate(          p0_old_fine(0:nr(finest_radial_level)-1))
    allocate(          p0_new_fine(0:nr(finest_radial_level)-1))
    allocate(          p0_nph_fine(0:nr(finest_radial_level)-1))
    allocate(        rho0_old_fine(0:nr(finest_radial_level)-1))
    allocate(        rho0_new_fine(0:nr(finest_radial_level)-1))
    allocate(        rho0_nph_fine(0:nr(finest_radial_level)-1))
    allocate(   gamma1bar_old_fine(0:nr(finest_radial_level)-1))
    allocate(   gamma1bar_new_fine(0:nr(finest_radial_level)-1))
    allocate(   gamma1bar_nph_fine(0:nr(finest_radial_level)-1))
    allocate(p0_minus_peosbar_fine(0:nr(finest_radial_level)-1))
    allocate(       etarho_cc_fine(0:nr(finest_radial_level)-1))
    allocate(         Sbar_in_fine(0:nr(finest_radial_level)-1))
    allocate(       grav_edge_fine(0:nr(finest_radial_level)))


    ! 2) copy the data into the temp, uniformly-gridded basestate arrays.
    call prolong_base_to_uniform(p0_old,p0_old_fine)
    call prolong_base_to_uniform(p0_new,p0_new_fine)
    call prolong_base_to_uniform(rho0_old,rho0_old_fine)
    call prolong_base_to_uniform(rho0_new,rho0_new_fine)
    call prolong_base_to_uniform(gamma1bar_old,gamma1bar_old_fine)
    call prolong_base_to_uniform(gamma1bar_new,gamma1bar_new_fine)
    call prolong_base_to_uniform(p0_minus_peosbar,p0_minus_peosbar_fine)
    call prolong_base_to_uniform(etarho_cc,etarho_cc_fine)
    call prolong_base_to_uniform(Sbar_in,Sbar_in_fine)

    ! create time-centered base-state quantities
    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr(finest_radial_level)-1
       p0_nph_fine(r)        = HALF*(p0_old_fine(r)        + p0_new_fine(r))
       rho0_nph_fine(r)      = HALF*(rho0_old_fine(r)      + rho0_new_fine(r))
       gamma1bar_nph_fine(r) = HALF*(gamma1bar_old_fine(r) + gamma1bar_new_fine(r))
    enddo
    !$OMP END PARALLEL DO

    ! 3) solve to w0bar -- here we just take into account the Sbar and
    !    volume discrepancy terms
    w0bar_fine(:) = ZERO

    ! lower boundary condition
    w0bar_fine(0) = ZERO

    do r=1,nr(finest_radial_level)
       gamma1bar_p0_avg = gamma1bar_nph_fine(r-1) * p0_nph_fine(r-1)

       if (r-1 .lt. base_cutoff_density_coord(finest_radial_level)) then
          volume_discrepancy = dpdt_factor * p0_minus_peosbar_fine(r-1)/dt
       else
          volume_discrepancy = ZERO
       end if

       w0bar_fine(r) =  w0bar_fine(r-1) + Sbar_in_fine(r-1) * dr(finest_radial_level) &
            - (volume_discrepancy / gamma1bar_p0_avg ) * dr(finest_radial_level)

    enddo

    ! 4) get the edge-centered gravity on the uniformly-gridded
    ! basestate arrays
    call amrex_error("make_w0.f90: need to write make_grav_edge_uniform")
!    call make_grav_edge_uniform(grav_edge_fine, rho0_nph_fine)


    ! 5) solve for delta w0
    deltaw0_fine(:) = ZERO

    ! this takes the form of a tri-diagonal matrix:
    ! A_j (dw_0)_{j-3/2} +
    ! B_j (dw_0)_{j-1/2} +
    ! C_j (dw_0)_{j+1/2} = F_j

    allocate(A(0:nr(finest_radial_level)))
    allocate(B(0:nr(finest_radial_level)))
    allocate(C(0:nr(finest_radial_level)))
    allocate(u(0:nr(finest_radial_level)))
    allocate(F(0:nr(finest_radial_level)))

    A   = ZERO
    B   = ZERO
    C   = ZERO
    F   = ZERO
    u   = ZERO

    !$OMP PARALLEL DO PRIVATE(r,dpdr)
    do r=1,base_cutoff_density_coord(finest_radial_level)
       A(r) = gamma1bar_nph_fine(r-1) * p0_nph_fine(r-1)
       A(r) = A(r) / dr(finest_radial_level)**2

       dpdr = (p0_nph_fine(r)-p0_nph_fine(r-1))/dr(finest_radial_level)

       B(r) = -(gamma1bar_nph_fine(r-1) * p0_nph_fine(r-1) + &
                gamma1bar_nph_fine(r  ) * p0_nph_fine(r  )) / dr(finest_radial_level)**2
       B(r) = B(r) - TWO * dpdr / (r_edge_loc(finest_radial_level,r))

       C(r) = gamma1bar_nph_fine(r) * p0_nph_fine(r)
       C(r) = C(r) / dr(finest_radial_level)**2

       F(r) = TWO * dpdr * w0bar_fine(r) / r_edge_loc(finest_radial_level,r) - &
              grav_edge_fine(r) * (etarho_cc_fine(r) - etarho_cc_fine(r-1)) / &
              dr(finest_radial_level)
    end do
    !$OMP END PARALLEL DO

    ! Lower boundary
    A(0) = zero
    B(0) = one
    C(0) = zero
    F(0) = zero

    ! Upper boundary
    A(base_cutoff_density_coord(finest_radial_level)+1) = -one
    B(base_cutoff_density_coord(finest_radial_level)+1) = one
    C(base_cutoff_density_coord(finest_radial_level)+1) = zero
    F(base_cutoff_density_coord(finest_radial_level)+1) = zero

    ! Call the tridiagonal solver
    call tridiag(A, B, C, F, u, base_cutoff_density_coord(finest_radial_level)+2)

    do r=1,base_cutoff_density_coord(finest_radial_level)+1
       deltaw0_fine(r) = u(r)
    end do

    do r=base_cutoff_density_coord(finest_radial_level)+2,nr(finest_radial_level)
       deltaw0_fine(r) = deltaw0_fine(base_cutoff_density_coord(finest_radial_level)+1)
    end do

    ! 6) compute w0 = w0bar + deltaw0
    w0_fine = w0bar_fine + deltaw0_fine

    ! 7) fill the multilevel w0 array from the uniformly-gridded w0 we
    ! just solved for.  Here, we make the coarse edge underneath equal
    ! to the fine edge value.
    w0(finest_radial_level,:) = w0_fine(:)
    do n = finest_radial_level, 1, -1
       do r = 0, nr(n), 2
          w0(n-1,r/2) = w0(n,r)
       enddo
    enddo

    ! 8) zero w0 where there is no corresponding full state array
    do n=1,finest_radial_level
       do j=1,numdisjointchunks(n)
          if (j .eq. numdisjointchunks(n)) then
             do r=r_end_coord(n,j)+2,nr(n)
                w0(n,r) = ZERO
             end do
          else
             do r=r_end_coord(n,j)+2,r_start_coord(n,j+1)-1
                w0(n,r) = ZERO
             end do
          end if
       end do
    end do

    call restrict_base(w0,0)
    call fill_ghost_base(w0,0)

    ! compute the forcing terms
    do n=0,finest_radial_level
       do j=1,numdisjointchunks(n)

          ! Compute the forcing term in the base state velocity
          ! equation, - 1/rho0 grad pi0
          dt_avg = HALF * (dt + dtold)
          do r=r_start_coord(n,j),r_end_coord(n,j)
             w0_old_cen(n,r) = HALF * (w0_old(n,r) + w0_old(n,r+1))
             w0_new_cen(n,r) = HALF * (w0(n,r) + w0(n,r+1))
             w0_avg = HALF * (dt * w0_old_cen(n,r) + dtold *  w0_new_cen(n,r)) / dt_avg
             div_avg = HALF * (dt * (w0_old(n,r+1)-w0_old(n,r)) + &
                  dtold * (w0(n,r+1)-w0(n,r))) / dt_avg
             w0_force(n,r) = (w0_new_cen(n,r)-w0_old_cen(n,r))/dt_avg + w0_avg*div_avg/dr(n)
          end do

       end do
    end do

    call restrict_base(w0_force,1)
    call fill_ghost_base(w0_force,1)

  end subroutine make_w0_planar_var_g

  subroutine make_w0_spherical(w0,w0_old,Sbar_in, &
                               rho0_old,rho0_new,p0_old,p0_new, &
                               gamma1bar_old,gamma1bar_new, &
                               p0_minus_peosbar, &
                               etarho_ec,etarho_cc,w0_force, &
                               r_cc_loc,r_edge_loc,dt,dtold)

    double precision, intent(  out) ::               w0(0:nr_fine  )
    double precision, intent(in   ) ::           w0_old(0:nr_fine  )
    double precision, intent(in   ) ::          Sbar_in(0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_old(0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_new(0:nr_fine-1)
    double precision, intent(in   ) ::           p0_old(0:nr_fine-1)
    double precision, intent(in   ) ::           p0_new(0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_old(0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_new(0:nr_fine-1)
    double precision, intent(in   ) :: p0_minus_peosbar(0:nr_fine-1)
    double precision, intent(in   ) ::        etarho_ec(0:nr_fine  )
    double precision, intent(in   ) ::        etarho_cc(0:nr_fine-1)
    double precision, intent(  out) ::         w0_force(0:nr_fine-1)
    double precision, intent(in   ) ::         r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::       r_edge_loc(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: dt,dtold

    ! Local variables
    integer                    :: r
    double precision            :: dpdr, volume_discrepancy, w0_avg, div_avg, dt_avg

    double precision ::    w0_old_cen(0:nr_fine-1)
    double precision ::    w0_new_cen(0:nr_fine-1)
    double precision :: gamma1bar_nph(0:nr_fine-1)
    double precision ::        p0_nph(0:nr_fine-1)
    double precision ::             A(0:nr_fine)
    double precision ::             B(0:nr_fine)
    double precision ::             C(0:nr_fine)
    double precision ::             u(0:nr_fine)
    double precision ::             F(0:nr_fine)
    double precision ::  w0_from_Sbar(0:nr_fine)

    ! These need the extra dimension so we can call make_grav_edge
    double precision ::  rho0_nph(0:0,0:nr_fine-1)
    double precision :: grav_edge(0:0,0:nr_fine  )

    ! create time-centered base-state quantities
    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       p0_nph(r)        = HALF*(p0_old(r)        + p0_new(r))
       rho0_nph(0,r)    = HALF*(rho0_old(r)      + rho0_new(r))
       gamma1bar_nph(r) = HALF*(gamma1bar_old(r) + gamma1bar_new(r))
    enddo
    !$OMP END PARALLEL DO

    ! NOTE: We first solve for the w0 resulting only from Sbar,
    !      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    !      (r^2 Sbar).  Then we will solve for the update, delta w0.

    w0_from_Sbar = ZERO

    do r=1,nr_fine

       if (rho0_old(r-1) .gt. base_cutoff_density) then
          volume_discrepancy = dpdt_factor * p0_minus_peosbar(r-1)/dt
       else
          volume_discrepancy = ZERO
       endif

       w0_from_Sbar(r) = w0_from_Sbar(r-1) + dr(0) * Sbar_in(r-1) * r_cc_loc(0,r-1)**2 - &
            dr(0)* volume_discrepancy * r_cc_loc(0,r-1)**2 / (gamma1bar_nph(r-1)*p0_nph(r-1))

    end do

    !$OMP PARALLEL DO PRIVATE(r)
    do r=1,nr_fine
       w0_from_Sbar(r) = w0_from_Sbar(r) / r_edge_loc(0,r)**2
    end do
    !$OMP END PARALLEL DO

    ! make the edge-centered gravity
    call make_grav_edge(grav_edge,rho0_nph,r_edge_loc)

    ! NOTE:  now we solve for the remainder, (r^2 * delta w0)
    ! this takes the form of a tri-diagonal matrix:
    ! A_j (r^2 dw_0)_{j-3/2} +
    ! B_j (r^2 dw_0)_{j-1/2} +
    ! C_j (r^2 dw_0)_{j+1/2} = F_j

    A   = ZERO
    B   = ZERO
    C   = ZERO
    F   = ZERO
    u   = ZERO

    ! Note that we are solving for (r^2 delta w0), not just w0.

    !$OMP PARALLEL DO PRIVATE(r,dpdr)
    do r=1,base_cutoff_density_coord(0)
       A(r) = gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(0,r-1)**2
       A(r) = A(r) / dr(0)**2

       B(r) = -( gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(0,r-1)**2 &
                +gamma1bar_nph(r  ) * p0_nph(r  ) / r_cc_loc(0,r  )**2 ) / dr(0)**2

       dpdr = (p0_nph(r)-p0_nph(r-1))/dr(0)

       B(r) = B(r) - four * dpdr / (r_edge_loc(0,r))**3

       C(r) = gamma1bar_nph(r) * p0_nph(r) / r_cc_loc(0,r)**2
       C(r) = C(r) / dr(0)**2

       F(r) = four * dpdr * w0_from_Sbar(r) / r_edge_loc(0,r) - &
              grav_edge(0,r) * (r_cc_loc(0,r  )**2 * etarho_cc(r  ) - &
              r_cc_loc(0,r-1)**2 * etarho_cc(r-1)) / &
              (dr(0) * r_edge_loc(0,r)**2) - &
              four * M_PI * Gconst * HALF * &
              (rho0_nph(0,r) + rho0_nph(0,r-1)) * etarho_ec(r)
    end do
    !$OMP END PARALLEL DO

    ! Lower boundary
    A(0) = zero
    B(0) = one
    C(0) = zero
    F(0) = zero

    ! Upper boundary
    A(base_cutoff_density_coord(0)+1) = -one
    B(base_cutoff_density_coord(0)+1) = one
    C(base_cutoff_density_coord(0)+1) = zero
    F(base_cutoff_density_coord(0)+1) = zero

    ! Call the tridiagonal solver
    call tridiag(A, B, C, F, u, base_cutoff_density_coord(0)+2)

    w0(0) = ZERO + w0_from_Sbar(0)

    !$OMP PARALLEL DO PRIVATE(r)
    do r=1,base_cutoff_density_coord(0)+1
       w0(r) = u(r) / r_edge_loc(0,r)**2 + w0_from_Sbar(r)
    end do
    !$OMP END PARALLEL DO

    do r=base_cutoff_density_coord(0)+2,nr_fine
       w0(r) = w0(base_cutoff_density_coord(0)+1)&
            *r_edge_loc(0,base_cutoff_density_coord(0)+1)**2/r_edge_loc(0,r)**2
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    dt_avg = HALF * (dt + dtold)

    !$OMP PARALLEL DO PRIVATE(r,w0_avg,div_avg)
    do r = 0,nr_fine-1
       w0_old_cen(r) = HALF * (w0_old(r) + w0_old(r+1))
       w0_new_cen(r) = HALF * (w0    (r) + w0    (r+1))
       w0_avg = HALF * (dt *  w0_old_cen(r)           + dtold *  w0_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (w0_old(r+1)-w0_old(r)) + dtold * (w0(r+1)-w0(r))) / dt_avg
       w0_force(r) = (w0_new_cen(r)-w0_old_cen(r)) / dt_avg + w0_avg * div_avg / dr(0)
    end do
    !$OMP END PARALLEL DO

  end subroutine make_w0_spherical

  subroutine make_w0_spherical_simple(w0,Sbar_in, &
                                      rho0_old,rho0_new,p0_old,p0_new, &
                                      gamma1bar_old,gamma1bar_new,etarho_cc, &
                                      r_cc_loc,r_edge_loc,dt)

    double precision, intent(  out) ::            w0(0:nr_fine  )
    double precision, intent(in   ) ::       Sbar_in(0:nr_fine-1)
    double precision, intent(in   ) ::      rho0_old(0:nr_fine-1)
    double precision, intent(in   ) ::      rho0_new(0:nr_fine-1)
    double precision, intent(in   ) ::        p0_old(0:nr_fine-1)
    double precision, intent(in   ) ::        p0_new(0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar_old(0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar_new(0:nr_fine-1)
    double precision, intent(in   ) ::     etarho_cc(0:nr_fine-1) ! store dp0dt here
    double precision, intent(in   ) ::      r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    r_edge_loc(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: dt

    ! Local variables
    integer :: r
    double precision :: LHS, RHS

    double precision ::      rho0_nph(0:0,0:nr_fine-1)
    double precision :: gamma1bar_nph(0:nr_fine-1)
    double precision ::        p0_nph(0:nr_fine-1)

    double precision ::     grav_cell(0:0,0:nr_fine-1)

    do r=0,nr_fine-1
       p0_nph(r)        = HALF*(p0_old(r)        + p0_new(r))
       rho0_nph(0,r)    = HALF*(rho0_old(r)      + rho0_new(r))
       gamma1bar_nph(r) = HALF*(gamma1bar_old(r) + gamma1bar_new(r))
    end do

    call make_grav_cell(grav_cell,rho0_nph,r_cc_loc,r_edge_loc)

    w0(0) = 0.d0

    do r=1,base_cutoff_density_coord(0)+1
       LHS = r_edge_loc(0,r)**2 / (dr(0)*r_cc_loc(0,r-1)**2) - &
            (rho0_nph(0,r-1)*grav_cell(0,r-1) / (2.d0*gamma1bar_nph(r-1)*p0_nph(r-1)) )

       RHS = r_edge_loc(0,r-1)**2 / (dr(0)*r_cc_loc(0,r-1)**2) + &
            (rho0_nph(0,r-1)*grav_cell(0,r-1) / (2.d0*gamma1bar_nph(r-1)*p0_nph(r-1)) )
       RHS = RHS * w0(r-1)
       RHS = RHS + Sbar_in(r-1) - 1.d0 / (gamma1bar_nph(r-1)*p0_nph(r-1)) * etarho_cc(r-1)

       w0(r) = RHS / LHS
    end do

    do r=base_cutoff_density_coord(0)+2,nr_fine
       w0(r) = w0(base_cutoff_density_coord(0)+1)&
            *r_edge_loc(0,base_cutoff_density_coord(0)+1)**2/r_edge_loc(0,r)**2
    end do

!!$    ! DEBUG - output to file
!!$    open(unit=1234, file="w0_simple.out")
!!$    do r=0,nr_fine
!!$       write(1234,*) r_edge_loc(0,r),w0(r)
!!$    end do
!!$    close(1234)
  end subroutine make_w0_spherical_simple

  subroutine prolong_base_to_uniform(base_ml, base_fine)

    double precision, intent(in   ) :: base_ml(0:max_radial_level,0:nr_fine)
    double precision, intent(inout) :: base_fine(0:nr_fine)

    ! local
    integer :: n, j, r
    logical, allocatable :: imask_fine(:)
    integer :: r1

    ! the mask array will keep track of whether we've filled in data
    ! in a corresponding radial bin.  .false. indicates that we've
    ! already output there.
    allocate(imask_fine(0:nr_fine-1))
    imask_fine(:) = .true.

    ! r1 is the factor between the current level grid spacing and the
    ! FINEST level
    r1 = 1

    do n = finest_radial_level, 0, -1
       do j = 1,numdisjointchunks(n)
          do r = r_start_coord(n,j), r_end_coord(n,j)

             if (any(imask_fine(r*r1:(r+1)*r1-1) ) ) then
                 base_fine(r*r1:(r+1)*r1-1) = base_ml(n,r)
                imask_fine(r*r1:(r+1)*r1-1) = .false.
             endif

          enddo
       enddo

       ! update r1 for the next coarsest level -- assume a jump by
       ! factor of 2
       r1 = r1*2

    enddo

    ! check to make sure that no mask values are still true
    if (any(imask_fine)) then
       call amrex_error("ERROR: unfilled cells in prolong_base_to_uniform")
    endif


  end subroutine prolong_base_to_uniform

  subroutine tridiag(a,b,c,r,u,n)

      integer           , intent(in   ) :: n
      double precision, intent(in   ) :: a(1:n), b(1:n), c(1:n), r(1:n)
      double precision, intent(inout) :: u(1:n)

      integer j
      double precision, allocatable :: gam(:)
      double precision :: bet

      allocate(gam(n))

      if (b(1) .eq. 0) call amrex_error('tridiag: CANT HAVE B(1) = ZERO')

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) call amrex_error('tridiag: TRIDIAG FAILED')
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

  end subroutine tridiag

  subroutine make_w0_sphr_irreg(w0,w0_old,Sbar_in, &
                                 rho0_old,rho0_new,p0_old,p0_new, &
                                 gamma1bar_old,gamma1bar_new, &
                                 p0_minus_peosbar, &
                                 etarho_ec,etarho_cc,w0_force, &
                                 r_cc_loc,r_edge_loc,dt,dtold)

    double precision, intent(  out) ::               w0(0:nr_fine  )
    double precision, intent(in   ) ::           w0_old(0:nr_fine  )
    double precision, intent(in   ) ::          Sbar_in(0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_old(0:nr_fine-1)
    double precision, intent(in   ) ::         rho0_new(0:nr_fine-1)
    double precision, intent(in   ) ::           p0_old(0:nr_fine-1)
    double precision, intent(in   ) ::           p0_new(0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_old(0:nr_fine-1)
    double precision, intent(in   ) ::    gamma1bar_new(0:nr_fine-1)
    double precision, intent(in   ) :: p0_minus_peosbar(0:nr_fine-1)
    double precision, intent(in   ) ::        etarho_ec(0:nr_fine  )
    double precision, intent(in   ) ::        etarho_cc(0:nr_fine-1)
    double precision, intent(  out) ::         w0_force(0:nr_fine-1)
    double precision, intent(in   ) ::         r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::       r_edge_loc(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: dt,dtold

    ! Local variables
    integer                    :: r
    double precision           :: dpdr, volume_discrepancy, w0_avg, div_avg, dt_avg
    double precision           :: dr1, dr2, dr3

    double precision ::    w0_old_cen(0:nr_fine-1)
    double precision ::    w0_new_cen(0:nr_fine-1)
    double precision :: gamma1bar_nph(0:nr_fine-1)
    double precision ::        p0_nph(0:nr_fine-1)
    double precision ::             A(0:nr_fine)
    double precision ::             B(0:nr_fine)
    double precision ::             C(0:nr_fine)
    double precision ::             u(0:nr_fine)
    double precision ::             F(0:nr_fine)
    double precision ::  w0_from_Sbar(0:nr_fine)

    ! These need the extra dimension so we can call make_grav_edge
    double precision ::  rho0_nph(0:0,0:nr_fine-1)
    double precision :: grav_edge(0:0,0:nr_fine  )

    ! create time-centered base-state quantities
    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_fine-1
       p0_nph(r)        = HALF*(p0_old(r)        + p0_new(r))
       rho0_nph(0,r)    = HALF*(rho0_old(r)      + rho0_new(r))
       gamma1bar_nph(r) = HALF*(gamma1bar_old(r) + gamma1bar_new(r))
    enddo
    !$OMP END PARALLEL DO

    ! NOTE: We first solve for the w0 resulting only from Sbar,
    !      w0_from_sbar by integrating d/dr (r^2 w0_from_sbar) =
    !      (r^2 Sbar).  Then we will solve for the update, delta w0.

    w0_from_Sbar = ZERO

    do r=1,nr_fine

       if (rho0_old(r-1) .gt. base_cutoff_density) then
          volume_discrepancy = dpdt_factor * p0_minus_peosbar(r-1)/dt
       else
          volume_discrepancy = ZERO
       endif

       dr1 = r_edge_loc(0,r) - r_edge_loc(0,r-1)
       w0_from_Sbar(r) = w0_from_Sbar(r-1) + dr1 * Sbar_in(r-1) * r_cc_loc(0,r-1)**2 - &
            dr1* volume_discrepancy * r_cc_loc(0,r-1)**2 / (gamma1bar_nph(r-1)*p0_nph(r-1))

    end do

    !$OMP PARALLEL DO PRIVATE(r)
    do r=1,nr_fine
       w0_from_Sbar(r) = w0_from_Sbar(r) / r_edge_loc(0,r)**2
    end do
    !$OMP END PARALLEL DO

    ! make the edge-centered gravity
    call make_grav_edge(grav_edge,rho0_nph,r_edge_loc)

    ! NOTE:  now we solve for the remainder, (r^2 * delta w0)
    ! this takes the form of a tri-diagonal matrix:
    ! A_j (r^2 dw_0)_{j-3/2} +
    ! B_j (r^2 dw_0)_{j-1/2} +
    ! C_j (r^2 dw_0)_{j+1/2} = F_j

    A   = ZERO
    B   = ZERO
    C   = ZERO
    F   = ZERO
    u   = ZERO

    ! Note that we are solving for (r^2 delta w0), not just w0.

    !$OMP PARALLEL DO PRIVATE(r,dpdr)
    do r=1,base_cutoff_density_coord(0)
       dr1 = r_edge_loc(0,r) - r_edge_loc(0,r-1)
       dr2 = r_edge_loc(0,r+1) - r_edge_loc(0,r)
       dr3 = r_cc_loc(0,r) - r_cc_loc(0,r-1)

       A(r) = gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(0,r-1)**2
       A(r) = A(r) / dr1 / dr3

       B(r) = -( gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(0,r-1)**2 / dr1 &
                +gamma1bar_nph(r  ) * p0_nph(r  ) / r_cc_loc(0,r  )**2 / dr2 ) / dr3

       dpdr = (p0_nph(r)-p0_nph(r-1))/dr3

       B(r) = B(r) - four * dpdr / (r_edge_loc(0,r))**3

       C(r) = gamma1bar_nph(r) * p0_nph(r) / r_cc_loc(0,r)**2
       C(r) = C(r) / dr2 / dr3

       F(r) = four * dpdr * w0_from_Sbar(r) / r_edge_loc(0,r) - &
              grav_edge(0,r) * (r_cc_loc(0,r  )**2 * etarho_cc(r  ) - &
              r_cc_loc(0,r-1)**2 * etarho_cc(r-1)) / &
              (dr3 * r_edge_loc(0,r)**2) - &
              four * M_PI * Gconst * HALF * &
              (rho0_nph(0,r) + rho0_nph(0,r-1)) * etarho_ec(r)
    end do
    !$OMP END PARALLEL DO

    ! Lower boundary
    A(0) = zero
    B(0) = one
    C(0) = zero
    F(0) = zero

    ! Upper boundary
    A(base_cutoff_density_coord(0)+1) = -one
    B(base_cutoff_density_coord(0)+1) = one
    C(base_cutoff_density_coord(0)+1) = zero
    F(base_cutoff_density_coord(0)+1) = zero

    ! Call the tridiagonal solver
    call tridiag(A, B, C, F, u, base_cutoff_density_coord(0)+2)

    w0(0) = ZERO + w0_from_Sbar(0)

    !$OMP PARALLEL DO PRIVATE(r)
    do r=1,base_cutoff_density_coord(0)+1
       w0(r) = u(r) / r_edge_loc(0,r)**2 + w0_from_Sbar(r)
    end do
    !$OMP END PARALLEL DO

    do r=base_cutoff_density_coord(0)+2,nr_fine
       w0(r) = w0(base_cutoff_density_coord(0)+1)&
            *r_edge_loc(0,base_cutoff_density_coord(0)+1)**2/r_edge_loc(0,r)**2
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0
    dt_avg = HALF * (dt + dtold)

    !$OMP PARALLEL DO PRIVATE(r,w0_avg,div_avg)
    do r = 0,nr_fine-1
       dr1 = r_edge_loc(0,r+1) - r_edge_loc(0,r)
       w0_old_cen(r) = HALF * (w0_old(r) + w0_old(r+1))
       w0_new_cen(r) = HALF * (w0    (r) + w0    (r+1))
       w0_avg = HALF * (dt *  w0_old_cen(r)           + dtold *  w0_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (w0_old(r+1)-w0_old(r)) + dtold * (w0(r+1)-w0(r))) / dt_avg
       w0_force(r) = (w0_new_cen(r)-w0_old_cen(r)) / dt_avg + w0_avg * div_avg / dr1
    end do
    !$OMP END PARALLEL DO

  end subroutine make_w0_sphr_irreg

end module make_w0_module
