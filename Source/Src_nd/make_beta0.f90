module make_beta0_module

    use amrex_constants_module
    use base_state_geometry_module, only: nr_fine, dr, anelastic_cutoff_density_coord, &
                                          r_start_coord, r_end_coord, &
                                          nr, numdisjointchunks, finest_radial_level, &
                                          max_radial_level, restrict_base, fill_ghost_base
    use meth_params_module, only: beta0_type, use_linear_grav_in_beta0

    implicit none

 contains

 subroutine make_beta0(beta0,rho0,p0,gamma1bar,grav_cell) bind(C, name="make_beta0")

    double precision, intent(  out) :: beta0    (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0     (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: grav_cell(0:max_radial_level,0:nr_fine-1)

    ! local
    integer :: r, n, i, refrat, j
    double precision :: integral
    double precision :: lambda, mu, nu, kappa
    double precision :: denom, coeff1, coeff2, coeff3
    double precision :: del,dpls,dmin,slim,sflag
    double precision :: offset

    double precision, allocatable :: beta0_edge(:,:)

    call bl_proffortfuncstart("Maestro::make_beta0")

    allocate(beta0_edge(0:finest_radial_level,0:nr_fine))

    beta0 = 0.d0

    if (beta0_type .eq. 1) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Compute beta0 on the edges and average to the center
       !
       ! Multilevel Outline:
       !
       ! First, compute beta0 on edges and centers at level 0 only
       ! Obtain the starting value from rho0 at the bottom of the domain.
       ! do n=1,finest_radial_level
       !   Compute beta0 on edges and centers at level n
       !   Obtain the starting value of beta0_edge_lo from the coarser grid
       !   if n>0, compare the difference between beta0 at the top of level n to the
       !           corresponding point on level n-1
       !   do i=n-1,0,-1
       !     Offset the centered beta on level i above this point so the total integral
       !      is consistent
       !     Redo the anelastic cutoff part
       !   end do
       ! end do
       ! call restrict_base and fill_ghost_base
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do n=0,finest_radial_level

          do j=1,numdisjointchunks(n)

             ! Compute beta0 on edges and centers at level n

             if (n .eq. 0) then
                beta0_edge(0,0) = rho0(0,0)
             else
                ! Obtain the starting value of beta0_edge_lo from the coarser grid
                beta0_edge(n,r_start_coord(n,j)) = beta0_edge(n-1,r_start_coord(n,j)/2)
             end if

             do r=r_start_coord(n,j),r_end_coord(n,j)

                if (r < anelastic_cutoff_density_coord(n)) then

                   if (r .eq. 0 .or. r .eq. nr(n)-1) then

                      lambda = ZERO
                      mu = ZERO
                      nu = ZERO

                   else

                      ! piecewise linear reconstruction of rho0,
                      ! gamma1bar, and p0 -- see paper III, appendix C
                      del    = HALF* (rho0(n,r+1) - rho0(n,r-1))/dr(n)
                      dpls   = TWO * (rho0(n,r+1) - rho0(n,r  ))/dr(n)
                      dmin   = TWO * (rho0(n,r  ) - rho0(n,r-1))/dr(n)
                      slim   = min(abs(dpls), abs(dmin))
                      slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag  = sign(ONE,del)
                      lambda = sflag*min(slim,abs(del))

                      del   = HALF* (gamma1bar(n,r+1) - gamma1bar(n,r-1))/dr(n)
                      dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/dr(n)
                      dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      mu    = sflag*min(slim,abs(del))

                      del   = HALF* (p0(n,r+1) - p0(n,r-1))/dr(n)
                      dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/dr(n)
                      dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/dr(n)
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      nu    = sflag*min(slim,abs(del))

                   end if

                   if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
                        (nu*gamma1bar(n,r) - mu*p0(n,r)) .eq. ZERO .or. &
                        ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
                        (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
                        ((p0(n,r) + HALF*nu*dr(n))/ &
                        (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

                      ! just do piecewise constant integration
                      integral = abs(grav_cell(n,r))*rho0(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))

                   else

                      if ( use_linear_grav_in_beta0 ) then

                         ! also do piecewise linear reconstruction of
                         ! gravity -- not documented in publication yet.
                         del   = HALF* (grav_cell(n,r+1) - grav_cell(n,r-1))/dr(n)
                         dpls  = TWO * (grav_cell(n,r+1) - grav_cell(n,r  ))/dr(n)
                         dmin  = TWO * (grav_cell(n,r  ) - grav_cell(n,r-1))/dr(n)
                         slim  = min(abs(dpls), abs(dmin))
                         slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                         sflag = sign(ONE,del)
                         kappa = sflag*min(slim,abs(del))

                         denom = nu*gamma1bar(n,r) - mu*p0(n,r)
                         coeff1 = (lambda*gamma1bar(n,r) - mu*rho0(n,r)) * &
                              (kappa *gamma1bar(n,r) + mu*abs(grav_cell(n,r))) / &
                              (mu*mu*denom)
                         coeff2 = (lambda*p0(n,r) - nu*rho0(n,r))* &
                              (-kappa*p0(n,r) - nu*abs(grav_cell(n,r))) / &
                              (nu*nu*denom)
                         coeff3 = kappa*lambda / (mu*nu)

                         integral =  &
                              coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
                                          (gamma1bar(n,r) - HALF*mu*dr(n)) ) + &
                              coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
                                          (p0(n,r) - HALF*nu*dr(n)) ) - &
                              coeff3*dr(n)

                      else

                         ! paper III, equation C2
                         denom = nu*gamma1bar(n,r) - mu*p0(n,r)
                         coeff1 = lambda*gamma1bar(n,r)/mu - rho0(n,r)
                         coeff2 = lambda*p0(n,r)/nu - rho0(n,r)

                         integral = (abs(grav_cell(n,r))/denom)* &
                              (coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
                                           (gamma1bar(n,r) - HALF*mu*dr(n))) - &
                               coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
                                           (p0(n,r) - HALF*nu*dr(n))) )

                      end if
                   endif

                   beta0_edge(n,r+1) = beta0_edge(n,r) * exp(-integral)
                   beta0(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))

                else ! r >= anelastic_cutoff_density

                   beta0(n,r) = beta0(n,r-1) * (rho0(n,r)/rho0(n,r-1))
                   beta0_edge(n,r+1) = 2.d0*beta0(n,r) - beta0_edge(n,r)
                endif

             end do

             if (n .gt. 0) then

                ! Compare the difference between beta0 at the top of level n to the
                ! corresponding point on level n-1
                offset = beta0_edge(n,r_end_coord(n,j)+1) &
                     - beta0_edge(n-1,(r_end_coord(n,j)+1)/2)

                do i=n-1,0,-1

                   refrat = 2**(n-i)

                   ! Offset the centered beta on level i above this point so the total
                   ! integral is consistent
                   do r=r_end_coord(n,j)/refrat+1,nr(i)
                      beta0(i,r) = beta0(i,r) + offset
                   end do

                   ! Redo the anelastic cutoff part
                   do r=anelastic_cutoff_density_coord(i),nr(i)
                      if (rho0(i,r-1) /= ZERO) then
                         beta0(i,r) = beta0(i,r-1) * (rho0(i,r)/rho0(i,r-1))
                      endif
                   end do

                   ! This next piece of coded is needed for the case when the anelastic
                   ! cutoff coordinate lives on level n.  We first average beta0 from
                   ! level i+1 to level i in the region between the anelastic cutoff and
                   ! the top of grid n.  Then recompute beta0 at level i above the top
                   ! of grid n.
                   if (r_end_coord(n,j) .ge. anelastic_cutoff_density_coord(n)) then

                      do r=anelastic_cutoff_density_coord(i),(r_end_coord(n,j)+1)/refrat-1
                         beta0(i,r) = HALF*(beta0(i+1,2*r)+beta0(i+1,2*r+1))
                      end do

                      do r=(r_end_coord(n,j)+1)/refrat,nr(i)
                         if (rho0(i,r-1) /= ZERO) then
                            beta0(i,r) = beta0(i,r-1) * (rho0(i,r)/rho0(i,r-1))
                         endif
                      end do

                   end if

                end do ! end loop over i=n-1,0,-1

             end if ! end if (n .gt. 0)

          end do ! end loop over disjoint chunks

       end do ! end loop over levels

       ! zero the beta0 where there is no corresponding full state array
       do n=1,finest_radial_level
          do j=1,numdisjointchunks(n)
             if (j .eq. numdisjointchunks(n)) then
                do r=r_end_coord(n,j)+1,nr(n)-1
                   beta0(n,r) = ZERO
                end do
             else
                do r=r_end_coord(n,j)+1,r_start_coord(n,j+1)-1
                   beta0(n,r) = ZERO
                end do
             end if
          end do
       end do

    else if (beta0_type .eq. 2) then

       ! beta_0 = rho_0
       do n=0,finest_radial_level
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                beta0(n,r) = rho0(n,r)
             end do
          end do
       end do

    else if (beta0_type .eq. 3) then

       ! beta_0 = 1.d0
       do n=0,finest_radial_level
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                beta0 = 1.d0
             end do
          end do
       end do

    end if

    call restrict_base(beta0,1)
    call fill_ghost_base(beta0,1)

    deallocate(beta0_edge)

    call bl_proffortfuncstop("Maestro::make_beta0")

  end subroutine make_beta0

  subroutine make_beta0_irreg(beta0,rho0,p0,gamma1bar,grav_cell, &
                               r_cc_loc,r_edge_loc) &
                               bind(C, name="make_beta0_irreg")

    double precision, intent(  out) :: beta0    (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0     (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0       (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: grav_cell(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

    ! local
    integer :: r, n, i, refrat, j
    double precision :: integral
    double precision :: lambda, mu, nu, kappa
    double precision :: denom, coeff1, coeff2, coeff3
    double precision :: del,dpls,dmin,slim,sflag
    double precision :: offset
    double precision :: drp, drm, drc

    double precision, allocatable :: beta0_edge(:,:)

    call bl_proffortfuncstart("Maestro::make_beta0_irreg")

    allocate(beta0_edge(0:finest_radial_level,0:nr_fine))

    beta0 = 0.d0

    if (beta0_type .eq. 1) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Compute beta0 on the edges and average to the center
       !
       ! Multilevel Outline:
       !
       ! First, compute beta0 on edges and centers at level 0 only
       ! Obtain the starting value from rho0 at the bottom of the domain.
       ! do n=1,finest_radial_level
       !   Compute beta0 on edges and centers at level n
       !   Obtain the starting value of beta0_edge_lo from the coarser grid
       !   if n>0, compare the difference between beta0 at the top of level n to the
       !           corresponding point on level n-1
       !   do i=n-1,0,-1
       !     Offset the centered beta on level i above this point so the total integral
       !      is consistent
       !     Redo the anelastic cutoff part
       !   end do
       ! end do
       ! call restrict_base and fill_ghost_base
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do n=0,finest_radial_level

          do j=1,numdisjointchunks(n)

             ! Compute beta0 on edges and centers at level n

             if (n .eq. 0) then
                beta0_edge(0,0) = rho0(0,0)
             else
                ! Obtain the starting value of beta0_edge_lo from the coarser grid
                beta0_edge(n,r_start_coord(n,j)) = beta0_edge(n-1,r_start_coord(n,j)/2)
             end if

             do r=r_start_coord(n,j),r_end_coord(n,j)

                if (r < anelastic_cutoff_density_coord(n)) then

                   drp = r_edge_loc(n,r+1) - r_edge_loc(n,r)
                   drm = r_edge_loc(n,r) - r_edge_loc(n,r-1)

                   if (r .eq. 0 .or. r .eq. nr(n)-1) then

                      lambda = ZERO
                      mu = ZERO
                      nu = ZERO

                   else

                      drc = r_cc_loc(n,r+1) - r_cc_loc(n,r-1)

                      ! piecewise linear reconstruction of rho0,
                      ! gamma1bar, and p0 -- see paper III, appendix C
                      del    =       (rho0(n,r+1) - rho0(n,r-1))/drc
                      dpls   = TWO * (rho0(n,r+1) - rho0(n,r  ))/drp
                      dmin   = TWO * (rho0(n,r  ) - rho0(n,r-1))/drm
                      slim   = min(abs(dpls), abs(dmin))
                      slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag  = sign(ONE,del)
                      lambda = sflag*min(slim,abs(del))

                      del   =       (gamma1bar(n,r+1) - gamma1bar(n,r-1))/drc
                      dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/drp
                      dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/drm
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      mu    = sflag*min(slim,abs(del))

                      del   =       (p0(n,r+1) - p0(n,r-1))/drc
                      dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/drp
                      dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/drm
                      slim  = min(abs(dpls), abs(dmin))
                      slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
                      sflag = sign(ONE,del)
                      nu    = sflag*min(slim,abs(del))

                   end if

                   ! edge-to-cell-center spacings (~ dr/2)
                   drp = r_edge_loc(n,r+1) - r_cc_loc(n,r)
                   drm = r_cc_loc(n,r) - r_edge_loc(n,r)

                   if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
                        (nu*gamma1bar(n,r) - mu*p0(n,r)) .eq. ZERO .or. &
                        ((gamma1bar(n,r) + mu*drp)/ &
                        (gamma1bar(n,r) - mu*drm)) .le. ZERO .or. &
                        ((p0(n,r) + nu*drp)/ &
                        (p0(n,r) - nu*drm)) .le. ZERO) then

                      ! just do piecewise constant integration
                      integral = abs(grav_cell(n,r))*rho0(n,r)*(drp+drm)/(p0(n,r)*gamma1bar(n,r))

                   else

                      ! paper III, equation C2
                      denom = nu*gamma1bar(n,r) - mu*p0(n,r)
                      coeff1 = lambda*gamma1bar(n,r)/mu - rho0(n,r)
                      coeff2 = lambda*p0(n,r)/nu - rho0(n,r)

                      integral = (abs(grav_cell(n,r))/denom)* &
                           (coeff1*log( (gamma1bar(n,r) + mu*drp)/ &
                                        (gamma1bar(n,r) - mu*drm)) - &
                            coeff2*log( (p0(n,r) + nu*drp)/ &
                                        (p0(n,r) - nu*drm)) )

                   endif

                   beta0_edge(n,r+1) = beta0_edge(n,r) * exp(-integral)
                   beta0(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))

                else ! r >= anelastic_cutoff_density
                   beta0(n,r) = beta0(n,r-1) * (rho0(n,r)/rho0(n,r-1))
                   beta0_edge(n,r+1) = 2.d0*beta0(n,r) - beta0_edge(n,r)
                endif

             end do

             if (n .gt. 0) then

                ! Compare the difference between beta0 at the top of level n to the
                ! corresponding point on level n-1
                offset = beta0_edge(n,r_end_coord(n,j)+1) &
                     - beta0_edge(n-1,(r_end_coord(n,j)+1)/2)

                do i=n-1,0,-1

                   refrat = 2**(n-i)

                   ! Offset the centered beta on level i above this point so the total
                   ! integral is consistent
                   do r=r_end_coord(n,j)/refrat+1,nr(i)
                      beta0(i,r) = beta0(i,r) + offset
                   end do

                   ! Redo the anelastic cutoff part
                   do r=anelastic_cutoff_density_coord(i),nr(i)
                      if (rho0(i,r-1) /= ZERO) then
                         beta0(i,r) = beta0(i,r-1) * (rho0(i,r)/rho0(i,r-1))
                      endif
                   end do

                   ! This next piece of coded is needed for the case when the anelastic
                   ! cutoff coordinate lives on level n.  We first average beta0 from
                   ! level i+1 to level i in the region between the anelastic cutoff and
                   ! the top of grid n.  Then recompute beta0 at level i above the top
                   ! of grid n.
                   if (r_end_coord(n,j) .ge. anelastic_cutoff_density_coord(n)) then

                      do r=anelastic_cutoff_density_coord(i),(r_end_coord(n,j)+1)/refrat-1
                         beta0(i,r) = HALF*(beta0(i+1,2*r)+beta0(i+1,2*r+1))
                      end do

                      do r=(r_end_coord(n,j)+1)/refrat,nr(i)
                         if (rho0(i,r-1) /= ZERO) then
                            beta0(i,r) = beta0(i,r-1) * (rho0(i,r)/rho0(i,r-1))
                         endif
                      end do

                   end if

                end do ! end loop over i=n-1,0,-1

             end if ! end if (n .gt. 0)

          end do ! end loop over disjoint chunks

       end do ! end loop over levels

       ! zero the beta0 where there is no corresponding full state array
       do n=1,finest_radial_level
          do j=1,numdisjointchunks(n)
             if (j .eq. numdisjointchunks(n)) then
                do r=r_end_coord(n,j)+1,nr(n)-1
                   beta0(n,r) = ZERO
                end do
             else
                do r=r_end_coord(n,j)+1,r_start_coord(n,j+1)-1
                   beta0(n,r) = ZERO
                end do
             end if
          end do
       end do

    else if (beta0_type .eq. 2) then

       ! beta_0 = rho_0
       do n=0,finest_radial_level
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                beta0(n,r) = rho0(n,r)
             end do
          end do
       end do

    else if (beta0_type .eq. 3) then

       ! beta_0 = 1.d0
       do n=0,finest_radial_level
          do j=1,numdisjointchunks(n)
             do r=r_start_coord(n,j),r_end_coord(n,j)
                beta0 = 1.d0
             end do
          end do
       end do

    end if

    call restrict_base(beta0,1)
    call fill_ghost_base(beta0,1)

    deallocate(beta0_edge)

    call bl_proffortfuncstop("Maestro::make_beta0_irreg")

  end subroutine make_beta0_irreg

end module make_beta0_module
