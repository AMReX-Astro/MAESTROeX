module make_edge_state_module
  ! make_edge_state is for reconstruction of the base state variables to
  ! predict the time-centered edge states.

  use amrex_constants_module
  use amrex_error_module
  use meth_params_module, only: spherical, rel_eps, slope_order, ppm_type
  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr, numdisjointchunks, &
                                            finest_radial_level, r_start_coord, r_end_coord

  implicit none

contains

  subroutine make_edge_state_1d(s,sedge,w0,force,dt)

    double precision, intent(in   ) ::     s(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: sedge(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) ::    w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: force(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dt

    call bl_proffortfuncstart("Maestro::make_edge_state_1d")

    if (spherical .eq. 1) then
       call make_edge_state_1d_sphr(s,sedge,w0,force,dt)
    else
       call make_edge_state_1d_planar(s,sedge,w0,force,dt)
    end if

    call bl_proffortfuncstop("Maestro::make_edge_state_1d")

  end subroutine make_edge_state_1d

  subroutine make_edge_state_1d_sphr(s,sedge,w0,force,dt)

    double precision, intent(in   ) ::     s(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: sedge(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) ::    w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: force(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dt

    double precision :: dmin,dpls,ds,del,slim,sflag,ubardth,dth,savg,u
    double precision :: sigmap,sigmam,s6,D2,D2L,D2R,D2C,D2LIM,alphap,alpham,sgn
    double precision :: dafacem,dafacep,dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp
    double precision :: amax,delam,delap,D2ABS

    integer :: r,lo,hi

    logical :: extremum, bigp, bigm

    integer        , parameter :: cen=1, lim=2, flag=3, fromm=4
    double precision, parameter :: FOURTHIRDS = FOUR/THREE

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! cell based indexing
    double precision :: dsscr(-1:nr_fine,4)
    double precision :: slope(0:nr_fine-1)
    double precision ::    sp(0:nr_fine-1)
    double precision ::    sm(0:nr_fine-1)
    double precision ::    Ip(0:nr_fine-1)
    double precision ::    Im(0:nr_fine-1)

    ! edge based indexing
    double precision :: sedgel(-1:nr_fine+1)
    double precision :: sedger(-1:nr_fine+1)

    double precision :: s_ghost(-3:nr_fine+2)

    ! copy valid data into array with ghost cells
    s_ghost(0:nr_fine-1) = s(0,0:nr_fine-1)

    ! symmetry boundary condition at center
    s_ghost(-1) = s(0,0)
    s_ghost(-2) = s(0,1)
    s_ghost(-3) = s(0,2)

    ! first-order extrapolation at top of star
    s_ghost(nr_fine  ) = s(0,nr_fine-1)
    s_ghost(nr_fine+1) = s(0,nr_fine-1)
    s_ghost(nr_fine+2) = s(0,nr_fine-1)

    ! Need to initialize this because it's not always set.
    dsscr = ZERO

     dth = HALF*dt

     lo = 0
     hi = nr_fine-1

     if (ppm_type .eq. 0) then

        ! compute slopes
        if (slope_order .eq. 0) then

           slope = ZERO

        else if (slope_order .eq. 2) then

           do r=lo,hi
              ! do standard limiting on interior cells
              del = half*(s_ghost(r+1) - s_ghost(r-1))
              dpls = two*(s_ghost(r+1) - s_ghost(r  ))
              dmin = two*(s_ghost(r  ) - s_ghost(r-1))
              slim = min(abs(dpls), abs(dmin))
              slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
              sflag = sign(ONE,del)
              slope(r)= sflag*min(slim,abs(del))
           end do

        else if (slope_order .eq. 4) then

           do r=lo-1,hi+1
              ! do standard limiting to compute temporary slopes
              dsscr(r,cen) = half*(s_ghost(r+1)-s_ghost(r-1))
              dpls = two*(s_ghost(r+1)-s_ghost(r  ))
              dmin = two*(s_ghost(r  )-s_ghost(r-1))
              dsscr(r,lim)= min(abs(dmin),abs(dpls))
              dsscr(r,lim) = merge(dsscr(r,lim),ZERO,dpls*dmin.gt.ZERO)
              dsscr(r,flag) = sign(ONE,dsscr(r,cen))
              dsscr(r,fromm)= dsscr(r,flag)*min(dsscr(r,lim),abs(dsscr(r,cen)))
           end do

           do r=lo,hi
              ! fourth-order limited slopes
              ds = FOURTHIRDS*dsscr(r,cen) - SIXTH*(dsscr(r+1,fromm)+dsscr(r-1,fromm))
              slope(r) = dsscr(r,flag)*min(abs(ds),dsscr(r,lim))
           end do

        end if ! which slope order

        ! compute sedgel and sedger
        do r=lo,hi
           u = HALF*(w0(0,r)+w0(0,r+1))
           ubardth = dth*u/dr(0)  ! NOTE: ubardth=0 for use_exact_base_state case
           sedgel(r+1)= s(0,r) + (HALF-ubardth)*slope(r) + dth*force(0,r)
           sedger(r  )= s(0,r) - (HALF+ubardth)*slope(r) + dth*force(0,r)
        end do

     else if (ppm_type .eq. 1) then

        ! interpolate s to radial edges, store these temporary values into sedgel
        !$OMP PARALLEL DO PRIVATE(r,del,dmin,dpls)
        do r=lo-1,hi+1
           ! compute van Leer slopes
           del  = HALF * (s_ghost(r+1) - s_ghost(r-1))
           dmin = TWO  * (s_ghost(r  ) - s_ghost(r-1))
           dpls = TWO  * (s_ghost(r+1) - s_ghost(r  ))
           if (dmin*dpls .gt. ZERO) then
              dsscr(r,1) = sign(ONE,del)*min(abs(del),abs(dmin),abs(dpls))
           end if
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO PRIVATE(r)
        do r=lo,hi+1
           ! 4th order interpolation of s to radial faces
           sedgel(r) = &
                HALF*(s_ghost(r)+s_ghost(r-1))-SIXTH*(dsscr(r,1)-dsscr(r-1,1))
           ! make sure sedgel lies in between adjacent cell-centered values
           sedgel(r) = max(sedgel(r),min(s_ghost(r),s_ghost(r-1)))
           sedgel(r) = min(sedgel(r),max(s_ghost(r),s_ghost(r-1)))
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO PRIVATE(r)
        do r=lo,hi

           ! first copy sedgel into sp and sm
           sp(r) = sedgel(r+1)
           sm(r) = sedgel(r  )

           ! modify using quadratic limiters
           if ((sp(r)-s(0,r))*(s(0,r)-sm(r)) .le. ZERO) then
              sp(r) = s(0,r)
              sm(r) = s(0,r)
           else if (abs(sp(r)-s(0,r)) .ge. TWO*abs(sm(r)-s(0,r))) then
              sp(r) = THREE*s(0,r) - TWO*sm(r)
           else if (abs(sm(r)-s(0,r)) .ge. TWO*abs(sp(r)-s(0,r))) then
              sm(r) = THREE*s(0,r) - TWO*sp(r)
           end if

        end do
        !$OMP END PARALLEL DO

     else if (ppm_type .eq. 2) then

        ! interpolate s to radial edges, store these temporary values into sedgel
        !$OMP PARALLEL DO PRIVATE(r,D2,D2L,D2R,sgn,D2LIM)
        do r=lo-1,hi+2

           ! fourth-order stencil
           sedgel(r) = (7.d0/12.d0)*(s_ghost(r-1)+s_ghost(r)) &
                - (1.d0/12.d0)*(s_ghost(r-2)+s_ghost(r+1))

           ! limit sedge
           if ((sedgel(r)-s_ghost(r-1))*(s_ghost(r)-sedgel(r)) .lt. ZERO) then
              D2  = THREE*(s_ghost(r-1)-TWO*sedgel(r)+s_ghost(r))
              D2L = s_ghost(r-2)-TWO*s_ghost(r-1)+s_ghost(r)
              D2R = s_ghost(r-1)-TWO*s_ghost(r)+s_ghost(r+1)
              sgn = sign(ONE,D2)
              D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
              sedgel(r) = HALF*(s_ghost(r-1)+s_ghost(r)) - SIXTH*D2LIM
           end if

        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO PRIVATE(r,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
        !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM) &
        !$OMP PRIVATE(amax,delam,delap,D2ABS)
        do r=lo,hi

           ! use Colella 2008 limiters
           ! This is a new version of the algorithm
           ! to eliminate sensitivity to roundoff.
           alphap = sedgel(r+1)-s_ghost(r)
           alpham = sedgel(r  )-s_ghost(r)
           bigp = abs(alphap).gt.TWO*abs(alpham)
           bigm = abs(alpham).gt.TWO*abs(alphap)
           extremum = .false.

           if (alpham*alphap .ge. ZERO) then
              extremum = .true.
           else if (bigp .or. bigm) then
              ! Possible extremum. We look at cell centered values and face
              ! centered values for a change in sign in the differences adjacent to
              ! the cell. We use the pair of differences whose minimum magnitude is
              ! the largest, and thus least susceptible to sensitivity to roundoff.
              dafacem = sedgel(r) - sedgel(r-1)
              dafacep = sedgel(r+2) - sedgel(r+1)
              dabarm = s_ghost(r) - s_ghost(r-1)
              dabarp = s_ghost(r+1) - s_ghost(r)
              dafacemin = min(abs(dafacem),abs(dafacep))
              dabarmin= min(abs(dabarm),abs(dabarp))
              if (dafacemin.ge.dabarmin) then
                 dachkm = dafacem
                 dachkp = dafacep
              else
                 dachkm = dabarm
                 dachkp = dabarp
              endif
              extremum = (dachkm*dachkp .le. 0.d0)
           end if

           if (extremum) then
              D2  = SIX*(alpham + alphap)
              D2L = s_ghost(r-2)-TWO*s_ghost(r-1)+s_ghost(r)
              D2R = s_ghost(r)-TWO*s_ghost(r+1)+s_ghost(r+2)
              D2C = s_ghost(r-1)-TWO*s_ghost(r)+s_ghost(r+1)
              sgn = sign(ONE,D2)
              D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
              D2ABS = max(abs(D2),1.d-10)
              alpham = alpham*D2LIM/D2ABS
              alphap = alphap*D2LIM/D2ABS
           else
              if (bigp) then
                 sgn = sign(ONE,alpham)
                 amax = -alphap**2 / (4*(alpham + alphap))
                 delam = s_ghost(r-1) - s_ghost(r)
                 if (sgn*amax .ge. sgn*delam) then
                    if (sgn*(delam - alpham).ge.1.d-10) then
                       alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                    else
                       alphap = -TWO*alpham
                    endif
                 endif
              end if
              if (bigm) then
                 sgn = sign(ONE,alphap)
                 amax = -alpham**2 / (4*(alpham + alphap))
                 delap = s_ghost(r+1) - s_ghost(r)
                 if (sgn*amax .ge. sgn*delap) then
                    if (sgn*(delap - alphap).ge.1.d-10) then
                       alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                    else
                       alpham = -TWO*alphap
                    endif
                 endif
              end if
           end if

           sm(r) = s_ghost(r) + alpham
           sp(r) = s_ghost(r) + alphap

        end do ! loop over r
        !$OMP END PARALLEL DO

     end if

     if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

        !$OMP PARALLEL DO PRIVATE(r,sigmap,sigmam,s6)
        do r=lo,hi
           ! compute Ip and Im
           sigmap = abs(w0(0,r+1))*dt/dr(0)  ! NOTE: sigmap=0 for use_exact_base_state case
           sigmam = abs(w0(0,r  ))*dt/dr(0)  ! NOTE: sigmam=0 for use_exact_base_state case
           s6 = SIX*s(0,r) - THREE*(sm(r)+sp(r))
           if (w0(0,r+1) .gt. rel_eps) then
              Ip(r) = sp(r) - (sigmap/TWO)*(sp(r)-sm(r)-(ONE-TWO3RD*sigmap)*s6)
           else
              Ip(r) = s(0,r)
           end if
           if (w0(0,r) .lt. -rel_eps) then
              Im(r) = sm(r) + (sigmam/TWO)*(sp(r)-sm(r)+(ONE-TWO3RD*sigmam)*s6)
           else
              Im(r) = s(0,r)
           end if

           ! compute sedgel and sedger
           sedgel(r+1) = Ip(r) + dth*force(0,r)
           sedger(r  ) = Im(r) + dth*force(0,r)
        end do
        !$OMP END PARALLEL DO

     end if

     ! Fix center and edge of star by reflecting the extrapolated state.
     ! An alternate way would be to compute these values using the entire algorithm,
     ! but that would require more ghost cells at several stages.
     ! By symmetry arguments, this would make no difference at the center of the star
     ! and the accuracy at the edge of the star is not important here
     sedgel(0)       = sedger(0)
     sedger(nr_fine) = sedgel(nr_fine)

     !$OMP PARALLEL DO PRIVATE(r,savg)
     do r=lo,hi+1
        ! solve Riemann problem to get final edge state
        sedge(0,r)=merge(sedgel(r),sedger(r),w0(0,r).gt.ZERO)
        savg = HALF*(sedger(r) + sedgel(r))
        sedge(0,r)=merge(savg,sedge(0,r),abs(w0(0,r)) .lt. rel_eps)
     end do
     !$OMP END PARALLEL DO

   end subroutine make_edge_state_1d_sphr

   subroutine make_edge_state_1d_planar(s,sedge,w0,force,dt)

     double precision, intent(in   ) ::     s(0:max_radial_level,0:nr_fine-1)
     double precision, intent(inout) :: sedge(0:max_radial_level,0:nr_fine)
     double precision, intent(in   ) ::    w0(0:max_radial_level,0:nr_fine)
     double precision, intent(in   ) :: force(0:max_radial_level,0:nr_fine-1)
     double precision, intent(in   ) :: dt

     double precision :: dmin,dpls,ds,del,slim,sflag,ubardth,dth,dtdr,savg,u
     double precision :: sigmap,sigmam,s6,D2,D2L,D2R,D2C,D2LIM,alphap,alpham,sgn
     double precision :: dafacem,dafacep,dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp
     double precision :: amax,delam,delap,D2ABS

     integer :: r,lo,hi,n,i

     logical :: extremum, bigp, bigm

     integer        , parameter :: cen=1, lim=2, flag=3, fromm=4
     double precision, parameter :: FOURTHIRDS = FOUR/THREE

     ! constant used in Colella 2008
     double precision, parameter :: C = 1.25d0

     ! cell based indexing
     double precision :: slope(0:finest_radial_level,0:nr_fine-1)
     double precision :: dxscr(0:finest_radial_level,0:nr_fine-1,4)
     double precision ::  dsvl(0:finest_radial_level,-1:nr_fine)
     double precision ::    sp(0:finest_radial_level,0:nr_fine-1)
     double precision ::    sm(0:finest_radial_level,0:nr_fine-1)
     double precision ::    Ip(0:finest_radial_level,0:nr_fine-1)
     double precision ::    Im(0:finest_radial_level,0:nr_fine-1)

     ! edge based indexing
     double precision :: sedgel(0:finest_radial_level,0:nr_fine)
     double precision :: sedger(0:finest_radial_level,0:nr_fine)

     dth = HALF*dt
     dtdr = dt/dr(0)

     dsvl = ZERO

     ! error checking to make sure that there is a 2 cell buffer at the top and bottom
     ! of the domain for finer levels in planar geometry.  This can be removed if
     ! blocking_factor is implemented at set > 1.
     if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
        do n=1,finest_radial_level
           do i=1,numdisjointchunks(n)
              if (r_start_coord(n,i) .eq. 2) then
                 call amrex_error("make_edge_state assumes blocking_factor > 1 at lo boundary")
              else if (r_end_coord(n,i) .eq. nr(n)-3) then
                 call amrex_error("make_edge_state assumes blocking_factor > 1 at hi boundary")
              end if
           end do
        end do
     end if

     if (ppm_type .eq. 0) then

        ! compute slopes
        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              if (slope_order .eq. 0) then

                 slope(n,:) = ZERO

              else if (slope_order .eq. 2) then

                 do r=lo,hi
                    if (r .eq. 0) then
                       ! one-sided difference
                       slope(n,r) = s(n,r+1)-s(n,r)
                    else if (r .eq. nr(n)-1) then
                       ! one-sided difference
                       slope(n,r) = s(n,r)-s(n,r-1)
                    else
                       ! do standard limiting on interior cells
                       del = half*(s(n,r+1) - s(n,r-1))
                       dpls = two*(s(n,r+1) - s(n,r  ))
                       dmin = two*(s(n,r  ) - s(n,r-1))
                       slim = min(abs(dpls), abs(dmin))
                       slim = merge(slim, ZERO, dpls*dmin.gt.ZERO)
                       sflag = sign(ONE,del)
                       slope(n,r)= sflag*min(slim,abs(del))
                    end if
                 end do

              else if (slope_order .eq. 4) then

                 do r=lo-1,hi+1
                    if (r .eq. 0) then
                       ! one-sided difference
                       dxscr(n,r,fromm) = s(n,r+1)-s(n,r)
                    else if (r .eq. nr(n)-1) then
                       ! one-sided difference
                       dxscr(n,r,fromm) = s(n,r)-s(n,r-1)
                    else if (r .gt. 0 .and. r .lt. nr(n)-1) then
                       ! do standard limiting to compute temporary slopes
                       dxscr(n,r,cen) = half*(s(n,r+1)-s(n,r-1))
                       dpls = two*(s(n,r+1)-s(n,r  ))
                       dmin = two*(s(n,r  )-s(n,r-1))
                       dxscr(n,r,lim)= min(abs(dmin),abs(dpls))
                       dxscr(n,r,lim) = merge(dxscr(n,r,lim),ZERO,dpls*dmin.gt.ZERO)
                       dxscr(n,r,flag) = sign(ONE,dxscr(n,r,cen))
                       dxscr(n,r,fromm)= dxscr(n,r,flag) &
                            *min(dxscr(n,r,lim),abs(dxscr(n,r,cen)))
                    end if
                 end do

                 do r=lo,hi
                    if (r .eq. 0) then
                       ! one-sided difference
                       slope(n,r) = s(n,r+1)-s(n,r)
                    else if (r .eq. nr(n)-1) then
                       ! one-sided difference
                       slope(n,r) = s(n,r)-s(n,r-1)
                    else
                       ! fourth-order limited slopes on interior
                       ds = FOURTHIRDS*dxscr(n,r,cen) &
                            - SIXTH*(dxscr(n,r+1,fromm) + dxscr(n,r-1,fromm))
                       slope(n,r) = dxscr(n,r,flag)*min(abs(ds),dxscr(n,r,lim))
                    end if
                 end do

              end if ! which slope order

           end do ! loop over disjointchunks
        end do ! loop over levels

        ! compute sedgel and sedger
        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              do r = lo,hi
                 u = HALF*(w0(n,r)+w0(n,r+1))
                 ubardth = dth*u/dr(n)
                 sedgel(n,r+1)= s(n,r) + (HALF-ubardth)*slope(n,r) + dth * force(n,r)
                 sedger(n,r  )= s(n,r) - (HALF+ubardth)*slope(n,r) + dth * force(n,r)
              end do

           end do ! loop over disjointchunks
        end do ! loop over levels

     else if (ppm_type .eq. 1) then

        ! interpolate s to radial edges, store these temporary values into sedgel
        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              ! compute van Leer slopes
              !$OMP PARALLEL DO PRIVATE(r,del,dmin,dpls)
              do r=lo-1,hi+1
                 if (r .eq. 0) then
                    ! one-sided difference
                    dsvl(n,r) = s(n,r+1)-s(n,r)
                 else if (r .eq. nr(n)-1) then
                    ! one-sided difference
                    dsvl(n,r) = s(n,r)-s(n,r-1)
                 else if (r .gt. 0 .and. r .lt. nr(n)-1) then
                    del  = HALF * (s(n,r+1) - s(n,r-1))
                    dmin = TWO  * (s(n,r  ) - s(n,r-1))
                    dpls = TWO  * (s(n,r+1) - s(n,r  ))
                    if (dmin*dpls .gt. ZERO) &
                         dsvl(n,r) = sign(ONE,del)*min(abs(del),abs(dmin),abs(dpls))
                 end if
              end do
              !$OMP END PARALLEL DO

              !$OMP PARALLEL DO PRIVATE(r)
              do r=lo,hi+1
                 if (r .eq. 0) then
                    ! 2nd order interpolation to boundary face
                    sedgel(n,r) = s(n,r) - half*dsvl(n,r)
                 else if (r .eq. nr(n)) then
                    ! 2nd order interpolation to boundary face
                    sedgel(n,r) = s(n,r-1) + half*dsvl(n,r)
                 else
                    ! 4th order interpolation of s to radial faces
                    sedgel(n,r) = HALF*(s(n,r)+s(n,r-1)) - SIXTH*(dsvl(n,r)-dsvl(n,r-1))
                    ! make sure sedgel lies in between adjacent cell-centered values
                    sedgel(n,r) = max(sedgel(n,r),min(s(n,r),s(n,r-1)))
                    sedgel(n,r) = min(sedgel(n,r),max(s(n,r),s(n,r-1)))
                 end if
              end do
              !$OMP END PARALLEL DO

           end do ! loop over disjointchunks
        end do ! loop over levels

        ! copy sedgel into sp and sm
        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              !$OMP PARALLEL DO PRIVATE(r)
              do r=lo,hi

                 sp(n,r) = sedgel(n,r+1)
                 sm(n,r) = sedgel(n,r  )

                 ! modify using quadratic limiters
                 if ((sp(n,r)-s(n,r))*(s(n,r)-sm(n,r)) .le. ZERO) then
                    sp(n,r) = s(n,r)
                    sm(n,r) = s(n,r)
                 else if (abs(sp(n,r)-s(n,r)) .ge. TWO*abs(sm(n,r)-s(n,r))) then
                    sp(n,r) = THREE*s(n,r) - TWO*sm(n,r)
                 else if (abs(sm(n,r)-s(n,r)) .ge. TWO*abs(sp(n,r)-s(n,r))) then
                    sm(n,r) = THREE*s(n,r) - TWO*sp(n,r)
                 end if
              end do
              !$OMP END PARALLEL DO

           end do ! loop over disjointchunks
        end do ! loop over levels

     else if (ppm_type .eq. 2) then

        ! interpolate s to radial edges, store these temporary values into sedgel
        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              ! store centered differences in dsvl
              !$OMP PARALLEL DO PRIVATE(r)
              do r=lo-3,hi+3
                 if (r .eq. 0) then
                    ! one-sided difference
                    dsvl(n,r) = s(n,r+1)-s(n,r)
                 else if (r .eq. nr(n)-1) then
                    ! one-sided difference
                    dsvl(n,r) = s(n,r)-s(n,r-1)
                 else if (r .gt. 0 .and. r .lt. nr(n)-1) then
                    ! centered difference
                    dsvl(n,r) = HALF * (s(n,r+1) - s(n,r-1))
                 end if
              end do
              !$OMP END PARALLEL DO

              !$OMP PARALLEL DO PRIVATE(r,D2,D2L,D2R,sgn,D2LIM)
              do r=lo-2,hi+3
                 if (r .eq. 0) then
                    ! 2nd order interpolation to boundary face
                    sedgel(n,r) = s(n,r) - half*dsvl(n,r)
                 else if (r .eq. nr(n)) then
                    ! 2nd order interpolation to boundary face
                    sedgel(n,r) = s(n,r-1) + half*dsvl(n,r)
                 else if (r .gt. 0 .and. r .lt. nr(n)) then
                    ! 4th order interpolation of s to radial faces
                    sedgel(n,r) = HALF*(s(n,r)+s(n,r-1)) - SIXTH*(dsvl(n,r)-dsvl(n,r-1))
                    if (r .ge. 2 .and. r .le. nr(n)-2) then
                       ! limit sedge
                       if ((sedgel(n,r)-s(n,r-1))*(s(n,r)-sedgel(n,r)) .lt. ZERO) then
                          D2  = THREE*(s(n,r-1)-TWO*sedgel(n,r)+s(n,r))
                          D2L = s(n,r-2)-TWO*s(n,r-1)+s(n,r)
                          D2R = s(n,r-1)-TWO*s(n,r)+s(n,r+1)
                          sgn = sign(ONE,D2)
                          D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                          sedgel(n,r) = HALF*(s(n,r-1)+s(n,r)) - SIXTH*D2LIM
                       end if
                    end if
                 end if
              end do
              !$OMP END PARALLEL DO

           end do ! loop over disjointchunks
        end do ! loop over levels

        ! use Colella 2008 limiters
        ! This is a new version of the algorithm
        ! to eliminate sensitivity to roundoff.
        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              !$OMP PARALLEL DO PRIVATE(r,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep) &
              !$OMP PRIVATE(dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM) &
              !$OMP PRIVATE(amax,delam,delap,D2ABS)
              do r=lo,hi

                 if (r .ge. 2 .and. r .le. nr(n)-3) then

                    alphap = sedgel(n,r+1)-s(n,r)
                    alpham = sedgel(n,r  )-s(n,r)
                    bigp = abs(alphap).gt.TWO*abs(alpham)
                    bigm = abs(alpham).gt.TWO*abs(alphap)
                    extremum = .false.

                    if (alpham*alphap .ge. ZERO) then
                       extremum = .true.
                    else if (bigp .or. bigm) then
                       ! Possible extremum. We look at cell centered values and face
                       ! centered values for a change in sign in the differences adjacent to
                       ! the cell. We use the pair of differences whose minimum magnitude is
                       ! the largest, and thus least susceptible to sensitivity to roundoff.
                       dafacem = sedgel(n,r) - sedgel(n,r-1)
                       dafacep = sedgel(n,r+2) - sedgel(n,r+1)
                       dabarm = s(n,r) - s(n,r-1)
                       dabarp = s(n,r+1) - s(n,r)
                       dafacemin = min(abs(dafacem),abs(dafacep))
                       dabarmin= min(abs(dabarm),abs(dabarp))
                       if (dafacemin.ge.dabarmin) then
                          dachkm = dafacem
                          dachkp = dafacep
                       else
                          dachkm = dabarm
                          dachkp = dabarp
                       endif
                       extremum = (dachkm*dachkp .le. 0.d0)
                    end if

                    if (extremum) then
                       D2  = SIX*(alpham + alphap)
                       D2L = s(n,r-2)-TWO*s(n,r-1)+s(n,r)
                       D2R = s(n,r)-TWO*s(n,r+1)+s(n,r+2)
                       D2C = s(n,r-1)-TWO*s(n,r)+s(n,r+1)
                       sgn = sign(ONE,D2)
                       D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                       D2ABS = max(abs(D2),1.d-10)
                       alpham = alpham*D2LIM/D2ABS
                       alphap = alphap*D2LIM/D2ABS
                    else
                       if (bigp) then
                          sgn = sign(ONE,alpham)
                          amax = -alphap**2 / (4*(alpham + alphap))
                          delam = s(n,r-1) - s(n,r)
                          if (sgn*amax .ge. sgn*delam) then
                             if (sgn*(delam - alpham).ge.1.d-10) then
                                alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                             else
                                alphap = -TWO*alpham
                             endif
                          endif
                       end if
                       if (bigm) then
                          sgn = sign(ONE,alphap)
                          amax = -alpham**2 / (4*(alpham + alphap))
                          delap = s(n,r+1) - s(n,r)
                          if (sgn*amax .ge. sgn*delap) then
                             if (sgn*(delap - alphap).ge.1.d-10) then
                                alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                             else
                                alpham = -TWO*alphap
                             endif
                          endif
                       end if
                    end if

                    sm(n,r) = s(n,r) + alpham
                    sp(n,r) = s(n,r) + alphap

                 else

                    sp(n,r) = sedgel(n,r+1)
                    sm(n,r) = sedgel(n,r  )

                 end if ! test (r .ge. 2 .and. r .le. nr(n)-3)

              end do ! loop over r
              !$OMP END PARALLEL DO

           end do ! loop over disjointchunks
        end do ! loop over levels

     end if

     ! compute Ip and Im
     if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then

        do n=0,finest_radial_level
           do i=1,numdisjointchunks(n)

              lo = r_start_coord(n,i)
              hi = r_end_coord(n,i)

              !$OMP PARALLEL DO PRIVATE(r,sigmap,sigmam,s6)
              do r=lo,hi
                 sigmap = abs(w0(n,r+1))*dtdr
                 sigmam = abs(w0(n,r  ))*dtdr
                 s6 = SIX*s(n,r) - THREE*(sm(n,r)+sp(n,r))
                 if (w0(n,r+1) .gt. rel_eps) then
                    Ip(n,r) = sp(n,r) - (sigmap/TWO)*(sp(n,r)-sm(n,r)-(ONE-TWO3RD*sigmap)*s6)
                 else
                    Ip(n,r) = s(n,r)
                 end if
                 if (w0(n,r) .lt. -rel_eps) then
                    Im(n,r) = sm(n,r) + (sigmam/TWO)*(sp(n,r)-sm(n,r)+(ONE-TWO3RD*sigmam)*s6)
                 else
                    Im(n,r) = s(n,r)
                 end if

                 ! compute sedgel and sedger
                 sedgel(n,r+1) = Ip(n,r) + dth * force(n,r)
                 sedger(n,r  ) = Im(n,r) + dth * force(n,r)
              end do
              !$OMP END PARALLEL DO

           end do ! loop over disjointchunks
        end do ! loop over levels

     end if

     ! sync up edge states at coarse-fine interface
     do n=0,finest_radial_level
        do i=1,numdisjointchunks(n)

           lo = r_start_coord(n,i)
           hi = r_end_coord(n,i)

           ! if we are not at the finest level, copy in the sedger and sedgel states
           ! from the next finer level at the c-f interface
           if (n .ne. finest_radial_level) then
              sedger(n,r_start_coord(n+1,i)/2) = sedger(n+1,r_start_coord(n+1,i))
              sedgel(n,(r_end_coord(n+1,i)+1)/2) = sedgel(n+1,r_end_coord(n+1,i)+1)
           end if

           ! if we are not at the coarsest level, copy in the sedgel and sedger states
           ! from the next coarser level at the c-f interface
           if (n .ne. 0) then
              sedgel(n,lo) = sedgel(n-1,lo/2)
              sedger(n,hi+1) = sedger(n-1,(hi+1)/2)
           end if

        end do ! loop over disjointchunks
     end do ! loop over levels

     ! solve Riemann problem to get final edge state
     do n=0,finest_radial_level
        do i=1,numdisjointchunks(n)

           lo = r_start_coord(n,i)
           hi = r_end_coord(n,i)

           !$OMP PARALLEL DO PRIVATE(r,savg)
           do r=lo,hi+1
              if (r .eq. 0) then
                 ! pick interior state at lo domain boundary
                 sedge(n,r) = sedger(n,r)
              else if (r .eq. nr(n)) then
                 ! pick interior state at hi domain boundary
                 sedge(n,r) = sedgel(n,r)
              else
                 ! upwind
                 sedge(n,r)=merge(sedgel(n,r),sedger(n,r),w0(n,r).gt.ZERO)
                 savg = HALF*(sedger(n,r) + sedgel(n,r))
                 sedge(n,r)=merge(savg,sedge(n,r),abs(w0(n,r)) .lt. rel_eps)
              end if
           end do
           !$OMP END PARALLEL DO

        end do  ! loop over disjointchunks
     end do ! loop over levels

   end subroutine make_edge_state_1d_planar

 end module make_edge_state_module
