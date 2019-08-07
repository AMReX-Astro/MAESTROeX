! diag module for MAESTROeX.
! currently, there are 3 output files:
!
!   diag_enuc.out:
!          peak nuc energy generation rate (erg / g / s)
!          x/y/z location of peak enuc
!          velocity components at location of peak enuc (including vr)
!          radius of peak enuc
!          total nuclear energy release (erg / s)
!
!   diag_temp.out:
!          peak temperature
!          x/y/z location of peak temperature
!          velocity components at location of peak temperature (including vr)
!          radius of peak temperature
!          central temperature
!
!   diag_vel.out:
!          peak total velocity
!          peak Mach number
!          total kinetic energy
!          gravitational potential energy
!          total internal energy
!          velocity components at the center of the star
!
module diag_module

  use meth_params_module, only: rho_comp, spec_comp, temp_comp, prob_lo, &
       sponge_start_factor, sponge_center_density, base_cutoff_density
  use network, only: nspec
  use base_state_geometry_module, only:  max_radial_level, nr_fine, center
  use eos_module
  use eos_type_module
  use amrex_constants_module
  use fundamental_constants_module, only: Gconst

  implicit none

  private

contains

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine diag(lev, lo, hi, &
                  scal, s_lo, s_hi, nc_s, &
                  rho_Hnuc, hn_lo, hn_hi, &
                  rho_Hext, he_lo, he_hi, &
                  rho0, p0, &
                  u, u_lo, u_hi, &
                  w0, dx, &
                  Mach_max,temp_max,enuc_max,Hext_max, &
                  mask,     m_lo, m_hi, use_mask) &
                  bind(C, name="diag")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    double precision, intent(in   ) ::     scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    integer         , intent(in   ) :: hn_lo(3), hn_hi(3)
    double precision, intent(in   ) :: rho_Hnuc(hn_lo(1):hn_hi(1),hn_lo(2):hn_hi(2),hn_lo(3):hn_hi(3))
    integer         , intent(in   ) :: he_lo(3), he_hi(3)
    double precision, intent(in   ) :: rho_Hext(he_lo(1):he_hi(1),he_lo(2):he_hi(2),he_lo(3):he_hi(3))
    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: p0(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) ::  u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(inout) :: Mach_max, temp_max, enuc_max, Hext_max
    integer         , intent(in   ) :: m_lo(3), m_hi(3)
    integer         , intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    integer         , intent(in   ) :: use_mask

    !     Local variables
    integer            :: i, j, k
    double precision   :: weight
    double precision   :: x, y, z
    double precision   :: vel
    logical            :: cell_valid

    type (eos_t) :: eos_state


    ! weight is the factor by which the volume of a cell at the current level
    ! relates to the volume of a cell at the coarsest level of refinement.
#if (AMREX_SPACEDIM == 2)
    weight = 1.d0 / 4.d0**(lev)
#elif (AMREX_SPACEDIM == 3)
    weight = 1.d0 / 8.d0**(lev)
#endif

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k) + 0.5d0) * dx(3)

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j) + 0.5d0) * dx(2)

          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + 0.5d0) * dx(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( use_mask .eq. 1) then
                if ( (mask(i,j,k).eq.1) ) cell_valid = .false.
             end if

             if (cell_valid) then

                ! vel is the magnitude of the velocity, including w0
#if (AMREX_SPACEDIM == 2)
                vel = sqrt(  u(i,j,k,1)**2 + &
                     ( u(i,j,k,2) + 0.5d0*(w0(lev,j) + w0(lev,j+1)) )**2 )
#elif (AMREX_SPACEDIM == 3)
                vel = sqrt(  u(i,j,k,1)**2 + &
                     u(i,j,k,2)**2 + &
                     ( u(i,j,k,3) + 0.5d0*(w0(lev,k) + w0(lev,k+1)) )**2 )
#endif

                ! call the EOS to get the sound speed and internal energy
                eos_state%T     = scal(i,j,k,temp_comp)
                eos_state%rho   = scal(i,j,k,rho_comp)
                eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                call eos(eos_input_rt, eos_state)

                ! max Mach number
                Mach_max = max(Mach_max,vel/eos_state%cs)

                ! max temp and enuc
                temp_max = max(temp_max,scal(i,j,k,temp_comp))
                enuc_max = max(enuc_max,rho_Hnuc(i,j,k)/scal(i,j,k,rho_comp))
                Hext_max = max(Hext_max,rho_Hext(i,j,k)/scal(i,j,k,rho_comp))

             endif
          enddo
       enddo
    enddo

  end subroutine diag


  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine diag_sphr(lev, lo, hi, &
                       scal, s_lo, s_hi, nc_s, &
                       rho_Hnuc, hn_lo, hn_hi, &
                       rho_Hext, he_lo, he_hi, &
                       u, u_lo, u_hi, &
                       w0macx, x_lo, x_hi, &
                       w0macy, y_lo, y_hi, &
                       w0macz, z_lo, z_hi, &
                       w0r, wr_lo, wr_hi, &
                       dx, &
                       normal, n_lo, n_hi, &
                       T_max, coord_Tmax, vel_Tmax, &
                       enuc_max, coord_enucmax, vel_enucmax, &
                       kin_ener, int_ener, nuc_ener, &
                       U_max, Mach_max, &
                       ncenter, T_center, vel_center, &
                       mask,     m_lo, m_hi, use_mask) &
                       bind(C, name="diag_sphr")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3), nc_s
    double precision, intent(in   ) ::     scal(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc_s)
    integer         , intent(in   ) :: hn_lo(3), hn_hi(3)
    double precision, intent(in   ) :: rho_Hnuc(hn_lo(1):hn_hi(1),hn_lo(2):hn_hi(2),hn_lo(3):hn_hi(3))
    integer         , intent(in   ) :: he_lo(3), he_hi(3)
    double precision, intent(in   ) :: rho_Hext(he_lo(1):he_hi(1),he_lo(2):he_hi(2),he_lo(3):he_hi(3))
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) ::    u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) ::   w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) ::   w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) ::   w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    integer         , intent(in   ) :: wr_lo(3), wr_hi(3)
    double precision, intent(in   ) :: w0r(wr_lo(1):wr_hi(1),wr_lo(2):wr_hi(2),wr_lo(3):wr_hi(3))
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: n_lo(3), n_hi(3)
    double precision, intent(in   ) ::   normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),3)
    double precision, intent(inout) :: T_max, coord_Tmax(3), vel_Tmax(3)
    double precision, intent(inout) :: enuc_max, coord_enucmax(3), vel_enucmax(3)
    double precision, intent(inout) :: kin_ener, int_ener, nuc_ener
    double precision, intent(inout) :: U_max, Mach_max
    integer         , intent(inout) :: ncenter
    double precision, intent(inout) :: T_center, vel_center(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    integer         , intent (in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    integer         , intent (in   ) :: use_mask

    !     Local variables
    integer          :: i, j, k
    double precision :: weight
    double precision :: x, y, z
    double precision :: vel, velr
    logical          :: cell_valid

    type (eos_t) :: eos_state


    ! weight is the factor by which the volume of a cell at the current level
    ! relates to the volume of a cell at the coarsest level of refinement.
    weight = 1.d0 / 8.d0**(lev)

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k) + 0.5d0) * dx(3)

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j) + 0.5d0) * dx(2)

          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + 0.5d0) * dx(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( use_mask .eq. 1 ) then
                if ( (mask(i,j,k).eq.1) ) cell_valid = .false.
             end if

             ! we only consider cells inside of where the sponging begins
             if ( cell_valid .and. &
                  scal(i,j,k,rho_comp) >= sponge_start_factor*sponge_center_density ) then

                ! is it one of the 8 zones surrounding the center?
                if ( abs(x - center(1)) < dx(1)  .and. &
                     abs(y - center(2)) < dx(2)  .and. &
                     abs(z - center(3)) < dx(3) ) then

                   ncenter = ncenter + 1

                   T_center = T_center + scal(i,j,k,temp_comp)

                   vel_center(1) = vel_center(1) + &
                        u(i,j,k,1) + 0.5d0*(w0macx(i,j,k)+w0macx(i+1,j,k))
                   vel_center(2) = vel_center(2) + &
                        u(i,j,k,2) + 0.5d0*(w0macy(i,j,k)+w0macy(i,j+1,k))
                   vel_center(3) = vel_center(3) + &
                        u(i,j,k,3) + 0.5d0*(w0macz(i,j,k)+w0macz(i,j,k+1))
                endif

                ! velr is the projection of the velocity (including w0) onto
                ! the radial unit vector
                velr = u(i,j,k,1)*normal(i,j,k,1) + &
                     u(i,j,k,2)*normal(i,j,k,2) + &
                     u(i,j,k,3)*normal(i,j,k,3) + w0r(i,j,k)

                ! vel is the magnitude of the velocity, including w0
                vel = sqrt( (u(i,j,k,1)+0.5d0*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                     (u(i,j,k,2)+0.5d0*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                     (u(i,j,k,3)+0.5d0*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)

                !
                ! max T, location, and velocity at that location (including w0)
                !
                if ( scal(i,j,k,temp_comp) > T_max ) then
                   T_max         = scal(i,j,k,temp_comp)
                   coord_Tmax(1) = x
                   coord_Tmax(2) = y
                   coord_Tmax(3) = z
                   vel_Tmax(1)   = u(i,j,k,1)+0.5d0*(w0macx(i,j,k)+w0macx(i+1,j,k))
                   vel_Tmax(2)   = u(i,j,k,2)+0.5d0*(w0macy(i,j,k)+w0macy(i,j+1,k))
                   vel_Tmax(3)   = u(i,j,k,3)+0.5d0*(w0macz(i,j,k)+w0macz(i,j,k+1))
                end if
                !
                ! max enuc
                !
                if ( rho_Hnuc(i,j,k)/scal(i,j,k,rho_comp) > enuc_max ) then
                   enuc_max         = rho_Hnuc(i,j,k)/scal(i,j,k,rho_comp)
                   coord_enucmax(1) = x
                   coord_enucmax(2) = y
                   coord_enucmax(3) = z
                   vel_enucmax(1)   = u(i,j,k,1)+0.5d0*(w0macx(i,j,k)+w0macx(i+1,j,k))
                   vel_enucmax(2)   = u(i,j,k,2)+0.5d0*(w0macy(i,j,k)+w0macy(i,j+1,k))
                   vel_enucmax(3)   = u(i,j,k,3)+0.5d0*(w0macz(i,j,k)+w0macz(i,j,k+1))
                end if

                ! call the EOS to get the sound speed and internal energy
                eos_state%T     = scal(i,j,k,temp_comp)
                eos_state%rho   = scal(i,j,k,rho_comp)
                eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                call eos(eos_input_rt, eos_state)

                ! kinetic, internal, and nuclear energies
                kin_ener = kin_ener + weight*scal(i,j,k,rho_comp)*vel**2
                int_ener = int_ener + weight*scal(i,j,k,rho_comp)*eos_state%e
                nuc_ener = nuc_ener + weight*rho_Hnuc(i,j,k)

                ! max vel and Mach number
                U_max = max(U_max,vel)
                Mach_max = max(Mach_max,vel/eos_state%cs)

             end if  ! end cell_valid and density check

          enddo
       enddo
    enddo

  end subroutine diag_sphr

  subroutine diag_grav_energy(grav_ener, rho0, &
                               r_cc_loc, r_edge_loc) &
       bind(C, name="diag_grav_energy")

    double precision, intent(inout) :: grav_ener
    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)

    !     Local variables
    integer            :: r
    double precision   :: m(0:nr_fine-1)
    double precision   :: term1, term2

    grav_ener = ZERO

    ! m(r) will contain mass enclosed by the center
    m(0) = FOUR3RD*M_PI*rho0(0,0)*r_cc_loc(0,0)**3

    ! dU = - G M dM / r;  dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
    grav_ener = -FOUR*M_PI*Gconst*m(0)*r_cc_loc(0,0)*rho0(0,0)*(r_edge_loc(0,1)-r_edge_loc(0,0))

    do r=1,nr_fine-1

       ! the mass is defined at the cell-centers, so to compute the
       ! mass at the current center, we need to add the contribution
       ! of the upper half of the zone below us and the lower half of
       ! the current zone.

       ! don't add any contributions from outside the star -- i.e.
       ! rho < base_cutoff_density
       if ( rho0(0,r-1) > base_cutoff_density ) then
          term1 = FOUR3RD*M_PI*rho0(0,r-1) * &
               (r_edge_loc(0,r) - r_cc_loc(0,r-1)) * &
               (r_edge_loc(0,r)**2 + &
                r_edge_loc(0,r)*r_cc_loc(0,r-1) + &
                r_cc_loc(0,r-1)**2)
       else
          term1 = ZERO
       endif

       if ( rho0(0,r) > base_cutoff_density ) then
          term2 = FOUR3RD*M_PI*rho0(0,r  )*&
               (r_cc_loc(0,r) - r_edge_loc(0,r  )) * &
               (r_cc_loc(0,r)**2 + &
                r_cc_loc(0,r)*r_edge_loc(0,r  ) + &
                r_edge_loc(0,r  )**2)
       else
          term2 = ZERO
       endif

       m(r) = m(r-1) + term1 + term2

       ! dU = - G M dM / r;
       ! dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
       grav_ener = grav_ener - &
            FOUR*M_PI*Gconst*m(r)*r_cc_loc(0,r)*rho0(0,r)*(r_edge_loc(0,r+1)-r_edge_loc(0,r))

    enddo

  end subroutine diag_grav_energy

end module diag_module
