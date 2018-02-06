! The make_flux routines take the predicted edges states of the scalars
! and the MAC velocities and compute the fluxes through the
! interfaces.

! For the species fluxes, the construction of the fluxes depends on
! what form the incoming edge states take.  This depends on
! species_pred_type:
!
! predict_rhoprime_and_X: 
!    We have rho' and X, and need a edge-centered base state to
!    make the final fluxes
!
! predict_rhoX:
!    We use the (rho X) edge state directly to compute the fluxes.
!    No base state input needed.
!
! predict_rho_and_X:
!   The fluxes are computed from the product of the rho and X 
!   edge states, again, no base state input needed.
!
!
! For enthalpy, there are a wide range of quantities that we predict,
! but they fall into 2 categories.  The enthalpy edge states either
! contain predictions of h or (rho h)'.  (There is limited support for
! h' prediction, but it is not well tested).  If we have h, then we
! construct a rho depending on the species states (i.e. species_pred_type).
! If we have (rho h)', then we use the base state to make (rho h)_0 on
! edges.

module make_flux_module

  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: enthalpy_pred_type

  implicit none

  private
  
    integer, parameter :: rho_comp  = 1
    integer, parameter :: rhoh_comp = 2

    integer, parameter :: predict_rhoh             = 0;
    integer, parameter :: predict_rhohprime        = 1;
    integer, parameter :: predict_h                = 2;
    integer, parameter :: predict_T_then_rhohprime = 3;
    integer, parameter :: predict_T_then_h         = 4;
    integer, parameter :: predict_hprime           = 5;
    integer, parameter :: predict_Tprime_then_h    = 6;

contains
  
#if (AMREX_SPACEDIM == 1)
  subroutine make_rhoX_flux_1d(lev, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 umac,   u_lo, u_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 w0, & 
                                 startcomp, endcomp) bind(C,name="make_rhoX_flux_1d")

    integer         , intent(in   ) :: lev, lo(1), hi(1)
    integer         , intent(in   ) :: fx_lo(1), fx_hi(1), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),nc_fx)
    integer         , intent(in   ) :: x_lo(1), x_hi(1), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),nc_x)
    integer         , intent(in   ) :: u_lo(1), u_hi(1)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_levl,0:nr_fine)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i
    double precision :: rho0_edge

    do comp = startcomp, endcomp

       ! create x-fluxes
       do i = lo(1),hi(1)+1

             ! edge states are rho' and X.  To make the (rho X) flux,
             ! we need the edge state of rho0
             rho0_edge = HALF*(rho0_edge_old(lev,i)+rho0_edge_new(lev,i))
             sfluxx(i,comp) = &
                  (umac(i)+w0(lev,i))*(rho0_edge+sedgex(i,rho_comp))*sedgex(i,comp)
       end do
    end do ! end comp loop

  end subroutine make_rhoX_flux_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine make_rhoX_flux_2d(lev, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sfluxy, fy_lo, fy_hi, nc_fy, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 sedgey, y_lo, y_hi, nc_y, &
                                 umac,   u_lo, u_hi, &
                                 vmac,   v_lo, v_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 w0, & 
                                 startcomp, endcomp) bind(C,name="make_rhoX_flux_2d")

    integer         , intent(in   ) :: lev, lo(2), hi(2)
    integer         , intent(in   ) :: fx_lo(2), fx_hi(2), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),nc_fx)
    integer         , intent(in   ) :: fy_lo(2), fy_hi(2), nc_fy
    double precision, intent(inout) :: sfluxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),nc_fy)
    integer         , intent(in   ) :: x_lo(2), x_hi(2), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),nc_x)
    integer         , intent(in   ) :: y_lo(2), y_hi(2), nc_y
    double precision, intent(inout) :: sedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),nc_y)
    integer         , intent(in   ) :: u_lo(2), u_hi(2)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2))
    integer         , intent(in   ) :: v_lo(2), v_hi(2)
    double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i,j
    double precision :: rho0_edge

    do comp = startcomp, endcomp 

       ! create x-fluxes
       do j=lo(2),hi(2)
          rho0_edge = HALF*(rho0_old(lev,j)+rho0_new(lev,j))
          do i=lo(1),hi(1)+1

                ! edge states are rho' and X.  To make the (rho X) flux,
                ! we need the edge state of rho0
                sfluxx(i,j,comp) = umac(i,j)* &
                     (rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,comp)
          end do
       end do
             
       ! create y-fluxes
       do j = lo(2),hi(2)+1
          rho0_edge = HALF*(rho0_edge_old(lev,j)+rho0_edge_new(lev,j))
          do i = lo(1),hi(1)

                ! edge states are rho' and X.  To make the (rho X) flux,
                ! we need the edge state of rho0
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(lev,j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,comp)
          end do
       end do
    end do

  end subroutine make_rhoX_flux_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine make_rhoX_flux_3d(lev, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sfluxy, fy_lo, fy_hi, nc_fy, &
                                 sfluxz, fz_lo, fz_hi, nc_fz, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 sedgey, y_lo, y_hi, nc_y, &
                                 sedgez, z_lo, z_hi, nc_z, &
                                 umac,   u_lo, u_hi, &
                                 vmac,   v_lo, v_hi, &
                                 wmac,   w_lo, w_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 w0, & 
                                 startcomp, endcomp) bind(C,name="make_rhoX_flux_3d")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),nc_fx)
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3), nc_fy
    double precision, intent(inout) :: sfluxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),nc_fy)
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3), nc_fz
    double precision, intent(inout) :: sfluxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),nc_fz)
    integer         , intent(in   ) :: x_lo(3), x_hi(3), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: y_lo(3), y_hi(3), nc_y
    double precision, intent(inout) :: sedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nc_y)
    integer         , intent(in   ) :: z_lo(3), z_hi(3), nc_z
    double precision, intent(inout) :: sedgez(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),nc_z)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i,j,k
    double precision :: rho0_edge
    
    do comp = startcomp, endcomp 
       ! create x-fluxes and y-fluxes

       !$OMP PARALLEL PRIVATE(i,j,k,rho0_edge)
       !$OMP DO 
       do k=lo(3),hi(3)
          rho0_edge = HALF*(rho0_old(lev,k)+rho0_new(lev,k))

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)
             end do
          end do
          
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       ! create z-fluxes
       !$OMP DO
       do k=lo(3),hi(3)+1
          rho0_edge = HALF*(rho0_edge_old(lev,k)+rho0_edge_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(lev,k))* &
                     (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    end do

  end subroutine make_rhoX_flux_3d
#endif


#if (AMREX_SPACEDIM == 1)
  subroutine make_rhoh_flux_1d(lev, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 umac,   u_lo, u_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 rhoh0_old, rhoh0_edge_old, &
                                 rhoh0_new, rhoh0_edge_new, &
                                 w0) bind(C,name="make_rhoh_flux_1d")

    integer         , intent(in   ) :: lev, lo(1), hi(1)
    integer         , intent(in   ) :: fx_lo(1), fx_hi(1), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),nc_fx)
    integer         , intent(in   ) :: x_lo(1), x_hi(1), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),nc_x)
    integer         , intent(in   ) :: u_lo(1), u_hi(1)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_levl,0:nr_fine)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rhoh0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rhoh0_edge_old(0:max_radial_levl,0:nr_fine)
    double precision, intent(in   ) :: rhoh0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rhoh0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)

    ! Local variables
    integer          :: i
    double precision :: rho0_edge, rhoh0_edge
    logical          :: have_h, have_hprime, have_rhoh

    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    ! create x-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       ! density edge state is rho'
       do i=lo(1),hi(1)+1
          rho0_edge = HALF*(rho0_edge_old(lev,i)+rho0_edge_new(lev,i))
          sfluxx(i,rhoh_comp) = &
               (umac(i)+w0(lev,i))*(rho0_edge+sedgex(i,rho_comp))*sedgex(i,rhoh_comp)
       end do

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("make_rhoh_flux_1d : predict_hprime not coded yet")

    else if (have_rhoh) then

       do i=lo(1),hi(1)+1
          sfluxx(i,rhoh_comp) = (umac(i)+w0(lev,i))*sedgex(i,rhoh_comp)
       end do

    else if (enthalpy_pred_type.eq.predict_rhohprime) then
       ! enthalpy edge state is (rho h)'
       do i=lo(1),hi(1)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(lev,i)+rhoh0_edge_new(lev,i))
          sfluxx(i,rhoh_comp) = (umac(i)+w0(lev,i))*(sedgex(i,rhoh_comp)+rhoh0_edge)
       end do

    else
       call bl_error("make_rhoh_flux_1d : enthalpy_pred_type not recognized.")
    end if

  end subroutine make_rhoh_flux_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine make_rhoh_flux_2d(lev, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sfluxy, fy_lo, fy_hi, nc_fy, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 sedgey, y_lo, y_hi, nc_y, &
                                 umac,   u_lo, u_hi, &
                                 vmac,   v_lo, v_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 rhoh0_old, rhoh0_edge_old, &
                                 rhoh0_new, rhoh0_edge_new, &
                                 w0) bind(C,name="make_rhoh_flux_2d")

    integer         , intent(in   ) :: lev, lo(2), hi(2)
    integer         , intent(in   ) :: fx_lo(2), fx_hi(2), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),nc_fx)
    integer         , intent(in   ) :: fy_lo(2), fy_hi(2), nc_fy
    double precision, intent(inout) :: sfluxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),nc_fy)
    integer         , intent(in   ) :: x_lo(2), x_hi(2), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),nc_x)
    integer         , intent(in   ) :: y_lo(2), y_hi(2), nc_y
    double precision, intent(inout) :: sedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),nc_y)
    integer         , intent(in   ) :: u_lo(2), u_hi(2)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2))
    integer         , intent(in   ) :: v_lo(2), v_hi(2)
    double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rhoh0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rhoh0_edge_old(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rhoh0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rhoh0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)

    ! Local variables
    integer          :: i,j
    double precision :: rho0_edge, rhoh0_edge
    logical          :: have_h, have_hprime, have_rhoh

    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh
    
    ! create x-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       ! density edge state is rho'          
       do j=lo(2),hi(2)
          rho0_edge = HALF*(rho0_old(lev,j)+rho0_new(lev,j))
          do i=lo(1),hi(1)+1
             sfluxx(i,j,rhoh_comp) = &
                  umac(i,j)*(rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,rhoh_comp)
          end do
       end do

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("make_rhoh_flux_2d : predict_hprime not coded yet")

    else if (have_rhoh) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             sfluxx(i,j,rhoh_comp) = umac(i,j)*sedgex(i,j,rhoh_comp)
          end do
       end do

    else if (enthalpy_pred_type.eq.predict_rhohprime) then 
       ! enthalpy edge state is (rho h)'
       do j=lo(2),hi(2)
          rhoh0_edge = HALF*(rhoh0_old(lev,j)+rhoh0_new(lev,j))
          do i=lo(1),hi(1)+1
             sfluxx(i,j,rhoh_comp) = umac(i,j)*(rhoh0_edge+sedgex(i,j,rhoh_comp))
          end do
       end do

    else
       call bl_error("make_rhoh_flux_2d : enthalpy_pred_type not recognized.")
    end if
    
    ! create y-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       ! density edge state is rho'
       do j=lo(2),hi(2)+1
          rho0_edge = HALF*(rho0_edge_old(lev,j)+rho0_edge_new(lev,j))
          do i=lo(1),hi(1)
             sfluxy(i,j,rhoh_comp) = &
                  (vmac(i,j)+w0(lev,j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,rhoh_comp)
          end do
       end do

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("make_rhoh_flux_2d : predict_hprime not coded yet")

    else if (have_rhoh) then

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             sfluxy(i,j,rhoh_comp) = (vmac(i,j)+w0(lev,j))*sedgey(i,j,rhoh_comp)
          end do
       end do

    else if (enthalpy_pred_type.eq.predict_rhohprime) then
       ! enthalpy edge state is (rho h)'
       do j=lo(2),hi(2)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(lev,j)+rhoh0_edge_new(lev,j))
          do i=lo(1),hi(1)
             sfluxy(i,j,rhoh_comp) = (vmac(i,j)+w0(lev,j))*(sedgey(i,j,rhoh_comp)+rhoh0_edge)
          end do
       end do

    else
       call bl_error("make_rhoh_flux_2d : enthalpy_pred_type not recognized.")
    end if

  end subroutine make_rhoh_flux_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine make_rhoh_flux_3d(lev, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sfluxy, fy_lo, fy_hi, nc_fy, &
                                 sfluxz, fz_lo, fz_hi, nc_fz, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 sedgey, y_lo, y_hi, nc_y, &
                                 sedgez, z_lo, z_hi, nc_z, &
                                 umac,   u_lo, u_hi, &
                                 vmac,   v_lo, v_hi, &
                                 wmac,   w_lo, w_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 rhoh0_old, rhoh0_edge_old, &
                                 rhoh0_new, rhoh0_edge_new, &
                                 w0) bind(C,name="make_rhoh_flux_3d")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),nc_fx)
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3), nc_fy
    double precision, intent(inout) :: sfluxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),nc_fy)
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3), nc_fz
    double precision, intent(inout) :: sfluxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),nc_fz)
    integer         , intent(in   ) :: x_lo(3), x_hi(3), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: y_lo(3), y_hi(3), nc_y
    double precision, intent(inout) :: sedgey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nc_y)
    integer         , intent(in   ) :: z_lo(3), z_hi(3), nc_z
    double precision, intent(inout) :: sedgez(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),nc_z)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac  (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac  (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rhoh0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rhoh0_edge_old(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: rhoh0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rhoh0_edge_new(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)

    ! Local variables
    integer          :: i,j,k
    double precision :: rho0_edge, rhoh0_edge
    logical         :: have_h, have_hprime, have_rhoh
    
    have_h = enthalpy_pred_type.eq.predict_h .or. &
             enthalpy_pred_type.eq.predict_T_then_h .or. &
             enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    ! create x-fluxes and y-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       ! density edge state is rho'

       !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
       do k=lo(3),hi(3)
          rho0_edge = HALF*(rho0_old(lev,k)+rho0_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                sfluxx(i,j,k,rhoh_comp) = &
                     umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp)
             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                sfluxy(i,j,k,rhoh_comp) = &
                     vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if (have_hprime) then
       
       ! enthalpy edge state is h'
       call bl_error("make_rhoh_flux_3d : predict_hprime not coded yet")

    else if (have_rhoh) then
       
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                sfluxx(i,j,k,rhoh_comp) = umac(i,j,k)*sedgex(i,j,k,rhoh_comp)
             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*sedgey(i,j,k,rhoh_comp)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if (enthalpy_pred_type.eq.predict_rhohprime) then
       ! enthalpy edge state is (rho h)'

       !$OMP PARALLEL DO PRIVATE(i,j,k,rhoh0_edge)
       do k=lo(3),hi(3)
          rhoh0_edge = HALF*(rhoh0_old(lev,k)+rhoh0_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                sfluxx(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edge+sedgex(i,j,k,rhoh_comp))
             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*(rhoh0_edge+sedgey(i,j,k,rhoh_comp))
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else
       call bl_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
    end if

    ! create z-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       ! density edge state is rho'
       
       !$OMP PARALLEL DO PRIVATE(i,j,k,rho0_edge)
       do k=lo(3),hi(3)+1
          rho0_edge = HALF*(rho0_edge_old(lev,k)+rho0_edge_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(lev,k))* &
                     (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call bl_error("make_rhoh_flux_3d : predict_hprime not coded yet")

    else if (have_rhoh) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(lev,k))*sedgez(i,j,k,rhoh_comp)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if (enthalpy_pred_type.eq.predict_rhohprime) then
       ! enthalpy edge state is (rho h)'

       !$OMP PARALLEL DO PRIVATE(i,j,k,rhoh0_edge)
       do k=lo(3),hi(3)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(lev,k)+rhoh0_edge_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sfluxz(i,j,k,rhoh_comp) = &
                     (wmac(i,j,k)+w0(lev,k))*(sedgez(i,j,k,rhoh_comp)+rhoh0_edge)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       
    else
       call bl_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
    end if
    

  end subroutine make_rhoh_flux_3d
#endif

end module make_flux_module
