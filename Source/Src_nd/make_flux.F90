module make_flux_module
  ! The make_flux routines take the predicted edges states of the scalars
  ! and the MAC velocities and compute the fluxes through the
  ! interfaces.
  !
  ! For the species fluxes, the construction of the fluxes depends on
  ! what form the incoming edge states take.  This depends on
  ! species_pred_type:
  !
  ! predict_rhoprime_and_X:
  ! We have rho' and X, and need a edge-centered base state to
  ! make the final fluxes
  !
  ! predict_rhoX:
  ! We use the (rho X) edge state directly to compute the fluxes.
  ! No base state input needed.
  !
  ! predict_rho_and_X:
  ! The fluxes are computed from the product of the rho and X
  ! edge states, again, no base state input needed.
  !
  !
  ! For enthalpy, there are a wide range of quantities that we predict,
  ! but they fall into 2 categories.  The enthalpy edge states either
  ! contain predictions of h or (rho h)'.  (There is limited support for
  ! h' prediction, but it is not well tested).  If we have h, then we
  ! construct a rho depending on the species states (i.e. species_pred_type).
  ! If we have (rho h)', then we use the base state to make (rho h)_0 on
  ! edges.

  use amrex_error_module
  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use network, only: nspec
  use meth_params_module, only: rho_comp, rhoh_comp, enthalpy_pred_type, species_pred_type, &
       evolve_base_state, use_exact_base_state, spec_comp

  implicit none

  private

  integer, parameter :: predict_rhoh             = 0;
  integer, parameter :: predict_rhohprime        = 1;
  integer, parameter :: predict_h                = 2;
  integer, parameter :: predict_T_then_rhohprime = 3;
  integer, parameter :: predict_T_then_h         = 4;
  integer, parameter :: predict_hprime           = 5;
  integer, parameter :: predict_Tprime_then_h    = 6;

  integer, parameter :: predict_rhoprime_and_X   = 1;
  integer, parameter :: predict_rhoX             = 2;
  integer, parameter :: predict_rho_and_X        = 3;

contains

#if (AMREX_SPACEDIM == 2)
  subroutine make_rhoX_flux_2d(lev, lo, hi, &
       sfluxx, fx_lo, fx_hi, nc_fx, &
       sfluxy, fy_lo, fy_hi, nc_fy, &
       etarhoflux, eta_lo, eta_hi, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, &
       rho0_old, rho0_edge_old, &
       rho0_new, rho0_edge_new, &
       rho0_predicted_edge, &
       w0, &
       startcomp, endcomp) bind(C,name="make_rhoX_flux_2d")
    ! Binds to C function ``make_rhoX_flux_2d``

    integer         , intent(in   ) :: lev, lo(2), hi(2)
    integer         , intent(in   ) :: fx_lo(2), fx_hi(2), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),nc_fx)
    integer         , intent(in   ) :: fy_lo(2), fy_hi(2), nc_fy
    double precision, intent(inout) :: sfluxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),nc_fy)
    integer         , intent(in   ) :: eta_lo(2), eta_hi(2)
    double precision, intent(inout) :: etarhoflux(eta_lo(1):eta_hi(1),eta_lo(2):eta_hi(2))
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
    double precision, intent(in   ) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine)
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

             if (species_pred_type == predict_rhoprime_and_X) then
                ! edge states are rho' and X.  To make the (rho X) flux,
                ! we need the edge state of rho0
                sfluxx(i,j,comp) = umac(i,j)* &
                     (rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,comp)

             else if (species_pred_type == predict_rhoX) then
                ! edge states are (rho X)
                sfluxx(i,j,comp) = umac(i,j)*sedgex(i,j,comp)

             else if (species_pred_type == predict_rho_and_X) then
                ! edge states are rho and X
                sfluxx(i,j,comp) = umac(i,j)* &
                     sedgex(i,j,rho_comp)*sedgex(i,j,comp)

             end if
          end do
       end do

       ! create y-fluxes
       do j = lo(2),hi(2)+1
          rho0_edge = HALF*(rho0_edge_old(lev,j)+rho0_edge_new(lev,j))
          do i = lo(1),hi(1)

             if (species_pred_type == predict_rhoprime_and_X) then
                ! edge states are rho' and X.  To make the (rho X) flux,
                ! we need the edge state of rho0
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(lev,j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,comp)

             else if (species_pred_type == predict_rhoX) then
                ! edge states are (rho X)
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(lev,j))*sedgey(i,j,comp)

             else if (species_pred_type == predict_rho_and_X) then
                ! edge state are rho and X
                sfluxy(i,j,comp) = &
                     (vmac(i,j)+w0(lev,j))*sedgey(i,j,rho_comp)*sedgey(i,j,comp)

             endif

             if (evolve_base_state .and. .not.use_exact_base_state) then
                if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) then
                   etarhoflux(i,j) = etarhoflux(i,j) + sfluxy(i,j,comp)
                end if

                if ( comp.eq.spec_comp+nspec-1) then
                   etarhoflux(i,j) = etarhoflux(i,j) - w0(lev,j)*rho0_predicted_edge(lev,j)
                end if
             endif  ! evolve_base_state
          end do
       end do
    end do

    ! compute the density fluxes by summing the species fluxes

    ! loop for x-fluxes
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          sfluxx(i,j,1) = sum(sfluxx(i,j,startcomp:endcomp))
       end do
    end do

    ! loop for y-fluxes
    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          sfluxy(i,j,1) = sum(sfluxy(i,j,startcomp:endcomp))
       end do
    end do

  end subroutine make_rhoX_flux_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine make_rhoX_flux_3d(lev, lo, hi, &
       sfluxx, fx_lo, fx_hi, nc_fx, &
       sfluxy, fy_lo, fy_hi, nc_fy, &
       sfluxz, fz_lo, fz_hi, nc_fz, &
       etarhoflux, eta_lo, eta_hi, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       sedgez, z_lo, z_hi, nc_z, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, &
       wmac,   w_lo, w_hi, &
       rho0_old, rho0_edge_old, &
       rho0_new, rho0_edge_new, &
       rho0_predicted_edge, &
       w0, &
       startcomp, endcomp) bind(C,name="make_rhoX_flux_3d")
    ! Binds to C function ``make_rhoX_flux_3d``

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),nc_fx)
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3), nc_fy
    double precision, intent(inout) :: sfluxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),nc_fy)
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3), nc_fz
    double precision, intent(inout) :: sfluxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),nc_fz)
    integer         , intent(in   ) :: eta_lo(3), eta_hi(3)
    double precision, intent(inout) :: etarhoflux(eta_lo(1):eta_hi(1),eta_lo(2):eta_hi(2),eta_lo(3):eta_hi(3))
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
    double precision, intent(in   ) :: rho0_predicted_edge(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i,j,k
    double precision :: rho0_edge

    do comp = startcomp, endcomp
       ! create x-fluxes and y-fluxes

       do k=lo(3),hi(3)
          rho0_edge = HALF*(rho0_old(lev,k)+rho0_new(lev,k))

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*sedgex(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxx(i,j,k,comp) = &
                        umac(i,j,k)*sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp)

                endif

             end do
          end do

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxy(i,j,k,comp) = &
                        vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp)

                endif

             end do
          end do
       end do

       ! create z-fluxes
       do k=lo(3),hi(3)+1
          rho0_edge = HALF*(rho0_edge_old(lev,k)+rho0_edge_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(lev,k))* &
                        (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(lev,k))*sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxz(i,j,k,comp) = (wmac(i,j,k)+w0(lev,k))* &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp)

                endif

                if (evolve_base_state .and. .not.use_exact_base_state) then
                   if (comp .ge. spec_comp .and. comp .le. spec_comp+nspec-1) then
                      etarhoflux(i,j,k) = etarhoflux(i,j,k) + sfluxz(i,j,k,comp)
                   end if

                   if ( comp.eq.spec_comp+nspec-1) then
                      etarhoflux(i,j,k) = etarhoflux(i,j,k) - w0(lev,k)*rho0_predicted_edge(lev,k)
                   end if
                endif ! evolve_base_state
             end do
          end do
       end do

    end do

    ! compute the density fluxes by summing the species fluxes

    ! loop for x-fluxes
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             sfluxx(i,j,k,1) = sum(sfluxx(i,j,k,startcomp:endcomp))
          end do
       end do
    end do

    ! loop for y-fluxes
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             sfluxy(i,j,k,1) = sum(sfluxy(i,j,k,startcomp:endcomp))
          end do
       end do
    end do

    ! loop for z-fluxes
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             sfluxz(i,j,k,1) = sum(sfluxz(i,j,k,startcomp:endcomp))
          end do
       end do
    end do

  end subroutine make_rhoX_flux_3d

  !----------------------------------------------------------------------------
  ! make_rhoX_flux_3d_sphr
  !----------------------------------------------------------------------------
  subroutine make_rhoX_flux_3d_sphr(lo, hi, &
       sfluxx, fx_lo, fx_hi, nc_fx, &
       sfluxy, fy_lo, fy_hi, nc_fy, &
       sfluxz, fz_lo, fz_hi, nc_fz, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       sedgez, z_lo, z_hi, nc_z, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, &
       wmac,   w_lo, w_hi, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       rho0_edgex, rx_lo, rx_hi, &
       rho0_edgey, ry_lo, ry_hi, &
       rho0_edgez, rz_lo, rz_hi, &
       startcomp, endcomp) bind(C,name="make_rhoX_flux_3d_sphr")
    ! Binds to C function ``make_rhoX_flux_3d_sphr``

    integer         , intent(in   ) :: lo(3), hi(3)
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
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    double precision, intent(inout) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    double precision, intent(inout) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    double precision, intent(inout) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    integer         , intent(in   ) :: rx_lo(3), rx_hi(3)
    double precision, intent(in   ) :: rho0_edgex(rx_lo(1):rx_hi(1),rx_lo(2):rx_hi(2),rx_lo(3):rx_hi(3))
    integer         , intent(in   ) :: ry_lo(3), ry_hi(3)
    double precision, intent(in   ) :: rho0_edgey(ry_lo(1):ry_hi(1),ry_lo(2):ry_hi(2),ry_lo(3):ry_hi(3))
    integer         , intent(in   ) :: rz_lo(3), rz_hi(3)
    double precision, intent(in   ) :: rho0_edgez(rz_lo(1):rz_hi(1),rz_lo(2):rz_hi(2),rz_lo(3):rz_hi(3))
    integer         , intent(in   ) :: startcomp, endcomp

    ! local
    integer          :: comp
    integer          :: i,j,k
    ! double precision :: rho0_edge

    do comp = startcomp, endcomp

       ! loop for x-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   ! rho0_edge = HALF*(rho0macx_old(i,j,k)+rho0macx_new(i,j,k))
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        (rho0_edgex(i,j,k) + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * sedgex(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxx(i,j,k,comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp)

                endif

             end do
          end do
       end do

       ! loop for y-fluxes
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   ! rho0_edge = HALF*(rho0macy_old(i,j,k)+rho0macy_new(i,j,k))
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        (rho0_edgey(i,j,k) + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * sedgey(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxy(i,j,k,comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp)

                endif

             end do
          end do
       end do

       ! loop for z-fluxes
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   ! rho0_edge = HALF*(rho0macz_old(i,j,k)+rho0macz_new(i,j,k))
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        (rho0_edgez(i,j,k) + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * sedgez(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sfluxz(i,j,k,comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp)

                endif

             end do
          end do
       end do

    end do ! end loop over components

    ! compute the density fluxes by summing the species fluxes

    ! loop for x-fluxes
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             sfluxx(i,j,k,1) = sum(sfluxx(i,j,k,startcomp:endcomp))
          end do
       end do
    end do

    ! loop for y-fluxes
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             sfluxy(i,j,k,1) = sum(sfluxy(i,j,k,startcomp:endcomp))
          end do
       end do
    end do

    ! loop for z-fluxes
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             sfluxz(i,j,k,1) = sum(sfluxz(i,j,k,startcomp:endcomp))
          end do
       end do
    end do

  end subroutine make_rhoX_flux_3d_sphr
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
    ! Binds to C function ``make_rhoh_flux_2d``

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

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'
          do j=lo(2),hi(2)
             rho0_edge = HALF*(rho0_old(lev,j)+rho0_new(lev,j))
             do i=lo(1),hi(1)+1
                sfluxx(i,j,rhoh_comp) = &
                     umac(i,j)*(rho0_edge+sedgex(i,j,rho_comp))*sedgex(i,j,rhoh_comp)
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                sfluxx(i,j,rhoh_comp) = &
                     umac(i,j)*sedgex(i,j,rho_comp)*sedgex(i,j,rhoh_comp)
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("make_rhoh_flux_2d : predict_hprime not coded yet")

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
       call amrex_error("make_rhoh_flux_2d : enthalpy_pred_type not recognized.")
    end if

    ! create y-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'
          do j=lo(2),hi(2)+1
             rho0_edge = HALF*(rho0_edge_old(lev,j)+rho0_edge_new(lev,j))
             do i=lo(1),hi(1)
                sfluxy(i,j,rhoh_comp) = &
                     (vmac(i,j)+w0(lev,j))*(rho0_edge+sedgey(i,j,rho_comp))*sedgey(i,j,rhoh_comp)
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                sfluxy(i,j,rhoh_comp) = &
                     (vmac(i,j)+w0(lev,j))*sedgey(i,j,rho_comp)*sedgey(i,j,rhoh_comp)
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("make_rhoh_flux_2d : predict_hprime not coded yet")

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
       call amrex_error("make_rhoh_flux_2d : enthalpy_pred_type not recognized.")
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
    ! Binds to C function ``make_rhoh_flux_3d``

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

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

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

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   sfluxx(i,j,k,rhoh_comp) = &
                        umac(i,j,k)*sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp)
                end do
             end do

             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   sfluxy(i,j,k,rhoh_comp) = &
                        vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp)
                end do
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("make_rhoh_flux_3d : predict_hprime not coded yet")

    else if (have_rhoh) then

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

    else if (enthalpy_pred_type.eq.predict_rhohprime) then
       ! enthalpy edge state is (rho h)'

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

    else
       call amrex_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
    end if

    ! create z-fluxes
    if (have_h) then
       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          do k=lo(3),hi(3)+1
             rho0_edge = HALF*(rho0_edge_old(lev,k)+rho0_edge_new(lev,k))
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(lev,k))* &
                        (rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(lev,k))* &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,rhoh_comp)
                end do
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("make_rhoh_flux_3d : predict_hprime not coded yet")

    else if (have_rhoh) then

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0(lev,k))*sedgez(i,j,k,rhoh_comp)
             end do
          end do
       end do

    else if (enthalpy_pred_type.eq.predict_rhohprime) then
       ! enthalpy edge state is (rho h)'

       do k=lo(3),hi(3)+1
          rhoh0_edge = HALF*(rhoh0_edge_old(lev,k)+rhoh0_edge_new(lev,k))
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sfluxz(i,j,k,rhoh_comp) = &
                     (wmac(i,j,k)+w0(lev,k))*(sedgez(i,j,k,rhoh_comp)+rhoh0_edge)
             end do
          end do
       end do

    else
       call amrex_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
    end if


  end subroutine make_rhoh_flux_3d

  !----------------------------------------------------------------------------
  ! mk_rhoh_flux_3d_sphr
  !----------------------------------------------------------------------------
  subroutine mk_rhoh_flux_3d_sphr(lo, hi, &
       sfluxx, fx_lo, fx_hi, nc_fx, &
       sfluxy, fy_lo, fy_hi, nc_fy, &
       sfluxz, fz_lo, fz_hi, nc_fz, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       sedgez, z_lo, z_hi, nc_z, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, &
       wmac,   w_lo, w_hi, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       rho0_edgex, rx_lo, rx_hi, &
       rho0_edgey, ry_lo, ry_hi, &
       rho0_edgez, rz_lo, rz_hi, &
       h0_edgex, hx_lo, hx_hi, &
       h0_edgey, hy_lo, hy_hi, &
       h0_edgez, hz_lo, hz_hi) bind(C,name="make_rhoh_flux_3d_sphr")
    ! Binds to C function ``make_rhoh_flux_3d_sphr``

    integer         , intent(in   ) :: lo(3), hi(3)
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
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    integer         , intent(in   ) :: rx_lo(3), rx_hi(3)
    double precision, intent(in   ) :: rho0_edgex(rx_lo(1):rx_hi(1),rx_lo(2):rx_hi(2),rx_lo(3):rx_hi(3))
    integer         , intent(in   ) :: ry_lo(3), ry_hi(3)
    double precision, intent(in   ) :: rho0_edgey(ry_lo(1):ry_hi(1),ry_lo(2):ry_hi(2),ry_lo(3):ry_hi(3))
    integer         , intent(in   ) :: rz_lo(3), rz_hi(3)
    double precision, intent(in   ) :: rho0_edgez(rz_lo(1):rz_hi(1),rz_lo(2):rz_hi(2),rz_lo(3):rz_hi(3))
    integer         , intent(in   ) :: hx_lo(3), hx_hi(3)
    double precision, intent(in   ) :: h0_edgex(hx_lo(1):hx_hi(1),hx_lo(2):hx_hi(2),hx_lo(3):hx_hi(3))
    integer         , intent(in   ) :: hy_lo(3), hy_hi(3)
    double precision, intent(in   ) :: h0_edgey(hy_lo(1):hy_hi(1),hy_lo(2):hy_hi(2),hy_lo(3):hy_hi(3))
    integer         , intent(in   ) :: hz_lo(3), hz_hi(3)
    double precision, intent(in   ) :: h0_edgez(hz_lo(1):hz_hi(1),hz_lo(2):hz_hi(2),hz_lo(3):hz_hi(3))

    ! local
    integer          :: i,j,k
    logical          :: have_h, have_hprime, have_rhoh

    have_h = enthalpy_pred_type.eq.predict_h .or. &
         enthalpy_pred_type.eq.predict_T_then_h .or. &
         enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    ! create x-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1

                   sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        (rho0_edgex(i,j,k) + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp)

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1

                   sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k) + w0macx(i,j,k)) * &
                        sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp)

                end do
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          ! (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
          ! computed from (rho h)_0 / rho_0
          ! sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1

                   sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k)+w0macx(i,j,k)) * &
                        (sedgex(i,j,k,rho_comp)+rho0_edgex(i,j,k)) * (sedgex(i,j,k,rhoh_comp)+h0_edgex(i,j,k))

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")

       endif

    else if (have_rhoh) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k) + w0macx(i,j,k))*sedgex(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                ! where h_0 is computed from (rho h)_0 / rho_0

                sfluxx(i,j,k,rhoh_comp) = &
                     (umac(i,j,k)+w0macx(i,j,k))*(rho0_edgex(i,j,k)*h0_edgex(i,j,k)+sedgex(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

    ! create y-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)

                   sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        (rho0_edgey(i,j,k) + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp)

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)

                   sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k) + w0macy(i,j,k)) * &
                        sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp)

                end do
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          ! (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
          ! computed from (rho h)_0 / rho_0
          ! sfluxy = (vmac(i,j,k)+w0macy(i,j,k)) * (rho h)_edge

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)

                   sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k)+w0macy(i,j,k)) * &
                        (sedgey(i,j,k,rho_comp)+rho0_edgey(i,j,k)) * (sedgey(i,j,k,rhoh_comp)+h0_edgey(i,j,k))

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")

       endif


    else if (have_rhoh) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k) + w0macy(i,j,k))*sedgey(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                ! where h_0 is computed from (rho h)_0 / rho_0

                sfluxy(i,j,k,rhoh_comp) = &
                     (vmac(i,j,k)+w0macy(i,j,k))*(rho0_edgey(i,j,k)*h0_edgey(i,j,k)+sedgey(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

    ! create z-fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        (rho0_edgez(i,j,k) + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp)

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k) + w0macz(i,j,k)) * &
                        sedgez(i,j,k,rho_comp)*sedgez(i,j,k,rhoh_comp)

                end do
             end do
          end do

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          ! (rho h)_edge = (h' + h_0) * (rho' + rho_0)
          ! where h0 is computed from (rho h)_0 / rho_0
          ! sfluxz = (wmac(i,j,k)+w0macz(i,j,k)) * (rho h)_edge

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k)+w0macz(i,j,k)) * &
                        (sedgez(i,j,k,rho_comp)+rho0_edgez(i,j,k)) * (sedgez(i,j,k,rhoh_comp)+h0_edgez(i,j,k))

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")

       endif

    else if (have_rhoh) then

       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k) + w0macz(i,j,k))*sedgez(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                ! where h_0 is computed from (rho h)_0 / rho_0

                sfluxz(i,j,k,rhoh_comp) = &
                     (wmac(i,j,k)+w0macz(i,j,k))*(rho0_edgez(i,j,k)*h0_edgez(i,j,k)+sedgez(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

  end subroutine mk_rhoh_flux_3d_sphr

  subroutine mk_rhoh_flux_3d_sphr_irreg(lo, hi, &
       sfluxx, fx_lo, fx_hi, nc_fx, &
       sfluxy, fy_lo, fy_hi, nc_fy, &
       sfluxz, fz_lo, fz_hi, nc_fz, &
       sedgex, x_lo, x_hi, nc_x, &
       sedgey, y_lo, y_hi, nc_y, &
       sedgez, z_lo, z_hi, nc_z, &
       umac,   u_lo, u_hi, &
       vmac,   v_lo, v_hi, &
       wmac,   w_lo, w_hi, &
       w0macx, wx_lo, wx_hi, &
       w0macy, wy_lo, wy_hi, &
       w0macz, wz_lo, wz_hi, &
       rhoh0_edgex, hx_lo, hx_hi, &
       rhoh0_edgey, hy_lo, hy_hi, &
       rhoh0_edgez, hz_lo, hz_hi) bind(C,name="make_rhoh_flux_3d_sphr_irreg")
    ! Binds to C function ``make_rhoh_flux_3d_sphr_irreg``

    integer         , intent(in   ) :: lo(3), hi(3)
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
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    double precision, intent(in   ) :: w0macx(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    integer         , intent(in   ) :: wy_lo(3), wy_hi(3)
    double precision, intent(in   ) :: w0macy(wy_lo(1):wy_hi(1),wy_lo(2):wy_hi(2),wy_lo(3):wy_hi(3))
    integer         , intent(in   ) :: wz_lo(3), wz_hi(3)
    double precision, intent(in   ) :: w0macz(wz_lo(1):wz_hi(1),wz_lo(2):wz_hi(2),wz_lo(3):wz_hi(3))
    integer         , intent(in   ) :: hx_lo(3), hx_hi(3)
    double precision, intent(in   ) :: rhoh0_edgex(hx_lo(1):hx_hi(1),hx_lo(2):hx_hi(2),hx_lo(3):hx_hi(3))
    integer         , intent(in   ) :: hy_lo(3), hy_hi(3)
    double precision, intent(in   ) :: rhoh0_edgey(hy_lo(1):hy_hi(1),hy_lo(2):hy_hi(2),hy_lo(3):hy_hi(3))
    integer         , intent(in   ) :: hz_lo(3), hz_hi(3)
    double precision, intent(in   ) :: rhoh0_edgez(hz_lo(1):hz_hi(1),hz_lo(2):hz_hi(2),hz_lo(3):hz_hi(3))

    ! local
    integer          :: i,j,k
    logical          :: have_h, have_hprime, have_rhoh

    have_h = enthalpy_pred_type.eq.predict_h .or. &
         enthalpy_pred_type.eq.predict_T_then_h .or. &
         enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    ! create x-fluxes
    if (have_h) then

       ! enthalpy edge state is h
       call amrex_error("have_h not supported on irregular-spaced base state")

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("have_hprime not supported on irregular-spaced base state")


       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")

       endif

    else if (have_rhoh) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                ! sfluxx(i,j,k,rhoh_comp) = (umac(i,j,k) + w0macx(i,j,k))*sedgex(i,j,k,rhoh_comp)
                sfluxx(i,j,k,rhoh_comp) = umac(i,j,k)*sedgex(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rhoh_0)

                ! sfluxx(i,j,k,rhoh_comp) = &
                !      (umac(i,j,k)+w0macx(i,j,k))*(rhoh0_edgex(i,j,k)+sedgex(i,j,k,rhoh_comp))
                sfluxx(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edgex(i,j,k)+sedgex(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

    ! create y-fluxes
    if (have_h) then

       ! enthalpy edge state is h
       call amrex_error("have_h not supported on irregular-spaced base state")

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("have_hprime not supported on irregular-spaced base state")

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")

       endif


    else if (have_rhoh) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                ! sfluxy(i,j,k,rhoh_comp) = (vmac(i,j,k) + w0macy(i,j,k))*sedgey(i,j,k,rhoh_comp)
                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*sedgey(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rhoh_0)

                ! sfluxy(i,j,k,rhoh_comp) = &
                !      (vmac(i,j,k)+w0macy(i,j,k))*(rhoh0_edgey(i,j,k)+sedgey(i,j,k,rhoh_comp))
                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*(rhoh0_edgey(i,j,k)+sedgey(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

    ! create z-fluxes
    if (have_h) then

       ! enthalpy edge state is h
       call amrex_error("have_h not supported on irregular-spaced base state")

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
       call amrex_error("have_hprime not supported on irregular-spaced base state")

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")

       endif

    else if (have_rhoh) then

       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! sfluxz(i,j,k,rhoh_comp) = (wmac(i,j,k) + w0macz(i,j,k))*sedgez(i,j,k,rhoh_comp)
                sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k)*sedgez(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rhoh_0)

                ! sfluxz(i,j,k,rhoh_comp) = &
                !      (wmac(i,j,k)+w0macz(i,j,k))*(rhoh0_edgez(i,j,k)+sedgez(i,j,k,rhoh_comp))
                sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k)*(rhoh0_edgez(i,j,k)+sedgez(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

  end subroutine mk_rhoh_flux_3d_sphr_irreg
#endif

end module make_flux_module
