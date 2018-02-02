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
  use meth_params_module, only: rel_eps
  use base_state_geometry_module, only: nr_fine, max_radial_level

  implicit none

  private
  
contains
  
#if (AMREX_SPACEDIM == 1)
  subroutine make_rhoX_flux_1d(domlo, domhi, lo, hi, &
                                 sfluxx, fx_lo, fx_hi, nc_fx, &
                                 sedgex, x_lo, x_hi, nc_x, &
                                 umac,   u_lo, u_hi, &
                                 rho0_old, rho0_edge_old, &
                                 rho0_new, rho0_edge_new, &
                                 w0, & 
                                 startcomp, endcomp) bind(C,name="make_rhoX_flux_1d")

    integer         , intent(in   ) :: domlo(1), domhi(1), lo(1), hi(1)
    integer         , intent(in   ) :: fx_lo(1), fx_hi(1), nc_fx
    double precision, intent(inout) :: sfluxx(fx_lo(1):fx_hi(1),nc_fx)
    integer         , intent(in   ) :: x_lo(1), x_hi(1), nc_x
    double precision, intent(inout) :: sedgex(x_lo(1):x_hi(1),nc_x)
    integer         , intent(in   ) :: u_lo(1), u_hi(1)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i
    double precision :: rho0_edge

    do comp = startcomp, endcomp

       ! create x-fluxes
       do i = lo(1),hi(1)+1

          if (species_pred_type == predict_rhoprime_and_X) then
             ! edge states are rho' and X.  To make the (rho X) flux,
             ! we need the edge state of rho0
             rho0_edge = HALF*(rho0_edge_old(i)+rho0_edge_new(i))
             sfluxx(i,comp) = &
                  (umac(i)+w0(i))*(rho0_edge+sedgex(i,rho_comp))*sedgex(i,comp)
          
          else if (species_pred_type == predict_rhoX) then
             ! edge states are (rho X)
             sfluxx(i,comp) = &
                  (umac(i)+w0(i))*sedgex(i,comp)             
             
          else if (species_pred_type == predict_rho_and_X) then
             ! edge states are rho and X
             sfluxx(i,comp) = &
                  (umac(i)+w0(i))*sedgex(i,rho_comp)*sedgex(i,comp)
          endif
          
       end do
    end do ! end comp loop

  end subroutine make_rhoX_flux_1d
#endif

#if (AMREX_SPACEDIM == 2)
  subroutine make_rhoX_flux_2d(domlo, domhi, lo, hi, &
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

    integer         , intent(in   ) :: domlo(2), domhi(2), lo(2), hi(2)
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
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i,j
    double precision :: rho0_edge

    do comp = startcomp, endcomp 

    end do

  end subroutine make_rhoX_flux_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine make_rhoX_flux_3d(domlo, domhi, lo, hi, &
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

    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
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
    double precision, intent(in   ) :: rho0_edge_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1) 
    double precision, intent(in   ) :: rho0_edge_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: startcomp, endcomp

    ! Local variables
    integer          :: comp
    integer          :: i,j,k
    double precision :: rho0_edge
    
    do comp = startcomp, endcomp 

    end do

  end subroutine make_rhoX_flux_3d
#endif

end module make_flux_module
