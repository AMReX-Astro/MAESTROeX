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

#if (AMREX_SPACEDIM == 3)
  
  !----------------------------------------------------------------------------
  ! make_rhoX_flux_3d_sphr
  !----------------------------------------------------------------------------
  subroutine make_rhoX_flux_3d_sphr(lo, hi, &
       sflux, f_lo, f_hi, nc_f, &
       sedge, x_lo, x_hi, nc_x, &
       umac,  u_lo, u_hi, &
       w0mac, w_lo, w_hi, &
       rho0_edge, r_lo, r_hi, &
       startcomp, endcomp) bind(C,name="make_rhoX_flux_3d_sphr")
    ! Binds to C function ``make_rhoX_flux_3d_sphr``

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer  , value, intent(in   ) :: nc_f
    double precision, intent(inout) :: sflux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer  , value, intent(in   ) :: nc_x
    double precision, intent(inout) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(inout) :: w0mac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    integer         , intent(in   ) :: r_lo(3), r_hi(3)
    double precision, intent(in   ) :: rho0_edge(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    integer  , value, intent(in   ) :: startcomp, endcomp

    ! local
    integer          :: comp
    integer          :: i,j,k

    !$gpu

    ! reset density flux
    sflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = ZERO

    ! loop for fluxes
    do comp = startcomp, endcomp
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (species_pred_type == predict_rhoprime_and_X) then
                   ! edge states are rho' and X.  To make the (rho X)
                   ! flux, we need the edge state of rho0
                   sflux(i,j,k,comp) = umac(i,j,k) * &
                        (rho0_edge(i,j,k) + sedge(i,j,k,rho_comp))*sedge(i,j,k,comp)

                else if (species_pred_type == predict_rhoX) then
                   ! edge states are (rho X)
                   sflux(i,j,k,comp) = umac(i,j,k) * sedge(i,j,k,comp)

                else if (species_pred_type == predict_rho_and_X) then
                   ! edge states are rho and X
                   sflux(i,j,k,comp) = umac(i,j,k) * &
                        sedge(i,j,k,rho_comp)*sedge(i,j,k,comp)

                endif

                ! compute the density fluxes by summing the species fluxes
                sflux(i,j,k,1) = sflux(i,j,k,1) + sflux(i,j,k,comp)

             end do
          end do
       end do
    end do ! end loop over components

  end subroutine make_rhoX_flux_3d_sphr
#endif


#if (AMREX_SPACEDIM == 2)
  subroutine make_rhoh_flux_2d(lo, hi, lev, idir, &
       sflux, f_lo, f_hi, nc_f, &
       sedge, x_lo, x_hi, nc_x, &
       umac,   u_lo, u_hi, &
       rho0_old, rho0_edge_old, &
       rho0_new, rho0_edge_new, &
       rhoh0_old, rhoh0_edge_old, &
       rhoh0_new, rhoh0_edge_new, &
       w0) bind(C,name="make_rhoh_flux_2d")
    ! Binds to C function ``make_rhoh_flux_2d``

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev, idir
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer  , value, intent(in   ) :: nc_f
    double precision, intent(inout) :: sflux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer  , value, intent(in   ) :: nc_x
    double precision, intent(in   ) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
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
    logical          :: have_h, have_hprime, have_rhoh

    !$gpu

    have_h = enthalpy_pred_type.eq.predict_h .or. &
         enthalpy_pred_type.eq.predict_T_then_h .or. &
         enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    if (idir == 1) then

       ! create x-fluxes
       if (have_h) then
          ! enthalpy edge state is h

          if (species_pred_type == predict_rhoprime_and_X) then
             ! density edge state is rho'
             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   rho0_edge = HALF*(rho0_old(lev,j)+rho0_new(lev,j))
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*(rho0_edge+sedge(i,j,k,rho_comp))*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          else if (species_pred_type == predict_rho_and_X .or. &
               species_pred_type == predict_rhoX) then
             ! density edge state is rho
             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*sedge(i,j,k,rho_comp)*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          endif

       else if (have_hprime) then

          ! enthalpy edge state is h'
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_2d : predict_hprime not coded yet")
#endif

       else if (have_rhoh) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (enthalpy_pred_type.eq.predict_rhohprime .or. &
        enthalpy_pred_type.eq.predict_T_then_rhohprime) then
          ! enthalpy edge state is (rho h)'
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                rhoh0_edge = HALF*(rhoh0_old(lev,j)+rhoh0_new(lev,j))
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edge+sedge(i,j,k,rhoh_comp))
                end do
             end do
          end do

       else
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_2d : enthalpy_pred_type not recognized.")
#endif
       end if

    else ! idir == 2

       ! create y-fluxes
       if (have_h) then
          ! enthalpy edge state is h

          if (species_pred_type == predict_rhoprime_and_X) then
             ! density edge state is rho'
             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   rho0_edge = HALF*(rho0_edge_old(lev,j)+rho0_edge_new(lev,j))
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*(rho0_edge+sedge(i,j,k,rho_comp))*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          else if (species_pred_type == predict_rho_and_X .or. &
               species_pred_type == predict_rhoX) then
             ! density edge state is rho
             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*sedge(i,j,k,rho_comp)*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          endif

       else if (have_hprime) then

          ! enthalpy edge state is h'
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_2d : predict_hprime not coded yet")
#endif

       else if (have_rhoh) then
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (enthalpy_pred_type.eq.predict_rhohprime .or. &
        enthalpy_pred_type.eq.predict_T_then_rhohprime) then
          ! enthalpy edge state is (rho h)'
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                rhoh0_edge = HALF*(rhoh0_edge_old(lev,j)+rhoh0_edge_new(lev,j))
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*(sedge(i,j,k,rhoh_comp)+rhoh0_edge)
                end do
             end do
          end do

       else
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_2d : enthalpy_pred_type not recognized.")
#endif
       end if

    endif

  end subroutine make_rhoh_flux_2d
#endif


#if (AMREX_SPACEDIM == 3)
  subroutine make_rhoh_flux_3d(lo, hi, lev, idir, &
       sflux, f_lo, f_hi, nc_f, &
       sedge, x_lo, x_hi, nc_x, &
       umac,   u_lo, u_hi, &
       rho0_old, rho0_edge_old, &
       rho0_new, rho0_edge_new, &
       rhoh0_old, rhoh0_edge_old, &
       rhoh0_new, rhoh0_edge_new, &
       w0) bind(C,name="make_rhoh_flux_3d")
    ! Binds to C function ``make_rhoh_flux_3d``

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev, idir
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer  , value, intent(in   ) :: nc_f
    double precision, intent(inout) :: sflux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer  , value, intent(in   ) :: nc_x
    double precision, intent(in   ) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
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

    !$gpu

    have_h = enthalpy_pred_type.eq.predict_h .or. &
         enthalpy_pred_type.eq.predict_T_then_h .or. &
         enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    if (idir == 1) then

       ! create x-fluxes
       if (have_h) then
          ! enthalpy edge state is h

          if (species_pred_type == predict_rhoprime_and_X) then
             ! density edge state is rho'

             do k=lo(3),hi(3)
                rho0_edge = HALF*(rho0_old(lev,k)+rho0_new(lev,k))
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*(rho0_edge+sedge(i,j,k,rho_comp))*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          else if (species_pred_type == predict_rho_and_X .or. &
               species_pred_type == predict_rhoX) then
             ! density edge state is rho

             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*sedge(i,j,k,rho_comp)*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          endif

       else if (have_hprime) then

          ! enthalpy edge state is h'
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_3d : predict_hprime not coded yet")
#endif

       else if (have_rhoh) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (enthalpy_pred_type.eq.predict_rhohprime .or. &
        enthalpy_pred_type.eq.predict_T_then_rhohprime) then
          ! enthalpy edge state is (rho h)'

          do k=lo(3),hi(3)
             rhoh0_edge = HALF*(rhoh0_old(lev,k)+rhoh0_new(lev,k))
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edge+sedge(i,j,k,rhoh_comp))
                end do
             end do
          end do

       else
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
#endif
       end if

    elseif (idir == 2) then

       ! create y-fluxes
       if (have_h) then
          ! enthalpy edge state is h

          if (species_pred_type == predict_rhoprime_and_X) then
             ! density edge state is rho'

             do k=lo(3),hi(3)
                rho0_edge = HALF*(rho0_old(lev,k)+rho0_new(lev,k))
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*(rho0_edge+sedge(i,j,k,rho_comp))*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          else if (species_pred_type == predict_rho_and_X .or. &
               species_pred_type == predict_rhoX) then
             ! density edge state is rho

             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = &
                           umac(i,j,k)*sedge(i,j,k,rho_comp)*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          endif

       else if (have_hprime) then

          ! enthalpy edge state is h'
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_3d : predict_hprime not coded yet")
#endif

       else if (have_rhoh) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (enthalpy_pred_type.eq.predict_rhohprime .or. &
        enthalpy_pred_type.eq.predict_T_then_rhohprime) then
          ! enthalpy edge state is (rho h)'

          do k=lo(3),hi(3)
             rhoh0_edge = HALF*(rhoh0_old(lev,k)+rhoh0_new(lev,k))
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edge+sedge(i,j,k,rhoh_comp))
                end do
             end do
          end do

       else
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
#endif
       end if


    else ! idir == 3

       ! create z-fluxes
       if (have_h) then
          ! enthalpy edge state is h

          if (species_pred_type == predict_rhoprime_and_X) then
             ! density edge state is rho'

             do k=lo(3),hi(3)
                rho0_edge = HALF*(rho0_edge_old(lev,k)+rho0_edge_new(lev,k))
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = umac(i,j,k)* &
                           (rho0_edge+sedge(i,j,k,rho_comp))*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          else if (species_pred_type == predict_rho_and_X .or. &
               species_pred_type == predict_rhoX) then
             ! density edge state is rho

             do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      sflux(i,j,k,rhoh_comp) = umac(i,j,k)* &
                           sedge(i,j,k,rho_comp)*sedge(i,j,k,rhoh_comp)
                   end do
                end do
             end do

          endif

       else if (have_hprime) then

          ! enthalpy edge state is h'
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_3d : predict_hprime not coded yet")
#endif

       else if (have_rhoh) then
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)
                end do
             end do
          end do

       else if (enthalpy_pred_type.eq.predict_rhohprime .or. &
        enthalpy_pred_type.eq.predict_T_then_rhohprime) then
          ! enthalpy edge state is (rho h)'

          do k=lo(3),hi(3)
             rhoh0_edge = HALF*(rhoh0_edge_old(lev,k)+rhoh0_edge_new(lev,k))
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   sflux(i,j,k,rhoh_comp) = &
                        umac(i,j,k)*(sedge(i,j,k,rhoh_comp)+rhoh0_edge)
                end do
             end do
          end do

       else
#ifndef AMREX_USE_CUDA
          call amrex_error("make_rhoh_flux_3d : enthalpy_pred_type not recognized.")
#endif
       end if
    end if


  end subroutine make_rhoh_flux_3d

  !----------------------------------------------------------------------------
  ! mk_rhoh_flux_3d_sphr
  !----------------------------------------------------------------------------
  subroutine make_rhoh_flux_3d_sphr(lo, hi, &
       sflux, fx_lo, fx_hi, nc_fx, &
       sedge, x_lo, x_hi, nc_x, &
       umac,   u_lo, u_hi, &
       w0mac, wx_lo, wx_hi, &
       rho0_edge, rx_lo, rx_hi, &
       h0_edge, hx_lo, hx_hi) bind(C,name="make_rhoh_flux_3d_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3)
    integer  , value, intent(in   ) :: nc_fx
    double precision, intent(inout) :: sflux(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),nc_fx)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer  , value, intent(in   ) :: nc_x
    double precision, intent(in   ) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    double precision, intent(in   ) :: w0mac(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    integer         , intent(in   ) :: rx_lo(3), rx_hi(3)
    double precision, intent(in   ) :: rho0_edge(rx_lo(1):rx_hi(1),rx_lo(2):rx_hi(2),rx_lo(3):rx_hi(3))
    integer         , intent(in   ) :: hx_lo(3), hx_hi(3)
    double precision, intent(in   ) :: h0_edge(hx_lo(1):hx_hi(1),hx_lo(2):hx_hi(2),hx_lo(3):hx_hi(3))

    ! local
    integer          :: i,j,k
    logical          :: have_h, have_hprime, have_rhoh

    !$gpu

    have_h = enthalpy_pred_type.eq.predict_h .or. &
         enthalpy_pred_type.eq.predict_T_then_h .or. &
         enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    ! create fluxes
    if (have_h) then

       ! enthalpy edge state is h

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   sflux(i,j,k,rhoh_comp) = umac(i,j,k) * &
                        (rho0_edge(i,j,k) + sedge(i,j,k,rho_comp))*sedge(i,j,k,rhoh_comp)

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   sflux(i,j,k,rhoh_comp) = umac(i,j,k) * &
                        sedge(i,j,k,rho_comp)*sedge(i,j,k,rhoh_comp)

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
                do i = lo(1), hi(1)

                   sflux(i,j,k,rhoh_comp) = umac(i,j,k) * &
                        (sedge(i,j,k,rho_comp)+rho0_edge(i,j,k)) * (sedge(i,j,k,rhoh_comp)+h0_edge(i,j,k))

                end do
             end do
          end do

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
#ifndef AMREX_USE_CUDA
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")
#endif

       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
#ifndef AMREX_USE_CUDA
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")
#endif
       endif

    else if (have_rhoh) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                ! where h_0 is computed from (rho h)_0 / rho_0

                sflux(i,j,k,rhoh_comp) = &
                     umac(i,j,k)*(rho0_edge(i,j,k)*h0_edge(i,j,k)+sedge(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

  end subroutine make_rhoh_flux_3d_sphr

  subroutine make_rhoh_flux_3d_sphr_irreg(lo, hi, &
       sflux, fx_lo, fx_hi, nc_fx, &
       sedge, x_lo, x_hi, nc_x, &
       umac,   u_lo, u_hi, &
       w0mac, wx_lo, wx_hi, &
       rhoh0_edge, hx_lo, hx_hi) bind(C,name="make_rhoh_flux_3d_sphr_irreg")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3)
    integer  , value, intent(in   ) :: nc_fx
    double precision, intent(inout) :: sflux(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),nc_fx)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer  , value, intent(in   ) :: nc_x
    double precision, intent(in   ) :: sedge(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),nc_x)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac  (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: wx_lo(3), wx_hi(3)
    double precision, intent(in   ) :: w0mac(wx_lo(1):wx_hi(1),wx_lo(2):wx_hi(2),wx_lo(3):wx_hi(3))
    integer         , intent(in   ) :: hx_lo(3), hx_hi(3)
    double precision, intent(in   ) :: rhoh0_edge(hx_lo(1):hx_hi(1),hx_lo(2):hx_hi(2),hx_lo(3):hx_hi(3))

    ! local
    integer          :: i,j,k
    logical          :: have_h, have_hprime, have_rhoh

    !$gpu

    have_h = enthalpy_pred_type.eq.predict_h .or. &
         enthalpy_pred_type.eq.predict_T_then_h .or. &
         enthalpy_pred_type.eq.predict_Tprime_then_h

    have_hprime = enthalpy_pred_type.eq.predict_hprime

    have_rhoh = enthalpy_pred_type.eq.predict_rhoh

    ! create fluxes
    if (have_h) then

       ! enthalpy edge state is h
#ifndef AMREX_USE_CUDA
       call amrex_error("have_h not supported on irregular-spaced base state")
#endif

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X .or. &
            species_pred_type == predict_rhoX) then
          ! density edge state is rho

       endif

    else if (have_hprime) then

       ! enthalpy edge state is h'
#ifndef AMREX_USE_CUDA
       call amrex_error("have_hprime not supported on irregular-spaced base state")
#endif

       if (species_pred_type == predict_rhoprime_and_X) then
          ! density edge state is rho'

       else if (species_pred_type == predict_rho_and_X) then
          ! density edge state is rho
#ifndef AMREX_USE_CUDA
          call amrex_error("ERROR: predict_rho_and_X and predict_hprime not supported together")
#endif
       else if (species_pred_type == predict_rhoX) then
          ! density edge state is rho
#ifndef AMREX_USE_CUDA
          call amrex_error("ERROR: predict_rhoX and predict_hprime not supported together")
#endif
       endif

    else if (have_rhoh) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                sflux(i,j,k,rhoh_comp) = umac(i,j,k)*sedge(i,j,k,rhoh_comp)

             end do
          end do
       end do

    else

       ! enthalpy edge state is (rho h)'

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Average (rho h) onto edges by averaging rho and h
                ! separately onto edges.
                !  (rho h)_edge = (rho h)' + (rhoh_0)

                sflux(i,j,k,rhoh_comp) = umac(i,j,k)*(rhoh0_edge(i,j,k)+sedge(i,j,k,rhoh_comp))

             end do
          end do
       end do

    endif

  end subroutine make_rhoh_flux_3d_sphr_irreg
#endif

end module make_flux_module
