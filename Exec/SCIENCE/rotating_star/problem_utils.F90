
module problem_utils

  use model_parser_module
  use amrex_constants_module
  use network, only: nspec
  use fill_3d_data_module, only: put_1d_array_on_cart, put_1d_array_on_cart_sphr
  use base_state_geometry_module, only: nr_fine, max_radial_level, dr, nr
  use meth_params_module, only: prob_lo, spherical, model_file

  implicit none

  private

contains

  subroutine set_enuc_omegadot(rho0, rho_omegadot, rho_Hnuc) &
    bind(C, name="set_enuc_omegadot")

    double precision, intent(inout) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: rho_omegadot(0:max_radial_level,0:nr_fine-1, 1:nspec)
    double precision, intent(inout) :: rho_Hnuc(0:max_radial_level,0:nr_fine-1)

    integer          :: n, r, spec 
    double precision :: model_dr, rmax, starting_rad, rloc

    ! put the model quantities on the correct grid

    model_dr = (model_r(npts_model) - model_r(1)) / dble(npts_model-1)
    rmax = model_r(npts_model)

    do n=0,max_radial_level

        if (spherical .eq. 0) then
            starting_rad = prob_lo(AMREX_SPACEDIM)
         else
            starting_rad = ZERO
         endif

        do r=0,nr(n)-1
            rloc = starting_rad + (dble(r) + HALF)*dr(n)

            ! here we account for r > rmax of the model.hse array, assuming
            ! that the state stays constant beyond rmax
            rloc = min(rloc, rmax)

            do spec=1,nspec
                rho_omegadot(n, r, spec) = interpolate(rloc, ispec_model+nspec+spec-1) * rho0(n, r)
            enddo 

            rho_Hnuc(n, r) = interpolate(rloc, ienuc_model) * rho0(n, r)
        enddo
    enddo

  end subroutine set_enuc_omegadot

  subroutine put_1d_spec_array_on_cart(lo, hi, lev, &
    s0_cart, s0_cart_lo, s0_cart_hi, &
    s0, is_input_edge_centered) &
    bind(C, name="put_1d_spec_array_on_cart")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer  , value, intent(in   ) :: lev
    integer         , intent(in   ) :: s0_cart_lo(3), s0_cart_hi(3)
    double precision, intent(inout) :: s0_cart(s0_cart_lo(1):s0_cart_hi(1), &
        s0_cart_lo(2):s0_cart_hi(2), &
        s0_cart_lo(3):s0_cart_hi(3), 1:nspec)
    double precision, intent(inout) :: s0(0:max_radial_level,0:nr_fine-1+is_input_edge_centered,1:nspec)
    integer  , value, intent(in   ) :: is_input_edge_centered

    ! local
    integer spec

    !$gpu

    do spec = 1, nspec 
        call put_1d_array_on_cart(lo, hi, lev, s0_cart(:,:,:,spec), s0_cart_lo, s0_cart_hi, 1, &
        s0(:,:,spec), is_input_edge_centered, 0)
    enddo 

  end subroutine put_1d_spec_array_on_cart

  subroutine put_1d_spec_array_on_cart_sphr(lo, hi, &
    s0_cart, s0_cart_lo, s0_cart_hi, &
    s0, dx, &
    is_input_edge_centered, &
    r_cc_loc, r_edge_loc, &
    cc_to_r, ccr_lo, ccr_hi) &
    bind(C, name="put_1d_spec_array_on_cart_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s0_cart_lo(3), s0_cart_hi(3)
    double precision, intent(inout) :: s0_cart(s0_cart_lo(1):s0_cart_hi(1), &
        s0_cart_lo(2):s0_cart_hi(2), &
        s0_cart_lo(3):s0_cart_hi(3), 1:nspec)
    double precision, intent(in   ) :: s0(0:max_radial_level,0:nr_fine-1+is_input_edge_centered,1:nspec)
    double precision, intent(in   ) :: dx(3)
    integer  , value, intent(in   ) :: is_input_edge_centered
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
        ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

        ! Local variables
    integer          :: spec

    !$gpu

    do spec = 1, nspec 
        call put_1d_array_on_cart_sphr(lo, hi, &
        s0_cart(:,:,:,spec), s0_cart_lo, s0_cart_hi, 1, &
        s0(:,:,spec), dx, &
        is_input_edge_centered, &
        0, &
        r_cc_loc, r_edge_loc, &
        cc_to_r, ccr_lo, ccr_hi)
    enddo 

  end subroutine put_1d_spec_array_on_cart_sphr

end module problem_utils
