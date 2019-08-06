module burner_loop_module

  use amrex_error_module
  use burner_module
  use burn_type_module, only: burn_t, copy_burn_t
  use network, only: nspec, network_species_index
  use meth_params_module, only: rho_comp, rhoh_comp, temp_comp, spec_comp, &
       pi_comp, nscal, burner_threshold_cutoff, burner_threshold_species, &
       burning_cutoff_density, reaction_sum_tol, &
       drive_initial_convection
  use base_state_geometry_module, only: max_radial_level, nr_fine

  implicit none

  public :: burner_loop_init

  private

  integer, allocatable, save :: ispec_threshold

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ispec_threshold
#endif

contains

  subroutine burner_loop_init()

      allocate(ispec_threshold)

      ispec_threshold = network_species_index(burner_threshold_species)

  end subroutine burner_loop_init

  subroutine burner_loop(lo, hi, &
       lev, &
       s_in,     i_lo, i_hi, &
       s_out,    o_lo, o_hi, &
       rho_Hext, e_lo, e_hi, &
       rho_odot, r_lo, r_hi, &
       rho_Hnuc, n_lo, n_hi, &
       tempbar_init_in, dt_in, &
       mask,     m_lo, m_hi, use_mask) &
       bind (C,name="burner_loop")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer, value  , intent (in   ) :: lev
    integer         , intent (in   ) :: i_lo(3), i_hi(3)
    integer         , intent (in   ) :: o_lo(3), o_hi(3)
    integer         , intent (in   ) :: e_lo(3), e_hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    double precision, intent (in   ) ::    s_in (i_lo(1):i_hi(1),i_lo(2):i_hi(2),i_lo(3):i_hi(3),nscal)
    double precision, intent (inout) ::    s_out(o_lo(1):o_hi(1),o_lo(2):o_hi(2),o_lo(3):o_hi(3),nscal)
    double precision, intent (in   ) :: rho_Hext(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    double precision, intent (inout) :: rho_odot(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec)
    double precision, intent (inout) :: rho_Hnuc(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3))
    double precision, intent (in   ) :: tempbar_init_in(0:max_radial_level,0:nr_fine-1)
    double precision, value, intent (in) :: dt_in
    integer         , intent (in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    integer, value  , intent (in   ) :: use_mask

    ! local
    integer          :: i, j, k, n, r
    double precision :: rho,T_in

    double precision :: x_in(nspec)
    double precision :: x_out(nspec)
    double precision :: rhowdot(nspec)
    double precision :: rhoH
    double precision :: x_test
    logical          :: cell_valid
    double precision :: sumX

    type (burn_t)    :: state_in, state_out

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( use_mask .eq. 1 ) then
                if ( (mask(i,j,k).eq.1) ) cell_valid = .false.
             endif

             if (cell_valid) then
                rho = s_in(i,j,k,rho_comp)
                do n = 1, nspec
                   x_in(n) = s_in(i,j,k,n+spec_comp-1) / rho
                enddo

                if (drive_initial_convection) then
#if (AMREX_SPACEDIM == 2)
                   r = j
#elif (AMREX_SPACEDIM == 3)
                   r = k
#endif
                   T_in = tempbar_init_in(lev,r)
                else
                   T_in = s_in(i,j,k,temp_comp)
                endif

                ! Fortran doesn't guarantee short-circuit evaluation of logicals
                ! so we need to test the value of ispec_threshold before using it
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = x_in(ispec_threshold)
                else
                   x_test = 0.d0
                endif

                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho > burning_cutoff_density .and.                &
                     ( ispec_threshold < 0 .or.                       &
                     (ispec_threshold > 0 .and. x_test > burner_threshold_cutoff) ) ) then
                   ! Initialize burn state_in and state_out
                   state_in % e   = 0.0d0
                   state_in % rho = rho
                   state_in % T   = T_in
                   do n = 1, nspec
                      state_in % xn(n) = x_in(n)
                   enddo
                   call copy_burn_t(state_out, state_in)
                   call burner(state_in, state_out, dt_in)
                   do n = 1, nspec
                      x_out(n) = state_out % xn(n)
                   enddo
                   do n = 1, nspec
                      rhowdot(n) = state_out % rho * &
                           (state_out % xn(n) - state_in % xn(n)) / dt_in
                   enddo
                   rhoH = state_out % rho * (state_out % e - state_in % e) / dt_in
                else
                   x_out = x_in
                   rhowdot = 0.d0
                   rhoH = 0.d0
                endif

                ! check if sum{X_k} = 1
                sumX = 0.d0
                do n = 1, nspec
                   sumX = sumX + x_out(n)
                enddo
                if (abs(sumX - 1.d0) > reaction_sum_tol) then
#ifndef AMREX_USE_GPU
                   call amrex_error("ERROR: abundances do not sum to 1", abs(sumX-1.d0))
#endif
                   do n = 1, nspec
                      state_out % xn(n) = state_out % xn(n)/sumX
                   enddo
                endif

                ! pass the density and pi through
                s_out(i,j,k,rho_comp) = s_in(i,j,k,rho_comp)
                s_out(i,j,k,pi_comp) = s_in(i,j,k,pi_comp)

                ! update the species
                do n = 1, nspec
                   s_out(i,j,k,n+spec_comp-1) = x_out(n) * rho
                enddo

                ! store the energy generation and species create quantities
                do n = 1, nspec
                   rho_odot(i,j,k,n) = rhowdot(n)
                enddo
                rho_Hnuc(i,j,k) = rhoH

                ! update the enthalpy -- include the change due to external heating
                s_out(i,j,k,rhoh_comp) = s_in(i,j,k,rhoh_comp) &
                     + dt_in*rho_Hnuc(i,j,k) + dt_in*rho_Hext(i,j,k)

             endif
          enddo
       enddo
    enddo

  end subroutine burner_loop

  subroutine burner_loop_sphr(lo, hi, &
       s_in,     i_lo, i_hi, &
       s_out,    o_lo, o_hi, &
       rho_Hext, e_lo, e_hi, &
       rho_odot, r_lo, r_hi, &
       rho_Hnuc, n_lo, n_hi, &
       tempbar_init_cart, t_lo, t_hi, dt_in, &
       mask,     m_lo, m_hi, use_mask) &
       bind (C,name="burner_loop_sphr")

    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: i_lo(3), i_hi(3)
    integer         , intent (in   ) :: o_lo(3), o_hi(3)
    integer         , intent (in   ) :: e_lo(3), e_hi(3)
    integer         , intent (in   ) :: r_lo(3), r_hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: m_lo(3), m_hi(3)
    double precision, intent (in   ) ::    s_in (i_lo(1):i_hi(1),i_lo(2):i_hi(2),i_lo(3):i_hi(3),nscal)
    double precision, intent (inout) ::    s_out(o_lo(1):o_hi(1),o_lo(2):o_hi(2),o_lo(3):o_hi(3),nscal)
    double precision, intent (in   ) :: rho_Hext(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    double precision, intent (inout) :: rho_odot(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec)
    double precision, intent (inout) :: rho_Hnuc(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3))
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    double precision, intent (in   ) :: tempbar_init_cart(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, value, intent (in) :: dt_in
    integer         , intent (in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    integer, value  , intent (in   ) :: use_mask

    ! local
    integer          :: i, j, k, n
    double precision :: rho,T_in

    double precision :: x_in(nspec)
    double precision :: x_out(nspec)
    double precision :: rhowdot(nspec)
    double precision :: rhoH
    double precision :: x_test
    logical          :: cell_valid
    double precision :: sumX

    type (burn_t)    :: state_in, state_out

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( use_mask .eq. 1 ) then
                if ( (mask(i,j,k).eq.1) ) cell_valid = .false.
             endif

             if (cell_valid) then
                rho = s_in(i,j,k,rho_comp)
                do n = 1, nspec
                   x_in(n) = s_in(i,j,k,n+spec_comp-1) / rho
                enddo

                if (drive_initial_convection) then
                   T_in = tempbar_init_cart(i,j,k)
                else
                   T_in = s_in(i,j,k,temp_comp)
                endif

                ! Fortran doesn't guarantee short-circuit evaluation of logicals
                ! so we need to test the value of ispec_threshold before using it
                ! as an index in x_in
                if (ispec_threshold > 0) then
                   x_test = x_in(ispec_threshold)
                else
                   x_test = 0.d0
                endif

                ! if the threshold species is not in the network, then we burn
                ! normally.  if it is in the network, make sure the mass
                ! fraction is above the cutoff.
                if (rho > burning_cutoff_density .and.                &
                     ( ispec_threshold < 0 .or.                       &
                     (ispec_threshold > 0 .and. x_test > burner_threshold_cutoff) ) ) then
                   ! Initialize burn state_in and state_out
                   state_in % e   = 0.0d0
                   state_in % rho = rho
                   state_in % T   = T_in
                   do n = 1, nspec
                      state_in % xn(n) = x_in(n)
                   enddo
                   call copy_burn_t(state_out, state_in)
                   call burner(state_in, state_out, dt_in)
                   do n = 1, nspec
                      x_out(n) = state_out % xn(n)
                   enddo
                   do n = 1, nspec
                      rhowdot(n) = state_out % rho * &
                           (state_out % xn(n) - state_in % xn(n)) / dt_in
                   enddo
                   rhoH = state_out % rho * (state_out % e - state_in % e) / dt_in
                else
                   x_out = x_in
                   rhowdot = 0.d0
                   rhoH = 0.d0

                   ! if we didn't burn, make sure that our abundances sum to
                   ! 1 -- this shouldn't normally be an issue, but some
                   ! combination of AMR + hitting the low density cutoff
                   ! can introduce a small error
                   sumX = 0.d0
                   do n = 1, nspec
                      sumX = sumX + x_out(n)
                   enddo
                   if (abs(sumX - 1.d0) > reaction_sum_tol) then
                      do n = 1, nspec
                         state_out % xn(n) = state_out % xn(n)/sumX
                      enddo
                   endif
                endif

                ! pass the density and pi through
                s_out(i,j,k,rho_comp) = s_in(i,j,k,rho_comp)
                s_out(i,j,k,pi_comp) = s_in(i,j,k,pi_comp)

                ! update the species
                do n = 1, nspec
                   s_out(i,j,k,n+spec_comp-1) = x_out(n) * rho
                enddo

                ! store the energy generation and species create quantities
                do n = 1, nspec
                   rho_odot(i,j,k,n) = rhowdot(n)
                enddo
                rho_Hnuc(i,j,k) = rhoH

                ! update the enthalpy -- include the change due to external heating
                s_out(i,j,k,rhoh_comp) = s_in(i,j,k,rhoh_comp) &
                     + dt_in*rho_Hnuc(i,j,k) + dt_in*rho_Hext(i,j,k)

             endif
          enddo
       enddo
    enddo

  end subroutine burner_loop_sphr

end module burner_loop_module
