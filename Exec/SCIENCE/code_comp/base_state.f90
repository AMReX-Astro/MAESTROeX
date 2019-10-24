module base_state_module
  ! init_base_state is used to initialize the base state arrays from the
  ! model file.  The actual reading of the model file is handled by the
  ! model_parser_module in Util/
  !
  ! Note: The initial base state quantities returned from this routine
  ! are only a temporary base state.  These quantities are mapped onto
  ! the full 2- or 3-d state in initscaldata.f90 and a new base state is
  ! created after initialization by averaging the density and calling
  ! enforce_HSE in initialize.f90.

  use model_parser_module
  use eos_type_module
  use eos_module
  use amrex_constants_module
  use simple_log_module
  use inlet_bc_module
  use fundamental_constants_module, only: Gconst
  use amrex_fort_module, only: amrex_spacedim
  use network, only: nspec
  use meth_params_module, only: nscal, model_file, spherical, base_cutoff_density, &
                                do_2d_planar_octant, do_planar_invsq_grav, rho_comp, &
                                rhoh_comp, spec_comp, temp_comp, grav_const, &
                                planar_invsq_mass, print_init_hse_diag, prob_lo
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use model_util_module, only : set_species, fv, dUdy
  use probin_module, only: rho_0, p_0

  implicit none

  private

contains

  subroutine init_base_state(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init) &
       bind(C, name="init_base_state")
    ! Binds to C function ``init_base_state``

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)

    ! local variables

    double precision :: x, y, z, fheat, rhopert, U_old(2), U_new(2), h
    double precision :: k1(2), k2(2), k3(2), k4(2)
    double precision, allocatable :: pres(:), dens(:)
    double precision :: xn(nspec)
    integer :: i, j, k, n, r, n_dy
    type (eos_t) :: eos_state

    do n=0,max_radial_level

        allocate(pres(0:nr(n)-1))
        allocate(dens(0:nr(n)-1))

        U_old(1) = log(rho_0)
        U_old(2) = log(p_0)

        do r=0,nr(n)-1

            y = prob_lo(amrex_spacedim) + (dble(r) + HALF)*dr(n)

            ! do HSE using RK2

            ! our integration starts at y - h
            if (r .eq. 0) then
                h = dr(n) * HALF
            else
                h = dr(n)
            endif

            k1(:) = dUdy(y - h, U_old)
            U_new(:) = U_old(:) + h * dUdy(y - HALF*h, U_old + HALF*h * k1)

            dens(r) = exp(U_new(1))
            pres(r) = exp(U_new(2))

            U_old(:) = U_new(:)

        end do

        do r=0,nr(n)-1

            y = prob_lo(amrex_spacedim) + (dble(r) + HALF)*dr(n)

            if (r < 0) then
                eos_state % rho = rho_0
                eos_state % p = p_0
            else
                eos_state%rho = dens(r)
                eos_state%p = pres(r)
            endif

            eos_state%xn(:) = set_species(y)

            call eos(eos_input_rp, eos_state)

            s0_init(n, r, rho_comp) = eos_state % rho
            s0_init(n, r, rhoh_comp) = eos_state % rho * eos_state % h
            s0_init(n, r, spec_comp:spec_comp+nspec-1) = eos_state % rho * eos_state % xn(:)
            p0_init(n, r) = eos_state % p
            s0_init(n, r, temp_comp) = eos_state % t

        end do

       deallocate(pres)
       deallocate(dens)

       ! initialize any inlet BC parameters
       call set_inlet_bcs()

    end do ! end loop over levels

    ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    rho0(:,:) = s0_init(:,:,rho_comp)
    rhoh0(:,:) = s0_init(:,:,rhoh_comp)
    tempbar(:,:) = s0_init(:,:,temp_comp)
    tempbar_init(:,:) = s0_init(:,:,temp_comp)
    p0(:,:) = p0_init(:,:)

  end subroutine init_base_state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_base_state_irreg(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init, &
                                     r_cc_loc, r_edge_loc) &
       bind(C, name="init_base_state_irreg")
    ! Binds to C function ``init_base_state_irreg``

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

  end subroutine init_base_state_irreg
end module base_state_module
