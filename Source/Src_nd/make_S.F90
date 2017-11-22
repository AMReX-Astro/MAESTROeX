module make_S_module

  use eos_type_module
  use eos_module
  use amrex_fort_module, only: amrex_spacedim
  use network, only: nspec
  use meth_params_module, only: rho_comp, temp_comp, spec_comp
  use base_state_geometry_module, only:  max_radial_level, nr_fine

  implicit none

  private

  public :: make_S_cc, make_ccrhs

contains

  subroutine make_S_cc(lo, hi, &
                       S_cc,  s_lo, s_hi, &
                       scal,  c_lo, c_hi, nc_c, &
                       rodot, r_lo, r_hi, nc_r, &
                       rHnuc, n_lo, n_hi, &
                       rHext, e_lo, e_hi, &
                       therm, t_lo, t_hi) bind (C,name="make_S_cc")
    
    integer         , intent (in   ) :: lo(3), hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3), nc_c
    integer         , intent (in   ) :: r_lo(3), r_hi(3), nc_r
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: e_lo(3), e_hi(3)
    integer         , intent (in   ) :: t_lo(3), t_hi(3)
    double precision, intent (inout) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) :: scal (c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nc_c)
    double precision, intent (in   ) :: rodot(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nc_r)
    double precision, intent (in   ) :: rHnuc(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3))
    double precision, intent (in   ) :: rHext(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    double precision, intent (in   ) :: therm(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    integer i,j,k
    integer pt_index(3)
    type(eos_t) :: eos_state

    integer comp
    double precision sigma, xi_term, pres_term

    ! loop over the data
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
       
       eos_state%rho   = scal(i,j,k,rho_comp)
       eos_state%T     = scal(i,j,k,temp_comp)
       eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

       pt_index(:) = (/i, j, k/)

       ! dens, temp, and xmass are inputs
       call eos(eos_input_rt, eos_state, pt_index)

       sigma = eos_state%dpdt / &
            (eos_state%rho * eos_state%cp * eos_state%dpdr)

       xi_term = 0.d0
       pres_term = 0.d0
       do comp = 1, nspec
          xi_term = xi_term - &
               eos_state%dhdX(comp)*rodot(i,j,k,comp)/eos_state%rho 

          pres_term = pres_term + &
               eos_state%dpdX(comp)*rodot(i,j,k,comp)/eos_state%rho
       enddo

       S_cc(i,j,k) = (sigma/eos_state%rho) * &
            ( rHext(i,j,k) + rHnuc(i,j,k) + therm(i,j,k) ) &
            + sigma*xi_term &
            + pres_term/(eos_state%rho*eos_state%dpdr)

    enddo
    enddo
    enddo

  end subroutine make_S_cc

  subroutine make_ccrhs(lev, lo, hi, &
                        ccrhs, c_lo, c_hi, &
                        S_cc,  s_lo, s_hi, &
                        Sbar, div_coeff) bind (C,name="make_ccrhs")
    
    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    integer         , intent (in   ) :: s_lo(3), s_hi(3)
    double precision, intent (inout) :: ccrhs(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent (in   ) :: S_cc (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent (in   ) ::      Sbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent (in   ) :: div_coeff(0:max_radial_level,0:nr_fine-1)

    integer i,j,k,r

    ! loop over the data
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

#if (AMREX_SPACEDIM == 1)
       r = i
#elif (AMREX_SPACEDIM == 2)
       r = j
#elif (AMREX_SPACEDIM == 3)
       r = k
#endif
       ccrhs(i,j,k) = div_coeff(lev,r) * (S_cc(i,j,k) - Sbar(lev,r))

    enddo
    enddo
    enddo

  end subroutine make_ccrhs

  subroutine make_nodalrhs(lev, lo, hi, &
                           nodalrhs, n_lo, n_hi, &
                           ccrhs,    c_lo, c_hi) bind (C,name="make_nodalrhs")
    
    integer         , intent (in   ) :: lev, lo(3), hi(3)
    integer         , intent (in   ) :: n_lo(3), n_hi(3)
    integer         , intent (in   ) :: c_lo(3), c_hi(3)
    double precision, intent (inout) :: nodalrhs(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3))
    double precision, intent (in   ) ::    ccrhs(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

    integer i,j,k
    integer joff,koff

    joff = 0
    koff = 0
    if (amrex_spacedim .ge. 2) joff=1
    if (amrex_spacedim .ge. 3) koff=1

    ! loop over the data
    do k = lo(3),hi(3)+koff
    do j = lo(2),hi(2)+joff
    do i = lo(1),hi(1)+1

#if (AMREX_SPACEDIM == 1)
       nodalrhs(i,j,k) = 0.5d0 * ( ccrhs(i,j,k) + ccrhs(i-1,j,k) )
#elif (AMREX_SPACEDIM == 2)
       nodalrhs(i,j,k) = 0.25d0 * ( ccrhs(i,j  ,k) + ccrhs(i-1,j  ,k) &
                                  + ccrhs(i,j-1,k) + ccrhs(i-1,j-1,k) )
#elif (AMREX_SPACEDIM == 3)
       nodalrhs(i,j,k) = 0.125d0 * ( ccrhs(i,j  ,k-1) + ccrhs(i-1,j  ,k-1) &
                                   + ccrhs(i,j-1,k-1) + ccrhs(i-1,j-1,k-1) &
                                   + ccrhs(i,j  ,k  ) + ccrhs(i-1,j  ,k  ) &
                                   + ccrhs(i,j-1,k  ) + ccrhs(i-1,j-1,k  ) )
#endif

    enddo
    enddo
    enddo

  end subroutine make_nodalrhs

end module make_S_module
