!
! Compute eta_rho = Avg { rho' U dot e_r }  (see paper III, Eq. 30)
!
! We keep make three quantities here:
!    etarho     is edge-centered
!    etarho_cc  is cell-centered
!
! For plane-parallel geometries, we compute etarho by averaging up
! interface fluxes (etarho_flux) created in mkflux.  We compute etarho_cc
! from etarho.
!
! For spherical geometries,
!      We construct a multifab containing {rho' (U dot e_r)} and
!      use the average routine to put it in cell-centers
!      on the base state to get etarho_cc.  We compute etarho from these
!      cell-centered quantites by averaging to the center.
!

module make_eta_module

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module
  use base_state_geometry_module, only: nr_fine, dr, &
       max_radial_level, numdisjointchunks, &
       r_start_coord, r_end_coord, finest_radial_level, &
       restrict_base, fill_ghost_base
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

  private

contains

  subroutine make_etarho_planar(etarho_ec, etarho_cc, &
       etarhosum, ncell) bind(C, name="make_etarho_planar")

    double precision, intent(  out) :: etarho_ec(0:max_radial_level,0:nr_fine)
    double precision, intent(  out) :: etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: etarhosum(0:nr_fine,0:max_radial_level) ! note swapped shaping
    double precision, intent(in   ) ::     ncell(0:nr_fine)

    ! Local variables
    integer :: r,i,n

    etarho_ec       = ZERO
    etarho_cc       = ZERO

    do n=0,finest_radial_level
       do i=1,numdisjointchunks(n)
          do r = r_start_coord(n,i), r_end_coord(n,i)+1
             etarho_ec(n,r) = etarhosum(r,n) / dble(ncell(n))
          end do
       end do
    end do

    ! These calls shouldn't be needed since the planar algorithm doesn't use
    ! these outside of this function, but this is just to be safe in case
    ! things change in the future.
    call restrict_base(etarho_ec,0)
    call fill_ghost_base(etarho_ec,0)

    ! make the cell-centered etarho_cc by averaging etarho to centers
    do n=0,finest_radial_level
       do i=1,numdisjointchunks(n)
          do r=r_start_coord(n,i),r_end_coord(n,i)
             etarho_cc(n,r) = HALF*(etarho_ec(n,r) + etarho_ec(n,r+1))
          enddo
       enddo
    enddo

    ! These calls shouldn't be needed since the planar algorithm only uses
    ! etarho_cc to make_psi, and then we fill ghost cells in make_psi, but
    ! this is just to be safe in case things change in the future
    call restrict_base(etarho_ec,1)
    call fill_ghost_base(etarho_ec,1)

  end subroutine make_etarho_planar

  subroutine sum_etarho(lev, domlo, domhi, lo, hi, &
       etarhoflux, x_lo, x_hi, &
       etarhosum) bind(C, name="sum_etarho")

    integer         , intent(in   ) :: lev, domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) :: etarhoflux(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    double precision, intent(inout) :: etarhosum(0:nr_fine,0:max_radial_level) ! note swapped shaping

    ! local
    integer :: i,j,k
    logical :: top_edge

    ! Sum etarho
#if (AMREX_SPACEDIM == 2)
    k = lo(3)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          etarhosum(j,lev) = etarhosum(j,lev) + etarhoflux(i,j,k)
       end do
    end do

#elif (AMREX_SPACEDIM == 3)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k,lev) = etarhosum(k,lev) + etarhoflux(i,j,k)
          end do
       end do
    end do

#endif


    ! we only add the contribution at the top edge if we are at the top of the domain
    ! this prevents double counting
    top_edge = .false.
#if (AMREX_SPACEDIM == 2)
    do i=1,numdisjointchunks(lev)
       if (hi(2) .eq. r_end_coord(lev,i)) then
          top_edge = .true.
       end if
    end do
    if(top_edge) then
       k = hi(3)
       j = hi(2)+1
       do i=lo(1),hi(1)
          etarhosum(j,lev) = etarhosum(j,lev) + etarhoflux(i,j,k)
       end do
    end if

#elif (AMREX_SPACEDIM == 3)
    do i=1,numdisjointchunks(lev)
       if (hi(3) .eq. r_end_coord(lev,i)) then
          top_edge = .true.
       end if
    end do
    if(top_edge) then
       k=hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             etarhosum(k,lev) = etarhosum(k,lev) + etarhoflux(i,j,k)
          end do
       end do
    end if

#endif

  end subroutine sum_etarho

  !---------------------------------------------------------------------------
  ! spherical routines
  !---------------------------------------------------------------------------

  subroutine construct_eta_cart(lo, hi, &
       rho_old, ro_lo, ro_hi, &
       rho_new, rn_lo, rn_hi, &
       umac,     u_lo,  u_hi, &
       vmac,     v_lo,  v_hi, &
       wmac,     w_lo,  w_hi, &
       w0macx,   x_lo,  x_hi, &
       w0macy,   y_lo,  y_hi, &
       w0macz,   z_lo,  z_hi, &
       normal,   n_lo,  n_hi, &
       eta_cart, e_lo,  e_hi, &
       rho0_old, rho0_new, dx, &
       r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) &
       bind(C, name="construct_eta_cart")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: ro_lo(3), ro_hi(3)
    integer         , intent(in   ) :: rn_lo(3), rn_hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    integer         , intent(in   ) :: n_lo(3), n_hi(3)
    integer         , intent(in   ) :: e_lo(3), e_hi(3)
    double precision, intent(in   ) ::  rho_old(ro_lo(1):ro_hi(1),ro_lo(2):ro_hi(2),ro_lo(3):ro_hi(3))
    double precision, intent(in   ) ::  rho_new(rn_lo(1):rn_hi(1),rn_lo(2):rn_hi(2),rn_lo(3):rn_hi(3))
    double precision, intent(in   ) ::     umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) ::     vmac(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent(in   ) ::     wmac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(in   ) ::   w0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    double precision, intent(in   ) ::   w0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    double precision, intent(in   ) ::   w0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    double precision, intent(in   ) ::   normal(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),3)
    double precision, intent(inout) :: eta_cart(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    double precision, intent(in   ) :: rho0_old(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0_new(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    ! Local
    double precision ::      rho0_nph(0:max_radial_level,0:nr_fine-1)

    double precision, pointer :: rho0_new_cart(:,:,:,:)
    double precision, pointer :: rho0_nph_cart(:,:,:,:)

    double precision :: U_dot_er
    integer :: i,j,k,r

    call bl_allocate(rho0_new_cart,lo,hi,1)
    call bl_allocate(rho0_nph_cart,lo,hi,1)

    ! put the time-centered base state density on a Cartesian patch.
    do r = 0, nr_fine-1
       rho0_nph(0,r) = HALF*(rho0_old(0,r) + rho0_new(0,r))
    enddo

    call put_1d_array_on_cart_sphr(lo,hi,rho0_new_cart,lo,hi,1,rho0_new,dx,0,0, &
         r_cc_loc,r_edge_loc, cc_to_r,ccr_lo,ccr_hi)
    call put_1d_array_on_cart_sphr(lo,hi,rho0_nph_cart,lo,hi,1,rho0_nph,dx,0,0, &
         r_cc_loc,r_edge_loc, cc_to_r,ccr_lo,ccr_hi)

    !$OMP PARALLEL DO PRIVATE(i,j,k,U_dot_er)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             U_dot_er = HALF*(    umac(i,j,k) +   umac(i+1,j,k) &
                  + w0macx(i,j,k) + w0macx(i+1,j,k)) * normal(i,j,k,1) + &
                  HALF*(    vmac(i,j,k) +   vmac(i,j+1,k) &
                  + w0macy(i,j,k) + w0macy(i,j+1,k)) * normal(i,j,k,2) + &
                  HALF*(    wmac(i,j,k) +   wmac(i,j,k+1) &
                  + w0macz(i,j,k) + w0macz(i,j,k+1)) * normal(i,j,k,3)

             ! construct time-centered [ rho' (U dot e_r) ]
             eta_cart(i,j,k) = (HALF*(rho_old(i,j,k) + rho_new(i,j,k)) - &
                  rho0_nph_cart(i,j,k,1)) * U_dot_er

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    call bl_deallocate(rho0_new_cart)
    call bl_deallocate(rho0_nph_cart)

  end subroutine construct_eta_cart

end module make_eta_module
