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

  use bl_constants_module
  use base_state_geometry_module, only: nr_fine, dr, &
                                        max_radial_level, numdisjointchunks, & 
                                        r_start_coord, r_end_coord, finest_radial_level, &
                                        restrict_base, fill_ghost_base

  implicit none

  private

contains

  subroutine make_etarho_planar(etarho_ec, etarho_cc, &
                                  etarhosum, ncell) bind(C, name="make_etarho_planar")

    double precision, intent(  out) :: etarho_ec(0:max_radial_level,0:nr_fine)
    double precision, intent(  out) :: etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: etarhosum(0:nr_fine,0:max_radial_level)
    double precision, intent(in   ) :: ncell(0:nr_fine,0:max_radial_level)
    
    ! Local variables
    integer :: r,i,n
   
    etarho_ec       = ZERO
    etarho_cc       = ZERO

    do n=0,finest_radial_level
       do i=1,numdisjointchunks(n)
          do r = r_start_coord(n,i), r_end_coord(n,i)
             etarho_ec(n,r) = etarhosum(r,n) / dble(ncell(r,n))
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
                         etarhosum, ncell) bind(C, name="sum_etarho")

    integer         , intent(in   ) :: lev, domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: x_lo(3), x_hi(3) 
    double precision, intent(in   ) :: etarhoflux(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    double precision, intent(inout) :: etarhosum(0:nr_fine,0:max_radial_level)
    double precision, intent(inout) :: ncell(0:nr_fine,0:max_radial_level)

    ! local
    integer :: i,j,k
    logical :: top_edge

#if (AMREX_SPACEDIM == 2) 
       ncell(:,lev) = domhi(1)-domlo(1)+1
#elif (AMREX_SPACEDIM == 3)
       ncell(:,lev) = (domhi(1)-domlo(1)+1)*(domhi(2)-domlo(2)+1)
#endif

       ! Sum etarho
#if (AMREX_SPACEDIM == 1)
       k = lo(3)
       j = lo(2)
       do i=lo(1),hi(1)
          etarhosum(i,lev) = etarhoflux(i,j,k)
       end do
       
#elif (AMREX_SPACEDIM == 2)
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
#if (AMREX_SPACEDIM == 1) 
       do k = 1,numdisjointchunks(lev)
          if (hi(1) .eq. r_end_coord(lev,k)) then
             top_edge = .true.
          end if
       end do
       if (top_edge) then
          k = hi(3)
          j = hi(2)
          i = hi(1)+1
          etarhosum(i,lev) = etarhoflux(i,j,k)
       end if

#elif (AMREX_SPACEDIM == 2)
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
  
end module make_eta_module
