! Create the psi term, where psi = D_0 p_0/Dt 

module make_psi_module

  use bl_constants_module
  use base_state_geometry_module, only: nr_fine, dr, &
                                        max_radial_level, numdisjointchunks, & 
                                        r_start_coord, r_end_coord, finest_radial_level, &
                                        restrict_base, fill_ghost_base, base_cutoff_density_coord
  use meth_params_module, only: grav_const

  implicit none

  private

contains

  subroutine make_psi_planar(etarho_cc,psi) bind(C, name="make_psi_planar")

    double precision, intent(in   ) :: etarho_cc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::       psi(0:max_radial_level,0:nr_fine-1)
    
    ! Local variables
    integer         :: r,i,n
   
    psi = ZERO

    do n=0,finest_radial_level
       do i=1,numdisjointchunks(n)
          do r = r_start_coord(n,i), r_end_coord(n,i)
             if (r .lt. base_cutoff_density_coord(n)) then
                psi(n,r) = etarho_cc(n,r) * abs(grav_const)
             end if
          end do
       end do
    end do

    call restrict_base(psi,1)
    call fill_ghost_base(psi,1)
    
  end subroutine make_psi_planar

  subroutine make_psi_spherical(psi,w0,gamma1bar,p0_avg,Sbar_in,r_cc_loc,r_edge_loc) bind(C, name="make_psi_spherical")

    double precision, intent(inout) ::       psi(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::        w0(0:max_radial_level,0:nr_fine  )
    double precision, intent(in   ) :: gamma1bar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::    p0_avg(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   Sbar_in(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::  r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine )
    
    ! local variables
    integer :: r
    real(kind=dp_t) :: div_w0_sph

    psi = ZERO

    !$OMP PARALLEL DO PRIVATE(r,div_w0_sph)
    do r=0,base_cutoff_density_coord(0)-1

       div_w0_sph = one/(r_cc_loc(0,r)**2)* &
            (r_edge_loc(0,r+1)**2 * w0(0,r+1) - &
             r_edge_loc(0,r  )**2 * w0(0,r  )) / dr(0)

       psi(0,r) = gamma1bar(0,r) * p0_avg(0,r) * (Sbar_in(0,r) - div_w0_sph)

    enddo
    !$OMP END PARALLEL DO

  end subroutine make_psi_spherical
  
end module make_psi_module
