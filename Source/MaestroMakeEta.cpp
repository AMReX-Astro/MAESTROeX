
#include <Maestro.H>

using namespace amrex;

// compute eta_rho at edge- and cell-centers
void
Maestro::MakeEtarho (Vector<Real>& etarho_edge,
		     Vector<Real>& etarho_cell,
		     const Vector<MultiFab>& etarho_flux)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEtarho()",MakeEtarho);

    // Local variables
    Vector<Real> etarhosum( (nr_fine+1)*(max_radial_level+1) ); 
    Vector<Real> ncell( (nr_fine+1)*(max_radial_level+1) );
    etarhosum.shrink_to_fit();
    ncell.shrink_to_fit();

    for (int i=0; i<(nr_fine+1)*(max_radial_level+1); ++i) {
	etarhosum[i] = 0.0;
	ncell[i] = 0.0;
    }

    if (spherical == 0) {
	for (int lev=0; lev<=finest_level; ++lev) {

	    // get references to the MultiFabs at level lev
	    const MultiFab& sold_mf = sold[lev];
	    const MultiFab& etarhoflux_mf = etarho_flux[lev];
	    
	    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
	    for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {
		
		// Get the index space of the valid region
		const Box& validBox = mfi.validbox();
		const Box& domainBox = geom[lev].Domain();

		// call fortran subroutine
		// use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
		// lo/hi coordinates (including ghost cells)
		// We will also pass "validBox", which specifies the "valid" region.
		sum_etarho(&lev, ARLIM_3D(domainBox.loVect()), ARLIM_3D(domainBox.hiVect()),
			   ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			   BL_TO_FORTRAN_3D(etarhoflux_mf[mfi]),
			   etarhosum.dataPtr(), ncell.dataPtr());
	    }	    

	}

	ParallelDescriptor::ReduceRealSum(etarhosum.dataPtr(),(nr_fine+1)*(max_radial_level+1));

	make_etarho_planar(etarho_edge.dataPtr(), etarho_cell.dataPtr(), 
			   etarhosum.dataPtr(), ncell.dataPtr());

    } else {
	// call make_etarho_spherical(s1,s2,umac,w0mac,rho0_old,rho0_new,dx,normal, &
	//                          etarho_ec,etarho_cc,mla,the_bc_tower%bc_tower_array)
	//
    }

}
