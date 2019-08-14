
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext
void
Maestro::MakeIntraCoeffs (const Vector<MultiFab>& s1,
			  const Vector<MultiFab>& s2,
			  Vector<MultiFab>& cp,
			  Vector<MultiFab>& xi)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeIntraCoeffs()",MakeIntraCoeffs);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
	const MultiFab& s1_mf = s1[lev];
	const MultiFab& s2_mf = s2[lev];
              MultiFab& cp_mf = cp[lev];
              MultiFab& xi_mf = xi[lev];


        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s1_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
	    
            // make_intra_coeffs(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
	    // 		      BL_TO_FORTRAN_3D(rho_Hext_mf[mfi]),
	    // 		      BL_TO_FORTRAN_FAB(scal_mf[mfi]), dx, &t_old);
        }
    }

    // average down
    AverageDown(cp,0,1);
}
