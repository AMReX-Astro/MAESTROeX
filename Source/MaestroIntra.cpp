
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute cp and xi
void
Maestro::MakeIntraCoeffs (const Vector<MultiFab>& scal1,
			  const Vector<MultiFab>& scal2,
			  Vector<MultiFab>& cp,
			  Vector<MultiFab>& xi)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeIntraCoeffs()",MakeIntraCoeffs);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
	const MultiFab& s1_mf = scal1[lev];
	const MultiFab& s2_mf = scal2[lev];
              MultiFab& cp_mf = cp[lev];
              MultiFab& xi_mf = xi[lev];

        Print() << "... Level " << lev << " create intra coeffs:" << std::endl;

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(s1_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& gtbx = mfi.growntilebox(1);
            const Real* dx = geom[lev].CellSize();
	    
            // call fortran subroutine
#pragma gpu box(gtbx)
            make_intra_coeffs(AMREX_INT_ANYD(gtbx.loVect()), AMREX_INT_ANYD(gtbx.hiVect()),
	     		      BL_TO_FORTRAN_ANYD(s1_mf[mfi]),
	     		      BL_TO_FORTRAN_ANYD(s2_mf[mfi]),
			      BL_TO_FORTRAN_ANYD(cp_mf[mfi]),
			      BL_TO_FORTRAN_ANYD(xi_mf[mfi]));
        }
    }
}
