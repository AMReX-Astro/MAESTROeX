
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeGamma1bar (const Vector<MultiFab>& scal,
                        Vector<Real>& gamma1bar,
                        const Vector<Real>& p0)
{
    Vector<MultiFab> gamma1(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        gamma1[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
              MultiFab& gamma1_mf = gamma1[lev];
        const MultiFab&   scal_mf =   scal[lev];
        
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(gamma1_mf); mfi.isValid(); ++mfi ) {
            
            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_gamma(lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                       BL_TO_FORTRAN_3D(gamma1_mf[mfi]),
                       BL_TO_FORTRAN_FAB(scal_mf[mfi]),
                       p0.dataPtr());
        }
    }

    // average fine data onto coarser cells
    AverageDown(gamma1,0,1);

    // call average to create gamma1bar
    Average(gamma1,gamma1bar,0);
}
