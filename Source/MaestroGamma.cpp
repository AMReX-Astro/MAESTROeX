
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
        const MultiFab& scal_mf = scal[lev];
        
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(gamma1_mf); mfi.isValid(); ++mfi ) {
            FArrayBox& gamma1_fab = gamma1_mf[mfi];
            const FArrayBox& scal_fab = scal_mf[mfi];

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const int* validLo = validBox.loVect();
            const int* validHi = validBox.hiVect();

            // Get the index space of the valid+ghost region for each FAB
            // Note each of these boxes may contain ghost cells, and thus are
            // larger than or equal to mfi.validbox().
            const Box& gamma1_box = gamma1_fab.box();
            const Box& scal_box = scal_fab.box();

            // We can now pass the information to a Fortran routine,
            // e.g. gamma1_fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "gamma1_box". We will also pass "box",
            // which specifies our "work" region .
            make_gamma(lev, ARLIM_3D(validLo), ARLIM_3D(validHi),
                       gamma1_fab.dataPtr(),
                       ARLIM_3D(gamma1_box.loVect()), ARLIM_3D(gamma1_box.hiVect()),
                       scal_fab.dataPtr(),
                       ARLIM_3D(scal_box.loVect()), ARLIM_3D(scal_box.hiVect()),
                       scal_fab.nComp(),
                       p0.dataPtr());
        }
    }

    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,gamma1); // average lev+1 down to lev
    }

    // call average to create gamma1bar
    Average(gamma1,gamma1bar,0);
}
