
#include <Maestro.H>

using namespace amrex;

void
Maestro::Make_S_cc (Vector<std::unique_ptr<MultiFab> >& S_cc)
{

    for (int lev=0; lev<=finest_level; ++lev) 
    {
        MultiFab& S_cc_mf = *S_cc[lev];

        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) // Loop over boxes
        {
            // Get the index space of this iteration
            const Box& box = mfi.validbox();

            // Get a reference to the FAB , which contains data and box
            FArrayBox& fab = S_cc_mf[mfi];

            // Get the index space for the data region in the FAB.
            // Note "abox" may have ghost cells, and is thus larger than
            // or equal to "box" obtained using mfi.validbox().
            const Box & abox = fab.box();

            // We can now pass the information to a Fortran routine,
            // fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "abox". We will also pass "box",
            // which specifies our "work" region .
            make_S_cc(ARLIM_3D(box.loVect()),ARLIM_3D(box.hiVect()),
                      fab.dataPtr(), ARLIM_3D(abox.loVect()),ARLIM_3D(abox.hiVect()),
                      S_cc_mf.nGrow(), fab.nComp());
        }

    }

    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,S_cc); // average lev+1 down to lev
    }

}
