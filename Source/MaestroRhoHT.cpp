
#include <Maestro.H>

using namespace amrex;

void
Maestro::TfromRhoH (Vector<MultiFab>& scal,
                    const Vector<Real>& p0)
{

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            makeTfromRhoH(lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                          BL_TO_FORTRAN_FAB(scal_mf[mfi]),
                          p0.dataPtr());
        }

    }

    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,scal,Temp,1); // average lev+1 down to lev
    }

}

void
Maestro::TfromRhoP (Vector<MultiFab>& scal,
                    const Vector<Real>& p0,
                    int updateRhoH)
{

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_mf = scal[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            makeTfromRhoP(lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                          BL_TO_FORTRAN_FAB(scal_mf[mfi]),
                          p0.dataPtr(),updateRhoH);
        }

    }

// average lev+1 down to lev
    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,scal,Temp,1); 
        if (updateRhoH == 1) {
            AverageDownTo(lev,scal,RhoH,1);
        }
    }

}

