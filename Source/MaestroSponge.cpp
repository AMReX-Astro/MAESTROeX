
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeSponge (Vector<MultiFab>& sponge)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeSponge()",MakeSponge);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& sponge_mf = sponge[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sponge_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "tileBox", which specifies the "valid" region.
            mk_sponge(ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                      BL_TO_FORTRAN_3D(sponge_mf[mfi]),ZFILL(dx),&dt);
        }

    }

    // average fine data onto coarser cells
    AverageDown(sponge,0,1);

}
