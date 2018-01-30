
#include <Maestro.H>

using namespace amrex;


// print out the contents of a cell-centered Vector of MultiFabs
void
Maestro::PrintCC (const Vector<MultiFab>& CC)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& CC_mf = CC[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(CC_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            print_cc(lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                     BL_TO_FORTRAN_FAB(CC_mf[mfi]), CC_mf.nGrow());

            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

void
Maestro::PrintEdge (const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& EDGE_mf = EDGE[lev][dir];
        const MultiFab& sold_mf = sold[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            print_edge(lev, dir, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                       BL_TO_FORTRAN_FAB(EDGE_mf[mfi]), EDGE_mf.nGrow());

            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

void
Maestro::PrintNodal (const Vector<MultiFab>& NODAL)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& NODAL_mf = NODAL[lev];
        const MultiFab& sold_mf = sold[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            print_nodal(lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                       BL_TO_FORTRAN_FAB(NODAL_mf[mfi]), NODAL_mf.nGrow());

            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

