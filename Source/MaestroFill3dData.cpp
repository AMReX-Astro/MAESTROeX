/*
A collection of routines for mapping 1D arrays onto multi-D cartesian MultiFabs
*/

#include <Maestro.H>

using namespace amrex;

void
Maestro::Put1dArrayOnCart (const Vector<Real>& s0,
                           Vector<MultiFab>& s0_cart,
                           Vector<BCRec>& bcs,
                           int is_input_edge_centered,
                           int is_output_a_vector)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& s0_cart_mf = s0_cart[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(s0_cart_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            put_1d_array_on_cart(lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                                 BL_TO_FORTRAN_FAB(s0_cart_mf[mfi]),
                                 s0.dataPtr(), is_input_edge_centered, is_output_a_vector);
        }
    }

    int ncomp = is_output_a_vector ? AMREX_SPACEDIM : 1;

    // set covered coarse cells to be the average of overlying fine cells
    AverageDown(s0_cart,0,ncomp);

    // fill ghost cells using first-order extrapolation
    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, t_new, s0_cart[lev], s0_cart, s0_cart, 0, 0, ncomp, bcs);
    }
}
