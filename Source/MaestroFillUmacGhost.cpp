


#include <Maestro.H>

using namespace amrex;

void
Maestro::FillUmacGhost (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac)
{
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab& sold_mf = sold[lev];  // need a cell-centered MF for the MFIter        
        MultiFab& umacx_mf = umac[lev][0];
        MultiFab& umacy_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& umacz_mf = umac[lev][2];
#endif

        for ( MFIter mfi(sold_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            fill_umac_ghost(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                            BL_TO_FORTRAN_3D(umacx_mf[mfi]),
                            BL_TO_FORTRAN_3D(umacy_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                            BL_TO_FORTRAN_3D(umacz_mf[mfi]),
#endif
                            lo_bc.dataPtr(),hi_bc.data());

        }
    }
}
