
#include <Maestro.H>

using namespace amrex;

void
Maestro::Make_S_cc (Vector<MultiFab>& S_cc,
                    const Vector<MultiFab>& scal,
                    const Vector<MultiFab>& rho_omegadot,
                    const Vector<MultiFab>& rho_Hnuc,
                    const Vector<MultiFab>& rho_Hext,
                    const Vector<MultiFab>& thermal)
{

    for (int lev=0; lev<=finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        MultiFab& S_cc_mf = S_cc[lev];
        const MultiFab& scal_mf = scal[lev];
        const MultiFab& rho_odot_mf = rho_omegadot[lev];
        const MultiFab& rho_Hnuc_mf = rho_Hnuc[lev];
        const MultiFab& rho_Hext_mf = rho_Hext[lev];
        const MultiFab& thermal_mf = thermal[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) {
            // get references to the FABs, each containing data and the valid+ghost box
            FArrayBox& S_cc_fab = S_cc_mf[mfi];
            const FArrayBox& scal_fab = scal_mf[mfi];
            const FArrayBox& rho_odot_fab = rho_odot_mf[mfi];
            const FArrayBox& rho_Hnuc_fab = rho_Hnuc_mf[mfi];
            const FArrayBox& rho_Hext_fab = rho_Hext_mf[mfi];
            const FArrayBox& thermal_fab = thermal_mf[mfi];

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const int* validLo = validBox.loVect();
            const int* validHi = validBox.hiVect();

            // Get the index space of the valid+ghost region for each FAB
            // Note each of these boxes may contain ghost cells, and thus are
            // larger than or equal to mfi.validbox().
            const Box& S_cc_box = S_cc_fab.box();
            const Box& scal_box = scal_fab.box();
            const Box& rho_odot_box = rho_odot_fab.box();
            const Box& rho_Hnuc_box = rho_Hnuc_fab.box();
            const Box& rho_Hext_box = rho_Hext_fab.box();
            const Box& thermal_box = thermal_fab.box();

            // We can now pass the information to a Fortran routine,
            // e.g. S_cc_fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "S_cc_box". We will also pass "box",
            // which specifies our "work" region .
            make_S_cc(ARLIM_3D(validLo), ARLIM_3D(validHi),
                      S_cc_fab.dataPtr(),
                      ARLIM_3D(S_cc_box.loVect()), ARLIM_3D(S_cc_box.hiVect()),
                      scal_fab.dataPtr(),
                      ARLIM_3D(scal_box.loVect()), ARLIM_3D(scal_box.hiVect()),
                      scal_fab.nComp(),
                      rho_odot_fab.dataPtr(),
                      ARLIM_3D(rho_odot_box.loVect()), ARLIM_3D(rho_odot_box.hiVect()),
                      rho_odot_fab.nComp(),
                      rho_Hnuc_fab.dataPtr(),
                      ARLIM_3D(rho_Hnuc_box.loVect()), ARLIM_3D(rho_Hnuc_box.hiVect()),
                      rho_Hext_fab.dataPtr(),
                      ARLIM_3D(rho_Hext_box.loVect()), ARLIM_3D(rho_Hext_box.hiVect()),
                      thermal_fab.dataPtr(),
                      ARLIM_3D(thermal_box.loVect()), ARLIM_3D(thermal_box.hiVect()));
        }

    }

    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,S_cc); // average lev+1 down to lev
    }

}


void
Maestro::Make_NodalRHS (const Vector<MultiFab>& S_cc,
                        Vector<MultiFab>& nodalrhs,
                        const Vector<Real>& Sbar,
                        const Vector<Real>& div_coeff)
{

    Vector<MultiFab> ccrhs(finest_level+1);

    for (int lev=0; lev<=finest_level; ++lev) {
        // build ccrhs
        ccrhs[lev].define(grids[lev], dmap[lev], 1, 1);

        // fill ccrhs
        // get references to the MultiFabs at level lev
        MultiFab& ccrhs_mf = ccrhs[lev];
        const MultiFab& S_cc_mf = S_cc[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) {
            // get references to the FABs, each containing data and the valid+ghost box
            FArrayBox& ccrhs_fab = ccrhs_mf[mfi];
            const FArrayBox& S_cc_fab = S_cc_mf[mfi];

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const int* validLo = validBox.loVect();
            const int* validHi = validBox.hiVect();

            // Get the index space of the valid+ghost region for each FAB
            // Note each of these boxes may contain ghost cells, and thus are
            // larger than or equal to mfi.validbox().
            const Box& ccrhs_box = ccrhs_fab.box();
            const Box& S_cc_box = S_cc_fab.box();

            // We can now pass the information to a Fortran routine,
            // e.g. S_cc_fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "S_cc_box". We will also pass "box",
            // which specifies our "work" region .
            make_ccrhs(lev, ARLIM_3D(validLo), ARLIM_3D(validHi),
                       ccrhs_fab.dataPtr(),
                       ARLIM_3D(ccrhs_box.loVect()), ARLIM_3D(ccrhs_box.hiVect()),
                       S_cc_fab.dataPtr(),
                       ARLIM_3D(S_cc_box.loVect()), ARLIM_3D(S_cc_box.hiVect()),
                       Sbar.dataPtr(), div_coeff.dataPtr());
        }
    }

    // fill ghost cells using first-order extrapolation
    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, t_old, ccrhs[lev], ccrhs, ccrhs, 0, 1, bcs_f);
    }

    // make_nodalrhs
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& nodalrhs_mf = nodalrhs[lev];
        const MultiFab& ccrhs_mf = ccrhs[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(ccrhs_mf); mfi.isValid(); ++mfi ) {
            // get references to the FABs, each containing data and the valid+ghost box
            FArrayBox& nodalrhs_fab = nodalrhs_mf[mfi];
            const FArrayBox& ccrhs_fab = ccrhs_mf[mfi];

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const int* validLo = validBox.loVect();
            const int* validHi = validBox.hiVect();

            // Get the index space of the valid+ghost region for each FAB
            // Note each of these boxes may contain ghost cells, and thus are
            // larger than or equal to mfi.validbox().
            const Box& nodalrhs_box = nodalrhs_fab.box();
            const Box& ccrhs_box = ccrhs_fab.box();

            // We can now pass the information to a Fortran routine,
            // e.g. ccrhs_fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "ccrhs_box". We will also pass "box",
            // which specifies our "work" region .
            make_nodalrhs(lev, ARLIM_3D(validLo), ARLIM_3D(validHi),
                          nodalrhs_fab.dataPtr(),
                          ARLIM_3D(nodalrhs_box.loVect()), ARLIM_3D(nodalrhs_box.hiVect()),
                          ccrhs_fab.dataPtr(),
                          ARLIM_3D(ccrhs_box.loVect()), ARLIM_3D(ccrhs_box.hiVect()));
        }
    }
}
