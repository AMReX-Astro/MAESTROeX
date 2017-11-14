
#include <Maestro.H>

using namespace amrex;

void
Maestro::Make_S_cc (Vector<MultiFab>& S_cc,
                    Vector<MultiFab>& scal,
                    Vector<MultiFab>& rho_omegadot,
                    Vector<MultiFab>& rho_Hnuc,
                    Vector<MultiFab>& rho_Hext,
                    Vector<MultiFab>& thermal)
{

    for (int lev=0; lev<=finest_level; ++lev) 
    {
        MultiFab& S_cc_mf = S_cc[lev];
        MultiFab& scal_mf = scal[lev];
        MultiFab& rho_odot_mf = rho_omegadot[lev];
        MultiFab& rho_Hnuc_mf = rho_Hnuc[lev];
        MultiFab& rho_Hext_mf = rho_Hext[lev];
        MultiFab& thermal_mf = thermal[lev];

        for ( MFIter mfi(S_cc_mf); mfi.isValid(); ++mfi ) // Loop over boxes
        {
            // get references to the FABs, each containing data and the valid+ghost box
            FArrayBox& S_cc_fab = S_cc_mf[mfi];
            FArrayBox& scal_fab = scal_mf[mfi];
            FArrayBox& rho_odot_fab = rho_odot_mf[mfi];
            FArrayBox& rho_Hnuc_fab = rho_Hnuc_mf[mfi];
            FArrayBox& rho_Hext_fab = rho_Hext_mf[mfi];
            FArrayBox& thermal_fab = thermal_mf[mfi];

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // Get the index space of the valid+ghost region for each FAB
            // Note each of these boxes may contain ghost cells, and thus are
            // larger than or equal to mfi.validbox().
            const Box & S_cc_box = S_cc_fab.box();
            const Box & scal_box = scal_fab.box();
            const Box & rho_odot_box = rho_odot_fab.box();
            const Box & rho_Hnuc_box = rho_Hnuc_fab.box();
            const Box & rho_Hext_box = rho_Hext_fab.box();
            const Box & thermal_box = thermal_fab.box();

            // We can now pass the information to a Fortran routine,
            // e.g. S_cc_fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "S_cc_box". We will also pass "box",
            // which specifies our "work" region .
            make_S_cc(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                      S_cc_fab.dataPtr(),
                      ARLIM_3D(S_cc_box.loVect()), ARLIM_3D(S_cc_box.hiVect()),
                      S_cc_fab.nComp(),
                      scal_fab.dataPtr(),
                      ARLIM_3D(scal_box.loVect()), ARLIM_3D(scal_box.hiVect()),
                      scal_fab.nComp(),
                      rho_odot_fab.dataPtr(),
                      ARLIM_3D(rho_odot_box.loVect()), ARLIM_3D(rho_odot_box.hiVect()),
                      rho_odot_fab.nComp(),
                      rho_Hnuc_fab.dataPtr(),
                      ARLIM_3D(rho_Hnuc_box.loVect()), ARLIM_3D(rho_Hnuc_box.hiVect()),
                      rho_Hnuc_fab.nComp(),
                      rho_Hext_fab.dataPtr(),
                      ARLIM_3D(rho_Hext_box.loVect()), ARLIM_3D(rho_Hext_box.hiVect()),
                      rho_Hext_fab.nComp(),
                      thermal_fab.dataPtr(),
                      ARLIM_3D(thermal_box.loVect()), ARLIM_3D(thermal_box.hiVect()),
                      thermal_fab.nComp());
        }

    }

    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,S_cc); // average lev+1 down to lev
    }

}

/*
void
Maestro::Make_S_nodal (Vector<<MultiFab>& S_cc,
                       Vector<MultiFab>& S_nodal,
                       Real Sbar,
                       const Vector<Real>& div_coeff)
{}
*/

void
Maestro::Make_S_nodal ()
{
}
