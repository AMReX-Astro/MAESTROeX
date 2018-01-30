
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeVelForce (Vector<MultiFab>& vel_force,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                       const Vector<MultiFab>& rho,
                       const Vector<Real>& rho0,
                       const Vector<Real>& grav_cell,
                       const Vector<Real>& w0_force,
                       int do_add_utilde_force)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& vel_force_mf = vel_force[lev];
        const MultiFab& gpi_mf = gpi[lev];
        const MultiFab& rho_mf = rho[lev];
        const MultiFab& uedge_mf = uedge[lev][0];
        const MultiFab& vedge_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wedge_mf = uedge[lev][2];
#endif

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(vel_force_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_vel_force(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                           BL_TO_FORTRAN_FAB(gpi_mf[mfi]),
                           BL_TO_FORTRAN_N_3D(rho_mf[mfi],Rho),
                           BL_TO_FORTRAN_3D(uedge_mf[mfi]),
                           BL_TO_FORTRAN_3D(vedge_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                           BL_TO_FORTRAN_3D(wedge_mf[mfi]),
#endif
                           w0.dataPtr(),
                           w0_force.dataPtr(),
                           rho0.dataPtr(),
                           grav_cell.dataPtr(),
                           &do_add_utilde_force);
        }
    }

    // average fine data onto coarser cells
    AverageDown(vel_force,0,AMREX_SPACEDIM);

    // note - we need to reconsider the bcs type here
    // it matches fortran MAESTRO but is that correct?
    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, t_old, vel_force[lev], vel_force, vel_force, 0, 0, AMREX_SPACEDIM, 0, bcs_u);
    }


}


void
Maestro::ModifyScalForce(Vector<MultiFab>& scal_force,
                         const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                         const Vector<Real>& s0,
                         const Vector<Real>& s0_edge,
                         int comp,
                         const Vector<BCRec>& bcs,
                         int fullform)
{

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_force_mf = scal_force[lev];
        const MultiFab& sold_mf = sold[lev];
        const MultiFab& umac_mf = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
        const MultiFab& vmac_mf = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wmac_mf = umac[lev][2];
#endif
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_force_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Real* dx = geom[lev].CellSize();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.

            // be careful to pass in "comp+1" here since fortran uses 1-based indexing for components
            modify_scal_force(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                              scal_force_mf[mfi].dataPtr(comp), 
                              ARLIM_3D(scal_force_mf[mfi].loVect()), ARLIM_3D(scal_force_mf[mfi].hiVect()),
                              sold_mf[mfi].dataPtr(comp), 
                              ARLIM_3D(sold_mf[mfi].loVect()), ARLIM_3D(sold_mf[mfi].hiVect()),
                              BL_TO_FORTRAN_3D(umac_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
                              BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                              BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
#endif
                              s0.dataPtr(), s0_edge.dataPtr(), w0.dataPtr(),
                              dx, &fullform);
        }


    }


    // average fine data onto coarser cells
    AverageDown(scal_force,comp,1);

    // fill ghost cells
    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, t_old, scal_force[lev], scal_force, scal_force, comp, comp, 1, 0, bcs_f);
    }

}
