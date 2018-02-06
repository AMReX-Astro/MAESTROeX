
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
    FillPatch(t_old, vel_force, vel_force, vel_force, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

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
    FillPatch(t_old, scal_force, scal_force, scal_force, comp, comp, 1, 0, bcs_f);

}

void
Maestro::MakeRhoHForce(Vector<MultiFab>& scal_force,
                       int is_prediction,
                       const Vector<MultiFab>& thermal,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       int add_thermal)

{
    // if we are doing the prediction, then it only makes sense to be in
    // this routine if the quantity we are predicting is rhoh', h, or rhoh
    if (is_prediction == 1 && !(enthalpy_pred_type == predict_rhohprime ||
                                enthalpy_pred_type == predict_h ||
                                enthalpy_pred_type == predict_rhoh) ) {
        Abort("ERROR: should only call mkrhohforce when predicting rhoh', h, or rhoh");
    }

    Vector<Real> rho0( (max_radial_level+1)*nr_fine );
    Vector<Real>   p0( (max_radial_level+1)*nr_fine );
    Vector<Real> grav( (max_radial_level+1)*nr_fine );
    rho0.shrink_to_fit();
      p0.shrink_to_fit();
    grav.shrink_to_fit();

    if (is_prediction == 1) {
        rho0 = rho0_old;
          p0 =   p0_old;
    }
    else {
        for(int i=0; i<rho0.size(); ++i) {
            rho0[i] = 0.5*(rho0_old[i]+rho0_new[i]);
              p0[i] = 0.5*(  p0_old[i]+  p0_new[i]);
        }
    }

    if (spherical == 1) {
    }

    make_grav_cell(grav.dataPtr(),
                   rho0.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());


    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& scal_force_mf = scal_force[lev];
#if (AMREX_SPACEDIM == 1)
        const MultiFab& umac_mf = umac[lev][0];
#elif (AMREX_SPACEDIM == 2)
        const MultiFab& vmac_mf = umac[lev][1];
#elif (AMREX_SPACEDIM == 3)
        const MultiFab& wmac_mf = umac[lev][2];
#endif
        const MultiFab& thermal_mf = thermal[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_force_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            mkrhohforce(&lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                        scal_force_mf[mfi].dataPtr(RhoH), 
                        ARLIM_3D(scal_force_mf[mfi].loVect()), ARLIM_3D(scal_force_mf[mfi].hiVect()),
#if (AMREX_SPACEDIM == 1)
                        BL_TO_FORTRAN_3D(umac_mf[mfi]),
#elif (AMREX_SPACEDIM == 2)
                        BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#elif (AMREX_SPACEDIM == 3)
                        BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                        BL_TO_FORTRAN_3D(thermal_mf[mfi]),
                        p0.dataPtr(), rho0.dataPtr(), grav.dataPtr(), psi.dataPtr(),
                        &is_prediction, &add_thermal);
        }
    }

    // average down and fill ghost cells
    AverageDown(scal_force,RhoH,1);
    FillPatch(t_old,scal_force,scal_force,scal_force,RhoH,RhoH,1,0,bcs_f);
}
