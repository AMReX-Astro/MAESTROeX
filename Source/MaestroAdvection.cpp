
#include <Maestro.H>

using namespace amrex;

// compute unprojected mac velocities
void
Maestro::AdvancePremac (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                        const Vector<Real>& w0_force)
{

    // create a uold with filled ghost cells
    Vector<MultiFab> utilde(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        utilde[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
        FillPatch(lev, t_new, utilde[lev], uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u);
    }

    // create a MultiFab to hold uold + w0
    Vector<MultiFab>      ufull(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
    }

    // create ufull = uold + w0
    Put1dArrayOnCart(w0,ufull,1,1,bcs_u,0);
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Add(ufull[lev],utilde[lev],0,0,AMREX_SPACEDIM,ng_adv);
    }

    // create a face-centered MultiFab to hold utrans
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > utrans(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        utrans[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        utrans[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
        utrans[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
    }

    // create utrans
    MakeUtrans(utilde,ufull,utrans);
    
    // create a MultiFab to hold the velocity forcing
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    int do_add_utilde_force = 1;
    MakeVelForce(vel_force,utrans,sold,rho0_old,grav_cell_old,w0_force,do_add_utilde_force);

    // add w0 to trans velocities
    Addw0 (utrans,1.);

    VelPred(utilde,ufull,utrans,umac,vel_force);
}

void
Maestro::MakeUtrans (const Vector<MultiFab>& utilde,
                     const Vector<MultiFab>& ufull,
                     Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& utrans_mf  = utrans[lev][0];
#if (AMREX_SPACEDIM >= 2)
              MultiFab& vtrans_mf  = utrans[lev][1];
#if (AMREX_SPACEDIM == 3)
              MultiFab& wtrans_mf  = utrans[lev][2];
#endif
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();
            const Real* dx = geom[lev].CellSize();

            const FArrayBox& utilde_fab  = utilde_mf[mfi];
            const FArrayBox& ufull_fab   = ufull_mf[mfi];
            FArrayBox& utrans_fab = utrans_mf[mfi];
#if (AMREX_SPACEDIM >= 2)
            FArrayBox& vtrans_fab = vtrans_mf[mfi];
#if (AMREX_SPACEDIM == 3)
            FArrayBox& wtrans_fab = wtrans_mf[mfi];
#endif
#endif

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
            mkutrans_1d(
#elif (AMREX_SPACEDIM == 2)
            mkutrans_2d(
#elif (AMREX_SPACEDIM == 3)
            mkutrans_3d(
#endif
                        lev, domainBox.loVect(), domainBox.hiVect(),
                        validBox.loVect(), validBox.hiVect(),
                        utilde_fab.dataPtr(), utilde_fab.loVect(), utilde_fab.hiVect(), utilde_fab.nCompPtr(),
                        utilde_mf.nGrow(),
                        ufull_fab.dataPtr(), ufull_fab.loVect(), ufull_fab.hiVect(), ufull_fab.nCompPtr(),
                        utrans_fab.dataPtr(), utrans_fab.loVect(), utrans_fab.hiVect(),
#if (AMREX_SPACEDIM >= 2)
                        vtrans_fab.dataPtr(), vtrans_fab.loVect(), vtrans_fab.hiVect(),
#if (AMREX_SPACEDIM == 3)
                        wtrans_fab.dataPtr(), wtrans_fab.loVect(), wtrans_fab.hiVect(),
#endif
#endif
                        w0.dataPtr(), dx, dt, bcs_u[0].data(), phys_bc.dataPtr());

        } // end MFIter loop
    } // end loop over levels
    
    // fill peroidic ghost cells
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            utrans[lev][d].FillBoundary(geom[lev].periodicity());
        }
    }

    // fill ghost cells behind physical boundaries
    FillUmacGhost(utrans);

    // FIXME need to add edge_restriction and create_umac_grown
    //
    //
}

void
Maestro::VelPred (const Vector<MultiFab>& utilde,
                  const Vector<MultiFab>& ufull,
                  const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans,
                  Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                  const Vector<MultiFab>& force)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
              MultiFab& umac_mf    = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
        const MultiFab& utrans_mf  = utrans[lev][0];
        const MultiFab& vtrans_mf  = utrans[lev][1];
              MultiFab& vmac_mf    = umac[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wtrans_mf  = utrans[lev][2];
              MultiFab& wmac_mf    = umac[lev][2];
#endif
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();
            const Real* dx = geom[lev].CellSize();

            const FArrayBox& utilde_fab = utilde_mf[mfi];
            const FArrayBox& ufull_fab  = ufull_mf[mfi];
                  FArrayBox& umac_fab   = umac_mf[mfi];
#if (AMREX_SPACEDIM >= 2)
            const FArrayBox& utrans_fab = utrans_mf[mfi];
            const FArrayBox& vtrans_fab = vtrans_mf[mfi];
                  FArrayBox& vmac_fab   = vmac_mf[mfi];
#if (AMREX_SPACEDIM == 3)
            const FArrayBox& wtrans_fab = wtrans_mf[mfi];
                  FArrayBox& wmac_fab   = wmac_mf[mfi];
#endif
#endif
            const FArrayBox& force_fab = force_mf[mfi];

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
            velpred_1d(
#elif (AMREX_SPACEDIM == 2)
            velpred_2d(
#elif (AMREX_SPACEDIM == 3)
            velpred_3d(
#endif
                        lev, domainBox.loVect(), domainBox.hiVect(),
                        validBox.loVect(), validBox.hiVect(),
                        utilde_fab.dataPtr(), utilde_fab.loVect(), utilde_fab.hiVect(), utilde_fab.nCompPtr(),
                        utilde_mf.nGrow(),
                        ufull_fab.dataPtr(), ufull_fab.loVect(), ufull_fab.hiVect(), ufull_fab.nCompPtr(),
#if (AMREX_SPACEDIM >= 2)
                        utrans_fab.dataPtr(), utrans_fab.loVect(), utrans_fab.hiVect(),
                        vtrans_fab.dataPtr(), vtrans_fab.loVect(), vtrans_fab.hiVect(),
#if (AMREX_SPACEDIM == 3)
                        wtrans_fab.dataPtr(), wtrans_fab.loVect(), wtrans_fab.hiVect(),
#endif
#endif
                        umac_fab.dataPtr(), umac_fab.loVect(), umac_fab.hiVect(),
#if (AMREX_SPACEDIM >= 2)
                        vmac_fab.dataPtr(), vmac_fab.loVect(), vmac_fab.hiVect(),
#if (AMREX_SPACEDIM == 3)
                        wmac_fab.dataPtr(), wmac_fab.loVect(), wmac_fab.hiVect(),
#endif
#endif
                        force_fab.dataPtr(), force_fab.loVect(), force_fab.hiVect(), force_fab.nCompPtr(),
                        w0.dataPtr(), dx, dt, bcs_u[0].data(), phys_bc.dataPtr());





        } // end MFIter loop
    } // end loop over levels

    // FIXME need to add edge_restriction
    //
    //

}

void
Maestro::MakeEdgeScal (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<MultiFab>& force,
                       int is_vel, const Vector<BCRec>& bcs,
                       int nbccomp, int comp, int bccomp, int is_conservative)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& scal_mf   = sold[lev];
              MultiFab& sedgex_mf = sedge[lev][0];
        const MultiFab& umac_mf   = umac[lev][0];
#if (AMREX_SPACEDIM >= 2)
              MultiFab& sedgey_mf = sedge[lev][1];
        const MultiFab& vmac_mf   = umac[lev][1];

#if (AMREX_SPACEDIM == 3)
              MultiFab& sedgez_mf = sedge[lev][2];
        const MultiFab& wmac_mf   = umac[lev][2];

#endif
#endif
        const MultiFab& force_mf = force[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
            const Box& domainBox = geom[lev].Domain();
            const Real* dx = geom[lev].CellSize();

            const FArrayBox& scal_fab   = scal_mf[mfi];
                  FArrayBox& sedgex_fab = sedgex_mf[mfi];
            const FArrayBox& umac_fab   = umac_mf[mfi];
#if (AMREX_SPACEDIM >= 2)
                  FArrayBox& sedgey_fab = sedgey_mf[mfi];
            const FArrayBox& vmac_fab   = vmac_mf[mfi];
#if (AMREX_SPACEDIM == 3)
                  FArrayBox& sedgez_fab = sedgez_mf[mfi];
            const FArrayBox& wmac_fab   = wmac_mf[mfi];
#endif
#endif
            const FArrayBox& force_fab = force_mf[mfi];

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
#if (AMREX_SPACEDIM == 1)
            make_edge_scal_1d(
#elif (AMREX_SPACEDIM == 2)
            make_edge_scal_2d(
#elif (AMREX_SPACEDIM == 3)
            make_edge_scal_3d(
#endif
                        domainBox.loVect(), domainBox.hiVect(),
                        validBox.loVect(), validBox.hiVect(),
                        scal_fab.dataPtr(), scal_fab.loVect(), scal_fab.hiVect(), scal_fab.nCompPtr(), scal_mf.nGrow(),
                        sedgex_fab.dataPtr(), sedgex_fab.loVect(), sedgex_fab.hiVect(), sedgex_fab.nCompPtr(),
#if (AMREX_SPACEDIM >= 2)
                        sedgey_fab.dataPtr(), sedgey_fab.loVect(), sedgey_fab.hiVect(), sedgey_fab.nCompPtr(),
#if (AMREX_SPACEDIM == 3)
                        sedgez_fab.dataPtr(), sedgez_fab.loVect(), sedgez_fab.hiVect(), sedgez_fab.nCompPtr(),
#endif
#endif
                        umac_fab.dataPtr(), umac_fab.loVect(), umac_fab.hiVect(),
#if (AMREX_SPACEDIM >= 2)
                        vmac_fab.dataPtr(), vmac_fab.loVect(), vmac_fab.hiVect(),
#if (AMREX_SPACEDIM == 3)
                        wmac_fab.dataPtr(), wmac_fab.loVect(), wmac_fab.hiVect(),
#endif
#endif
                        force_fab.dataPtr(), force_fab.loVect(), force_fab.hiVect(), force_fab.nCompPtr(),
                        dx, dt, is_vel, bcs_u[0].data(),
                        nbccomp, comp, bccomp, is_conservative);





        } // end MFIter loop
    } // end loop over levels

    // FIXME need to add edge_restriction
    //
    //

}
