
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
        FillPatch(lev, t_new, utilde[lev], uold, uold, 0, 0, AMREX_SPACEDIM, bcs_u);
    }

    // create a MultiFab to hold uold + w0
    Vector<MultiFab>      ufull(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
    }

    // create ufull = uold + w0
    Put1dArrayOnCart(w0,ufull,bcs_u,1,1);
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

    // VelPred();


}

void
Maestro::MakeUtrans (const Vector<MultiFab>& utilde,
                     const Vector<MultiFab>& ufull,
                     Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& sold_mf    = sold[lev]; // only used for the MFIter since it's cell-centered
        const MultiFab& utilde_mf  = utilde[lev];
        const MultiFab& ufull_mf   = ufull[lev];
        MultiFab& utrans_mf        = utrans[lev][0];
#if (AMREX_SPACEDIM >= 2)
        MultiFab& vtrans_mf        = utrans[lev][1];
#if (AMREX_SPACEDIM == 3)
        MultiFab& wtrans_mf        = utrans[lev][2];
#endif
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(utilde_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();
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
                        lev, validBox.loVect(), validBox.hiVect(),
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
