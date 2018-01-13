
#include <Maestro.H>

using namespace amrex;

// compute unprojected mac velocities
void
Maestro::AdvancePremac (Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac)
{

    // create a uold with filled ghost cells
    Vector<MultiFab> uold_ghost(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        uold_ghost[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
        FillPatch(lev, t_new, uold_ghost[lev], uold, uold, 0, 0, AMREX_SPACEDIM, bcs_u);
    }

    // create a MultiFab to hold uold + w0
    Vector<MultiFab>      ufull(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        ufull[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_adv);
    }

    // create ufull = uold + w0
    Put1dArrayOnCart(w0,ufull,bcs_u,1,1);
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Add(ufull[lev],uold_ghost[lev],0,0,AMREX_SPACEDIM,ng_adv);
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

    // MakeUtrans();
    FillUmacGhost(utrans);

    
    // create a MultiFab to hold the velocity forcing
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    int is_final_update = 0;
    int do_add_utilde_force = 1;
//    MakeVelForce(vel_force,utran,w0_force,is_final_update,do_add_utilde_force);

    // add w0 to trans velocities

    // VelPred();


}

void
Maestro::MakeUtrans (const Vector<MultiFab>& utilde,
                     const Vector<MultiFab>& ufill,
                     Vector<std::array< MultiFab, AMREX_SPACEDIM > >& utrans)
{



    // fill ghost cells

}
