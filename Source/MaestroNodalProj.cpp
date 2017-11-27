
#include <Maestro.H>

using namespace amrex;


// perform a nodal projection
void
Maestro::NodalProj (int proj_type,
                    Vector<MultiFab>& rhohalf)
{
    // modify unew depending on proj_type
    CreateUVecForProjection(proj_type,rhohalf);

    // create a unew with a filled ghost cell
    Vector<MultiFab> unew_ghost(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        unew_ghost[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        FillPatch(lev, t_new, unew_ghost[lev], unew, unew, 0, AMREX_SPACEDIM, bcs_u);
    }



}


// modify unew depending on proj_type
// initial_projection_comp: leave unew alone
// divu_iters_comp:         leave unew alone
// pressure_iters_comp:     unew = (unew-uold)/dt
// regular_timestep_comp:   unew = unew + dt*gpi/rhohalf
void
Maestro::CreateUVecForProjection (int proj_type,
                                  Vector<MultiFab>& rhohalf) {

    if (proj_type == initial_projection_comp) {
        // leave unew alone
    }
    else if (proj_type == divu_iters_comp) {
        // leave unew alone
    }
    else if (proj_type == pressure_iters_comp) {
        // unew = (unew-uold)/dt
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Saxpy(unew[lev],-1.0,uold[lev],0,0,AMREX_SPACEDIM,0);
            unew[lev].mult(1/dt);
        }

    }
    else if (proj_type == regular_timestep_comp) {
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                // gpi = gpi/rhohalf
                MultiFab::Divide(gpi[lev],rhohalf[lev],0,dir,1,0);
            }
            // unew = unew + dt*gpi/rhohalf
            MultiFab::Saxpy(unew[lev],dt,gpi[lev],0,0,AMREX_SPACEDIM,0);
            for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                // revert gpi
                MultiFab::Multiply(gpi[lev],rhohalf[lev],0,dir,1,0);
            }
        }
    }
    else {
        amrex::Abort("MaestroNodalProj: invalid proj_type");
    }

}
