
#include <Maestro.H>

using namespace amrex;


// perform a nodal projection
// given a velocity field unew^* (unew is assembled in CreateUvecForProj) 
// We want to project the velocity to satisfy div(beta_0*unew) = beta_0*S
// Since unew^* = unew + grad(phi) (where unew satisfies the constraint)
// we can apply div beta to this equation to get:
// div(beta_0 unew^*) = div(beta_0 unew) + div(beta_0 grad phi)
// OR div(beta_0 grad phi) = div(beta_0 unew^*) - beta_0*S
// We do a density weighted projection, so
// div((beta_0/rho) grad phi) = div(beta_0 unew*) - beta_0*S
// solve for phi, then set unew = unew^* - grad(phi)/rho
// Then depending on proj_type modify unew, pi, and grad(pi)
void
Maestro::NodalProj (int proj_type,
                    Vector<MultiFab>& nodalrhs,
                    Vector<MultiFab>& rhohalf)
{
    // modify unew to create uvec depending on proj_type
    // initial_projection_comp: (leave unew alone)
    // divu_iters_comp:         (leave unew alone)
    // pressure_iters_comp:     unew = (unew-uold)/dt
    // regular_timestep_comp:   unew = unew + dt*gpi/rhohalf
    CreateUvecForProj(proj_type,rhohalf);

    // convert beta_0 to multi-D MultiFab with 1 ghost cell


    // convert unew to beta_0*unew in valid region
    

    // convert rhohalf to rhohalf/beta_0 in valid+ghost region


    // create a unew with a filled ghost cell
    Vector<MultiFab> unew_ghost(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        unew_ghost[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        FillPatch(lev, t_new, unew_ghost[lev], unew, unew, 0, 0, AMREX_SPACEDIM, bcs_u);
    }


    // create RHS = div((beta_0/rho)*unew) - (beta_0/rho)*S



    // solve for phi




    // convert beta_0*unew back to unew



    // convert (rhohalf/beta_0) back to rhohalf


    // make grad phi



    // update velocity


}


// modify unew to create uvec depending on proj_type
// initial_projection_comp: (leave unew alone)
// divu_iters_comp:         (leave unew alone)
// pressure_iters_comp:     unew = (unew-uold)/dt
// regular_timestep_comp:   unew = unew + dt*gpi/rhohalf
void
Maestro::CreateUvecForProj (int proj_type,
                            Vector<MultiFab>& rhohalf) {

    if (proj_type == initial_projection_comp) {
        // unew = unew (leave unew alone)
    }
    else if (proj_type == divu_iters_comp) {
        // unew = unew (leave unew alone)
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
                // convert gpi to gpi/rhohalf
                MultiFab::Divide(gpi[lev],rhohalf[lev],0,dir,1,0);
            }
            // unew = unew + dt*gpi/rhohalf
            MultiFab::Saxpy(unew[lev],dt,gpi[lev],0,0,AMREX_SPACEDIM,0);
            for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                // revert gpi/rhohalf back to gpi
                MultiFab::Multiply(gpi[lev],rhohalf[lev],0,dir,1,0);
            }
        }
    }
    else {
        amrex::Abort("MaestroNodalProj: invalid proj_type");
    }

}
