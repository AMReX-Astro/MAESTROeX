
#include <Maestro.H>

using namespace amrex;


// Perform a nodal projection.
// Given a cell-centered velocity field Vproj (assembled in CreateUvecForProj),
// Vproj can decomposed into Vproj = Utilde + (1/sig) grad phi,
// where Utilde satisfies the constraint div(beta0*Utilde) = beta0*(S-Sbar)
// Depending on proj_type we use different values of Vproj, sig, and beta0
// to solve for phi we use div (beta0/sig) grad phi = div(beta*Vproj) - beta0*(S-Sbar)
// then solve for Utilde based on proj_type.

// Since unew^* = unew + grad(phi) (where unew satisfies the constraint)
// we can apply div beta to this equation to get:
// div(beta0 unew^*) = div(beta0 unew) + div(beta0 grad phi)
// OR div(beta0 grad phi) = div(beta0 unew^*) - beta0*S
// We do a density weighted projection, so
// div((beta0/rho) grad phi) = div(beta0 unew*) - beta0*S
// solve for phi, then set unew = unew^* - grad(phi)/rho
// Then depending on proj_type modify unew, pi, and grad(pi)
void
Maestro::NodalProj (int proj_type,
                    Vector<MultiFab>& rhcc)
{

    // build a multifab "sig" div (1/sig) grad phi = RHS with 1 ghost cell
    Vector<MultiFab> sig(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        sig[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    // fill in sig
    // initial_projection_comp: 1
    // divu_iters_comp:         1
    // pressure_iters_comp:     rho^1/2   -- (rhoold+rhonew)/2
    // regular_timestep_comp:   rho^n+1/2 -- (rhoold+rhonew)/2
    // fixme



    // build a multifab Vproj (the quantity that will be projected) with 1 ghost cell
    Vector<MultiFab> Vproj(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        Vproj[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    // fill in Vproj
    // initial_projection_comp: Utilde^0                        -- uold
    // divu_iters_comp:         Utilde^0                        -- uold
    // pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
    // regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
    CreateUvecForProj(proj_type,Vproj);





    // build a multifab to store a Cartesian version of beta0 with 1 ghost cell
    Vector<MultiFab> beta0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        beta0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }


    // convert beta0 to multi-D MultiFab
    // initial_projection_comp: beta0_old
    // divu_iters_comp:         beta0_old
    // pressure_iters_comp:     (beta0_old+beta0_new)/2
    // regular_timestep_comp:   (beta0_old+beta0_new)/2
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        Put1dArrayOnCart(beta0_old,beta0_cart,bcs_f,0,0);
    }
    else {
        // fixme
    }

    // convert Vproj to beta0*Vproj in valid region
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            MultiFab::Multiply(Vproj[lev],beta0_cart[lev],0,dir,1,1);
        }
    }

    // divide sig by beta0
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Divide(sig[lev],beta0_cart[lev],0,0,1,1);
    }

    // create a unew with a filled ghost cell
    Vector<MultiFab> unew_ghost(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        unew_ghost[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        FillPatch(lev, t_new, unew_ghost[lev], unew, unew, 0, 0, AMREX_SPACEDIM, bcs_u);
    }


    // Assemble the nodal RHS as the sum of the cell-centered RHS averaged to nodes 
    // plus div (beta0*Vproj) on nodes
    // fix up the boundary condition ghost cells for Vproj
    // solve for phi


    // make grad phi



    // convert sig/beta0 back to sig
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Multiply(sig[lev],beta0_cart[lev],0,0,1,1);
    }

    // update velocity
    // initial_projection_comp: Utilde^0   = V - (1/sig)*grad(phi)
    // divu_iters_comp:         Utilde^0   = V - (1/sig)*grad(phi)
    // pressure_iters_comp:     Utilde^n+1 = Utilde^n + dt(V-(1/sig)*grad(phi))
    // regular_timestep_comp:   Utilde^n+1 = V - (1/sig)*grad(phi)




}

// fill in Vproj
// initial_projection_comp: Utilde^0                        -- uold
// divu_iters_comp:         Utilde^0                        -- uold
// pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
// regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
// sig contains rhohalf if proj_type == regular_timestep_comp
void
Maestro::CreateUvecForProj (int proj_type,
                            Vector<MultiFab>& sig) {

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
                MultiFab::Divide(gpi[lev],sig[lev],0,dir,1,0);
            }
            // unew = unew + dt*gpi/rhohalf
            MultiFab::Saxpy(unew[lev],dt,gpi[lev],0,0,AMREX_SPACEDIM,0);
            for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                // revert gpi/rhohalf back to gpi
                MultiFab::Multiply(gpi[lev],sig[lev],0,dir,1,0);
            }
        }
    }
    else {
        amrex::Abort("MaestroNodalProj: invalid proj_type");
    }

    Real time;
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        time = t_old;
    }
    else  {
        time = t_new;
    }

    // replace with a fillpatch so we get fine ghost cells overlying coarse filled
/*
    for (int lev=0; lev<=finest_level; ++lev) {
        Vproj[lev].FillBoundary(0,AMREX_SPACEDIM,geom[lev].periodicity());
        bcs_u.FillBoundary(Vproj[lev],0,AMREX_SPACEDIM,time);
    }
*/


}
