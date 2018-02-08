
#include <Maestro.H>

using namespace amrex;

void
Maestro::VelocityAdvance (const Vector<MultiFab>& rhohalf,
                          Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                          const Vector<Real>& w0_force,
                          const Vector<Real>& rho0_nph,
                          const Vector<Real>& grav_cell_nph, 
			  const Vector<MultiFab>& sponge)
{

    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    Vector<std::array< MultiFab, AMREX_SPACEDIM > > uedge(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        AMREX_D_TERM(uedge[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], AMREX_SPACEDIM, 0);,
                     uedge[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], AMREX_SPACEDIM, 0);,
                     uedge[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], AMREX_SPACEDIM, 0););
    }

    //////////////////////////////////
    // Create the velocity forcing term at time n using rho 
    //////////////////////////////////

    MakeVelForce(vel_force,umac,sold,rho0_old,grav_cell_old,w0_force,1);
    
    //////////////////////////////////
    // Add w0 to MAC velocities
    //////////////////////////////////

    Addw0(umac,1.);

    //////////////////////////////////
    // Create the edge states of velocity using the MAC velocity plus w0 on edges.
    //////////////////////////////////

    MakeEdgeScal(uold,uedge,umac,vel_force,1,bcs_u,AMREX_SPACEDIM,0,0,AMREX_SPACEDIM,0);
    
    //////////////////////////////////
    // Subtract w0 from MAC velocities.
    //////////////////////////////////

    Addw0(umac,-1.);

    //////////////////////////////////
    // Now create the force at half-time using rhohalf 
    //////////////////////////////////

    MakeVelForce(vel_force,umac,rhohalf,rho0_nph,grav_cell_nph,w0_force,1);

    //////////////////////////////////
    // Update the velocity with convective differencing
    //////////////////////////////////
    /*
    call update_velocity(uold,unew,umac,uedge,force,w0,w0mac, &
                         dx,dt,sponge,mla,the_bc_level)
    */
    UpdateVel(umac, uedge, vel_force, sponge);

}
