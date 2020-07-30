
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::VelocityAdvance(
    const Vector<MultiFab>& rhohalf,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const Vector<MultiFab>& w0_force_cart, const BaseState<Real>& rho0_nph,
    const BaseState<Real>& grav_cell_nph, const Vector<MultiFab>& sponge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelocityAdvance()", VelocityAdvance);

    Vector<MultiFab> vel_force(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (ppm_trace_forces == 0) {
            vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        } else {
            // tracing needs more ghost cells
            vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
        }
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    Vector<std::array<MultiFab, AMREX_SPACEDIM> > uedge(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_D_TERM(uedge[lev][0].define(convert(grids[lev], nodal_flag_x),
                                          dmap[lev], AMREX_SPACEDIM, 0);
                     , uedge[lev][1].define(convert(grids[lev], nodal_flag_y),
                                            dmap[lev], AMREX_SPACEDIM, 0);
                     , uedge[lev][2].define(convert(grids[lev], nodal_flag_z),
                                            dmap[lev], AMREX_SPACEDIM, 0););
    }

    //////////////////////////////////
    // Create the velocity forcing term at time n using rho
    //////////////////////////////////

    MakeVelForce(vel_force, umac, sold, rho0_old, grav_cell_old, w0_force_cart,
#ifdef ROTATION
                 w0mac, false,
#endif
                 1);

    //////////////////////////////////
    // Add w0 to MAC velocities
    //////////////////////////////////

    Addw0(umac, w0mac, 1.);

    //////////////////////////////////
    // Create the edge states of velocity using the MAC velocity plus w0 on edges.
    //////////////////////////////////

    MakeEdgeScal(uold, uedge, umac, vel_force, true, bcs_u, AMREX_SPACEDIM, 0,
                 0, AMREX_SPACEDIM, false);

    //////////////////////////////////
    // Subtract w0 from MAC velocities.
    //////////////////////////////////

    Addw0(umac, w0mac, -1.);

    //////////////////////////////////
    // Now create the force at half-time using rhohalf
    //////////////////////////////////

    MakeVelForce(vel_force, umac, rhohalf, rho0_nph, grav_cell_nph,
                 w0_force_cart,
#ifdef ROTATION
                 w0mac, true,
#endif
                 1);

    //////////////////////////////////
    // Update the velocity with convective differencing
    //////////////////////////////////

    UpdateVel(umac, uedge, vel_force, sponge, w0mac);
}
