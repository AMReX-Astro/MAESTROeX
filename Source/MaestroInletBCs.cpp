#include <Maestro.H>
#include <Maestro_F.H>

#include <MaestroInletBCs.H>

using namespace amrex;

/*
inlet_bc_module is a simple container module that holds the parameters
that are used by physbc to implement the inlet boundary conditions.
As these are problem-specific, any problem needing inlet boundary
conditions should create its own version of this module, using this
outline.
*/

// Inlet BCs defined in header file should go here

void Maestro::SetInletBCs() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SetInletBCs()", SetInletBCs);

    // here we would initialize the parameters that are module variables.
    // this routine is called when the base state is defined initially,
    // and upon restart, just after the base state is read in.
}