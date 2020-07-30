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

Real InletBCs::INLET_RHO;
Real InletBCs::INLET_RHOH;
Real InletBCs::INLET_TEMP;
Real InletBCs::INLET_CS;

void Maestro::SetInletBCs() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SetInletBCs()", SetInletBCs);

    // here we initialize the parameters that are module variables.
    // this routine is called when the base state is defined initially,
    // and upon restart, just after the base state is read in.

    eos_t eos_state;

    eos_state.T = 10.0;
    eos_state.rho = 1.e-3;
    eos_state.p = 1.e6;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = 1.0 / NumSpec;
    }

    eos(eos_input_rp, eos_state);

    InletBCs::INLET_CS = eos_state.cs;

    eos_state.T = 10.e0;
    eos_state.rho = 5.e-4;
    eos_state.p = 1.e6;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = 1.0 / NumSpec;
    }

    eos(eos_input_rp, eos_state);

    InletBCs::INLET_RHO = eos_state.rho;
    InletBCs::INLET_RHOH = eos_state.rho * eos_state.h;
    InletBCs::INLET_TEMP = eos_state.T;
}
