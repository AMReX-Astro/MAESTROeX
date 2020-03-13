#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

/*
inlet_bc_module is a simple container module that holds the parameters
that are used by physbc to implement the inlet boundary conditions.
As these are problem-specific, any problem needing inlet boundary
conditions should create its own version of this module, using this
outline.
*/

void 
Maestro::SetInletBCs()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SetInletBCs()", SetInletBCs); 

    // here we initialize the parameters that are module variables.
    // this routine is called when the base state is defined initially,
    // and upon restart, just after the base state is read in.

    // figure out the indices for different species
    const auto ic12 = network_spec_index("carbon-12");
    const auto io16 = network_spec_index("oxygen-16");

    if (ic12 < 0 || io16 < 0) {
        Abort("ERROR: species indices undefined in inlet_bc");
    }

    eos_t eos_state;

    eos_state.rho = dens_fuel;
    eos_state.T   = temp_fuel;

    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp]    = 0.0;
    }
    eos_state.xn[ic12] = xc12_fuel;
    eos_state.xn[io16] = 1.e0 - xc12_fuel;

    eos(eos_input_rt, eos_state);

    inlet_bcs.INLET_RHO     = dens_fuel;
    inlet_bcs.INLET_RHOH    = dens_fuel*eos_state.h;
    inlet_bcs.INLET_TEMP    = temp_fuel;
    inlet_bcs.INLET_RHOHX.resize(NumSpec);
    for (auto comp = 0; comp < NumSpec; ++comp) {
        inlet_bcs.INLET_RHOHX[comp] = dens_fuel*eos_state.xn[comp];
    }
    inlet_bcs.INLET_VEL     = vel_fuel;
}