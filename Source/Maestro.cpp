
#include <Maestro.H>

using namespace amrex;

// default values for class data members listed in _cpp_parameters
// that are declared in maestro_params.H
#include <maestro_defaults.H>

// default values are overridden in VariableSetup()
int Maestro::Rho       = -1;
int Maestro::RhoH      = -1;
int Maestro::FirstSpec = -1;
int Maestro::NumSpec   = -1;
int Maestro::Temp      = -1;
int Maestro::Pi        = -1;
int Maestro::Nscal     = -1;

int Maestro::initial_projection_comp = 1;
int Maestro::divu_iters_comp         = 2;
int Maestro::pressure_iters_comp     = 3;
int Maestro::regular_timestep_comp   = 4;

// helper IntVects used to define face/nodal MultiFabs
#if (AMREX_SPACEDIM == 2)
IntVect Maestro::nodal_flag(1,1);
IntVect Maestro::nodal_flag_x(1,0);
IntVect Maestro::nodal_flag_y(0,1);
#elif (AMREX_SPACEDIM == 3)
IntVect Maestro::nodal_flag(1,1,1);
IntVect Maestro::nodal_flag_x(1,0,0);
IntVect Maestro::nodal_flag_y(0,1,0);
IntVect Maestro::nodal_flag_z(0,0,1);
#endif

Maestro::Maestro ()
{}

Maestro::~Maestro ()
{}
