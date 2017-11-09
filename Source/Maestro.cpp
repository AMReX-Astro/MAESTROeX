
#include <Maestro.H>
#include <maestro_defaults.H>

using namespace amrex;

int Maestro::Rho       = -1;
int Maestro::RhoH      = -1;
int Maestro::FirstSpec = -1;
int Maestro::NumSpec   = -1;
int Maestro::Temp      = -1;
int Maestro::Pi        = -1;
int Maestro::NSCAL     = -1;

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
