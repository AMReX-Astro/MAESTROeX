
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

//////////////////////////////////
// solver tolerances

// tolerances for the initial projection
Real Maestro::eps_init_proj_cart = 1.e-12;
Real Maestro::eps_init_proj_sph  = 1.e-10;
// tolerances for the divu iterations
Real Maestro::eps_divu_cart = 1.e-12;
Real Maestro::eps_divu_sph  = 1.e-10;
Real Maestro::divu_iter_factor = 100.;
Real Maestro::divu_level_factor = 10.;
// tolerances for the MAC projection
Real Maestro::eps_mac = 1.e-10;
Real Maestro::eps_mac_max = 1.e-8;
Real Maestro::mac_level_factor = 10.;
Real Maestro::eps_mac_bottom = 1.e-3;
// tolerances for the nodal projection
Real Maestro::eps_hg = 1.e-12;
Real Maestro::eps_hg_max = 1.e-10;
Real Maestro::hg_level_factor = 10.;
Real Maestro::eps_hg_bottom = 1.e-4;

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
