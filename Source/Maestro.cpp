
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// default values for class data members listed in _cpp_parameters
// that are declared in maestro_params.H
#include <maestro_defaults.H>

// default values are overwritten in VariableSetup()
int Maestro::Rho       = -1;
int Maestro::RhoH      = -1;
int Maestro::FirstSpec = -1;
int Maestro::NumSpec   = -1;
int Maestro::Temp      = -1;
int Maestro::Pi        = -1;
int Maestro::Nscal     = -1;

// number of ghost cells for sold/new and uold/new
// overwritten in VariableSetup()
int Maestro::ng_s = -1;

int Maestro::initial_projection_comp = 1;
int Maestro::divu_iters_comp         = 2;
int Maestro::pressure_iters_comp     = 3;
int Maestro::regular_timestep_comp   = 4;

// species prediction
int Maestro::predict_rhoprime_and_X    = 1;
int Maestro::predict_rhoX              = 2;
int Maestro::predict_rho_and_X         = 3;
// enthalpy prediction
int Maestro::predict_rhoh             = 0;
int Maestro::predict_rhohprime        = 1;
int Maestro::predict_h                = 2;
int Maestro::predict_T_then_rhohprime = 3;
int Maestro::predict_T_then_h         = 4;
int Maestro::predict_hprime           = 5;
int Maestro::predict_Tprime_then_h    = 6;

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

#ifdef AMREX_USE_CUDA
int Maestro::numBCThreadsMin[3] = {1, 1, 1};
#endif

// this will be reset upon restart
Real Maestro::previousCPUTimeUsed = 0.0;

Real Maestro::startCPUTime = 0.0;

Maestro::Maestro ()
{
}

Maestro::~Maestro ()
{
}

Real
Maestro::getCPUTime()
{

		int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
		numCores = numCores*omp_get_max_threads();
#endif

		Real T = numCores*(ParallelDescriptor::second() - startCPUTime) +
		         previousCPUTimeUsed;

		return T;
}
