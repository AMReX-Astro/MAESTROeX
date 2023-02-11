
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// number of ghost cells for sold/new and uold/new
// overwritten in VariableSetup()
int Maestro::ng_s = -1;

// helper IntVects used to define face/nodal MultiFabs
#if (AMREX_SPACEDIM == 2)
IntVect Maestro::nodal_flag(1, 1);
IntVect Maestro::nodal_flag_x(1, 0);
IntVect Maestro::nodal_flag_y(0, 1);
#elif (AMREX_SPACEDIM == 3)
IntVect Maestro::nodal_flag(1, 1, 1);
IntVect Maestro::nodal_flag_x(1, 0, 0);
IntVect Maestro::nodal_flag_y(0, 1, 0);
IntVect Maestro::nodal_flag_z(0, 0, 1);
#endif

// this will be reset upon restart
Real Maestro::previousCPUTimeUsed = 0.0;

Real Maestro::startCPUTime = 0.0;

Maestro::Maestro() = default;

Maestro::~Maestro() = default;

Real Maestro::getCPUTime() {
    int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
    numCores = numCores * omp_get_max_threads();
#endif

    Real T = numCores * (ParallelDescriptor::second() - startCPUTime) +
             previousCPUTimeUsed;

    return T;
}
