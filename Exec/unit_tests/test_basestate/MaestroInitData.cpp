
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, const Array4<Real> scal, const Array4<Real> vel, 
                       const Real* s0_p, 
                       const Real* p0_p)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });  

    // set scalars to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, Nscal, i, j, k, n, {
        scal(i,j,k,n) = 0.0;
    });  
}