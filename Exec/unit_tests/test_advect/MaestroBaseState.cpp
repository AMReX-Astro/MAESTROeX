#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& s0_init, RealVector& p0_init, 
                       RealVector& rho0, RealVector& rhoh0, 
                       RealVector& p0, RealVector& tempbar, 
                       RealVector& tempbar_init,
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    if (spherical) {
        Abort("ERROR: test_advect InitBaseState is not valid for spherical");
    }

    const int max_lev = max_radial_level + 1;
    
    for (auto i = 0; i < nr_fine; ++i) {
        for (auto n = 0; n < Nscal; ++n) {
            s0_init[lev+max_lev*(i+nr_fine*n)] = 0.0;
        }
        rho0[lev+max_lev*i] = 0.0;
        rhoh0[lev+max_lev*i] = 0.0;
        tempbar[lev+max_lev*i] = 0.0;
        tempbar_init[lev+max_lev*i] = 0.0;
        p0[lev+max_lev*i] = 0.0;
        p0_init[lev+max_lev*i] = 0.0;
    }

}
