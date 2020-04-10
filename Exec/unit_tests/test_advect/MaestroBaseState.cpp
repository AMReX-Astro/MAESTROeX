#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& rho0, BaseState<Real>& rhoh0, 
                       BaseState<Real>& p0, 
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
        rhoh0(lev,i) = 0.0;
        tempbar[lev+max_lev*i] = 0.0;
        tempbar_init[lev+max_lev*i] = 0.0;
        p0(lev,i) = 0.0;
        p0_init[lev+max_lev*i] = 0.0;
    }

}
