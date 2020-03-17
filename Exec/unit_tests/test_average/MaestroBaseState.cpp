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

    const int max_lev = max_radial_level + 1;
    const int n = lev;

    for (auto r = 0; r< nr[n]; ++r) {
        // height above the bottom of the domain
        Real dist = (Real(r) + 0.5) * dr[n];

        s0_init[n+max_lev*(r+nr_fine*Rho)] = exp(-dist*dist/0.1);

        p0_init[n+max_lev*r] = s0_init[n+max_lev*(r+nr_fine*Rho)];
    }
}