#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(BaseState<Real>& rho0, BaseState<Real>& rhoh0, 
                       BaseState<Real>& p0, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    const int n = lev;

    for (auto r = 0; r < nr[n]; ++r) {
        // height above the bottom of the domain
        Real dist = (Real(r) + 0.5) * dr[n];

        s0_init(n,r,Rho) = exp(-dist*dist/0.1);

        p0_init(n,r) = s0_init(n,r,Rho);
    }
}
