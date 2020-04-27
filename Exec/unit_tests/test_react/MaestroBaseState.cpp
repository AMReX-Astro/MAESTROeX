#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(BaseState<Real>& rho0_s, BaseState<Real>& rhoh0_s, 
                       BaseState<Real>& p0_s, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    s0_init.setVal(0.0);
    rho0.setVal(0.0);
    rhoh0.setVal(0.0);
    tempbar.setVal(0.0);
    tempbar_init.setVal(0.0);
    p0.setVal(0.0);
    p0_init.setVal(0.0);

    // initialize any inlet BC parameters
    SetInletBCs();
}
