#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0, BaseState<Real>& rhoh0,
                            BaseState<Real>& p0, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    if (spherical) {
        Abort("ERROR: test_advect InitBaseState is not valid for spherical");
    }

    s0_init.setVal(0.0);
    rho0.setVal(0.0);
    rhoh0.setVal(0.0);
    tempbar.setVal(0.0);
    tempbar_init.setVal(0.0);
    p0.setVal(0.0);
    p0_init.setVal(0.0);
}
