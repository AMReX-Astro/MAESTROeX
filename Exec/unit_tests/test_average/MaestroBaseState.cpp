#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0, BaseState<Real>& rhoh0,
                            BaseState<Real>& p0, const int n) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    auto p0_init_arr = p0_init.array();
    auto s0_init_arr = s0_init.array();

    for (auto r = 0; r < base_geom.nr(n); ++r) {
        // height above the bottom of the domain
        Real dist = (Real(r) + 0.5) * base_geom.dr(n);

        s0_init_arr(n, r, Rho) = exp(-dist * dist / 0.1);

        p0_init_arr(n, r) = s0_init_arr(n, r, Rho);
    }
}
