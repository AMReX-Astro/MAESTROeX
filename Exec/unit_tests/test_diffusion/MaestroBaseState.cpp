#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0_s, BaseState<Real>& rhoh0_s,
                            BaseState<Real>& p0_s, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    const auto n = lev;
    auto rho0 = rho0_s.array();
    auto rhoh0 = rhoh0_s.array();
    auto p0 = p0_s.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();
    auto s0_init_arr = s0_init.array();

    // get species indices
    const int ihe4 = network_spec_index("helium-4");
    const int ic12 = network_spec_index("carbon-12");
    const int ife56 = network_spec_index("iron-56");

    if (ihe4 < 0 || ic12 < 0 || ife56 < 0) {
        Print() << ihe4 << ", " << ic12 << ", " << ife56 << std::endl;
        Abort("Invalid species in InitBaseState.");
    }

    eos_t eos_state;

    eos_state.h = ambient_h;
    eos_state.rho = ambient_dens;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = 0.0;
    }
    eos_state.xn[ihe4] = ambient_he4;
    eos_state.xn[ic12] = ambient_c12;
    eos_state.xn[ife56] = ambient_fe56;

    eos(eos_input_rh, eos_state);

    Real diffusion_coefficient =
        const_conductivity / (eos_state.cp * ambient_dens);

    for (auto r = 0; r < base_geom.nr(n); ++r) {
        s0_init_arr(n, r, Rho) = eos_state.rho;
        s0_init_arr(n, r, RhoH) = eos_state.rho * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init_arr(n, r, FirstSpec + comp) =
                eos_state.rho * eos_state.xn[comp];
        }
        p0_init_arr(n, r) = eos_state.p;
        s0_init_arr(n, r, Temp) = eos_state.T;
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        rho0(lev, r) = s0_init_arr(lev, r, Rho);
        rhoh0(lev, r) = s0_init_arr(lev, r, RhoH);
        tempbar_arr(lev, r) = s0_init_arr(lev, r, Temp);
        tempbar_init_arr(lev, r) = s0_init_arr(lev, r, Temp);
        p0(lev, r) = p0_init_arr(lev, r);
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
