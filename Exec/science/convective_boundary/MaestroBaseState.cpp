#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

auto get_rho0(Real z);
auto get_p0(Real rho_0, Real X0);

void Maestro::InitBaseState(BaseState<Real>& rho0, BaseState<Real>& rhoh0,
                            BaseState<Real>& p0, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    auto rho0_arr = rho0.array();
    auto rhoh0_arr = rhoh0.array();
    auto p0_arr = p0.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();
    auto s0_init_arr = s0_init.array();

    const auto iH = network_spec_index("helium");
    const auto iH20 = network_spec_index("H20");

    const auto alpha = -N_rho / (std::log(9.0_rt / (9.0_rt - 8.0_rt * X_b)));
    const auto rho_b = rho_t * std::exp(N_rho);
    const auto R = k_B / m_p;
    const auto D = (1.0_rt + alpha) * 4.0_rt * X_b * R * T_0 / (9.0_rt * g0);


    // define some helper functions with lambdas
    auto get_rho0 = [=](Real z) {
        // return the background density 
        return rho_b * std::pow(1.0_rt + 8.0_rt * X_b * z / (9.0_rt - 8.0_rt * X_b) / D, alpha);
    };

    auto get_p0 = [=](Real rho_0, Real X0) {
        // return the background pressure
        return R * rho_0 * T_0 * (9.0_rt - 8.0_rt * X0) / 18.0_rt;
    };

    const int n = lev;
    eos_t eos_state;

    for (auto r = 0; r < base_geom.nr(n); ++r) {
        Real z = geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                 (Real(r) + 0.5) * base_geom.dr[n];

        eos_state.rho = get_rho0(z);
        // Print() << "z/D = " << z/D << ",  rho/rho_t = " << eos_state.rho/rho_t << std::endl;
        auto X0 = X_b * (1.0_rt - z / D);
        eos_state.p = get_p0(eos_state.rho, X0);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 0.0_rt;
        }
        eos_state.xn[iH20] = X0;
        eos_state.xn[iH] = 1.0_rt - X0;

        eos(eos_input_rp, eos_state);

        s0_init_arr(n, r, Rho) = eos_state.rho;
        s0_init_arr(n, r, RhoH) = eos_state.rho * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init_arr(n, r, FirstSpec + comp) =
                eos_state.rho * eos_state.xn[comp];
        }
        p0_init_arr(n, r) = get_p0(eos_state.rho, X0);
        s0_init_arr(n, r, Temp) = T_0;
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto i = 0; i < base_geom.nr_fine; ++i) {
        rho0_arr(lev, i) = s0_init_arr(lev, i, Rho);
        rhoh0_arr(lev, i) = s0_init_arr(lev, i, RhoH);
        tempbar_arr(lev, i) = s0_init_arr(lev, i, Temp);
        tempbar_init_arr(lev, i) = s0_init_arr(lev, i, Temp);
        p0_arr(lev, i) = p0_init_arr(lev, i);

        // Print() << "p0 = " << p0_arr(lev,i) << std::endl;
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
