#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0, BaseState<Real>& rhoh0,
                            BaseState<Real>& p0, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    if (spherical) {
        Abort("ERROR: double_bubble InitBaseState is not valid for spherical");
    }

    auto rho0_arr = rho0.array();
    auto rhoh0_arr = rhoh0.array();
    auto p0_arr = p0.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();
    auto s0_init_arr = s0_init.array();

    RealVector xn_zone(NumSpec);

    // only initialize the first species
    for (auto comp = 0; comp < NumSpec; ++comp) {
        xn_zone[comp] = 0.0;
    }
    xn_zone[0] = 1.0;

    // compute the pressure scale height (for an isothermal, ideal-gas
    // atmosphere)
    Real H = pres_base / dens_base / amrex::Math::abs(grav_const);

    eos_t eos_state;

    // for isentropic, we satisfy p ~ rho^gamma, but we'll need to get gamma
    eos_state.rho = dens_base;
    eos_state.p = pres_base;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = xn_zone[comp];
    }

    // initial guess
    eos_state.T = 1000.0;

    eos(eos_input_rp, eos_state);

    Real gamma_const = pres_base / (dens_base * eos_state.e) + 1.0;

    p0_init_arr(lev, 0) = pres_base;
    s0_init_arr(lev, 0, Rho) = dens_base;
    s0_init_arr(lev, 0, RhoH) = dens_base * eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        s0_init_arr(lev, 0, FirstSpec + comp) = dens_base * xn_zone[comp];
    }
    s0_init_arr(lev, 0, Temp) = eos_state.T;

    Real z0 = 0.5 * base_geom.dr(lev);

    // set an initial guess for the temperature -- this will be reset
    // by the EOS
    Real temp_zone = 1000.0;

    for (auto r = 1; r < base_geom.nr(lev); ++r) {
        // height above the bottom of the domain
        Real z = (Real(r) + 0.5) * base_geom.dr(lev);

        Real dens_zone = 0.0;
        if (do_isentropic) {
            // we can integrate HSE with p = K rho^gamma analytically
            dens_zone =
                dens_base * pow(grav_const * dens_base * (gamma_const - 1.0) *
                                        (z - z0) / (gamma_const * pres_base) +
                                    1.0,
                                1.0 / (gamma_const - 1.0));
        } else {
            // the density of an isothermal gamma-law atm is exponential
            dens_zone = dens_base * exp(-z / H);
        }

        s0_init_arr(lev, r, Rho) = dens_zone;

        // compute the pressure by discretizing HSE
        p0_init_arr(lev, r) =
            p0_init_arr(lev, r - 1) -
            base_geom.dr(lev) * 0.5 *
                (s0_init_arr(lev, r, Rho) + s0_init_arr(lev, r - 1, Rho)) *
                amrex::Math::abs(grav_const);

        // use the EOS to make the state consistent
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        eos_state.p = p0_init_arr(lev, r);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn_zone[comp];
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        s0_init_arr(lev, r, Rho) = dens_zone;
        s0_init_arr(lev, r, RhoH) = dens_zone * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init_arr(lev, r, FirstSpec + comp) = dens_zone * xn_zone[comp];
        }
        s0_init_arr(lev, r, Temp) = eos_state.T;
    }

    // copy s0_init and p0_init_arr into rho0, rhoh0, p0, and tempbar
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        rho0_arr(lev, r) = s0_init_arr(lev, r, Rho);
        rhoh0_arr(lev, r) = s0_init_arr(lev, r, RhoH);
        tempbar_arr(lev, r) = s0_init_arr(lev, r, Temp);
        tempbar_init_arr(lev, r) = s0_init_arr(lev, r, Temp);
        p0_arr(lev, r) = p0_init_arr(lev, r);
    }

    Real min_temp = 1.e99;
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        min_temp = amrex::min(min_temp, s0_init_arr(lev, r, Temp));
    }

    if (min_temp < small_temp) {
        if (lev == 1) {
            Print() << " " << std::endl;
            Print() << "WARNING: minimum model temperature is lower than the "
                       "EOS cutoff"
                    << std::endl;
            Print() << "         temperature, small_temp" << std::endl;
        }
    }

    Real max_hse_error = -1.e30;

    for (auto r = 1; r < base_geom.nr(lev); ++r) {
        Real dpdr =
            (p0_init_arr(lev, r) - p0_init_arr(lev, r - 1)) / base_geom.dr(lev);
        Real rhog = 0.5 *
                    (s0_init_arr(lev, r, Rho) + s0_init_arr(lev, r - 1, Rho)) *
                    grav_const;

        max_hse_error =
            amrex::max(max_hse_error,
                       amrex::Math::abs(dpdr - rhog) / amrex::Math::abs(dpdr));
    }

    Print() << " " << std::endl;
    Print() << "Maximum HSE Error = " << max_hse_error << std::endl;
    Print() << "   (after putting initial model into base state arrays, and"
            << std::endl;
    Print() << "    for density < base_cutoff_density)" << std::endl;
    Print() << " " << std::endl;

    // initialize any inlet BC parameters
    SetInletBCs();
}
