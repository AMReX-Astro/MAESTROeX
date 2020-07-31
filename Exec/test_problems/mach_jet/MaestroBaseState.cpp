#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0_s, BaseState<Real>& rhoh0_s,
                            BaseState<Real>& p0_s, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    if (spherical) {
        Abort("ERROR: mach_jet base_state is not valid for spherical");
    }

    const int n = lev;

    auto rho0 = rho0_s.array();
    auto rhoh0 = rhoh0_s.array();
    auto p0 = p0_s.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();
    auto s0_init_arr = s0_init.array();

    Print() << "cutoff densities:" << std::endl;
    Print() << "    low density cutoff (for mapping the model) =      "
            << base_cutoff_density << std::endl;
    Print() << "    buoyancy cutoff density                           "
            << std::endl;
    Print() << "        (for zeroing rho - rho_0, centrifugal term) = "
            << buoyancy_cutoff_factor * base_cutoff_density << std::endl;
    Print() << "    anelastic cutoff =                                "
            << anelastic_cutoff_density << std::endl;
    Print() << " " << std::endl;

    if (do_stratified) {
        // use the EOS to make the state consistent
        Real dens_zone = 1.e-3;
        Real pres_zone = 1.e6;

        // only initialize the first species
        RealVector xn_zone(NumSpec, 0.0);
        xn_zone[0] = 1.0;

        p0_init_arr(n, 0) = pres_zone;

        // H = pres_base / dens_base / abs(grav_const)
        Real H = 1.e6 / 1.e-3 / amrex::Math::abs(grav_const);

        // set an initial guess for the temperature -- this will be reset
        // ! by the EOS
        Real temp_zone = 10.0;

        eos_t eos_state;

        Real eos_gamma = eos_state.gam1;

        for (auto r = 0; r < base_geom.nr(n); ++r) {
            Real z = (Real(r) + 0.5) * base_geom.dr(n);

            if (do_isentropic) {
                dens_zone = 1.e-3 * pow(grav_const * 1.e-3 * (eos_gamma - 1.0) *
                                                z / (eos_gamma * 1.e6) +
                                            1.0,
                                        1.0 / (eos_gamma - 1.0));
            } else {
                dens_zone = 1.e-3 * exp(-z / H);
            }

            s0_init_arr(n, r, Rho) = dens_zone;

            if (r == 0) {
                p0_init_arr(n, r) -= base_geom.dr(n) * 0.5 *
                                     s0_init_arr(n, r, Rho) *
                                     amrex::Math::abs(grav_const);
            } else {
                p0_init_arr(n, r) =
                    p0_init_arr(n, r - 1) -
                    base_geom.dr(n) * 0.5 *
                        (s0_init_arr(n, r, Rho) + s0_init_arr(n, r - 1, Rho)) *
                        amrex::Math::abs(grav_const);
            }

            pres_zone = p0_init_arr(n, r);

            // use the EOS to make the state consistent
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            eos_state.p = pres_zone;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = xn_zone[comp];
            }

            // (rho,p) --> T, h
            eos(eos_input_rp, eos_state);

            s0_init_arr(n, r, Rho) = dens_zone;
            s0_init_arr(n, r, RhoH) = dens_zone * eos_state.h;
            for (auto comp = 1; comp < NumSpec; ++comp) {
                s0_init_arr(n, r, FirstSpec + comp) = 0.0;
            }
            s0_init_arr(n, r, FirstSpec) = dens_zone;
            s0_init_arr(n, r, Temp) = eos_state.T;
        }
    } else {
        eos_t eos_state;

        eos_state.T = 10.0;
        eos_state.rho = 1.e-3;
        eos_state.p = 1.e6;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 1.0 / NumSpec;
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        for (auto r = 0; r < base_geom.nr(n); ++r) {
            s0_init_arr(lev, r, Rho) = eos_state.rho;
            s0_init_arr(lev, r, RhoH) = eos_state.rho * eos_state.h;
            s0_init_arr(n, r, FirstSpec) = eos_state.rho;
            s0_init_arr(n, r, Temp) = eos_state.T;

            p0_init_arr(n, r) = eos_state.p;
        }
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
