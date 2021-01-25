#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0_s, BaseState<Real>& rhoh0_s,
                            BaseState<Real>& p0_s, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    if (spherical) {
        Abort("ERROR: base_state is not valid for spherical");
    }

    const int max_lev = base_geom.max_radial_level + 1;
    const int nr_fine = base_geom.nr_fine;
    const int n = lev;

    auto rho0 = rho0_s.array();
    auto rhoh0 = rhoh0_s.array();
    auto p0 = p0_s.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();
    auto s0_init_arr = s0_init.array();

    // calculate H or dens_base, as necessary
    Real H = 10.0;

    // strictly single species
    RealVector xn_zone(NumSpec, 0.0);
    xn_zone[0] = 1.0;

    // Set the state in the first cell in the vertical direction
    // ("the bottom")
    const Real z0 =
        0.5 * base_geom.dr(n);  // generic "z", actually 2d/3d agnostic

    const Real a = 2.0;
    const Real b = 2.0;
    const Real Ts = 1.0;
    const Real pi = 3.14159265358979323844;
    const Real c = H / pi;

    // set p0 and rho0 based on analytical value at the cc coord
    p0_init_arr(n, 0) = a + b * H;     // 22
    s0_init_arr(n, 0, Rho) = b + 1.0;  // 3

    eos_t eos_state;

    // use eos call to be consistent
    eos_state.p = p0_init_arr(n, 0);
    eos_state.rho = s0_init_arr(n, 0, Rho);
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = xn_zone[comp];
    }
    eos_state.T = 1000.0;          // just a guess
    eos(eos_input_rp, eos_state);  // (rho, p) --> T, h

    // set other base vars
    s0_init_arr(n, 0, RhoH) = eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        s0_init_arr(n, 0, FirstSpec + comp) = xn_zone[comp];
    }
    s0_init_arr(n, 0, Temp) = eos_state.T;

    for (auto r = 1; r < base_geom.nr(n); ++r) {
        // height above the bottom of the domain
        Real z = (Real(r) + 0.5) * base_geom.dr(n);

        // set rho analytically
        Real dens_zone = cos(z / c) + b;
        // needs to be set before pressure
        s0_init_arr(n, r, Rho) = dens_zone;

        // compute the pressure by discretizing HSE
        p0_init_arr(n, r) =
            p0_init_arr(n, r - 1) -
            base_geom.dr(n) * 0.5 *
                (s0_init_arr(n, r, Rho) + s0_init_arr(n, r - 1, Rho)) *
                fabs(grav_const);

        // use the EOS to make the state consistent
        eos_state.rho = dens_zone;
        eos_state.p = p0_init_arr(n, r);
        eos_state.T = 1000.0;  // guess
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn_zone[comp];
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        s0_init_arr(n, r, RhoH) = dens_zone * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init_arr(n, r, FirstSpec + comp) = dens_zone * xn_zone[comp];
        }
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
