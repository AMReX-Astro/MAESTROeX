#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

auto dUdy(Real y, RealVector U);
auto fv(Real y);
auto set_species(Real y);
auto grav_zone(Real y);

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

    // define some helper functions with lambdas
    auto fv = [](Real y) {
        if (y < 1.9375 * 4.e8) {
            return 0.0;
        } else if (y > 2.0625 * 4.e8) {
            return 1.0;
        } else {
            return 0.5 * (1.0 + sin(8.0 * M_PI * (y / 4.e8 - 2.0)));
        }
    };

    auto set_species = [this, &fv](Real y) {
        RealVector xn(NumSpec, 0.0);

        xn[0] = 1.0 - fv(y);
        xn[1] = fv(y);

        return xn;
    };

    auto grav_zone = [](Real y) {
        Real fg = 1.0;

        if (y < 1.0625 * 4.e8) {
            fg = 0.5 * (1.0 + sin(16.0 * M_PI * (y / 4.e8 - 1.03125)));
        } else if (y > 2.9375 * 4.e8) {
            fg = 0.5 * (1.0 - sin(16.0 * M_PI * (y / 4.e8 - 2.96875)));
        }

        return fg * g0 / pow(y / 4.e8, 1.25);
    };

    auto dUdy = [this, &fv, &set_species, &grav_zone](Real y, RealVector U) {
        eos_t eos_state;
        RealVector xn = set_species(y);

        eos_state.rho = U[0];
        eos_state.p = U[1];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn[comp];
        }

        eos(eos_input_rp, eos_state);

        Real gamma0 = eos_state.gam1;
        Real gamma = gamma0 + fv(y) * (gamma1 - gamma0);

        RealVector dU(2);

        dU[1] = exp(U[0]) * grav_zone(y) / exp(U[1]);
        dU[0] = dU[1] / gamma;

        return dU;
    };

    const int n = lev;

    // allocate arrays
    RealVector pres(base_geom.nr(n));
    RealVector dens(base_geom.nr(n));

    RealVector U_old(2);
    RealVector U_new(2);

    Vector<RealVector> k(2);

    eos_t eos_state;

    // initialize U_old
    U_old[0] = log(rho_0);
    U_old[1] = log(p_0);

    for (auto r = 0; r < base_geom.nr(n); ++r) {
        // height above the bottom of the domain
        Real y = geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                 (Real(r) + 0.5) * base_geom.dr(n);

        // do HSE using RK2

        // out intergration starts at y - h
        Real h = r == 0 ? base_geom.dr(n) * 0.5 : base_geom.dr(n);

        auto k = dUdy(y - h, U_old);

        RealVector yprime(2);
        for (int i = 0; i < 2; ++i) {
            yprime[i] = U_old[i] + 0.5 * h * k[i];
        }

        auto kprime = dUdy(y - 0.5 * h, yprime);

        for (int i = 0; i < 2; ++i) {
            U_new[i] = U_old[i] + h * kprime[i];
        }

        dens[r] = exp(U_new[0]);
        pres[r] = exp(U_new[1]);

        for (int i = 0; i < 2; ++i) {
            U_old[i] = U_new[i];
        }
    }

    for (auto r = 0; r < base_geom.nr(n); ++r) {
        Real y = geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                 (Real(r) + 0.5) * base_geom.dr[n];
        RealVector xn = set_species(y);

        eos_state.rho = dens[r];
        eos_state.p = pres[r];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn[comp];
        }

        eos(eos_input_rp, eos_state);

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
    for (auto i = 0; i < base_geom.nr_fine; ++i) {
        rho0_arr(lev, i) = s0_init_arr(lev, i, Rho);
        rhoh0_arr(lev, i) = s0_init_arr(lev, i, RhoH);
        tempbar_arr(lev, i) = s0_init_arr(lev, i, Temp);
        tempbar_init_arr(lev, i) = s0_init_arr(lev, i, Temp);
        p0_arr(lev, i) = p0_init_arr(lev, i);
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
