#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

Real dUdy(Real y, RealVector U);
Real fv(Real y);
RealVector set_species(Real y);
Real grav_zone(Real y);

void 
Maestro::InitBaseState(RealVector& s0_init, RealVector& p0_init, 
                       RealVector& rho0, RealVector& rhoh0, 
                       RealVector& p0, RealVector& tempbar, 
                       RealVector& tempbar_init,
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    // define some helper functions with lambdas
    auto set_species = [](Real y)
    {
        RealVector xn(NumSpec, 0.0);

        xn[0] = 1.0 - fv(y);
        xn[1] = fv[y]; 

        return xn;   
    };

    auto fv = [](Real y)
    {
        if (y < 1.9375 * 4.e8) {
            return 0.0;
        } else if (y > 2.0625 * 4.e8) {
            return 1.0;
        } else {
            return 0.5 * (1.0 + sin(8.0 * M_PI * (y/4.e8 - 2.0)));
        }
    };

    auto dUdy = [](Real y, RealVector U) 
    {
        eos_t eos_state;

        eos_state.rho = U[0];
        eos_state.p = U[1];
        eos_state.xn = set_species(y);

        eos(eos_input_rp, eos_state);

        Real gamma0 = eos_state.gam1;
        Real gamma = gamma0 + fv(y) * (gamma1 - gamma0);

        RealVector dU(2);

        dU[1] = exp(U[0]) * grav_zone(y) / exp(U[1]);
        dU[0] = dU[1] / gamma;

        return dU;
    };

    auto grav_zone = [](Real y)
    {
        Real fg = 1.0;

        if (y < 1.0625 * 4.e8) {
            fg = 0.5 * (1.0 + sin(16.0 * M_PI * (y/4.e8 - 1.03125)));
        } else if (y > 2.9375 * 4.e8) {
            fg = 0.5 * (1.0 - sin(16.0 * M_PI * (y/4.e8 - 2.96875)));
        }

        return fg * g0 / pow(y / 4.e8, 1.25);
    };

    const int max_lev = max_radial_level + 1;
    const int n = lev;

    // allocate arrays
    RealVector pres(nr[n]);
    RealVector dens(nr[n]);

    RealVector U_old(2);
    RealVector U_new(2);

    Vector<RealVector> k(2);

    eos_t eos_state;

    // initialize U_old
    U_old[0] = log(rho_0);
    U_old[1] = log(p_0);

    for (auto r = 0; r < nr[n]; ++r) {

        // height above the bottom of the domain
        Real y = geom[lev].ProbLo(AMREX_SPACEDIM-1) + (Real(r) + 0.5) * dr[n];

        // do HSE using RK2

        // out intergration starts at y - h
        Real h = r == 0 ? dr[n] * 0.5 : dr[n];

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

    for (auto r = 0; r < nr[n]; ++r) {

        Real y = geom[lev].ProbLo(AMREX_SPACEDIM-1) + (Real(r) + 0.5) * dr[n];

        eos_state.rho = dens[r];
        eos_state.p = pres[r];

        eos_state.xn = set_species(y);

        eos(eos_input_rp, eos_state);

        s0_init[n+max_lev*(r+nr_fine*Rho)] = eos_state.rho;
        s0_init[n+max_lev*(r+nr_fine*RhoH)] = eos_state.rho * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                eos_state.rho * eos_state.xn[comp];
        }
        p0_init[n+max_lev*r] = eos_state.p;
        s0_init[n+max_lev*(r+nr_fine*Temp)] = eos_state.T;
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto i = 0; i < nr_fine; ++i) {
        rho0[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Rho)];
        rhoh0[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*RhoH)];
        tempbar[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        tempbar_init[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        p0[lev+max_lev*i] = p0_init[lev+max_lev*i];
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}