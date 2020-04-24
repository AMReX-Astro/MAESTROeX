#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& rho0, RealVector& rhoh0, 
                       RealVector& p0, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    if (spherical) {
        Abort("ERROR: double_bubble InitBaseState is not valid for spherical");
    }

    const int max_lev = base_geom.max_radial_level + 1;
    const auto nr_fine = base_geom.nr_fine;
    const int n = lev;

    RealVector xn_zone(NumSpec);

    // only initialize the first species
    for (auto comp = 0; comp < NumSpec; ++comp) {
        xn_zone[comp] = 0.0;
    }
    xn_zone[0] = 1.0;

    // compute the pressure scale height (for an isothermal, ideal-gas
    // atmosphere)
    Real H = pres_base / dens_base / fabs(grav_const);

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

    p0_init[lev] = pres_base;
    s0_init[lev+max_lev*nr_fine*Rho] = dens_base;
    s0_init[lev+max_lev*nr_fine*RhoH] = dens_base * eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        s0_init[lev+max_lev*nr_fine*(FirstSpec+comp)] = dens_base * xn_zone[comp];
    }
    s0_init[lev+max_lev*nr_fine*Temp] = eos_state.T;

    Real z0 = 0.5 * base_geom.dr(n);

    // set an initial guess for the temperature -- this will be reset
    // by the EOS
    Real temp_zone = 1000.0;

    for (auto r = 1; r < base_geom.nr(n); ++r) {

        // height above the bottom of the domain
        Real z = (Real(r) + 0.5) * base_geom.dr(n);

        Real dens_zone = 0.0;
        if (do_isentropic) {
            // we can integrate HSE with p = K rho^gamma analytically
            dens_zone = dens_base * pow(grav_const * dens_base * 
                                        (gamma_const - 1.0) * (z-z0) / 
                                        (gamma_const * pres_base) + 1.0, 
                                        1.0 / (gamma_const-1.0));
        } else {
            // the density of an isothermal gamma-law atm is exponential
            dens_zone = dens_base * exp(-z / H);
        }

        s0_init[lev+max_lev*(r+nr_fine*Rho)] = dens_zone;

        // compute the pressure by discretizing HSE
        p0_init[lev+max_lev*r] = p0_init[lev+max_lev*(r-1)] - base_geom.dr(n) * 0.5 *
            (s0_init[lev+max_lev*(r+nr_fine*Rho)] + 
             s0_init[lev+max_lev*(r-1+nr_fine*Rho)]) * fabs(grav_const);

        // use the EOS to make the state consistent
        eos_state.T     = temp_zone;
        eos_state.rho   = dens_zone;
        eos_state.p     = p0_init[lev+max_lev*r];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn_zone[comp];
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        s0_init[n+max_lev*(r+nr_fine*Rho)] = dens_zone;
        s0_init[n+max_lev*(r+nr_fine*RhoH)] = dens_zone * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                dens_zone * xn_zone[comp];
        }
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

    Real min_temp = 1.e99;
    for (auto i = 0; i < nr_fine; ++i) {
        min_temp = min(min_temp, s0_init[lev+max_lev*(i+nr_fine*Temp)]);
    }

    if (min_temp < small_temp) {
        if (n == 1) {
            Print() << " " << std::endl;
            Print() << "WARNING: minimum model temperature is lower than the EOS cutoff" << std::endl;
            Print() << "         temperature, small_temp" << std::endl;
        }
    }

    Real max_hse_error = -1.e30;

    for (auto r = 1; r < base_geom.nr(n); ++r) {

        Real rloc = geom[lev].ProbLo(AMREX_SPACEDIM-1) + (Real(r) + 0.5)*base_geom.dr(n);

        Real dpdr = (p0_init[n+max_lev*r] - p0_init[n+max_lev*(r-1)]) / base_geom.dr(n);
        Real rhog = 0.5*(s0_init[n+max_lev*(r+nr_fine*Rho)] + 
                         s0_init[n+max_lev*(r-1+nr_fine*Rho)])*grav_const;

        max_hse_error = max(max_hse_error, fabs(dpdr - rhog)/fabs(dpdr));
    }

    Print() << " " << std::endl;
    Print() << "Maximum HSE Error = " << max_hse_error << std::endl;
    Print() << "   (after putting initial model into base state arrays, and" << std::endl;
    Print() << "    for density < base_cutoff_density)" << std::endl;
    Print() << " " << std::endl;

    // initialize any inlet BC parameters
    SetInletBCs();
}
