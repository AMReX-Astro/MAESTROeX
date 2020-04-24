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
        Abort("ERROR: mach_jet base_state is not valid for spherical");
    }

    const int max_lev = base_geom.max_radial_level + 1;
    const auto nr_fine = base_geom.nr_fine;
    const int n = lev;

    Print() << "cutoff densities:" << std::endl;
    Print() << "    low density cutoff (for mapping the model) =      " << 
            base_cutoff_density << std::endl;
    Print() <<  "    buoyancy cutoff density                           " << std::endl;
    Print() <<  "        (for zeroing rho - rho_0, centrifugal term) = " << 
            buoyancy_cutoff_factor*base_cutoff_density << std::endl;
    Print() <<  "    anelastic cutoff =                                " << 
            anelastic_cutoff_density << std::endl;
    Print() <<  " " << std::endl;
    
    if (do_stratified) {

        // use the EOS to make the state consistent
        Real dens_zone = 1.e-3;
        Real pres_zone = 1.e6;

        // only initialize the first species
        RealVector xn_zone(NumSpec, 0.0);
        xn_zone[0] = 1.0;

        p0_init[n] = pres_zone;

        // H = pres_base / dens_base / abs(grav_const)
        Real H = 1.e6 / 1.e-3 / fabs(grav_const);
        
        // set an initial guess for the temperature -- this will be reset
        // ! by the EOS
        Real temp_zone = 10.0;

        eos_t eos_state;

        Real eos_gamma = eos_state.gam1;

        for (auto r = 0; r < base_geom.nr(n); ++r) {

            Real z = (Real(r) + 0.5) * base_geom.dr(n);

            if (do_isentropic) {
                dens_zone = 1.e-3 * pow(grav_const*1.e-3*(eos_gamma - 1.0)*z / (eos_gamma*1.e6) + 1.0, 1.0/(eos_gamma-1.0));
            } else {
                dens_zone = 1.e-3 * exp(-z/H);
            }

            s0_init[n+max_lev*(r+nr_fine*Rho)] = dens_zone;

            if (r == 0) {
                p0_init[n+max_lev*r] -= base_geom.dr(n) * 0.5 * 
                    s0_init[n+max_lev*(r+nr_fine*Rho)] * fabs(grav_const);
            } else {
                p0_init[n+max_lev*r] = p0_init[n+max_lev*(r-1)] - base_geom.dr(n) * 
                    0.5 * (s0_init[n+max_lev*(r+nr_fine*Rho)] + 
                           s0_init[n+max_lev*(r-1+nr_fine*Rho)]) * 
                    fabs(grav_const);
            }

            pres_zone = p0_init[n+max_lev*r];

            // use the EOS to make the state consistent
            eos_state.T     = temp_zone;
            eos_state.rho   = dens_zone;
            eos_state.p     = pres_zone;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = xn_zone[comp];
            }

            // (rho,p) --> T, h
            eos(eos_input_rp, eos_state);

            s0_init[n+max_lev*(r+nr_fine*Rho)] = dens_zone;
            s0_init[n+max_lev*(r+nr_fine*RhoH)] = dens_zone * eos_state.h;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 0.0;
            }
            s0_init[n+max_lev*(r+nr_fine*FirstSpec)] = dens_zone;
            s0_init[n+max_lev*(r+nr_fine*Temp)] = eos_state.T;
        }

    } else {
        eos_t eos_state;
        
        eos_state.T     = 10.0;
        eos_state.rho   = 1.e-3;
        eos_state.p     = 1.e6;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 1.0 / NumSpec;
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        for (auto r = 0; r < base_geom.nr(n); ++r) {

            s0_init[n+max_lev*(r+nr_fine*Rho)] = eos_state.rho;
            s0_init[n+max_lev*(r+nr_fine*RhoH)] = eos_state.rho * eos_state.h;
	    s0_init[n+max_lev*(r+nr_fine*FirstSpec)] = eos_state.rho;
            s0_init[n+max_lev*(r+nr_fine*Temp)] = eos_state.T;

            p0_init[n+max_lev*r] = eos_state.p;
        }
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
