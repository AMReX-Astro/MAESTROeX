#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& rho0, BaseState<Real>& rhoh0, 
                       BaseState<Real>& p0, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    if (spherical) {
        Abort("ERROR: rt InitBaseState is not valid for spherical");
    }

    const Real TINY = 1.e-10;
    const Real SMALL = 1.e-12;
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

    const Real min_dens = min(rho_1, rho_2);

    if (anelastic_cutoff_density > min_dens || base_cutoff_density > min_dens) {
       Abort("ERROR: for the RT problem, the anelastic and base cutoff densities > min(rho)");
    }

    if (min_dens < small_dens) {
        Print() << " " << std::endl;
        Print() << "WARNING: minimum model density is lower than the EOS cutoff" << std::endl;
        Print() << "         density, small_dens" << std::endl;
    }

    // rmid is the middle of the domain
    Real rmid = 0.5*(geom[lev].ProbLo(AMREX_SPACEDIM-1) + geom[lev].ProbHi(AMREX_SPACEDIM-1));

    // p0_light is the pressure at the base of the light fluid (lower
    // half of the domain)
    Real p0_light = p0_base;

    // p0heavy is the pressure at the base of the heavy fluid (top half
    // of the domain).  To find this, we integrate dp/dr = rho g from the
    // bottom of the domain (p = p0light) to the interface (y = rmid).  The
    // density is constant in that region, rho = rho_1
    Real p0_heavy = p0_base + rho_1*grav_const*(rmid - geom[lev].ProbLo(AMREX_SPACEDIM-1));

    // set the compositions
    const int ia = network_spec_index("A");
    const int ib = network_spec_index("B");

    RealVector xn_light(NumSpec);
    RealVector xn_heavy(NumSpec);
    for (auto comp = 0; comp < NumSpec; ++comp) {
        xn_light[comp] = SMALL;
        xn_heavy[comp] = SMALL;
    }

    xn_light[ia] = 1.0 - (NumSpec-1)*SMALL;
    xn_heavy[ib] = 1.0 - (NumSpec-1)*SMALL;

    // set a guess for the temperature for the EOS calls
    Real t_guess = 1.e-8;

    for (auto r = 0; r < base_geom.nr(n); ++r) {

        // height above the bottom of the domain
        Real rloc = (Real(r) + 0.5) * base_geom.dr(n);

        Real d_ambient = 0.0;
        Real p_ambient = 0.0;
        Real t_ambient = 0.0;
        RealVector xn_ambient(NumSpec);

        if (rloc > rmid) {
            // top half -- heavy fluid
            d_ambient = rho_2;
            p_ambient = p0_heavy + rho_2*grav_const*(rloc - rmid);
            t_ambient = t_guess;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = xn_heavy[comp];
            }
        } else {
            // lower half -- light fluid
            d_ambient = rho_1;
            p_ambient = p0_light + rho_1*grav_const*(rloc - geom[lev].ProbLo(AMREX_SPACEDIM-1));
            t_ambient = t_guess;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = xn_light[comp];
            }
        }

        eos_t eos_state;

        // use the EOS to make the state consistent
        eos_state.T     = t_ambient;
        eos_state.rho   = d_ambient;
        eos_state.p     = p_ambient;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn_ambient[comp];
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        s0_init[n+max_lev*(r+nr_fine*Rho)] = d_ambient;
        s0_init[n+max_lev*(r+nr_fine*RhoH)] = d_ambient * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                d_ambient * xn_ambient[comp];
        }
        p0_init[n+max_lev*r] = eos_state.p; // p_ambient !
        s0_init[n+max_lev*(r+nr_fine*Temp)] = t_ambient;
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto i = 0; i < nr_fine; ++i) {
        rho0[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Rho)];
        rhoh0(lev,i) = s0_init[lev+max_lev*(i+nr_fine*RhoH)];
        tempbar[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        tempbar_init[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        p0(lev,i) = p0_init[lev+max_lev*i];
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
