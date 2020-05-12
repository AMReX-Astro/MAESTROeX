#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(BaseState<Real>& rho0_s, BaseState<Real>& rhoh0_s, 
                       BaseState<Real>& p0_s, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    // sanity check
    if (dens_fuel < base_cutoff_density || dens_fuel < anelastic_cutoff_density) {
        Abort("ERROR: ERROR: fuel density < (base_cutoff_density or anelastic_cutoff_density)");
    }

    const Real TINY = 1.e-10;
    const Real SMALL = 1.e-12;
    const int max_lev = base_geom.max_radial_level + 1;
    const auto nr_fine = base_geom.nr_fine;
    const int n = lev;
    auto rho0 = rho0_s.array();
    auto rhoh0 = rhoh0_s.array();
    auto p0 = p0_s.array();
    auto p0_init_arr = p0_init.array();
    auto tempbar_arr = tempbar.array();
    auto tempbar_init_arr = tempbar_init.array();

    Print() << "cutoff densities:" << std::endl;
    Print() << "    low density cutoff (for mapping the model) =      " << 
            base_cutoff_density << std::endl;
    Print() <<  "    buoyancy cutoff density                           " << std::endl;
    Print() <<  "        (for zeroing rho - rho_0, centrifugal term) = " << 
            buoyancy_cutoff_factor*base_cutoff_density << std::endl;
    Print() <<  "    anelastic cutoff =                                " << 
            anelastic_cutoff_density << std::endl;
    Print() <<  " " << std::endl;

    // const Real min_dens = min(rho_1, rho_2);

    // if (anelastic_cutoff_density > min_dens || base_cutoff_density > min_dens) {
    //    Abort("ERROR: for the RT problem, the anelastic and base cutoff densities > min(rho)");
    // }

    // if (min_dens < small_dens) {
    //     Print() << " " << std::endl;
    //     Print() << "WARNING: minimum model density is lower than the EOS cutoff" << std::endl;
    //     Print() << "         density, small_dens" << std::endl;
    // }

    // // rmid is the middle of the domain
    // Real rmid = 0.5*(geom[lev].ProbLo(AMREX_SPACEDIM-1) + geom[lev].ProbHi(AMREX_SPACEDIM-1));

    // // p0_light is the pressure at the base of the light fluid (lower
    // // half of the domain)
    // Real p0_light = p0_base;

    // // p0heavy is the pressure at the base of the heavy fluid (top half
    // // of the domain).  To find this, we integrate dp/dr = rho g from the
    // // bottom of the domain (p = p0light) to the interface (y = rmid).  The
    // // density is constant in that region, rho = rho_1
    // Real p0_heavy = p0_base + rho_1*grav_const*(rmid - geom[lev].ProbLo(AMREX_SPACEDIM-1));

    // figure out the indices for different species
    const int ic12 = network_spec_index("carbon-12");
    const int io16 = network_spec_index("oxygen-16");
    const int img24 = network_spec_index("magnesium-24");

    if (ic12 < 0 || io16 < 0 || img24 < 0) {
        Abort("ERROR: species indices not defined");
    }

    // length of the domain 
    Real rlen = geom[lev].ProbHi(AMREX_SPACEDIM-1) - geom[lev].ProbLo(AMREX_SPACEDIM-1);

    // figure out the thermodynamics of the fuel and ash state
    RealVector xn_fuel(NumSpec, 0.0);
    RealVector xn_ash(NumSpec, 0.0);

    eos_t eos_state;

    // fuel
    xn_fuel[ic12] = xc12_fuel;
    xn_fuel[io16] = 1.0 - xc12_fuel;

    eos_state.rho = dens_fuel;
    eos_state.T = temp_fuel;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = xn_fuel[comp];
    }

    eos(eos_input_rt, eos_state);

    // note: p_ambient should be = p0_init
    Real p_ambient = eos_state.p;
    Real rhoh_fuel = dens_fuel * eos_state.h;

    // ash
    xn_ash[io16] = 1.0 - xc12_fuel;
    xn_ash[img24] = xc12_fuel;

    eos_state.rho = dens_fuel; // initial guess
    eos_state.T = temp_ash;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = xn_ash[comp];
    }
    eos_state.p = p_ambient;

    eos(eos_input_tp, eos_state);

    Real dens_ash = eos_state.rho;
    Real rhoh_ash = dens_ash * eos_state.h;

    for (auto r = 0; r < base_geom.nr(n); ++r) {

        // height above the bottom of the domain
        Real rloc = geom[lev].ProbLo(AMREX_SPACEDIM-1) + (Real(r) + 0.5) * base_geom.dr(n);

        if (rloc < geom[lev].ProbLo(AMREX_SPACEDIM-1) + interface_pos_frac*rlen) {
            // fuel 
            s0_init[n+max_lev*(r+nr_fine*Rho)] = dens_fuel;
            s0_init[n+max_lev*(r+nr_fine*RhoH)] = rhoh_fuel;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                    dens_fuel * xn_fuel[comp];
            }

        } else {
            // ash
            s0_init[n+max_lev*(r+nr_fine*Rho)] = dens_ash;
            s0_init[n+max_lev*(r+nr_fine*RhoH)] = rhoh_ash;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                    dens_fuel * xn_ash[comp];
            }
        }

        // give the temperature a smooth profile
        s0_init[n+max_lev*(r+nr_fine*Temp)] = temp_fuel + 
            (temp_ash - temp_fuel) * 0.5 * (1.0 + 
            tanh((rloc - (geom[lev].ProbLo(AMREX_SPACEDIM-1) + 
            interface_pos_frac*rlen)) / (smooth_len_frac*rlen)));

        // give the carbon mass fraction a smooth profile too
        RealVector xn_smooth(NumSpec, 0.0);

        xn_smooth[ic12] = xn_fuel[ic12] + (xn_ash[ic12] - xn_fuel[ic12]) * 
            0.5 * (1.0 + tanh((rloc - (geom[lev].ProbLo(AMREX_SPACEDIM-1) + 
            interface_pos_frac*rlen)) / (smooth_len_frac*rlen)));

        xn_smooth[io16] = xn_fuel[io16];
        xn_smooth[img24] = 1.0 - xn_smooth[ic12] - xn_smooth[io16];

        // get the new density and enthalpy
        eos_state.rho   = s0_init[n+max_lev*(r+nr_fine*Rho)];
        eos_state.T     = s0_init[n+max_lev*(r+nr_fine*Temp)];
        eos_state.p     = p_ambient;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn_smooth[comp];
        }

        // (T,p) --> rho, h
        eos(eos_input_tp, eos_state);

        s0_init[n+max_lev*(r+nr_fine*Rho)] = eos_state.rho;
        s0_init[n+max_lev*(r+nr_fine*RhoH)] = eos_state.rho * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                eos_state.rho * xn_smooth[comp];
        }
        p0_init_arr(n,r) = p_ambient;
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto i = 0; i < nr_fine; ++i) {
        rho0(lev,i) = s0_init[lev+max_lev*(i+nr_fine*Rho)];
        rhoh0(lev,i) = s0_init[lev+max_lev*(i+nr_fine*RhoH)];
        tempbar_arr(lev,i) = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        tempbar_init_arr(lev,i) = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        p0(lev,i) = p0_init_arr(lev,i);
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
