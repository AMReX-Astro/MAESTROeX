#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::InitBaseState(BaseState<Real>& rho0_s, BaseState<Real>& rhoh0_s,
                            BaseState<Real>& p0_s, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState);

    if (spherical) {
        Abort(
            "ERROR: Incompressible shear jet base_state is not valid for "
            "spherical");
    }

    const Real SMALL = 1.e-12;
    const int n = lev;

    auto rho0 = rho0_s.array();
    auto rhoh0 = rhoh0_s.array();
    auto p0 = p0_s.array();
    auto tempbar_arr = tempbar.array();
    auto p0_init_arr = p0_init.array();
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

    // Check min density
    const Real min_dens = rho_base;

    if (min_dens < small_dens) {
        Print() << " " << std::endl;
        Print() << "WARNING: minimum model density is lower than the EOS cutoff"
                << std::endl;
        Print() << "         density, small_dens" << std::endl;
    }

    // Location of the two thin shear layers
    const Real rshr1 = 0.25 * (geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                               geom[lev].ProbHi(AMREX_SPACEDIM - 1));
    const Real rshr2 = 0.75 * (geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                               geom[lev].ProbHi(AMREX_SPACEDIM - 1));

    // set the compositions
    // A and B act as tags.
    // A tags the part of the fluid initially in the "jet zone," i.e. the part of
    // the fluid with positive initial velocity.
    // B tags the part of the fluid with negative initial velocity
    const int ia = network_spec_index("A");
    const int ib = network_spec_index("B");

    // xn is an array where the ith element is the mass fraction (percent) of the
    // ith species.  Initially, the middle half of the domain is composed entirely of
    // "A" particles and negligible "B" particles, while the reverse is true of
    // the outer regions of the domain.

    // As we evolve in time these arrays track how the two fluids mix.
    RealVector xn_jet(NumSpec);
    RealVector xn_still(NumSpec);
    for (auto comp = 0; comp < NumSpec; ++comp) {
        xn_jet[comp] = SMALL;
        xn_still[comp] = SMALL;
    }

    xn_jet[ia] = 1.0 - (NumSpec - 1) * SMALL;
    xn_still[ib] = 1.0 - (NumSpec - 1) * SMALL;

    // set a guess for the temperature for the EOS calls
    Real t_guess = 1.e-8;

    // Iterate through each base state cell and initialize
    //   -The components of the fluid state 's':
    //       density, enthalpy, species mass fractions, and temperature
    //   -The pressure (note, pressure is NOT a component of the 's' multifab)
    for (auto r = 0; r < base_geom.nr(n); ++r) {
        // height above the bottom of the domain
        Real rloc = (Real(r) + 0.5) * base_geom.dr(n);

        // Init density, pressure, and temp
        Real d_ambient = rho_base;
        Real p_ambient = p_base;
        Real t_ambient = t_guess;
        RealVector xn_ambient(NumSpec);

        // Depending on location, initialize the mass fraction
        if (rloc > rshr1 && rloc < rshr2) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = xn_jet[comp];
            }
        } else {
            // Outer regions -- still fluid
            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = xn_still[comp];
            }
        }

        eos_t eos_state;

        // use the EOS to make the state consistent
        // We set density and pressure, and from this the EoS yields many
        // thermodynamic quantities (temperature and enthalpy being the two we
        // care about in this problem).
        eos_state.T = t_ambient;
        eos_state.rho = d_ambient;
        eos_state.p = p_ambient;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = xn_ambient[comp];
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        s0_init_arr(lev, r, Rho) = d_ambient;
        s0_init_arr(lev, r, RhoH) = d_ambient * eos_state.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            s0_init_arr(n, r, FirstSpec + comp) = d_ambient * xn_ambient[comp];
        }
        p0_init_arr(n, r) = p_base;
        s0_init_arr(lev, r, Temp) = eos_state.T;
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        rho0(lev, r) = s0_init_arr(lev, r, Rho);
        rhoh0(lev, r) = s0_init_arr(lev, r, RhoH);
        tempbar_arr(lev, r) = s0_init_arr(lev, r, Temp);
        tempbar_init_arr(lev, r) = s0_init_arr(lev, r, Temp);
        p0(lev, r) = p0_init_arr(lev, r);
    }

    // Check that the temperature is consistent with the EoS
    Real min_temp = 1.e99;
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
        min_temp = amrex::min(min_temp, s0_init_arr(lev, r, Temp));
    }

    if (min_temp < small_temp) {
        if (n == 1) {
            Print() << " " << std::endl;
            Print() << "WARNING: minimum model temperature is lower than the "
                       "EOS cutoff"
                    << std::endl;
            Print() << "         temperature, small_temp" << std::endl;
        }
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
