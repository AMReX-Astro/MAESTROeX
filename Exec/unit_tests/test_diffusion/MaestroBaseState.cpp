#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& rho0, BaseState<Real>& rhoh0_s, 
                       BaseState<Real>& p0_s, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    const int max_lev = base_geom.max_radial_level + 1;
    const auto nr_fine = base_geom.nr_fine;
    const int n = lev;
    auto rhoh0 = rhoh0_s.array();
    auto p0 = p0_s.array();
    
    // get species indices
    const int ihe4 = network_spec_index("helium-4");
    const int ic12 = network_spec_index("carbon-12");
    const int ife56 = network_spec_index("iron-56");

    if (ihe4 < 0 || ic12 < 0 || ife56 < 0) {
        Print() << ihe4 << ", " << ic12 << ", " << ife56 << std::endl;
        Abort("Invalid species in InitBaseState.");
    }

    eos_t eos_state;

    eos_state.h = ambient_h;
    eos_state.rho = ambient_dens;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = 0.0;
    }
    eos_state.xn[ihe4] = ambient_he4;
    eos_state.xn[ic12] = ambient_c12;
    eos_state.xn[ife56] = ambient_fe56;

    eos(eos_input_rh, eos_state);

    Real diffusion_coefficient = const_conductivity / (eos_state.cp * ambient_dens);
   
    for (auto r = 0; r < base_geom.nr(n); ++r) {

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
        rhoh0(lev,i) = s0_init[lev+max_lev*(i+nr_fine*RhoH)];
        tempbar[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        tempbar_init[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        p0(lev,i) = p0_init[lev+max_lev*i];
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
