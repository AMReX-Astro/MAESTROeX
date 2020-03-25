
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, const Array4<Real> scal, const Array4<Real> vel, 
                       const Real* s0_init, 
                       const Real* p0_init)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();
    const int max_lev = max_radial_level + 1;
    const auto nrf = nr_fine;

    if (AMREX_SPACEDIM != 2) {
        Abort("Error: InitLevelData only implemented for 2d.");
    }

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i,j,k,Rho) = s0_init[lev+max_lev*(r+nrf*Rho)];
        scal(i,j,k,RhoH) = s0_init[lev+max_lev*(r+nrf*RhoH)];
        scal(i,j,k,Temp) = s0_init[lev+max_lev*(r+nrf*Temp)];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i,j,k,FirstSpec+comp) = s0_init[lev+max_lev*(r+nrf*(FirstSpec+comp))];
        }
        // initialize pi to zero for now
        scal(i,j,k,Pi) = 0.0;
    });    

    // introduce density fluctuations 
    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    const auto pert_amp_loc = pert_amp;
    const auto scale_height_loc = scale_height;
    const auto pres_base_loc = pres_base;
    const auto k_hoz_loc = k_hoz;
    const auto k_vert_loc = k_vert;
    const auto grav_const_loc = fabs(grav_const);

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        const auto r = j;
        const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

        const Real rho_base = pres_base_loc / scale_height_loc / grav_const_loc;

        const Real rho0 = s0_init[lev+max_lev*(r+nrf*Rho)];

        Real rho_local = rho_base * pert_amp_loc * \
            std::exp(-y / (2.0 * scale_height_loc)) * \
            std::cos(x * k_hoz_loc * M_PI / (prob_hi[0] - prob_lo[0]));

        // if k_vert is 0, dont multiply by sin(0) = 0
        if (k_vert_loc != 0.0) {
            rho_local *= std::sin(y * k_vert_loc * M_PI / (prob_hi[1] - prob_lo[1]));
        } 
        rho_local += rho_base;

        eos_t eos_state;

        eos_state.rho = rho_local;
        eos_state.p = p0_init[lev+max_lev*r];
        eos_state.T = s0_init[lev+max_lev*(r+nrf*Temp)];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 1.0; // single fluid
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        scal(i,j,k,Rho) = eos_state.rho;
        scal(i,j,k,RhoH) = eos_state.rho * eos_state.h;
        scal(i,j,k,Temp) = eos_state.T;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i,j,k,FirstSpec+comp) = eos_state.rho * eos_state.xn[comp];
        }
    });
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
			   MultiFab& scal, MultiFab& vel)
{
    Abort("InitLevelDataSphr is not implemented.");
}
