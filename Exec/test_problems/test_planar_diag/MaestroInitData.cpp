
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    if (AMREX_SPACEDIM != 2) {
        Abort("Error: InitLevelData only implemented for 2d.");
    }

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    const auto s0_arr = s0_init.const_array();
    const auto p0_arr = p0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
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
    const auto grav_const_loc = amrex::Math::abs(grav_const);

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const auto r = j;
        const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
        const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

        const Real rho0 = s0_arr(lev, r, Rho);

        // This seems to work ok with sealed box?
        // Real rho_local = rho0 * (1.0 +
        // 			    pert_amp_loc * std::exp(-y/(2.0 * scale_height_loc)) *
        // 			    std::cos(x * k_hoz_loc * M_PI / (prob_hi[0] - prob_lo[0]) ) *
        // 			    sin(y * k_vert_loc * M_PI / (prob_hi[1] - prob_lo[1]) ) );

        // What about this ?

        const Real rho_base = pres_base_loc / scale_height_loc / grav_const_loc;

        Real rho_local =
            rho_base * pert_amp_loc * std::exp(-y / (2.0 * scale_height_loc)) *
            std::cos(x * k_hoz_loc * M_PI / (prob_hi[0] - prob_lo[0]));

        // if k_vert is 0, dont multiply by sin(0) = 0
        if (k_vert_loc != 0.0) {
            rho_local *=
                std::sin(y * k_vert_loc * M_PI / (prob_hi[1] - prob_lo[1]));
        }
        rho_local += rho0;

        eos_t eos_state;

        eos_state.rho = rho_local;
        eos_state.p = p0_arr(lev, r);
        eos_state.T = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 1.0;  // single fluid
        }

        // (rho,p) --> T, h
        eos(eos_input_rp, eos_state);

        scal(i, j, k, Rho) = eos_state.rho;
        scal(i, j, k, RhoH) = eos_state.rho * eos_state.h;
        scal(i, j, k, Temp) = eos_state.T;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) =
                eos_state.rho * eos_state.xn[comp];
        }
    });
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("InitLevelDataSphr is not implemented.");
}
