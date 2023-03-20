
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, [[maybe_unused]] const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    const auto s0_arr = s0_init.const_array();
    const auto p0_arr = p0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
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

    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    if (perturb_model) {
        const auto xcen = center[0];
        const auto ycen = AMREX_SPACEDIM == 2 ? xrb_pert_height : center[1];
#if AMREX_SPACEDIM == 3
        const auto zcen = xrb_pert_height;
#endif

        const auto rad_pert =
            -xrb_pert_size * xrb_pert_size / (4.0 * std::log(0.5));

        const bool perturb_temp_true = xrb_pert_type == 1;

        const auto xrb_pert_factor_loc = xrb_pert_factor;
        const auto rad_pert_loc = rad_pert;

        ParallelFor(tileBox, [=] (int i, int j, int k) {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            const auto x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - xcen;
            const auto y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - ycen;
#if AMREX_SPACEDIM == 3
            const auto z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - zcen;
#else
            const Real z = 0.0;
#endif

            const auto dist = std::sqrt(x * x + y * y + z * z);

            Real temp = 0.0;
            Real dens = 0.0;
            auto eos_input_flag = eos_input_tp;

            if (perturb_temp_true) {
                Real t0 = s0_arr(lev, r, Temp);
                temp = t0 * (1.0 + xrb_pert_factor_loc *
                                       std::exp(-dist * dist / rad_pert_loc));
                dens = s0_arr(lev, r, Rho);
                eos_input_flag = eos_input_tp;
            } else {
                Real d0 = s0_arr(lev, r, Rho);
                dens = d0 * (1.0 + xrb_pert_factor_loc *
                                       std::exp(-dist * dist / rad_pert_loc));
                temp = s0_arr(lev, r, Temp);
                eos_input_flag = eos_input_rp;
            }

            eos_t eos_state;

            eos_state.T = temp;
            eos_state.p = p0_arr(lev, r);
            eos_state.rho = dens;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] =
                    s0_arr(lev, r, FirstSpec + comp) / s0_arr(lev, r, Rho);
            }

            eos(eos_input_flag, eos_state);

            scal(i, j, k, Rho) = eos_state.rho;
            scal(i, j, k, RhoH) = eos_state.rho * eos_state.h;
            scal(i, j, k, Temp) = eos_state.T;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) =
                    eos_state.rho * eos_state.xn[comp];
            }
        });
    }

    if (apply_vel_field) {
        const auto velpert_amplitude_loc = velpert_amplitude;
        const auto velpert_scale_loc = velpert_scale;
        const auto num_vortices_loc = num_vortices;
        const auto velpert_height = velpert_height_loc;

        const Real offset = (prob_hi[0] - prob_lo[0]) / (num_vortices);

        // vortex x-coords
        RealVector vortices_xloc(num_vortices);
        for (auto i = 0; i < num_vortices; ++i) {
            vortices_xloc[i] = (static_cast<Real>(i) + 0.5_rt) * offset;
        }

        const Real* vortices_xloc_p = vortices_xloc.dataPtr();

        ParallelFor(tileBox, [=] (int i, int j, int k) {
            const auto x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
            const auto y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

            Real ydist = y - velpert_height;

            Real upert = 0.0;
            Real vpert = 0.0;

            for (auto vortex = 0; vortex < num_vortices_loc; ++vortex) {
                Real xdist = x - vortices_xloc_p[vortex];

                Real rad = std::sqrt(xdist * xdist + ydist * ydist);

                // e.g. Calder et al. ApJSS 143, 201-229 (2002)
                // we set things up so that every other vortex has the same
                // orientation
                upert -=
                    ydist * velpert_amplitude_loc *
                    std::exp(-rad * rad /
                             (2.0 * velpert_scale_loc * velpert_scale_loc)) *
                    pow(-1.0, vortex + 1);

                vpert +=
                    xdist * velpert_amplitude_loc *
                    std::exp(-rad * rad /
                             (2.0 * velpert_scale_loc * velpert_scale_loc)) *
                    pow(-1.0, vortex + 1);
            }
            vel(i, j, k, 0) += upert;
            vel(i, j, k, 1) += vpert;
        });
    }
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {

    amrex::ignore_unused(lev);
    amrex::ignore_unused(time);
    amrex::ignore_unused(scal);
    amrex::ignore_unused(vel);

    Abort("InitLevelDataSphr not implemented.");
}
