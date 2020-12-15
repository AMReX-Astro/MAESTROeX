
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    if (AMREX_SPACEDIM != 2 && apply_vel_field) {
        Abort("Error: apply_vel_field only supported for 2d.");
    }

    const auto tileBox = mfi.tilebox();

    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    const auto offset = (prob_hi[0] - prob_lo[0]) / num_vortices;

    // vortex x-coords
    RealVector vortices_xloc(num_vortices);
    for (auto i = 0; i < num_vortices; ++i) {
        vortices_xloc[i] = (Real(i) + 0.5) * offset + prob_lo[0];
    }

    // set velocity to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });

    if (apply_vel_field) {
        const auto velpert_height = velpert_height_loc;
        const auto num_vortices_loc = num_vortices;
        const auto velpert_scale_loc = velpert_scale;
        const auto velpert_amplitude_loc = velpert_amplitude;

        const Real* vortices_xloc_p = vortices_xloc.dataPtr();

        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
            const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1];
            const Real ydist = y - velpert_height;

            Real upert = 0.0;
            Real vpert = 0.0;

            // loop over each vortex
            for (auto vortex = 0; vortex < num_vortices_loc; ++vortex) {
                Real xdist = x - vortices_xloc_p[vortex];

                Real rad = std::sqrt(x * x + y * y);

                // e.g. Calder et al. ApJSS 143, 201-229 (2002)
                // we set things up so that every other vortex has the same
                // orientation
                upert -=
                    ydist / velpert_scale_loc * velpert_amplitude_loc *
                    std::exp(-rad * rad /
                             (2.0 * velpert_scale_loc * velpert_scale_loc)) *
                    pow(-1.0, vortex + 1);

                vpert +=
                    xdist / velpert_scale_loc * velpert_amplitude_loc *
                    std::exp(-rad * rad /
                             (2.0 * velpert_scale_loc * velpert_scale_loc)) *
                    pow(-1.0, vortex + 1);
            }

            vel(i, j, k, 0) += upert;
            vel(i, j, k, 1) += vpert;
        });
    }

    const auto s0_arr = s0_init.const_array();

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
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("InitLevelDataSphr not implemented");
}
