
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    // set velocity and scalars to zero
    ParallelFor(tileBox, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    vel(i, j, k, n) = 0.0;
                });
    ParallelFor(tileBox, Nscal,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    scal(i, j, k, n) = 0.0;
                });

    const auto center_p = center;
    const auto prob_lo = geom[lev].ProbLoArray();
    const auto dx = geom[lev].CellSizeArray();

    const auto base_cutoff_density_l = base_cutoff_density;
    // width of the gaussian
    const auto W = 0.05;

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const auto x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
        const auto y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
#if AMREX_SPACEDIM == 3
        const auto z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];
#else
        Real z = 0.0;
#endif

        Real dist = std::sqrt(x * x + y * y + z * z);

        scal(i, j, k, Rho) =
            amrex::max(std::exp(-dist * dist / (W * W)), base_cutoff_density_l);
        scal(i, j, k, FirstSpec) = scal(i, j, k, Rho);
    });
}
