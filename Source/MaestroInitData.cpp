
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

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
#if NAUX_NET > 0
        for (auto comp = 0; comp < NumAux; ++comp) {
            scal(i, j, k, FirstAux + comp) = s0_arr(lev, r, FirstAux + comp);
        }
#endif
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
    });
}

void Maestro::InitLevelDataSphr(const int lev, [[maybe_unused]] const Real time, MultiFab& scal,
                                MultiFab& vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        // set velocity to zero
        ParallelFor(tileBox, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        vel_arr(i, j, k, n) = 0.0;
                    });

        ParallelFor(tileBox, Nscal,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        scal_arr(i, j, k, n) = 0.0;
                    });
    }

    // if we are spherical, we want to make sure that p0 is good, since that is
    // what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    // then initialize h from rho, X, and p0.
    MultiFab p0_cart(scal.boxArray(), scal.DistributionMap(), 1, 0);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal
    MultiFab temp_mf(scal.boxArray(), scal.DistributionMap(), 1, 0);

    BaseState<Real> temp_vec(base_geom.max_radial_level + 1, base_geom.nr_fine);
    auto temp_arr = temp_vec.array();
    const auto s0_init_arr = s0_init.const_array();

    // initialize temperature
    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto r = 0; r < base_geom.nr_fine; ++r) {
            temp_arr(l, r) = s0_init_arr(l, r, Temp);
        }
    }

    Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_f, 0);
    MultiFab::Copy(scal, temp_mf, 0, Temp, 1, 0);

    // initialize p0_cart
    Put1dArrayOnCart(lev, p0_init, p0_cart, 0, 0);

    // initialize species
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            for (auto r = 0; r < base_geom.nr_fine; ++r) {
                temp_arr(l, r) = s0_init_arr(l, r, FirstSpec + comp);
            }
        }
        Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_s, FirstSpec + comp);
        MultiFab::Copy(scal, temp_mf, 0, FirstSpec + comp, 1, 0);
    }

    // initialize aux
#if NAUX_NET > 0
    for (auto comp = 0; comp < NumAux; ++comp) {
        for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
            for (auto r = 0; r < base_geom.nr_fine; ++r) {
                temp_arr(l, r) = s0_init_arr(l, r, FirstAux + comp);
            }
        }
        Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_s, FirstAux + comp);
        MultiFab::Copy(scal, temp_mf, 0, FirstAux + comp, 1, 0);
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        const Array4<const Real> p0_arr = p0_cart.array(mfi);

        // initialize rho as sum of partial densities rho*X_i
        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal_arr(i, j, k, Rho) += scal_arr(i, j, k, FirstSpec + comp);
            }

            // initialize (rho h) and T using the EOS
            eos_t eos_state;
            eos_state.T = scal_arr(i, j, k, Temp);
            eos_state.p = p0_arr(i, j, k);
            eos_state.rho = scal_arr(i, j, k, Rho);
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] =
                    scal_arr(i, j, k, FirstSpec + comp) / eos_state.rho;
            }
#if NAUX_NET > 0
            for (auto comp = 0; comp < NumAux; ++comp) {
                eos_state.aux[comp] =
                    scal_arr(i, j, k, FirstAux + comp) / eos_state.rho;
            }
#endif

            eos(eos_input_rp, eos_state);

            scal_arr(i, j, k, RhoH) = eos_state.rho * eos_state.h;
            scal_arr(i, j, k, Temp) = eos_state.T;
        });
    }
}
