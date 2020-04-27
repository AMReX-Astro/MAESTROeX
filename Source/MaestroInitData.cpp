
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, const Array4<Real> scal, const Array4<Real> vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();
    const int max_lev = base_geom.max_radial_level + 1;
    const auto nrf = base_geom.nr_fine;

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });

    const Real * AMREX_RESTRICT s0_p = s0_init.dataPtr();

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i,j,k,Rho) = s0_p[lev+max_lev*(r+nrf*Rho)];
        scal(i,j,k,RhoH) = s0_p[lev+max_lev*(r+nrf*RhoH)];
        scal(i,j,k,Temp) = s0_p[lev+max_lev*(r+nrf*Temp)];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i,j,k,FirstSpec+comp) = s0_p[lev+max_lev*(r+nrf*(FirstSpec+comp))];
        }
        // initialize pi to zero for now
        scal(i,j,k,Pi) = 0.0;
    });    
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
                           MultiFab& scal, MultiFab& vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);
    const int max_lev = base_geom.max_radial_level + 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        // set velocity to zero 
        AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
            vel_arr(i,j,k,n) = 0.0;
        });

        AMREX_PARALLEL_FOR_4D(tileBox, Nscal, i, j, k, n, {
            scal_arr(i,j,k,n) = 0.0;
        });
    }

    // if we are spherical, we want to make sure that p0 is good, since that is
    // what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    // then initialize h from rho, X, and p0.
    MultiFab p0_cart(scal.boxArray(), scal.DistributionMap(), 1, 0);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal 
    MultiFab temp_mf(scal.boxArray(), scal.DistributionMap(), 1, 0);

    RealVector temp_vec(max_lev*base_geom.nr_fine);

    // initialize temperature 
    for (auto l = 0; l < max_lev; ++l) {
        for (auto n = 0; n < base_geom.nr_fine; ++n) {
            temp_vec[l+max_lev*n] = s0_init[l+max_lev*(n+base_geom.nr_fine*Temp)];
        }
    }

    Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_f, 0);
    MultiFab::Copy(scal, temp_mf, 0, Temp, 1, 0);
    
    // initialize p0_cart
    Put1dArrayOnCart(lev, p0_init, p0_cart, 0, 0);

    // initialize species 
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto l = 0; l < max_lev; ++l) {
            for (auto n = 0; n < base_geom.nr_fine; ++n) {
                temp_vec[l+max_lev*n] = s0_init[l+max_lev*(n+base_geom.nr_fine*(FirstSpec+comp))];
            }
        }
        Put1dArrayOnCart(lev, temp_vec, temp_mf, 0, 0, bcs_s, FirstSpec+comp);
        MultiFab::Copy(scal, temp_mf, 0, FirstSpec+comp, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        const Array4<const Real> p0_arr = p0_cart.array(mfi);

        // initialize rho as sum of partial densities rho*X_i
        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal_arr(i,j,k,Rho) += scal_arr(i,j,k,FirstSpec+comp);
            }

            // initialize (rho h) and T using the EOS
            eos_t eos_state;
            eos_state.T     = scal_arr(i,j,k,Temp);
            eos_state.p     = p0_arr(i,j,k);
            eos_state.rho   = scal_arr(i,j,k,Rho);
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = scal_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
            }

            eos(eos_input_rp, eos_state);

            scal_arr(i,j,k,RhoH) = eos_state.rho*eos_state.h;
            scal_arr(i,j,k,Temp) = eos_state.T;
        });
    }
}
