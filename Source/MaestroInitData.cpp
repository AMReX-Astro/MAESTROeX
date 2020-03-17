
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
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
                       const MFIter& mfi, MultiFab& scal, MultiFab& vel, 
                       const RealVector& s0_init, 
                       const RealVector& p0_init)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);

    const auto tileBox = mfi.tilebox();
    const int max_lev = max_radial_level + 1;

    const Array4<Real> vel_arr = vel.array(mfi);
    const Array4<Real> scal_arr = scal.array(mfi);

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel_arr(i,j,k,n) = 0.0;
    });

    AMREX_PARALLEL_FOR_4D(tileBox, Nscal, i, j, k, n, {
        scal_arr(i,j,k,n) = 0.0;
    });

    // if we are spherical, we want to make sure that p0 is good, since that is
    // what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    // then initialize h from rho, X, and p0.
    Vector<MultiFab> p0_cart(finest_level);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal 
    Vector<MultiFab> temp_mf(finest_level);

    for (auto i = 0; i <= finest_level; ++i) {
        p0_cart[i].define(grids[lev], dmap[lev], 1, ng_s);
        temp_mf[i].define(grids[lev], dmap[lev], 1, ng_s);
    }

    RealVector temp_vec(max_lev*nr_fine);

    // initialize temperature 
    for (auto i = 0; i < max_lev*nr_fine; ++i) {
        temp_vec[i] = s0_init[i*Temp];
    }
    Put1dArrayOnCart(temp_vec, temp_mf, 0, 0, bcs_s, Temp);
    MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, scal.nGrow());
    
    // initialize p0_cart
    Put1dArrayOnCart(p0_init, p0_cart, 0, 0, bcs_f, 0);

    // initialize species 
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto i = 0; i < max_lev*nr_fine; ++i) {
            temp_vec[i] = s0_init[i*(FirstSpec+comp)];
        }
        Put1dArrayOnCart(temp_vec, temp_mf, 0, 0, bcs_s, FirstSpec+comp);
        MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, scal.nGrow());
    }

    const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

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
