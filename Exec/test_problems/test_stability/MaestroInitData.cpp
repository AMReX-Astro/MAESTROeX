
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
    Abort("InitLevelDataSphr not implemented;");
}
