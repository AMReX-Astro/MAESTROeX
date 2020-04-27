
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

    const auto vel_fuel_loc = vel_fuel;

    // initialize velocity
    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        for (auto n = 0; n < AMREX_SPACEDIM-1; ++n) {
            vel(i,j,k,n) = 0.0;
        }
        vel(i,j,k,AMREX_SPACEDIM-1) = vel_fuel_loc;
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
    Abort("InitLevelDataSphr not implemented");
}
