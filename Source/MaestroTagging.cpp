
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// 
void
Maestro::RetagArray(const int set, const int clear,
                    const Box& bx,
                    const int lev, RealVector& tag_array)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RetagArray()", RetagArray);

    Abort("Create a local copy of MaestroTagging.cpp in your build directory.\nHere is a sample that tags the refined region including buffer zone");

   // Tag on regions including buffer cells
   auto lo = bx.loVect3d()[AMREX_SPACEDIM-1];
   auto hi = bx.hiVect3d()[AMREX_SPACEDIM-1];
   const auto max_lev = max_radial_level + 1;

   AMREX_PARALLEL_FOR_1D((hi-lo)/2+1, i, {
       int r = i + lo/2;
       tag_array[lev-1+max_lev*r] = set;
   });
}


void
Maestro::TagBoxes(MultiFab& tag_mf, 
                  const int set, const int clear,
                  const MFIter mfi,
                  const int lev, RealVector& tag_array,
                  const Real time)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TagBoxes()", TagBoxes);

    Abort("Create a local copy of MaestroTagging.cpp in your build directory.\nHere is a sample that tags the full grid using a predetermined tag array");

    const Array4<Real> tag = tag_mf.array(mfi);
    const Real * AMREX_RESTRICT tag_array_p = tag_array.dataPtr();
    const int max_lev = max_radial_level + 1;

    const Box& tilebox  = mfi.tilebox();
    const auto dx = geom[lev].CellSizeArray();

    AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        if (tag_array_p[lev+max_lev*r] > 0) {
            tag(i,j,k) = set;
        }
    });
}

// 
void
Maestro::StateError(MultiFab& tag_mf, const MultiFab& state_mf, 
                   const int set, const int clear,
                   const MFIter mfi,
                   const int lev, RealVector& tag_array,
                   const Real time)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::StateError()", StateError);

    Abort("Create a local copy of MaestroTagging.cpp in your build directory.\nHere is a sample that tags the temperature above 6.5e8");

    const Array4<Real> tag = tag_mf.array(mfi);
    const Array4<const Real> state = state_mf.array(mfi);
    const Real * AMREX_RESTRICT tag_array_p = tag_array.dataPtr();
    const int max_lev = max_radial_level + 1;

    const Box& tilebox  = mfi.tilebox();
    const auto dx = geom[lev].CellSizeArray();

    // Tag on regions of high temperature
    AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
        if (state(i,j,k,Temp) >= 6.5e8) {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            tag(i,j,k) = set;
            tag_array[lev+max_lev*r] = set;
        }
    });
}