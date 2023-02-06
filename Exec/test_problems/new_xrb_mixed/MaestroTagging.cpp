
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::RetagArray(const Box& bx, const int lev) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RetagArray()", RetagArray);

    // re-compute tag_array since the actual grid structure changed due to buffering
    // this is required in order to compute numdisjointchunks, r_start_coord, r_end_coord

    // Tag on regions including buffer cells
    auto lo = bx.loVect3d()[AMREX_SPACEDIM - 1];
    auto hi = bx.hiVect3d()[AMREX_SPACEDIM - 1];
    const auto max_lev = base_geom.max_radial_level + 1;

    for (auto r = lo; r <= hi; ++r) {
        tag_array[lev - 1 + max_lev * (r / 2)] = TagBox::SET;
    }
}

void Maestro::TagBoxes(TagBoxArray& tags, const MFIter& mfi, const int lev,
                       const Real time) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TagBoxes()", TagBoxes);

    // tag all cells at a given height if any cells at that height were tagged

    const Array4<TagBox::TagType> tag = tags.array(mfi);
    const int* AMREX_RESTRICT tag_array_p = tag_array.dataPtr();
    const int max_lev = base_geom.max_radial_level + 1;

    const Box& tilebox = mfi.tilebox();

    ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        if (tag_array_p[lev + max_lev * r] > 0) {
            tag(i, j, k) = TagBox::SET;
        }
    });
}

void Maestro::StateError(TagBoxArray& tags, const MultiFab& state_mf,
                         const MFIter& mfi, const int lev, const Real time) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::StateError()", StateError);

    const Box& tileBox = mfi.tilebox();

    const auto tag = tags.array(mfi);
    const Array4<const Real> state = state_mf.array(mfi);
    int* AMREX_RESTRICT tag_array_p = tag_array.dataPtr();
    const int max_lev = base_geom.max_radial_level + 1;

    const Real dr_lev = base_geom.dr(lev);

    if (use_tpert_in_tagging) {
        // Tag on regions with largest temperature perturbation
        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            int r = AMREX_SPACEDIM == 2 ? j : k;
            Real height = (Real(r) + 0.5) * dr_lev;

            if ((height > 5.4e7 && height < 1.8e8) &&
                amrex::Math::abs(state(i, j, k, Temp)) > 3.e7) {
                tag(i, j, k) = TagBox::SET;
                tag_array_p[lev + max_lev * r] = TagBox::SET;
            }
        });
    } else {
        // Tag on regions of high temperature
        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (state(i, j, k, Temp) >= 6.5e8) {
                int r = AMREX_SPACEDIM == 2 ? j : k;

                tag(i, j, k) = TagBox::SET;
                tag_array_p[lev + max_lev * r] = TagBox::SET;
            }
        });
    }
}
