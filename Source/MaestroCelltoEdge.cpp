#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::CelltoEdge(const BaseState<Real>& s0_cell_s,
                         BaseState<Real>& s0_edge_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::CelltoEdge()", CelltoEdge);

    if (spherical) {
        Abort("Calling CelltoEdge with spherical == true");
    }

    const auto s0_cell = s0_cell_s.const_array();
    auto s0_edge = s0_edge_s.array();

    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(n); ++i) {
            Real nr_lev = base_geom.nr(n);
            const int lo = base_geom.r_start_coord(n, i);
            const int hi = base_geom.r_end_coord(n, i) + 1;
            ParallelFor(hi - lo + 1, [=] AMREX_GPU_DEVICE(int j) {
                int r = j + lo;

                if (r == 0) {
                    // if we are at lower domain boundary
                    s0_edge(n, r) = s0_cell(n, r);
                } else if (r == 1 || r == lo) {
                    // if we are at lower domain boundary+1 OR
                    // if we are at bottom of coarse-fine interface that is not a domain boundary
                    s0_edge(n, r) = 0.5 * (s0_cell(n, r - 1) + s0_cell(n, r));
                } else if (r == nr_lev) {
                    // if we are at upper domain boundary
                    s0_edge(n, r) = s0_cell(n, r - 1);
                } else if (r == nr_lev - 1 || r == hi) {
                    // if we are at upper domain boundary-1 OR
                    // if we are at top of coarse-fine interface that is not a domain boundary
                    s0_edge(n, r) = 0.5 * (s0_cell(n, r) + s0_cell(n, r - 1));
                } else {
                    // fourth order
                    Real tmp =
                        7.0 / 12.0 * (s0_cell(n, r) + s0_cell(n, r - 1)) -
                        1.0 / 12.0 * (s0_cell(n, r + 1) + s0_cell(n, r - 2));
                    Real s0min = amrex::min(s0_cell(n, r), s0_cell(n, r - 1));
                    Real s0max = amrex::max(s0_cell(n, r), s0_cell(n, r - 1));
                    s0_edge(n, r) = amrex::min(amrex::max(tmp, s0min), s0max);
                }
            });
            Gpu::synchronize();
        }
    }

    // make the edge values synchronous across levels
    RestrictBase(s0_edge_s, false);
}
