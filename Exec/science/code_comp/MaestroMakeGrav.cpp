#include <Maestro.H>

using namespace amrex;

// Define helper function
Real grav_zone(Real y) {
    Real fg;

    if (y < 1.0625 * 4.e8) {
        fg = 0.5 * (1.0 + std::sin(16.0 * M_PI * (y / 4.e8 - 1.03125)));
    } else if (y > 2.9375 * 4.e8) {
        fg = 0.5 * (1.0 - std::sin(16.0 * M_PI * (y / 4.e8 - 2.96875)));
    } else {
        fg = 1.0;
    }

    Real g = fg * g0 / std::pow((y / 4.e8), 1.25);

    return g;
}

void Maestro::MakeGravCell(BaseState<Real>& grav_cell,
                           const BaseState<Real>& rho0_s) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravCell()", MakeGravCell);

    const auto& r_cc_loc = base_geom.r_cc_loc;
    auto grav_cell_arr = grav_cell.array();

    if (spherical) {
        Abort("Error: MakeGravCell not defined for spherical!");
    }

    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        const int nr_lev = base_geom.nr(n);
        ParallelFor(nr_lev, [=] AMREX_GPU_DEVICE(long r) {
            grav_cell_arr(n, r) = grav_zone(r_cc_loc(n, r));
        });
        Gpu::synchronize();
    }

    RestrictBase(grav_cell, true);
    FillGhostBase(grav_cell, true);
}

void Maestro::MakeGravEdge(BaseState<Real>& grav_edge_state,
                           const BaseState<Real>& rho0_state) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeGravEdge()", MakeGravEdge);

    const auto& r_edge_loc = base_geom.r_edge_loc;
    auto grav_edge = grav_edge_state.array();

    if (spherical) {
        Abort("Error: MakeGravEdge not defined for spherical!");
    }

    for (auto n = 0; n <= base_geom.finest_radial_level; ++n) {
        const int nr_lev = base_geom.nr(n);
        ParallelFor(nr_lev, [=] AMREX_GPU_DEVICE(long r) {
            grav_edge(n, r) = grav_zone(r_edge_loc(n, r));
        });
        Gpu::synchronize();
    }

    RestrictBase(grav_edge, false);
    FillGhostBase(grav_edge, false);
}
