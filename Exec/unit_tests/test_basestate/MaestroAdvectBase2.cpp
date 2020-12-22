#include <Maestro.H>

using namespace amrex;

void Maestro::UpdateSpecies(const BaseState<Real>& rho0,
                            const BaseState<Real>& rho0_predicted_edge,
                            const BaseState<Real>& rhoX0_old,
                            BaseState<Real>& rhoX0_new) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::UpdateSpecies()", UpdateSpecies);

    BaseState<Real> force(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> X0(base_geom.max_radial_level + 1, base_geom.nr_fine);
    BaseState<Real> rhoX0_edge(base_geom.max_radial_level + 1,
                               base_geom.nr_fine + 1);

    auto rho0_predicted_edge_arr = rho0_predicted_edge.array();
    auto force_arr = force.array();
    auto X0_arr = X0.array();
    auto rhoX0_edge_arr = rhoX0_edge.array();
    const auto rho0_arr = rho0.const_array();
    const auto rhoX0_old_arr = rhoX0_old.const_array();
    auto rhoX0_new_arr = rhoX0_new.array();
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto& r_edge_loc = base_geom.r_edge_loc;
    const auto w0_arr = w0.const_array();

    // Update (rho X)_0
    for (int comp = 0; comp < NumSpec; ++comp) {
        for (int n = 0; n <= base_geom.max_radial_level; ++n) {
            const auto nr = base_geom.nr(n);
            ParallelFor(nr, [=] AMREX_GPU_DEVICE(int r) {
                X0_arr(n, r) =
                    amrex::max(rhoX0_old_arr(n, r) / rho0_arr(n, r), 0.0);
            });
            Gpu::synchronize();
        }

        force.setVal(0.0);

        // here we predict (rho X)_0 on the edges
        MakeEdgeState1d(X0, rhoX0_edge, force);

        for (int n = 0; n <= base_geom.max_radial_level; ++n) {
            const auto nr = base_geom.nr(n);
            ParallelFor(nr, [=] AMREX_GPU_DEVICE(int r) {
                rhoX0_edge_arr(n, r) *= rho0_predicted_edge_arr(n, r);
            });
            Gpu::synchronize();
        }

        for (int n = 0; n <= base_geom.max_radial_level; ++n) {
            const Real dr = base_geom.dr(n);
            const Real dt_loc = dt;

            const auto nr = base_geom.nr(n);

            // update (rho X)_0
            ParallelFor(nr, [=] AMREX_GPU_DEVICE(int r) {
                if (!spherical) {
                    rhoX0_new_arr(n, r, comp) =
                        rhoX0_old_arr(n, r, comp) -
                        dt_loc / dr *
                            (rhoX0_edge_arr(n, r + 1) * w0_arr(n, r + 1) -
                             rhoX0_edge_arr(n, r) * w0_arr(n, r));
                } else {
                    rhoX0_new_arr(n, r, comp) =
                        rhoX0_old_arr(n, r, comp) -
                        dt_loc / dr / (r_cc_loc(n, r) * r_cc_loc(n, r)) *
                            (r_edge_loc(n, r + 1) * r_edge_loc(n, r + 1) *
                                 rhoX0_edge_arr(n, r + 1) * w0_arr(n, r + 1) -
                             r_edge_loc(n, r) * r_edge_loc(n, r) *
                                 rhoX0_edge_arr(n, r) * w0_arr(n, r));
                }
            });
            Gpu::synchronize();
        }
    }

    // HACK: for some reason the left edge never works so
    // for now shall just copy it
    for (int n = 0; n <= base_geom.max_radial_level; ++n) {
        ParallelFor(NumSpec, [=] AMREX_GPU_DEVICE(int comp) {
            rhoX0_new_arr(n, 0, comp) = rhoX0_old_arr(n, 0, comp);
        });
        Gpu::synchronize();
    }

    // don't let the species leave here negative
    for (int comp = 0; comp < NumSpec; ++comp) {
        for (int n = 0; n <= base_geom.max_radial_level; ++n) {
            const auto nr = base_geom.nr(n);
            ParallelFor(nr, [=] AMREX_GPU_DEVICE(int r) {
                if (rhoX0_new_arr(n, r, comp) < 0.0) {
                    Real delta = rhoX0_new_arr(n, r, comp);
                    Real sumX = 0.0;
                    for (int comp2 = 0; comp2 < NumSpec; ++comp2) {
                        if (comp2 != comp &&
                            rhoX0_new_arr(n, r, comp2) >= 0.0) {
                            sumX += rhoX0_new_arr(n, r, comp2);
                        }
                    }
                    for (int comp2 = 0; comp2 < NumSpec; ++comp2) {
                        if (comp2 != comp &&
                            rhoX0_new_arr(n, r, comp2) >= 0.0) {
                            Real frac = rhoX0_new_arr(n, r, comp2) / sumX;
                            rhoX0_new_arr(n, r, comp2) -= frac * delta;
                        }
                    }
                    rhoX0_new_arr(n, r, comp) = 0.0;
                }
            });
            Gpu::synchronize();
        }
    }
}
