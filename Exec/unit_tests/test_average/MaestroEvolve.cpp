
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    for (auto lev = 0; lev <= base_geom.max_radial_level; ++lev) {
        InitBaseState(rho0_old, rhoh0_old, p0_old, lev);
    }
    const auto nr_fine = base_geom.nr_fine;

    InitFromScratch(0.0);

    Vector<MultiFab> phi(finest_level + 1);
    BaseState<Real> phi_exact(base_geom.max_radial_level + 1, nr_fine);
    BaseState<Real> phi_avg(base_geom.max_radial_level + 1, nr_fine);

    BaseState<Real> error(nr_fine);

    for (int lev = 0; lev <= finest_level; ++lev) {
        phi[lev].define(grids[lev], dmap[lev], 1, 0);
        phi[lev].setVal(0.);
    }

    phi_exact.setVal(0.);
    phi_avg.setVal(0.);
    error.setVal(0.);

    phi_exact.copy(p0_old);

    Print() << "Putting 1d array on Cartesian" << std::endl;

    Put1dArrayOnCart(phi_exact, phi, 0, 0);

    Average(phi, phi_avg, 0);

    auto error_arr = error.array();
    auto phi_exact_arr = phi_exact.array();
    auto phi_avg_arr = phi_avg.array();

    // Compare the initial and final phi
    for (int lev = 0; lev <= base_geom.max_radial_level; ++lev) {
        for (int r = 0; r < nr_fine; ++r) {
            error_arr(r) = phi_exact_arr(lev, r) - phi_avg_arr(lev, r);

            Real abs_norm = error_arr(r);

            Real rel_norm = (phi_exact_arr(lev, r) == 0.0)
                                ? 0.0
                                : error_arr(r) / phi_exact_arr(lev, r);

            Print() << "\tPhi = " << phi_exact_arr(lev, r)
                    << ",   Abs norm = " << abs_norm
                    << ",  Rel norm = " << rel_norm << std::endl;
        }
    }
}
