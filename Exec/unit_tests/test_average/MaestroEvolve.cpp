
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    for (auto lev = 0; lev <= base_geom.max_radial_level; ++lev) {
        InitBaseState(rho0_old, rhoh0_old,
                p0_old, lev);
	}
    const auto nr_fine = base_geom.nr_fine;
	
	InitFromScratch(0.0);

	Vector<MultiFab> phi(finest_level+1);
	Vector<Real> phi_exact((base_geom.max_radial_level+1)*nr_fine);
	Vector<Real> phi_avg((base_geom.max_radial_level+1)*nr_fine);

	Vector<Real> error(nr_fine);

	for (int lev=0; lev<=finest_level; ++lev) {
		phi[lev].define(grids[lev], dmap[lev], 1, 0);
		phi[lev].setVal(0.);
	}

	std::fill(phi_exact.begin(), phi_exact.end(), 0.);
	std::fill(phi_avg.begin(), phi_avg.end(), 0.);
	std::fill(error.begin(), error.end(), 0.);

	p0_old.setVal(0.);

	Print() << "Putting 1d array on Cartesian" << std::endl;

	Put1dArrayOnCart(phi_exact, phi, 0, 0);

	Average(phi, phi_avg, 0);

	// Compare the initial and final phi
	for (int lev=0; lev<=base_geom.max_radial_level; ++lev) {
		for (int r=0; r < nr_fine; ++r) {
			int idx = lev*nr_fine + r;
			error[r] = phi_exact[idx] - phi_avg[idx];

			Real abs_norm = error[r];

			Real rel_norm = (phi_exact[idx] == 0.0) ? 0.0 : error[r] / phi_exact[idx];

			Print() << "\tPhi = " << phi_exact[idx] << ",   Abs norm = " << abs_norm << ",  Rel norm = " << rel_norm << std::endl;
		}
	}
}
