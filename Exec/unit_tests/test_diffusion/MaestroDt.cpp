
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::EstDt() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()", EstDt);

    dt = 1.e20;

    const auto nr_fine = base_geom.nr_fine;
    const auto max_radial_level = base_geom.max_radial_level;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy((max_radial_level + 1) * nr_fine);
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(), w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_dummy(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev],
                                  1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev],
                                  1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev],
                                  1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Real dt_lev = 1.e50;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min : dt_lev)
#endif
        {
            dt_lev = 1.e50;
            Real dt_grid = 1.e50;

            for (MFIter mfi(uold[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                // Get the spacing of the valid region
                const auto dx = geom[lev].CellSizeArray();

                // calculate the timestep
                for (auto i = 0; i < AMREX_SPACEDIM; ++i) {
                    dt_grid = amrex::min(dt_grid,
                                         dx[i] * dx[i] / diffusion_coefficient);
                }
            }

            dt_lev = amrex::min(dt_lev, dt_mult_factor * dt_grid);
        }  //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to estdt for level " << lev
                    << " gives dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = amrex::min(dt, dt_lev);
    }  // end loop over levels

    if (maestro_verbose > 0) {
        Print() << "Minimum estdt over all levels = " << dt << std::endl;
    }

    if (dt < small_dt) {
        Abort("EstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << std::endl;
        }
        dt = max_dt;
    }

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }
}

void Maestro::FirstDt() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FirstDt()", FirstDt);

    dt = 1.e20;

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
        return;
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev)
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old, p0_cart, 0, 0, bcs_f, 0);
    Put1dArrayOnCart(gamma1bar_old, gamma1bar_cart, 0, 0, bcs_f, 0);

    for (int lev = 0; lev <= finest_level; ++lev) {
        Real dt_lev = 1.e99;

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min : dt_lev)
#endif
        {
            Real dt_grid = 1.e99;

            for (MFIter mfi(sold[lev], true); mfi.isValid(); ++mfi) {
                // Get the spacing of the valid region
                const Real* dx = geom[lev].CellSize();

                // local variables
                for (auto i = 0; i < AMREX_SPACEDIM; ++i) {
                    dt_grid = amrex::min(dt_grid,
                                         dx[i] * dx[i] / diffusion_coefficient);
                }
            }

            dt_lev = amrex::min(dt_lev, dt_mult_factor * dt_grid);
        }  // end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to firstdt for level " << lev
                    << " gives dt_lev = " << dt_lev << std::endl;
        }

        // multiply by init_shrink
        dt_lev *= init_shrink;

        if (maestro_verbose > 0) {
            Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev
                    << std::endl;
        }

        // update dt over all levels
        dt = amrex::min(dt, dt_lev);

    }  // end loop over levels

    if (maestro_verbose > 0) {
        Print() << "Minimum firstdt over all levels = " << dt << std::endl;
    }

    if (dt < small_dt) {
        Abort("FirstDt: dt < small_dt");
    }

    if (dt > max_dt) {
        if (maestro_verbose > 0) {
            Print() << "max_dt limits the new dt = " << max_dt << std::endl;
        }
        dt = max_dt;
    }
}
