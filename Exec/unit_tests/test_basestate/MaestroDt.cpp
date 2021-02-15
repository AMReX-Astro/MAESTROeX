
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::EstDt() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()", EstDt);

    // allocate a dummy w0_force and set equal to zero
    Vector<Real> w0_force_dummy((base_geom.max_radial_level + 1) *
                                base_geom.nr_fine);
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
#if (AMREX_SPACEDIM >= 2)
        umac_dummy[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev],
                                  1, 1);
        umac_dummy[lev][1].setVal(0.);
#endif
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev],
                                  1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // face-centered
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > w0mac(finest_level + 1);

    if (spherical == 1) {
        // initialize
        for (int lev = 0; lev <= finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev],
                                 1, 1);
#if (AMREX_SPACEDIM >= 2)
            w0mac[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev],
                                 1, 1);
#endif
#if (AMREX_SPACEDIM == 3)
            w0mac[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev],
                                 1, 1);
#endif
        }
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
#endif

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {
        Real dt_lev = dt;
        Real umax_lev = 0.;

        // get references to the MultiFabs at level lev
        const MultiFab& w0_mf = w0_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min : dt_lev) reduction(max : umax_lev)
#endif
        {
            Real dt_grid = dt;
            Real umax_grid = 0.;
            const Real SMALL = 1.0e-12;

            for (MFIter mfi(uold[lev], true); mfi.isValid(); ++mfi) {
                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const Real* dx = geom[lev].CellSize();

                const Real prob_lo = geom[lev].ProbLo()[0];
                const Real prob_hi = geom[lev].ProbHi()[0];
                Real dr_fine_loc;

                if (spherical == 0) {
                    const Real maxw0 = w0_mf[mfi].maxabs<RunOn::Device>(
                        tileBox, AMREX_SPACEDIM - 1);

                    dr_fine_loc = (prob_hi - prob_lo) / base_geom.nr_fine;
                    dt_grid = amrex::min(1.1 * dt_grid,
                                         cfl * dr_fine_loc / (maxw0 + SMALL));

                } else {
#if (AMREX_SPACEDIM == 3)
                    const Real maxw0 =
                        w0_mf[mfi].maxabs<RunOn::Device>(tileBox, 0);

                    if (prob_type == 1 || prob_type == 3 || prob_type == 4) {
                        dr_fine_loc = (prob_hi - prob_lo) / base_geom.nr_fine;
                    } else {
                        // need to compute it this way to agree with how the initial model
                        // was computed
                        dr_fine_loc = prob_hi * dx[0] / drdxfac;
                    }

                    dt_grid = amrex::min(1.1 * dt_grid,
                                         cfl * dr_fine_loc / (maxw0 + SMALL));
#endif
                }

                dt_lev = amrex::min(dt_lev, dt_grid);
                umax_lev = amrex::max(umax_lev, umax_grid);
            }  // end openmp
        }

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = amrex::max(umax, umax_lev);

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

    // set rel_eps in fortran module
    umax *= 1.e-8;
    rel_eps = umax;
}

void Maestro::FirstDt() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FirstDt()", FirstDt);

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    Vector<Real> w0_force_dummy((base_geom.max_radial_level + 1) *
                                base_geom.nr_fine);
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
#if (AMREX_SPACEDIM >= 2)
        umac_dummy[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev],
                                  1, 1);
        umac_dummy[lev][1].setVal(0.);
#endif
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

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force, umac_dummy, sold, rho0_old, grav_cell_old,
                 w0_force_cart_dummy, do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
#endif

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old, p0_cart, 0, 0, bcs_f, 0);
    Put1dArrayOnCart(gamma1bar_old, gamma1bar_cart, 0, 0, bcs_f, 0);

    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {
        Real dt_lev = 1.0e20;
        Real umax_lev = 0.;

        // get references to the MultiFabs at level lev
        MultiFab& uold_mf = uold[lev];
        MultiFab& sold_mf = sold[lev];
        MultiFab& vel_force_mf = vel_force[lev];
        MultiFab& S_cc_old_mf = S_cc_old[lev];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& gp0_cart_mf = gp0_cart[lev];
#endif
        const MultiFab& p0_mf = p0_cart[lev];
        const MultiFab& gamma1bar_mf = gamma1bar_cart[lev];

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min : dt_lev) reduction(max : umax_lev)
#endif
        for (MFIter mfi(sold_mf); mfi.isValid(); ++mfi) {
            Real dt_grid = initial_dt;
            Real umax_grid = 0.;

            dt_lev = amrex::min(dt_lev, dt_grid);
            umax_lev = amrex::max(umax_lev, umax_grid);
        }

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = amrex::max(umax, umax_lev);

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

    if (fixed_dt != -1.0) {
        // fixed dt
        dt = fixed_dt;
        if (maestro_verbose > 0) {
            Print() << "Setting fixed dt = " << dt << std::endl;
        }
    }

    // set rel_eps in fortran module
    umax *= 1.e-8;
    rel_eps = umax;
}
