
#include <Maestro.H>
#include <Maestro_F.H>

#ifdef AMREX_USE_CUDA
#include <AMReX_Arena.H>
#include <cuda_runtime_api.h>
#endif

using namespace amrex;

void Maestro::EstDt() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()", EstDt);

    dt = 1.e20;

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

    // face-centered
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > w0mac(finest_level + 1);

#if (AMREX_SPACEDIM == 3)
    if (spherical) {
        // initialize
        for (int lev = 0; lev <= finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev],
                                 1, 1);
            w0mac[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev],
                                 1, 1);
            w0mac[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev],
                                 1, 1);
        }

        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }
        }

        if (evolve_base_state &&
            (!use_exact_base_state && !average_base_state)) {
            MakeW0mac(w0mac);
        }
    }
#endif

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force, umac_dummy, sold, rho0_old, grav_cell_old,
                 w0_force_cart_dummy,
#ifdef ROTATION
                 w0mac, false,
#endif
                 do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    BaseState<Real> gp0(base_geom.max_radial_level + 1, base_geom.nr_fine + 1);
    gp0.setVal(0.);

    // divU constraint
    if (spherical) {
        EstDt_Divu(gp0, p0_old, gamma1bar_old);
    }

    Put1dArrayOnCart(gp0, gp0_cart, true, true, bcs_f, 0);
#endif

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old, p0_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(gamma1bar_old, gamma1bar_cart, false, false, bcs_f, 0);

    Real umax = 0.;

    Real dt_lev = 1.e50;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // create a MultiFab which will hold values for reduction
        // over
        MultiFab tmp(grids[lev], dmap[lev], 5, 0);

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min : dt_lev) reduction(max : umax_lev)
#endif
        {
            dt_lev = 1.e50;
            Real dt_grid = 1.e50;
            Real umax_grid = 0.;

            for (MFIter mfi(uold[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();

                const Array4<const Real> scal_arr = sold[lev].array(mfi);
                const Array4<const Real> u = uold[lev].array(mfi);
                const Array4<const Real> S_cc_arr = S_cc_old[lev].array(mfi);
                const Array4<const Real> dSdt_arr = dSdt[lev].array(mfi);
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
                const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
                const Array4<const Real> gamma1bar_arr =
                    gamma1bar_cart[lev].array(mfi);

                const Array4<Real> spd = tmp.array(mfi);

                if (!spherical) {
                    const Real rho_min = 1.e-20;
                    Real dt_temp = 1.e99;
                    const Real eps = 1.e-8;

                    Real spdx =
                        uold[lev][mfi].maxabs<RunOn::Device>(tileBox, 0);
                    tmp[mfi].setVal<RunOn::Device>(0.0, tileBox, 0, 1);
#if (AMREX_SPACEDIM == 2)
                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            spd(i, j, k, 0) =
                                u(i, j, k, 1) + 0.5 * (w0_arr(i, j, k, 1) +
                                                       w0_arr(i, j + 1, k, 1));
                        });
                    Real spdy = tmp[mfi].maxabs<RunOn::Device>(tileBox, 0);
                    Real spdz = 0.0;
#else
                    Real spdy =
                        uold[lev][mfi].maxabs<RunOn::Device>(tileBox, 1);

                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            spd(i, j, k, 0) =
                                u(i, j, k, 2) + 0.5 * (w0_arr(i, j, k, 2) +
                                                       w0_arr(i, j, k + 1, 2));
                        });
                    Real spdz = tmp[mfi].maxabs<RunOn::Device>(tileBox, 0);
#endif
                    Real spdr = w0_cart[lev][mfi].maxabs<RunOn::Device>(
                        tileBox, AMREX_SPACEDIM - 1);

                    umax_grid = amrex::max(
                        umax_grid,
                        amrex::max(spdx,
                                   amrex::max(spdy, amrex::max(spdz, spdr))));

                    if (spdx > eps) {
                        dt_temp = amrex::min(dt_temp, dx[0] / spdx);
                    }
                    if (spdy > eps) {
                        dt_temp = amrex::min(dt_temp, dx[1] / spdy);
                    }
#if AMREX_SPACEDIM == 3
                    if (spdz > eps) {
                        dt_temp = amrex::min(dt_temp, dx[2] / spdz);
                    }
#endif
                    if (spdr > eps) {
                        dt_temp =
                            amrex::min(dt_temp, dx[AMREX_SPACEDIM - 1] / spdr);
                    }

                    dt_temp *= cfl;

                    // Limit dt based on forcing terms
                    Real fx =
                        vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 0);
                    Real fy =
                        vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 1);
#if (AMREX_SPACEDIM == 3)
                    Real fz =
                        vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 2);
#endif

                    if (fx > eps) {
                        dt_temp =
                            amrex::min(dt_temp, std::sqrt(2.0 * dx[0] / fx));
                    }
                    if (fy > eps) {
                        dt_temp =
                            amrex::min(dt_temp, std::sqrt(2.0 * dx[1] / fy));
                    }
#if (AMREX_SPACEDIM == 3)
                    if (fz > eps) {
                        dt_temp =
                            amrex::min(dt_temp, std::sqrt(2.0 * dx[2] / fz));
                    }
#endif

                    dt_grid = amrex::min(dt_grid, dt_temp);

                    const auto nr_lev = base_geom.nr(lev);

                    tmp[mfi].setVal<RunOn::Device>(1.e99, tileBox, 1, 1);

                    // divU constraint
                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        Real gradp0 = 0.0;
#if (AMREX_SPACEDIM == 2)
                        if (j == 0) {
                            gradp0 =
                                (p0_arr(i, j + 1, k) - p0_arr(i, j, k)) / dx[1];
                        } else if (j == nr_lev - 1) {
                            gradp0 =
                                (p0_arr(i, j, k) - p0_arr(i, j - 1, k)) / dx[1];
                        } else {
                            gradp0 =
                                0.5 *
                                (p0_arr(i, j + 1, k) - p0_arr(i, j - 1, k)) /
                                dx[1];
                        }
#else 
                        if (k == 0) {
                            gradp0 = (p0_arr(i,j,k+1) - p0_arr(i,j,k)) / dx[2];
                        } else if (k == nr_lev-1) {
                            gradp0 = (p0_arr(i,j,k) - p0_arr(i,j,k-1)) / dx[2];
                        } else {
                            gradp0 = 0.5 * (p0_arr(i,j,k+1) - p0_arr(i,j,k-1)) / dx[2];
                        }
#endif
                        Real denom =
                            S_cc_arr(i, j, k) -
                            u(i, j, k, AMREX_SPACEDIM - 1) * gradp0 /
                                (gamma1bar_arr(i, j, k) * p0_arr(i, j, k));

                        if (denom > 0.0 &&
                            rho_min / scal_arr(i, j, k, Rho) < 1.0) {
                            spd(i, j, k, 1) =
                                0.4 * (1.0 - rho_min / scal_arr(i, j, k, Rho)) /
                                denom;
                        }
                    });

                    dt_temp = tmp[mfi].min<RunOn::Device>(tileBox, 1);

                    dt_grid = amrex::min(dt_grid, dt_temp);

                    tmp[mfi].setVal<RunOn::Device>(1.e99, tileBox, 2, 1);

                    // An additional dS/dt timestep constraint originally
                    // used in nova
                    // solve the quadratic equation
                    // (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
                    // which is equivalent to
                    // (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
                    // which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        if (dSdt_arr(i, j, k) > 1.e-20) {
                            auto a = 0.5 * scal_arr(i, j, k, Rho) *
                                     dSdt_arr(i, j, k);
                            auto b = scal_arr(i, j, k, Rho) * S_cc_arr(i, j, k);
                            auto c = rho_min - scal_arr(i, j, k, Rho);
                            spd(i, j, k, 2) =
                                0.4 * 2.0 * c /
                                (-b - std::sqrt(b * b - 4.0 * a * c));
                        }
                    });

                    dt_temp = tmp[mfi].min<RunOn::Device>(tileBox, 2);

                    dt_grid = amrex::min(dt_grid, dt_temp);
                } else {
#if (AMREX_SPACEDIM == 3)

                    const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                    const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                    const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
                    const Array4<const Real> gp0_arr = gp0_cart[lev].array(mfi);

                    const Real rho_min = 1.e-20;
                    Real dt_temp = 1.e99;
                    const Real eps = 1.e-8;

                    tmp[mfi].setVal<RunOn::Device>(0.0, tileBox, 0, 3);

                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            spd(i, j, k, 0) =
                                u(i, j, k, 0) +
                                0.5 * (w0macx(i, j, k) + w0macx(i + 1, j, k));
                        });
                    Real spdx = tmp[mfi].maxabs<RunOn::Device>(tileBox, 0);

                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            spd(i, j, k, 1) =
                                u(i, j, k, 1) +
                                0.5 * (w0macy(i, j, k) + w0macy(i, j + 1, k));
                        });
                    Real spdy = tmp[mfi].maxabs<RunOn::Device>(tileBox, 1);

                    ParallelFor(
                        tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            spd(i, j, k, 2) =
                                u(i, j, k, 2) +
                                0.5 * (w0macz(i, j, k) + w0macz(i, j, k + 1));
                        });
                    Real spdz = tmp[mfi].maxabs<RunOn::Device>(tileBox, 2);
                    Real spdr =
                        w0_cart[lev][mfi].maxabs<RunOn::Device>(tileBox, 0);

                    umax_grid = amrex::max(
                        umax_grid,
                        amrex::max(spdx,
                                   amrex::max(spdy, amrex::max(spdz, spdr))));

                    if (spdx > eps) {
                        dt_temp = amrex::min(dt_temp, dx[0] / spdx);
                    }
                    if (spdy > eps) {
                        dt_temp = amrex::min(dt_temp, dx[1] / spdy);
                    }
                    if (spdz > eps) {
                        dt_temp = amrex::min(dt_temp, dx[2] / spdz);
                    }
                    if (spdr > eps) {
                        dt_temp = amrex::min(dt_temp, base_geom.dr(0) / spdr);
                    }

                    dt_temp *= cfl;

                    // Limit dt based on forcing terms
                    Real fx =
                        vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 0);
                    Real fy =
                        vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 1);
                    Real fz =
                        vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 2);

                    if (fx > eps) {
                        dt_temp =
                            amrex::min(dt_temp, std::sqrt(2.0 * dx[0] / fx));
                    }
                    if (fy > eps) {
                        dt_temp =
                            amrex::min(dt_temp, std::sqrt(2.0 * dx[1] / fy));
                    }
                    if (fz > eps) {
                        dt_temp =
                            amrex::min(dt_temp, std::sqrt(2.0 * dx[2] / fz));
                    }

                    dt_grid = amrex::min(dt_grid, dt_temp);

                    tmp[mfi].setVal<RunOn::Device>(1.e50, tileBox, 3, 2);

                    // divU constraint
                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        Real gp_dot_u = 0.0;
                        for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                            gp_dot_u += u(i, j, k, n) * gp0_arr(i, j, k, n);
                        }

                        Real denom = S_cc_arr(i, j, k) - gp_dot_u;

                        if (denom > 0.0) {
                            spd(i, j, k, 3) =
                                0.4 * (1.0 - rho_min / scal_arr(i, j, k, Rho)) /
                                denom;
                        }
                    });

                    dt_temp = tmp[mfi].min<RunOn::Device>(tileBox, 3);

                    dt_grid = amrex::min(dt_grid, dt_temp);

                    // An additional dS/dt timestep constraint originally
                    // used in nova
                    // solve the quadratic equation
                    // (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
                    // which is equivalent to
                    // (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
                    // which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
                    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
                        if (dSdt_arr(i, j, k) > 1.e-20) {
                            auto a = 0.5 * scal_arr(i, j, k, Rho) *
                                     dSdt_arr(i, j, k);
                            auto b = scal_arr(i, j, k, Rho) * S_cc_arr(i, j, k);
                            auto c = rho_min - scal_arr(i, j, k, Rho);
                            spd(i, j, k, 4) =
                                0.4 * 2.0 * c /
                                (-b - std::sqrt(b * b - 4.0 * a * c));
                        }
                    });

                    dt_temp = tmp[mfi].min<RunOn::Device>(tileBox, 4);

                    dt_grid = amrex::min(dt_grid, dt_temp);
#else
                    Abort("EstDt: Spherical is not valid for DIM < 3");
#endif
                }
            }
            dt_lev = amrex::min(dt_lev, dt_grid);
            umax_lev = amrex::max(umax_lev, umax_grid);
        }  //end openmp

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

#ifdef ROTATION
    // face-centered
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > w0mac(finest_level + 1);
    // initialize
    for (int lev = 0; lev <= finest_level; ++lev) {
        w0mac[lev][0].define(convert(grids[lev], nodal_flag_x), dmap[lev], 1,
                             1);
        w0mac[lev][1].define(convert(grids[lev], nodal_flag_y), dmap[lev], 1,
                             1);
        w0mac[lev][2].define(convert(grids[lev], nodal_flag_z), dmap[lev], 1,
                             1);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            w0mac[lev][idim].setVal(0.);
        }
    }

    if (evolve_base_state &&
        (use_exact_base_state == 0 && average_base_state == 0)) {
        MakeW0mac(w0mac);
    }
#endif

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force, umac_dummy, sold, rho0_old, grav_cell_old,
                 w0_force_cart_dummy,
#ifdef ROTATION
                 w0mac, false,
#endif
                 do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    BaseState<Real> gp0(base_geom.max_radial_level + 1, base_geom.nr_fine + 1);
    gp0.setVal(0.);

    // divU constraint
    if (use_divu_firstdt && spherical) {
        EstDt_Divu(gp0, p0_old, gamma1bar_old);
    }

    Put1dArrayOnCart(gp0, gp0_cart, true, true, bcs_f, 0);
#endif

    Vector<MultiFab> p0_cart(finest_level + 1);
    Vector<MultiFab> gamma1bar_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old, p0_cart, false, false, bcs_f, 0);
    Put1dArrayOnCart(gamma1bar_old, gamma1bar_cart, false, false, bcs_f, 0);

    Real umax = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {
        Real dt_lev = 1.e99;
        Real umax_lev = 0.;

        // create a MultiFab which will hold values for reduction
        // over
        MultiFab tmp(grids[lev], dmap[lev], 1, 0);

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min : dt_lev) reduction(max : umax_lev)
#endif
        {
            Real dt_grid = 1.e50;
            Real umax_grid = 0.;

            for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();

                const Array4<const Real> scal_arr = sold[lev].array(mfi);
                const Array4<const Real> u = uold[lev].array(mfi);
                const Array4<const Real> S_cc_arr = S_cc_old[lev].array(mfi);
                const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
                const Array4<const Real> gamma1bar_arr =
                    gamma1bar_cart[lev].array(mfi);

                const Real eps = 1.e-8;
                const Real rho_min = 1.e-20;

                FArrayBox spd(tileBox);
                Elixir e_s = spd.elixir();

                const Array4<Real> spd_arr = spd.array();

                spd.setVal<RunOn::Device>(0.0, tileBox, 0, 1);

                ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    eos_t eos_state;

                    // compute the sound speed from rho and temp
                    eos_state.rho = scal_arr(i, j, k, Rho);
                    eos_state.T = scal_arr(i, j, k, Temp);
                    for (auto comp = 0; comp < NumSpec; ++comp) {
                        eos_state.xn[comp] =
                            scal_arr(i, j, k, FirstSpec + comp) /
                            scal_arr(i, j, k, Rho);
                    }
#if NAUX_NET > 0
                    for (auto comp = 0; comp < NumAux; ++comp) {
                        eos_state.aux[comp] =
                            scal_arr(i, j, k, FirstAux + comp) /
                            scal_arr(i, j, k, Rho);
                    }
#endif

                    // dens, temp, and xmass are inputs
                    eos(eos_input_rt, eos_state);

                    spd_arr(i, j, k) = eos_state.cs;
                });

                Real ux = uold[lev][mfi].maxabs<RunOn::Device>(tileBox, 0);
                Real uy = uold[lev][mfi].maxabs<RunOn::Device>(tileBox, 1);
                Real uz =
                    AMREX_SPACEDIM == 2
                        ? 0.1 * uy
                        : uold[lev][mfi].maxabs<RunOn::Device>(tileBox, 2);

                umax_grid =
                    amrex::max(umax_grid, amrex::max(ux, amrex::max(uy, uz)));

                ux /= dx[0];
                Real spdx = spd.max<RunOn::Device>(tileBox, 0) / dx[0];
                Real pforcex =
                    vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 0);
                uy /= dx[1];
                Real spdy = spd.max<RunOn::Device>(tileBox, 0) / dx[1];
                Real pforcey =
                    vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 1);
                Real spdz = spdy * 0.1;  // for 2d make sure this is < spdy
                uz /= AMREX_SPACEDIM == 2 ? dx[1] : dx[2];
#if (AMREX_SPACEDIM == 3)
                spdz = spd.max<RunOn::Device>(tileBox, 0) / dx[2];
                Real pforcez =
                    vel_force[lev][mfi].maxabs<RunOn::Device>(tileBox, 2);
#endif

                // use advective constraint unless velocities are zero everywhere
                // in which case we use the sound speed
                if (ux != 0.0 || uy != 0.0 || uz != 0.0) {
                    dt_grid = amrex::min(
                        dt_grid, cfl / amrex::max(ux, amrex::max(uy, uz)));
                } else if (spdx != 0.0 && spdy != 0.0 && spdz != 0.0) {
                    dt_grid = amrex::min(
                        dt_grid,
                        cfl / amrex::max(spdx, amrex::max(spdy, spdz)));
                }

                // sound speed constraint
                if (use_soundspeed_firstdt) {
                    Real dt_sound = 1.e99;
                    if (spdx == 0.0 && spdy == 0.0 && spdz == 0.0) {
                        dt_sound = 1.e99;
                    } else {
                        dt_sound =
                            cfl / amrex::max(spdx, amrex::max(spdy, spdz));
                    }
                    dt_grid = amrex::min(dt_grid, dt_sound);
                }

                // force constraints
                if (pforcex > eps) {
                    dt_grid =
                        amrex::min(dt_grid, std::sqrt(2.0 * dx[0] / pforcex));
                }
                if (pforcey > eps) {
                    dt_grid =
                        amrex::min(dt_grid, std::sqrt(2.0 * dx[1] / pforcey));
                }
#if (AMREX_SPACEDIM == 3)
                if (pforcez > eps) {
                    dt_grid =
                        amrex::min(dt_grid, std::sqrt(2.0 * dx[2] / pforcez));
                }
#endif

                // divU constraint
                if (use_divu_firstdt) {
                    Real dt_divu = 1.e99;

                    const Array4<Real> tmp_arr = tmp.array(mfi);
                    tmp[mfi].setVal<RunOn::Device>(1.e50, tileBox, 0, 1);

                    if (!spherical) {
                        const auto nr_lev = base_geom.nr(lev);

                        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                                  int k) {
                            Real gradp0 = 0.0;
#if (AMREX_SPACEDIM == 2)
                            if (j == 0) {
                                gradp0 =
                                    (p0_arr(i, j + 1, k) - p0_arr(i, j, k)) /
                                    dx[1];
                            } else if (j == nr_lev - 1) {
                                gradp0 =
                                    (p0_arr(i, j, k) - p0_arr(i, j - 1, k)) /
                                    dx[1];
                            } else {
                                gradp0 = 0.5 *
                                         (p0_arr(i, j + 1, k) -
                                          p0_arr(i, j - 1, k)) /
                                         dx[1];
                            }
#else 
                            if (k == 0) {
                                gradp0 = (p0_arr(i,j,k+1) - p0_arr(i,j,k)) / dx[2];
                            } else if (k == nr_lev-1) {
                                gradp0 = (p0_arr(i,j,k) - p0_arr(i,j,k-1)) / dx[2];
                            } else {
                                gradp0 = 0.5*(p0_arr(i,j,k+1) - p0_arr(i,j,k-1)) / dx[2];
                            }
#endif

                            Real denom =
                                S_cc_arr(i, j, k) -
                                u(i, j, k, AMREX_SPACEDIM - 1) * gradp0 /
                                    (gamma1bar_arr(i, j, k) * p0_arr(i, j, k));

                            if (denom > 0.0) {
                                tmp_arr(i, j, k) =
                                    0.4 *
                                    (1.0 - rho_min / scal_arr(i, j, k, Rho)) /
                                    denom;
                            }
                        });

                        dt_divu = tmp[mfi].min<RunOn::Device>(tileBox, 0);
                    } else {
#if (AMREX_SPACEDIM == 3)
                        const Array4<const Real> gp0_arr =
                            gp0_cart[lev].array(mfi);

                        ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                                  int k) {
                            Real gp_dot_u = 0.0;

                            for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                                gp_dot_u += gp0_arr(i, j, k, n) * u(i, j, k, n);
                            }

                            Real denom = S_cc_arr(i, j, k) - gp_dot_u;

                            if (denom > 0.0 && scal_arr(i, j, k, Rho) > 0.0) {
                                if (rho_min / scal_arr(i, j, k, Rho) < 1.0) {
                                    tmp_arr(i, j, k) =
                                        0.4 *
                                        (1.0 -
                                         rho_min / scal_arr(i, j, k, Rho)) /
                                        denom;
                                }
                            }
                        });
                        dt_divu = tmp[mfi].min<RunOn::Device>(tileBox, 0);
#endif
                    }
                    dt_grid = amrex::min(dt_grid, dt_divu);
                }
            }
            dt_lev = amrex::min(dt_lev, dt_grid);
            umax_lev = amrex::max(umax_lev, umax_grid);
        }  //end openmp

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

void Maestro::EstDt_Divu(BaseState<Real>& gp0, const BaseState<Real>& p0,
                         const BaseState<Real>& gamma1bar) {
    auto gp0_arr = gp0.array();
    const auto p0_arr = p0.const_array();
    const auto gamma1bar_arr = gamma1bar.const_array();
    const auto& r_cc_loc = base_geom.r_cc_loc;

    // spherical divU constraint
    if (use_exact_base_state) {
        ParallelFor(base_geom.nr_fine - 2, [=] AMREX_GPU_DEVICE(int i) {
            int r = i + 1;

            Real gamma1bar_p_avg =
                0.5 * (gamma1bar_arr(0, r) * p0_arr(0, r) +
                       gamma1bar_arr(0, r - 1) * p0_arr(0, r - 1));

            gp0_arr(0, r) = (p0_arr(0, r) - p0_arr(0, r - 1)) /
                            (r_cc_loc(0, r) - r_cc_loc(0, r - 1)) /
                            gamma1bar_p_avg;
        });
    } else {
        const auto dr0 = base_geom.dr(0);
        ParallelFor(base_geom.nr_fine - 2, [=] AMREX_GPU_DEVICE(int i) {
            int r = i + 1;

            Real gamma1bar_p_avg =
                0.5 * (gamma1bar_arr(0, r) * p0_arr(0, r) +
                       gamma1bar_arr(0, r - 1) * p0_arr(0, r - 1));

            gp0_arr(0, r) =
                (p0_arr(0, r) - p0_arr(0, r - 1)) / dr0 / gamma1bar_p_avg;
        });
    }
    Gpu::synchronize();

    gp0_arr(0, base_geom.nr_fine) = gp0_arr(0, base_geom.nr_fine - 1);
    gp0_arr(0, 0) = gp0_arr(0, 1);
}
