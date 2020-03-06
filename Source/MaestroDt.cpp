
#include <Maestro.H>
#include <Maestro_F.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <AMReX_Arena.H>
#endif

using namespace amrex;

void
Maestro::EstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EstDt()", EstDt);

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

    if (spherical == 1) {
        // initialize
        for (int lev=0; lev<=finest_level; ++lev) {
            w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
#if (AMREX_SPACEDIM == 3)
            w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
#endif
        }

        for (int lev=0; lev<=finest_level; ++lev) {
            for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }
        }

        if (evolve_base_state && (use_exact_base_state == 0 && average_base_state == 0)) {
            MakeW0mac(w0mac);
        }
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_cart_dummy,do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    RealVector gp0( (max_radial_level+1)*(nr_fine+1) );
    gp0.shrink_to_fit();
    std::fill(gp0.begin(),gp0.end(), 0.);

    // divU constraint
    if (spherical == 1) {
        EstDt_Divu(gp0, p0_old, gamma1bar_old);
    }
      
    Put1dArrayOnCart (gp0,gp0_cart,1,1,bcs_f,0);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old,p0_cart,0,0,bcs_f,0);
    Put1dArrayOnCart(gamma1bar_old,gamma1bar_cart,0,0,bcs_f,0);

    Real umax = 0.;

    Real dt_lev = 1.e99;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
        {

            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            for (MFIter mfi(uold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();

                const Array4<const Real> scal_arr = sold[lev].array(mfi);
                const Array4<const Real> u = uold[lev].array(mfi);
                const Array4<const Real> force = vel_force[lev].array(mfi);
                const Array4<const Real> S_cc_arr = S_cc_old[lev].array(mfi);
                const Array4<const Real> dSdt_arr = dSdt[lev].array(mfi);
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
                const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
                const Array4<const Real> gamma1bar_arr = gamma1bar_cart[lev].array(mfi);

                if (spherical == 0) {

                    const Real rho_min = 1.e-20;
                    Real dt_temp = 1.e99;
                    const Real eps = 1.e-8;

                    Real spdx = 0.0;
                    Real spdy = 0.0;
                    Real spdz = 0.0;
                    Real spdr = 0.0;

                    // Limit dt based on velocity terms.
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                        spdx = amrex::max(spdx, fabs(u(i,j,k,0)));
#if (AMREX_SPACEDIM == 2)
                        spdy = amrex::max(spdy, fabs(u(i,j,k,1) + 0.5 * (w0_arr(i,j,k,1)+w0_arr(i,j+1,k,1))));
#else   
                        spdy = amrex::max(spdy, fabs(u(i,j,k,1)));
                        spdz = amrex::max(spdz, fabs(u(i,j,k,2) + 0.5 * (w0_arr(i,j,k,2)+w0_arr(i,j,k+1,2))));
#endif
                        spdr = amrex::max(spdr, fabs(w0_arr(i,j,k,AMREX_SPACEDIM-1)));
                    });

                    umax_grid = amrex::max(umax_grid, std::max(spdx,std::max(spdy,std::max(spdz,spdr))));

                    if (spdx > eps) dt_temp = amrex::min(dt_temp, dx[0] / spdx);
                    if (spdy > eps) dt_temp = amrex::min(dt_temp, dx[1] / spdy);
                    if (spdz > eps) dt_temp = amrex::min(dt_temp, dx[2] / spdz);
                    if (spdr > eps) dt_temp = amrex::min(dt_temp, dx[AMREX_SPACEDIM-1] / spdr);

                    dt_temp *= cfl;

                    // Limit dt based on forcing terms
                    Real fx = 0.0;
                    Real fy = 0.0;
                    Real fz = 0.0;

                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                        fx = amrex::max(fx, fabs(force(i,j,k,0)));
                        fy = amrex::max(fy, fabs(force(i,j,k,1)));
#if (AMREX_SPACEDIM == 3)
                        fz = amrex::max(fz, fabs(force(i,j,k,2)));
#endif
                    });

                    if (fx > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[0] / fx));
                    if (fy > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[1] / fy));
#if (AMREX_SPACEDIM == 3)
                    if (fz > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[2] / fz));
#endif

                    const auto nr_lev = nr[lev];

                    // divU constraint
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                        Real gradp0 = 0.0;
#if (AMREX_SPACEDIM == 2)
                        if (j == 0) {
                            gradp0 = (p0_arr(i,j+1,k) - p0_arr(i,j,k)) / dx[1];
                        } else if (j == nr_lev-1) {
                            gradp0 = (p0_arr(i,j,k) - p0_arr(i,j-1,k)) / dx[1];
                        } else {
                            gradp0 = 0.5 * (p0_arr(i,j+1,k) - p0_arr(i,j-1,k)) / dx[1];
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
                        Real denom = S_cc_arr(i,j,k) - u(i,j,k,AMREX_SPACEDIM-1) * gradp0 / (gamma1bar_arr(i,j,k)*p0_arr(i,j,k));

                        if (denom > 0.0) {
                            dt_temp = amrex::min(dt_temp, 0.4*(1.0 - rho_min / scal_arr(i,j,k,Rho)) / denom);
                        }
                    });

                    // An additional dS/dt timestep constraint originally
                    // used in nova
                    // solve the quadratic equation
                    // (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
                    // which is equivalent to
                    // (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
                    // which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        if (dSdt_arr(i,j,k) > 1.e-20) {
                            auto a = 0.5 * scal_arr(i,j,k,Rho) * dSdt_arr(i,j,k);
                            auto b = scal_arr(i,j,k,Rho) * S_cc_arr(i,j,k);
                            auto c = rho_min - scal_arr(i,j,k,Rho);
                            dt_temp = amrex::min(dt_temp, 0.4*2.0*c / (-b-std::sqrt(b*b-4.0*a*c)));
                        }
                    });

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

                    Real spdx = 0.0;
                    Real spdy = 0.0;
                    Real spdz = 0.0;
                    Real spdr = 0.0;

                    // Limit dt based on velocity terms.
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                        spdx = amrex::max(spdx, fabs(u(i,j,k,0) + 0.5*(w0macx(i,j,k)+w0macx(i+1,j,k))));
                        spdy = amrex::max(spdy, fabs(u(i,j,k,1) + 0.5*(w0macy(i,j,k)+w0macy(i,j+1,k))));
                        spdz = amrex::max(spdz, fabs(u(i,j,k,2) + 0.5*(w0macz(i,j,k)+w0macz(i,j,k+1))));
                        spdr = amrex::max(spdr, fabs(w0_arr(i,j,k)));
                    });

                    umax_grid = amrex::max(umax_grid, std::max(spdx,std::max(spdy,std::max(spdz,spdr))));

                    if (spdx > eps) dt_temp = amrex::min(dt_temp, dx[0] / spdx);
                    if (spdy > eps) dt_temp = amrex::min(dt_temp, dx[1] / spdy);
                    if (spdz > eps) dt_temp = amrex::min(dt_temp, dx[2] / spdz);
                    if (spdr > eps) dt_temp = amrex::min(dt_temp, dr[0] / spdr);

                    dt_temp *= cfl;

                    // Limit dt based on forcing terms
                    Real fx = 0.0;
                    Real fy = 0.0;
                    Real fz = 0.0;

                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                        fx = amrex::max(fx, fabs(force(i,j,k,0)));
                        fy = amrex::max(fy, fabs(force(i,j,k,1)));
                        fz = amrex::max(fz, fabs(force(i,j,k,2)));
                    });

                    if (fx > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[0] / fx));
                    if (fy > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[1] / fy));
                    if (fz > eps) dt_temp = amrex::min(dt_temp, std::sqrt(2.0 * dx[2] / fz));
                    // divU constraint
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {

                        Real gp_dot_u = 0.0;
                        for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                            gp_dot_u += u(i,j,k,n) * gp0_arr(i,j,k,n);
                        }
                        
                        Real denom = S_cc_arr(i,j,k) - gp_dot_u;

                        if (denom > 0.0) {
                            dt_temp = amrex::min(dt_temp, 0.4*(1.0 - rho_min / scal_arr(i,j,k,Rho)) / denom);
                        }
                    });

                    // An additional dS/dt timestep constraint originally
                    // used in nova
                    // solve the quadratic equation
                    // (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
                    // which is equivalent to
                    // (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
                    // which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
                    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                        if (dSdt_arr(i,j,k) > 1.e-20) {
                            auto a = 0.5 * scal_arr(i,j,k,Rho) * dSdt_arr(i,j,k);
                            auto b = scal_arr(i,j,k,Rho) * S_cc_arr(i,j,k);
                            auto c = rho_min - scal_arr(i,j,k,Rho);
                            dt_temp = amrex::min(dt_temp, 0.4*2.0*c / (-b-std::sqrt(b*b-4.0*a*c)));
                        }
                    });

                    dt_grid = amrex::min(dt_grid, dt_temp);
#else
                    Abort("EstDt: Spherical is not valid for DIM < 3");
#endif
                }
            }

            dt_lev = std::min(dt_lev,dt_grid);
            umax_lev = std::max(umax_lev,umax_grid);

        } //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = std::max(umax,umax_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to estdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }     // end loop over levels

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
    c_rel_eps = umax;
    set_rel_eps(&umax);
}


void
Maestro::FirstDt ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FirstDt()", FirstDt);

    dt = 1.e20;

    // allocate a dummy w0_force and set equal to zero
    RealVector w0_force_dummy( (max_radial_level+1)*nr_fine );
    w0_force_dummy.shrink_to_fit();
    std::fill(w0_force_dummy.begin(),w0_force_dummy.end(), 0.);

    // build dummy w0_force_cart and set equal to zero
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], 3, 1);
        w0_force_cart_dummy[lev].setVal(0.);
    }

    // build a dummy umac and set equal to zero
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac_dummy(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        umac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        umac_dummy[lev][0].setVal(0.);
        umac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        umac_dummy[lev][1].setVal(0.);
#if (AMREX_SPACEDIM == 3)
        umac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
        umac_dummy[lev][2].setVal(0.);
#endif
    }

    // build and compute vel_force
    Vector<MultiFab> vel_force(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        vel_force[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        vel_force[lev].setVal(0.);
    }

    int do_add_utilde_force = 0;
    MakeVelForce(vel_force,umac_dummy,sold,rho0_old,grav_cell_old,
                 w0_force_cart_dummy,do_add_utilde_force);

#if (AMREX_SPACEDIM == 3)
    // build and initialize grad_p0 for spherical case
    Vector<MultiFab> gp0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        gp0_cart[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        gp0_cart[lev].setVal(0.);
    }
    RealVector gp0( (max_radial_level+1)*(nr_fine+1) );
    gp0.shrink_to_fit();
    std::fill(gp0.begin(),gp0.end(), 0.);

    // divU constraint
    if (use_divu_firstdt && spherical) {
        EstDt_Divu(gp0, p0_old, gamma1bar_old);
    }

    Put1dArrayOnCart (gp0,gp0_cart,1,1,bcs_f,0);
#endif

    Vector<MultiFab> p0_cart(finest_level+1);
    Vector<MultiFab> gamma1bar_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_old,p0_cart,0,0,bcs_f,0);
    Put1dArrayOnCart(gamma1bar_old,gamma1bar_cart,0,0,bcs_f,0);

    Real umax = 0.;

    Real dt_lev = 1.e99;
    Real umax_lev = 0.;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_lev) reduction(max:umax_lev)
#endif
        {
            Real dt_grid = 1.e99;
            Real umax_grid = 0.;

            for (MFIter mfi(sold[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();

                const Array4<const Real> scal_arr = sold[lev].array(mfi);
                const Array4<const Real> u = uold[lev].array(mfi);
                const Array4<const Real> force = vel_force[lev].array(mfi);
                const Array4<const Real> S_cc_arr = S_cc_old[lev].array(mfi);
                const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);
                const Array4<const Real> gamma1bar_arr = gamma1bar_cart[lev].array(mfi);

                const Real eps = 1.e-8;
                const Real rho_min = 1.e-20;

                Real spdx = 0.0;
                Real spdy = 0.0;
                Real spdz = 0.0;
                Real pforcex = 0.0;
                Real pforcey = 0.0;
                Real pforcez = 0.0;
                Real ux = 0.0;
                Real uy = 0.0;
                Real uz = 0.0;

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    eos_t eos_state;

                    // compute the sound speed from rho and temp
                    eos_state.rho = scal_arr(i,j,k,Rho);
                    eos_state.T = scal_arr(i,j,k,Temp);
                    for (auto comp = 0; comp < NumSpec; ++comp) {
                        eos_state.xn[comp] = scal_arr(i,j,k,FirstSpec+comp) / scal_arr(i,j,k,Rho);
                    }

                    // dens, temp, and xmass are inputs
                    eos(eos_input_rt, eos_state);

                    spdx = amrex::max(spdx, eos_state.cs);
                    ux = amrex::max(ux, fabs(u(i,j,k,0)));
                    pforcex = amrex::max(pforcex, fabs(force(i,j,k,0)));

                    spdy = amrex::max(spdy, eos_state.cs);
                    uy = amrex::max(uy, fabs(u(i,j,k,1)));
                    pforcey = amrex::max(pforcey, fabs(force(i,j,k,1)));
#if (AMREX_SPACEDIM == 3)
                    spdz = amrex::max(spdz, eos_state.cs);
                    uz = amrex::max(uz, fabs(u(i,j,k,2)));
                    pforcez = amrex::max(pforcez, fabs(force(i,j,k,2)));
#endif
                });

                umax_grid = amrex::max(umax_grid, std::max(ux,std::max(uy,uz)));

                ux /= dx[0];
                spdx /= dx[0];
                uy /= dx[1];
                spdy /= dx[1];
#if (AMREX_SPACEDIM == 3)
                uz /= dx[2];
                spdz /= dx[2];
#endif

                // use advective constraint unless velocities are zero everywhere
                // in which case we use the sound speed
                if (ux != 0.0 || uy != 0.0 || uz != 0.0) {
                    dt_grid = amrex::min(dt_grid, cfl / std::max(ux,std::max(uy,uz)));
                } else if (spdx != 0.0 && spdy != 0.0 && spdz != 0.0) {
                    dt_grid = amrex::min(dt_grid, cfl / std::max(spdx,std::max(spdy,spdz)));
                }

                // sound speed constraint
                if (use_soundspeed_firstdt) {
                    Real dt_sound = 0.0;

                    if (spdx == 0.0 && spdy == 0.0 && spdz == 0.0) {
                        dt_sound = 1.e99;
                    } else {
                        dt_sound = cfl / std::max(spdx,std::max(spdy,spdz));
                    }
                    dt_grid = amrex::min(dt_grid, dt_sound);
                }

                // force constraints
                if (pforcex > eps) dt_grid = amrex::min(dt_grid, std::sqrt(2.0 * dx[0] / pforcex));
                if (pforcey > eps) dt_grid = amrex::min(dt_grid, std::sqrt(2.0 * dx[1] / pforcey));
#if (AMREX_SPACEDIM == 3)
                if (pforcez > eps) dt_grid = amrex::min(dt_grid, std::sqrt(2.0 * dx[2] / pforcez));
#endif

                // divU constraint
                if (use_divu_firstdt) {
                    Real dt_divu = 1.e99;

                    if (!spherical) {
                        const auto nr_lev = nr[lev];

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                            Real gradp0 = 0.0;
#if (AMREX_SPACEDIM == 2)
                            if (j == 0) {
                                gradp0 = (p0_arr(i,j+1,k) - p0_arr(i,j,k)) / dx[1];
                            } else if (j == nr_lev-1) {
                                gradp0 = (p0_arr(i,j,k) - p0_arr(i,j-1,k)) / dx[1];
                            } else {
                                gradp0 = 0.5*(p0_arr(i,j+1,k) - p0_arr(i,j-1,k)) / dx[1];
                            }
#else 
                            if (k == 0) {
                                gradp0 = (p0_arr(i,j,k+1) - p0_arr(i,j,k)) / dx[2];
                            } else if (j == nr_lev-1) {
                                gradp0 = (p0_arr(i,j,k) - p0_arr(i,j,k-1)) / dx[2];
                            } else {
                                gradp0 = 0.5*(p0_arr(i,j,k+1) - p0_arr(i,j,k-1)) / dx[2];
                            }
#endif                        

                            Real denom = S_cc_arr(i,j,k) - u(i,j,k,AMREX_SPACEDIM-1) * gradp0 / (gamma1bar_arr(i,j,k) * p0_arr(i,j,k));

                            if (denom > 0.0) {
                                dt_divu = amrex::min(dt_divu, 0.4 * (1.0 - rho_min / scal_arr(i,j,k,Rho)) / denom);
                            }
                        });
                    } else {

                        const Array4<const Real> gp0_arr = gp0_cart[lev].array(mfi);

                        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                            Real gp_dot_u = 0.0;

                            for (auto n = 0; n < AMREX_SPACEDIM; ++n)    {
                                gp_dot_u += gp0_arr(i,j,k,n) * u(i,j,k,n);
                            }

                            Real denom = S_cc_arr(i,j,k) - gp_dot_u;

                            if (denom > 0.0) {
                                dt_divu = amrex::min(dt_divu, 0.4 * (1.0 - rho_min / scal_arr(i,j,k,Rho)) / denom);
                            }
                        });
                    }
                    dt_grid = amrex::min(dt_grid, dt_divu);
                }
            }

            dt_lev = std::min(dt_lev,dt_grid);
            umax_lev = std::max(umax_lev,umax_grid);
        } //end openmp

        // find the smallest dt over all processors
        ParallelDescriptor::ReduceRealMin(dt_lev);

        // find the largest umax over all processors
        ParallelDescriptor::ReduceRealMax(umax_lev);

        // update umax over all levels
        umax = std::max(umax,umax_lev);

        if (maestro_verbose > 0) {
            Print() << "Call to firstdt for level " << lev << " gives dt_lev = " << dt_lev << std::endl;
        }

        // multiply by init_shrink
        dt_lev *= init_shrink;

        if (maestro_verbose > 0) {
            Print() << "Multiplying dt_lev by init_shrink; dt_lev = " << dt_lev << std::endl;
        }

        // update dt over all levels
        dt = std::min(dt,dt_lev);

    }     // end loop over levels

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
    c_rel_eps = umax;
    set_rel_eps(&umax);
}


void
Maestro::EstDt_Divu(RealVector& gp0_vec, const RealVector& p0_vec, 
                    const RealVector& gamma1bar_vec)
{
    const int max_lev = max_radial_level + 1;

    Real * AMREX_RESTRICT gp0 = gp0_vec.dataPtr();
    const Real * AMREX_RESTRICT p0 = p0_vec.dataPtr();
    const Real * AMREX_RESTRICT gamma1bar = gamma1bar_vec.dataPtr();
    const Real * AMREX_RESTRICT r_cc_loc_p = r_cc_loc.dataPtr();
    const Real * AMREX_RESTRICT r_edge_loc_p = r_edge_loc.dataPtr();

    // spherical divU constraint
    if (use_exact_base_state) {
        AMREX_PARALLEL_FOR_1D(nr_fine-1, i, {
            int r = i + 1;

            Real gamma1bar_p_avg = 0.5 * (gamma1bar[max_lev*r]*p0[max_lev*r] + gamma1bar[max_lev*(r-1)]*p0[max_lev*(r-1)]);

            gp0[max_lev*r] = (p0[max_lev*r] - p0[max_lev*(r-1)]) / (r_cc_loc[max_lev*r] - r_cc_loc[max_lev*(r-1)])  / gamma1bar_p_avg;
        });
    } else {
        const auto dr0 = dr[0];
        AMREX_PARALLEL_FOR_1D(nr_fine-1, i, {
            int r = i + 1;

            Real gamma1bar_p_avg = 0.5 * (gamma1bar[max_lev*r]*p0[max_lev*r] + gamma1bar[max_lev*(r-1)]*p0[max_lev*(r-1)]);

            gp0[max_lev*r] = (p0[max_lev*r] - p0[max_lev*(r-1)]) / dr0 / gamma1bar_p_avg;
        });
    }
}